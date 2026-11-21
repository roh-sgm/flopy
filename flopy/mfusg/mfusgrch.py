"""
mfusgrch module.  Contains the MfUsgRch class. Note that the user can access
the MfUsgRch class as `flopy.mfusg.MfUsgRch`.

"""

import copy
import io

import numpy as np

from ..modflow.mfparbc import ModflowParBc as mfparbc
from ..modflow.mfrch import ModflowRch
from ..utils import Transient2d, Transient3d, Util2d, read1d
from ..utils.flopy_io import line_parse
from ..utils.utils_def import (
    get_pak_vals_shape,
    get_unitnumber_from_ext_unit_dict,
)
from .mfusg import MfUsg


def _is_int(value):
    return isinstance(value, (int, np.integer)) and not isinstance(value, bool)


def _period_flags(value, nper, name):
    """INRCHZONES/INSELEV values as {kper: int}; a list or array has one
    value per stress period."""
    if value is None:
        return {}
    if isinstance(value, dict):
        flags = dict(value)
    elif isinstance(value, (list, tuple, np.ndarray)):
        if len(value) != nper:
            raise ValueError(
                f"{name} has {len(value)} values for {nper} stress periods"
            )
        flags = dict(enumerate(value))
    else:
        raise ValueError(f"{name} must be a dict {{kper: int}} or one int per period")
    for kper, flag in flags.items():
        if not _is_int(kper) or not 0 <= kper < nper:
            raise ValueError(f"{name}: stress period {kper!r} not in 0..{nper - 1}")
        if not _is_int(flag):
            raise ValueError(f"{name}[{kper}] = {flag!r} is not an integer")
    return {int(k): int(v) for k, v in flags.items()}


def _array_keys(value, nper, name):
    """Check the stress-period keys of selev/iznrch as given: a dict or a
    Transient2d needs integer keys in 0..nper-1 (no float, bool or text keys,
    which Transient2d would convert); other values are period 1."""
    if isinstance(value, Transient2d):
        keys = value.transient_2ds
    elif isinstance(value, dict):
        keys = value
    else:
        return
    for kper in keys:
        if not _is_int(kper) or not 0 <= kper < nper:
            raise ValueError(
                f"{name}: stress period key {kper!r} is not an integer in 0..{nper - 1}"
            )


def _as_transient2d(model, value, dtype, name):
    """value as a new Transient2d; a given Transient2d is not changed
    (Transient2d(value=Transient2d) rebuilds the Util2d of its dict)."""
    if value is None:
        return None
    if isinstance(value, Transient2d):
        value = copy.copy(value)
        value.transient_2ds = dict(value.transient_2ds)
        shape = value.shape
    else:
        shape = get_pak_vals_shape(model, value)
    return Transient2d(model, shape, dtype, value, name=name)


def _rts_seepelev_options(nper, mxrtzones, seepelev, iznrch, selev, iniznrch, inselev):
    """Checks that need no array conversion; returns the INRCHZONES and
    INSELEV values as {kper: int}."""
    if not _is_int(mxrtzones) or mxrtzones < 0:
        raise ValueError(f"mxrtzones = {mxrtzones!r} must be an integer >= 0")
    if not _is_int(seepelev) or seepelev not in (0, 1):
        raise ValueError(f"seepelev = {seepelev!r} must be 0 or 1")
    izn = _period_flags(iniznrch, nper, "iniznrch")
    sel = _period_flags(inselev, nper, "inselev")
    if mxrtzones == 0 and (iznrch is not None or izn):
        raise ValueError("iznrch/iniznrch need the RTS option (mxrtzones > 0)")
    if not seepelev and (selev is not None or sel):
        raise ValueError("selev/inselev need the SEEPELEV option (seepelev=1)")
    _array_keys(iznrch, nper, "iznrch")
    _array_keys(selev, nper, "selev")
    return izn, sel


def _rts_seepelev_plan(
    nper, shape, mxrtzones, seepelev, iznrch, selev, iniznrch, inselev
):
    """Check the RTS zones (INRCHZONES, IZNRCH) and seepage elevations
    (INSELEV, SELEV) and return the [INRCHZONES, INSELEV] values of every
    stress period, None for a keyword that is not written.

    USG-T RP (gwf2rch8u1.f): INRCHZONES n > 0 reads IZNRCH and sets the
    number of RTS zones to n (the RTS header count is not used then); the
    RTS array is allocated with the period 1 count, so period 1 must read
    zones and later positive values must keep that count. n <= 0 or no
    keyword reuses the zones. INSELEV n > 0 reads SELEV and turns seepage
    on; 0 or no keyword turns it off; n < 0 keeps the previous state and
    SELEV. An array is read only with a value > 0, so a period with an
    array needs a value > 0 and the other way round. Without explicit
    values: a period with an array writes INRCHZONES with the period 1
    count (mxrtzones in period 1) / INSELEV 1, any other period writes no
    INRCHZONES / INSELEV -1."""
    izn, sel = _rts_seepelev_options(
        nper, mxrtzones, seepelev, iznrch, selev, iniznrch, inselev
    )
    defined = {}
    for name, t2d in (("iznrch", iznrch), ("selev", selev)):
        defined[name] = set()
        if t2d is None:
            continue
        for kper, u2d in t2d.transient_2ds.items():
            if tuple(u2d.shape) != tuple(shape):
                raise ValueError(
                    f"{name} period {kper + 1} has shape {u2d.shape}, RECH {shape}"
                )
            defined[name].add(kper)
    plan = []
    count = None  # number of RTS zones set in period 1
    for kper in range(nper):
        n = s = None
        if mxrtzones > 0:
            new = kper in defined["iznrch"]
            default = (mxrtzones if count is None else count) if new else None
            n = izn.get(kper, default)
            if n is not None and n > 0:
                if not new:
                    raise ValueError(
                        f"INRCHZONES {n} in period {kper + 1} needs iznrch data"
                    )
                if count is None:
                    count = n
                elif n != count:
                    raise ValueError(
                        f"INRCHZONES {n} in period {kper + 1}: the number of RTS "
                        f"zones is fixed in period 1 to {count}"
                    )
            elif new:
                raise ValueError(
                    f"iznrch is given for period {kper + 1} but INRCHZONES {n} "
                    "does not read it"
                )
            elif kper == 0:
                raise ValueError("RTS needs iznrch and INRCHZONES > 0 in period 1")
        if seepelev and (selev is not None or sel):
            new = kper in defined["selev"]
            s = sel.get(kper, 1 if new else -1)
            if s > 0 and not new:
                raise ValueError(f"INSELEV {s} in period {kper + 1} needs selev data")
            if s <= 0 and new:
                raise ValueError(
                    f"selev is given for period {kper + 1} but INSELEV {s} "
                    "does not read it"
                )
        plan.append([n, s])
    return plan


class MfUsgRch(ModflowRch):
    """
    MFUSG TRANSPORT Recharge Package Class.

    Parameters
    ----------
    model : model object
        The model object (of type :class:`flopy.modflow.mf.Modflow`) to which
        this package will be added.
    ipakcb : int, optional
        Toggles whether cell-by-cell budget data should be saved. If None or zero,
        budget data will not be saved (default is None).
    nrchop : int
        is the recharge option code.
        1: Recharge to top grid layer only
        2: Recharge to layer defined in irch
        3: Recharge to highest active cell (default is 3).
    rech : float or filename or ndarray or dict keyed on kper (zero-based)
        Recharge flux (default is 1.e-3, which is used for all stress periods)
    irch : int or filename or ndarray or dict keyed on kper (zero-based)
        Layer (for an unstructured grid) or node (for an unstructured grid) to
        which recharge is applied in each vertical column (only used when
        nrchop=2). Default is 0, which is used for all stress periods.
    extension : string
        Filename extension (default is 'rch')
    unitnumber : int
        File unit number (default is None).
    filenames : str or list of str
        Filenames to use for the package and the output files. If
        filenames=None the package name will be created using the model name
        and package extension and the cbc output name will be created using
        the model name and .cbc extension (for example, modflowtest.cbc),
        if ipakcb is a number greater than zero. If a single string is passed
        the package will be set to the string and cbc output names will be
        created using the model name and .cbc extension, if ipakcb is a
        number greater than zero. To define the names for all package files
        (input and output) the length of the list of strings should be 2.
        Default is None.
    inselev : dict {kper: int} or one int per stress period, optional
        INSELEV of each stress period (SEEPELEV option): > 0 reads selev of
        the period and turns seepage on, 0 turns it off, < 0 keeps the
        previous state and SELEV. Periods without a value write 1 when selev
        has an array for the period and -1 otherwise (default is None).
    iniznrch : dict {kper: int} or one int per stress period, optional
        INRCHZONES of each stress period (RTS option): > 0 reads iznrch of
        the period; the period 1 value is the number of RTS zones (it may
        differ from mxrtzones, which USG-T does not use then) and later
        values > 0 must keep it; <= 0 reuses the zones. Periods without a
        value write that number (mxrtzones in period 1) when iznrch has an
        array for the period and no keyword otherwise. Period 1 must read
        zones (default is None).

    Attributes
    ----------

    Methods
    -------

    See Also
    --------

    Notes
    -----
    Parameters are supported in Flopy only when reading in existing models.
    Parameter values are converted to native values in Flopy and the
    connection to "parameters" is thus nonexistent.

    Examples
    --------

    """

    def __init__(
        self,
        model,
        nrchop=3,
        ipakcb=None,
        rech=1e-3,
        irch=0,
        seepelev=0,
        mxrtzones=0,
        iconc=0,
        selev=None,
        iznrch=None,
        rchconc=None,
        unitnumber=None,
        filenames=None,
        inselev=None,
        iniznrch=None,
    ):
        """Package constructor."""
        msg = (
            "Model object must be of type flopy.mfusg.MfUsg\n"
            f"but received type: {type(model)}."
        )
        assert isinstance(model, MfUsg), msg

        # SELEV/IZNRCH are checked before the package is added to the model:
        # options and keys first, then the arrays converted without changing
        # the given objects
        nper = model.nrow_ncol_nlay_nper[3]
        _rts_seepelev_options(
            nper, mxrtzones, seepelev, iznrch, selev, iniznrch, inselev
        )
        selev = _as_transient2d(model, selev, np.float32, "rech_selev")
        iznrch = _as_transient2d(model, iznrch, np.int32, "rech_izn")
        _rts_seepelev_plan(
            nper,
            get_pak_vals_shape(model, rech),
            mxrtzones,
            seepelev,
            iznrch,
            selev,
            iniznrch,
            inselev,
        )

        # call base package constructor
        super().__init__(
            model,
            nrchop=nrchop,
            ipakcb=ipakcb,
            rech=rech,
            irch=irch,
            unitnumber=unitnumber,
            filenames=filenames,
        )

        self.mxrtzones = mxrtzones
        self.seepelev = seepelev
        self.iconc = iconc
        self.selev = selev
        self.iznrch = iznrch
        self.inselev = inselev
        self.iniznrch = iniznrch

        self.rchconc = None
        if rchconc is not None:
            mcomp = model.mcomp
            if model.iheat > 0:
                mcomp = model.mcomp + 1
            self.irchconc = [1] * mcomp
            rchconc_u3d_shape = (mcomp,) + get_pak_vals_shape(model, rech)
            self.rchconc = Transient3d(
                model, rchconc_u3d_shape, np.float32, rchconc, name="rech_conc"
            )

    def write_file(self, check=True, f=None):
        """
        Write the package file.

        Parameters
        ----------
        check : boolean
            Check package data for common errors. (default True)

        Returns
        -------
        None

        """
        nrow, ncol, nlay, nper = self.parent.nrow_ncol_nlay_nper
        # checked again (attributes may have changed) before any output
        plan = _rts_seepelev_plan(
            nper,
            self.rech.shape,
            self.mxrtzones,
            self.seepelev,
            self.iznrch,
            self.selev,
            self.iniznrch,
            self.inselev,
        )
        # allows turning off package checks when writing files at model level
        if check:
            self.check(
                f=f"{self.name[0]}.chk",
                verbose=self.parent.verbose,
                level=1,
            )
        # the text is built in memory; the file is opened at the end
        f_rch = io.StringIO()
        f_rch.write(f"{self.heading}\n")
        f_rch.write(f"{self.nrchop:10d}{self.ipakcb:10d}")
        if self.seepelev:
            f_rch.write(" SEEPELEV")
        if self.iconc:
            f_rch.write(" CONC")
        if self.mxrtzones:
            f_rch.write(f" RTS {self.mxrtzones:4d}")
        f_rch.write("\n")

        mcomp = self.parent.mcomp
        if self.parent.iheat > 0:
            mcomp = mcomp + 1

        if self.nrchop == 2:
            irch = {}
            for kper, u2d in self.irch.transient_2ds.items():
                irch[kper] = u2d.array + 1
            irch = Transient2d(
                self.parent,
                self.irch.shape,
                self.irch.dtype,
                irch,
                self.irch.name,
            )
            if not self.parent.structured:
                mxndrch = np.max(
                    [u2d.array.size for kper, u2d in self.irch.transient_2ds.items()]
                )
                f_rch.write(f"{mxndrch:10d}\n")

        if self.iconc:
            for icomp in range(mcomp):
                f_rch.write(f"{self.irchconc[icomp]:10.0f}")
            f_rch.write("\n")

        for kper in range(nper):
            inrech, file_entry_rech = self.rech.get_kper_entry(kper)
            if self.nrchop == 2:
                inirch, file_entry_irch = irch.get_kper_entry(kper)
                if not self.parent.structured:
                    inirch = self.rech[kper].array.size
            else:
                inirch = -1

            if self.nrchop == 2:
                f_rch.write(f"{inrech:10d}{inirch:10d} ")
            else:
                f_rch.write(f"{inrech:10d} ")

            inrchzones, inselev = plan[kper]
            if inrchzones is not None:
                f_rch.write(f" INRCHZONES {inrchzones}")
            if inselev is not None:
                f_rch.write(f" INSELEV {inselev}")
            if self.rchconc is not None:
                f_rch.write(" INCONC")
            f_rch.write(f" # Stress period {kper + 1}\n")

            if inrech >= 0:
                f_rch.write(file_entry_rech)
            if self.nrchop == 2:
                if inirch >= 0:
                    f_rch.write(file_entry_irch)

            # arrays only for INRCHZONES/INSELEV > 0 (RP reads them then)
            if inrchzones is not None and inrchzones > 0:
                f_rch.write(self.iznrch.transient_2ds[kper].get_file_entry())
            if inselev is not None and inselev > 0:
                f_rch.write(self.selev.transient_2ds[kper].get_file_entry())
            if self.rchconc is not None:
                inrech, file_entry_rchconc = self.rchconc.get_kper_entry(kper)
                f_rch.write(file_entry_rchconc)
        # a handle given by the caller is written and left open
        if f is not None:
            f.write(f_rch.getvalue())
        else:
            with open(self.fn_path, "w") as f_out:
                f_out.write(f_rch.getvalue())

    @classmethod
    def load(cls, f, model, nper=None, ext_unit_dict=None, check=True):
        """
        Load an existing package.

        Parameters
        ----------
        f : filename or file handle
            File to load.
        model : model object
            The model object (of type :class:`flopy.modflow.mf.Modflow`) to
            which this package will be added.
        nper : int
            The number of stress periods.  If nper is None, then nper will be
            obtained from the model object. (default is None).
        ext_unit_dict : dictionary, optional
            If the arrays in the file are specified using EXTERNAL,
            or older style array control records, then `f` should be a file
            handle.  In this case ext_unit_dict is required, which can be
            constructed using the function
            :class:`flopy.utils.mfreadnam.parsenamefile`.
        check : boolean
            Check package data for common errors. (default True)

        Returns
        -------
        rch : MfUsgRch object
            MfUsgRch object.

        Examples
        --------

        >>> import flopy
        >>> m = flopy.modflow.Modflow()
        >>> rch = flopy.modflow.MfUsgRch.load('test.rch', m)

        """
        if model.verbose:
            print("loading rch package file...")

        openfile = not hasattr(f, "read")
        if openfile:
            filename = f
            f = open(filename, "r")

        # dataset 0 -- header
        while True:
            line = f.readline()
            if line[0] != "#":
                break
        npar = 0

        if "parameter" in line.lower():
            raw = line.strip().split()
            npar = int(raw[1])
            if npar > 0:
                if model.verbose:
                    print(f"   Parameters detected. Number of parameters = {npar}")
            line = f.readline()
        # dataset 2
        t = line_parse(line)
        nrchop = int(t[0])
        ipakcb = int(t[1])

        # item 2 - options
        seepelev = 0
        if "SEEPELEV" in t:
            seepelev = 1

        iconc = 0
        if "CONC" in line:
            iconc = 1

        mxrtzones = 0
        if "RTS" in t:
            idx = t.index("RTS")
            mxrtzones = int(t[idx + 1])

        # dataset 2b for mfusg
        if not model.structured and nrchop == 2:
            line = f.readline()
            t = line_parse(line)
            mxndrch = int(t[0])

        mcomp = model.mcomp
        if model.iheat > 0:
            mcomp = model.mcomp + 1

        irchconc = np.empty((mcomp), dtype=np.int32)
        if iconc:
            if model.verbose:
                print(f"   loading IRCHCONC[{mcomp}] ...")
            irchconc = read1d(f, irchconc)

        if model.verbose:
            print(f"   irchconc: {irchconc}")

        # dataset 3 and 4 - parameters data
        pak_parms = None
        if npar > 0:
            pak_parms = mfparbc.loadarray(f, npar, model.verbose)

        if nper is None:
            nrow, ncol, nlay, nper = model.get_nrow_ncol_nlay_nper()
        else:
            nrow, ncol, nlay, _ = model.get_nrow_ncol_nlay_nper()

        u2d_shape = (nrow, ncol)

        # read data for every stress period
        current_rech = []
        rech = {}
        irch = None
        if nrchop == 2:
            irch = {}
        current_irch = []

        # arrays only of the periods that read them (flag > 0); flags as read,
        # INSELEV 0 for a period without the keyword (RP sets ISELEV = 0)
        iznrch = {}
        iniznrch = {}
        selev = {}
        inselev = {}

        rchconc = None
        if iconc:
            rchconc = {}
        inconc = [0] * nper
        current_rchconc = []

        for iper in range(nper):
            line = f.readline()
            t = line_parse(line)
            inrech = int(t[0])

            if nrchop == 2:
                inirch = int(t[1])
                if (not model.structured) and (inirch >= 0):
                    u2d_shape = (1, inirch)
            elif not model.structured:
                u2d_shape = (1, ncol[0])

            inselev[iper] = 0
            if "INSELEV" in t:
                idx = t.index("INSELEV")
                inselev[iper] = int(t[idx + 1])

            if "INRCHZONES" in t:
                idx = t.index("INRCHZONES")
                iniznrch[iper] = int(t[idx + 1])

            if "INCONC" in t:
                inconc[iper] = 1

            if inrech >= 0:
                if npar == 0:
                    if model.verbose:
                        print(f"   loading rech stress period {iper + 1:3d}...")
                    t = Util2d.load(
                        f,
                        model,
                        u2d_shape,
                        np.float32,
                        "rech",
                        ext_unit_dict,
                    )
                else:
                    parm_dict = {}
                    for ipar in range(inrech):
                        line = f.readline()
                        t = line.strip().split()
                        pname = t[0].lower()
                        try:
                            c = t[1].lower()
                            instance_dict = pak_parms.bc_parms[pname][1]
                            if c in instance_dict:
                                iname = c
                            else:
                                iname = "static"
                        except:
                            iname = "static"
                        parm_dict[pname] = iname
                    t = mfparbc.parameter_bcfill(model, u2d_shape, parm_dict, pak_parms)

                current_rech = t
            rech[iper] = current_rech

            if nrchop == 2:
                if inirch >= 0:
                    if model.verbose:
                        print(f"   loading irch stress period {iper + 1:3d}...")
                    t = Util2d.load(
                        f, model, u2d_shape, np.int32, "irch", ext_unit_dict
                    )
                    current_irch = Util2d(
                        model, u2d_shape, np.int32, t.array - 1, "irch"
                    )
                irch[iper] = current_irch

            # Item 9 for mfusg transport
            if mxrtzones:
                if iniznrch.get(iper, 0) > 0:
                    if model.verbose:
                        print(f"   loading iznrch stress period {iper + 1:3d}...")
                    iznrch[iper] = Util2d.load(
                        f, model, u2d_shape, np.int32, "iznrch", ext_unit_dict
                    )

            # Item 10 for mfusg transport
            if seepelev:
                if inselev[iper] > 0:
                    if model.verbose:
                        print(f"   loading selev stress period {iper + 1:3d}...")
                    selev[iper] = Util2d.load(
                        f, model, u2d_shape, np.float32, "selev", ext_unit_dict
                    )

            # Item 11 for mfusg transport
            if iconc:
                if inconc[iper]:
                    if model.verbose:
                        print(f"   loading rch conc stress period {iper + 1:3d}...")
                    if model.iheat > 0:
                        mcomp = model.mcomp + 1

                    current_rchconc = [0] * mcomp
                    for icomp in range(mcomp):
                        if irchconc[icomp]:
                            if model.verbose:
                                print(
                                    f"   loading rch conc stress period {iper + 1:3d}"
                                    f" component {icomp + 1}..."
                                )
                            current_rchconc[icomp] = Util2d.load(
                                f,
                                model,
                                u2d_shape,
                                np.float32,
                                "rchconc",
                                ext_unit_dict,
                            )
                rchconc[iper] = current_rchconc

        if openfile:
            f.close()

        unitnumber, filenames = get_unitnumber_from_ext_unit_dict(
            model, cls, ext_unit_dict, ipakcb
        )

        # create recharge package instance
        rch = cls(
            model,
            nrchop=nrchop,
            ipakcb=ipakcb,
            rech=rech,
            irch=irch,
            seepelev=seepelev,
            mxrtzones=mxrtzones,
            iconc=iconc,
            selev=selev or None,
            iznrch=iznrch or None,
            rchconc=rchconc,
            unitnumber=unitnumber,
            filenames=filenames,
            # without SEEPELEV/RTS, RP ignores INSELEV/INRCHZONES
            inselev=inselev if seepelev else None,
            iniznrch=iniznrch if mxrtzones else None,
        )
        if check:
            rch.check(
                f=f"{rch.name[0]}.chk",
                verbose=rch.parent.verbose,
                level=0,
            )
        return rch
