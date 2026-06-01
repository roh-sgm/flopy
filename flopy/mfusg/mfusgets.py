"""
mfusgets module.  Contains the MfUsgEts class. Note that the user can access
the MfUsgEts class as `flopy.mfusg.MfUsgEts`.

Implements the Segmented Evapotranspiration (ETS) Package for MODFLOW-USG-T,
as described in the GWT v2.7.0 IO documentation (section: Segmented
Evapotranspiration (ETS) Package Input Instructions).

This is a distinct package from the simpler EVT package (MfUsgEvt).  Key
differences from EVT:

  * NAM file type token is "ETS" (not "EVT").
  * Header format (item 2a): NETSOP IETSCB NPETS NETSEG IESFACTOR
      NETSEG    – number of piecewise-linear segments (>1 adds PXDP/PETM arrays)
      NPETS     – number of named ETS parameters
      IESFACTOR – transport solute removal flag
  * Per-stress-period header (item 5a, when NETSOP==2 or NETSEG>1):
        INETSS INETSR INETSX [INIETS [INSGDF]]
  * Per stress period, (NETSEG-1) pairs of segment arrays are read:
        PXDP – proportional extinction depth at segment intersection
        PETM – ET rate multiplier at segment intersection

References
----------
GWT v2.7.0 IO Doc, pp. 187+, "Segmented Evapotranspiration (ETS) Package".
"""

import numpy as np

from ..modflow.mfparbc import ModflowParBc as mfparbc
from ..pakbase import Package
from ..utils import Transient2d, Util2d
from ..utils.utils_def import get_pak_vals_shape
from ._usgt_parameters import (
    read_active_array_parameters,
    write_active_array_parameters,
    write_array_parameter_defs,
)


class MfUsgEts(Package):
    """
    MODFLOW-USG-T Segmented Evapotranspiration (ETS) Package Class.

    Parameters
    ----------
    model : MfUsg
        The model object to which this package will be added.
    netsop : int
        ET option code (NETSOP).
        1 = ET applied to top grid layer only (default).
        2 = ET applied to layer specified by IETS.
        3 = ET applied to highest active cell.
    ipakcb : int or None
        Unit number for cell-by-cell budget output (IETSCB).
        0 or None = no output.
    surf : float or array or dict keyed on kper
        ET surface elevation (ETSS). Default 0.0.
    evtr : float or array or dict keyed on kper
        Maximum ET flux (ETSR). Default 1e-3.
    exdp : float or array or dict keyed on kper
        ET extinction depth (ETSX). Default 1.0.
    ievt : int or array or dict keyed on kper
        Layer indicator (IETS). Used only when netsop==2. Default 1.
    netseg : int
        Number of piecewise-linear segments for the ET rate vs head
        function (NETSEG). 1 = simple (no PXDP/PETM). Default 1.
    npets : int
        Number of named ETS parameters (NPETS). Default 0.
    iesfactor : int
        Transport solute removal flag (IESFACTOR).
        1 = read ESFACTOR per component; 0 = solutes left behind. Default 0.
    esfactor : float or array of shape (mcomp,) or None
        Solute removal fraction per component (ESFACTOR(MCOMP)). Read only
        when transport is active (model.itrnsp != 0) and iesfactor == 1.
    pxdp : list of (float or array or dict keyed on kper) or None
        Proportional extinction depth at each of the (netseg-1) segment
        intersections. Length must equal netseg-1. Default 0.5 for each.
    petm : list of (float or array or dict keyed on kper) or None
        ET rate multiplier at each of the (netseg-1) segment intersections.
        Length must equal netseg-1. Default 0.5 for each.
    extension : str
        Filename extension. Default 'ets'.
    unitnumber : int or None
        File unit number. Default uses _defaultunit().
    filenames : str or list of str or None
        Filenames for package and CBC output files.

    parameters : flopy.modflow.ModflowParBc or None
        Parsed ETS array-parameter definitions (set by :meth:`load` when
        ``npets > 0`` and parameters are preserved). When present, the writer
        emits the parameter definition blocks and per-period active-parameter
        records instead of expanded ETSR arrays. Authoring this from scratch is
        not supported; see Notes.
    evtr_parm : dict or None
        Per-stress-period active-parameter records for ETSR, keyed by 0-based
        stress period: ``{kper: [(name, instance_or_None), ...]}``. Set by
        :meth:`load`. A period absent from the dict reuses the previous
        period's ETSR (``INETSR < 0``).

    Notes
    -----
    Only the ETSR (max ET rate) array can be parameterized in USG-T
    (``gwf2ets8u1.f``); ETSS/ETSX/IETS/PXDP/PETM are always plain arrays, so a
    parameterized ETS file naturally mixes a parameterized ETSR with
    non-parameterized arrays.

    USG-T reads ``NPETS`` from item 2a and calls ``UPARARRAL`` with ``IN=-1``
    (``parutl7.f``), so -- unlike MODFLOW-2005 ETS -- there is **no** separate
    ``PARAMETER N`` line. The writer puts ``NPETS`` in item 2a and follows it
    with the parameter definition blocks. The loader still tolerates a leading
    ``PARAMETER`` line for back-compatibility.

    Parameter *preservation* (load -> write -> reload with parameter syntax
    intact) is supported for ETSR. Parameter *authoring* from scratch
    (``npets > 0`` with no loaded definitions) raises ``NotImplementedError``.
    Loading with ``expand_parameters=True`` keeps the older "Expanded valid
    write" behavior (parameters expanded to arrays, ``NPETS=0`` on output).

    Examples
    --------
    >>> import flopy
    >>> m = flopy.mfusg.MfUsg(structured=False)
    >>> ets = flopy.mfusg.MfUsgEts(m, netsop=1, evtr=1.2e-4, netseg=2,
    ...                             pxdp=[0.5], petm=[0.5])
    """

    def __init__(
        self,
        model,
        netsop=1,
        ipakcb=None,
        surf=0.0,
        evtr=1e-3,
        exdp=1.0,
        ievt=1,
        netseg=1,
        npets=0,
        iesfactor=0,
        esfactor=None,
        pxdp=None,
        petm=None,
        parameters=None,
        evtr_parm=None,
        extension="ets",
        unitnumber=None,
        filenames=None,
    ):
        if unitnumber is None:
            unitnumber = MfUsgEts._defaultunit()

        filenames = self._prepare_filenames(filenames, 2)
        self.set_cbc_output_file(ipakcb, model, filenames[1])

        super().__init__(
            model,
            extension=extension,
            name=self._ftype(),
            unit_number=unitnumber,
            filenames=filenames[0],
        )

        self._generate_heading()
        self.url = "ets.html"

        self.netsop = netsop
        self.netseg = netseg
        self.npets = npets
        self.iesfactor = iesfactor
        self.esfactor = esfactor

        self.surf = Transient2d(
            model, get_pak_vals_shape(model, surf), np.float32, surf, name="surf"
        )
        self.evtr = Transient2d(
            model, get_pak_vals_shape(model, evtr), np.float32, evtr, name="evtr"
        )
        self.exdp = Transient2d(
            model, get_pak_vals_shape(model, exdp), np.float32, exdp, name="exdp"
        )
        self.ievt = Transient2d(
            model, get_pak_vals_shape(model, ievt), np.int32, ievt, name="ievt"
        )

        # Segment intersection arrays: one Transient2d per intersection
        nseg_int = max(0, netseg - 1)
        if pxdp is None:
            pxdp = [0.5] * nseg_int
        if petm is None:
            petm = [0.5] * nseg_int
        if not isinstance(pxdp, list):
            pxdp = [pxdp] * nseg_int
        if not isinstance(petm, list):
            petm = [petm] * nseg_int

        self.pxdp = [
            Transient2d(
                model,
                get_pak_vals_shape(model, pxdp[i]),
                np.float32,
                pxdp[i],
                name=f"pxdp{i}",
            )
            for i in range(nseg_int)
        ]
        self.petm = [
            Transient2d(
                model,
                get_pak_vals_shape(model, petm[i]),
                np.float32,
                petm[i],
                name=f"petm{i}",
            )
            for i in range(nseg_int)
        ]

        # Preserved MODFLOW array-parameter definitions for ETSR (set by load).
        self.parameters = parameters
        self.evtr_parm = evtr_parm if evtr_parm is not None else {}

        self.np = 0
        self.parent.add_package(self)

    # ------------------------------------------------------------------
    # internal helpers
    # ------------------------------------------------------------------

    def _ncells(self):
        nrow, ncol, nlay, nper = self.parent.nrow_ncol_nlay_nper
        return nrow * ncol

    # ------------------------------------------------------------------
    # write
    # ------------------------------------------------------------------

    def write_file(self, f=None):
        """Write the ETS package file."""
        preserve = self.npets > 0 and self.parameters is not None
        if self.npets > 0 and self.parameters is None:
            raise NotImplementedError(
                "MfUsgEts.write_file cannot author ETS parameter definitions "
                "from scratch (npets>0 without loaded parameter data). "
                "Parameter preservation is supported for files read by "
                "MfUsgEts.load; for from-scratch input use npets=0 (expanded "
                "arrays)."
            )
        nrow, ncol, nlay, nper = self.parent.nrow_ncol_nlay_nper
        close_on_exit = f is None
        if f is None:
            f = open(self.fn_path, "w")

        # Item 1 – heading
        f.write(f"{self.heading}\n")

        # Item 2a – NETSOP IETSCB NPETS NETSEG IESFACTOR
        f.write(
            f"{self.netsop:10d}{self.ipakcb:10d}"
            f"{self.npets:10d}{self.netseg:10d}{self.iesfactor:10d}\n"
        )

        # Item 2b – MXNDETS (unstructured + NETSOP==2 only)
        if not self.parent.structured and self.netsop == 2:
            mxndets = max(u2d.array.size for _, u2d in self.ievt.transient_2ds.items())
            f.write(f"{mxndets:10d}\n")

        # Item 2c – ESFACTOR(MCOMP) (transport active and IESFACTOR==1)
        if self.parent.itrnsp and self.iesfactor == 1:
            mcomp = self.parent.mcomp
            factors = (
                np.atleast_1d(self.esfactor)
                if self.esfactor is not None
                else np.ones(mcomp)
            )
            for icomp in range(mcomp):
                val = float(factors[icomp]) if icomp < len(factors) else 1.0
                f.write(f"{val:10.6f}")
            f.write("\n")

        # Items 3-4: ETS parameter definitions (preserved from load)
        if preserve:
            write_array_parameter_defs(f, self.parameters)

        nseg_int = max(0, self.netseg - 1)
        use_5a = (self.netsop == 2) or (self.netseg > 1)

        # Prepare 1-based layer indices for IEVT output
        if self.netsop == 2:
            ievt_out = {
                kper: u2d.array + 1 for kper, u2d in self.ievt.transient_2ds.items()
            }
            ievt_t2d = Transient2d(
                self.parent,
                self.ievt.shape,
                self.ievt.dtype,
                ievt_out,
                self.ievt.name,
            )

        for n in range(nper):
            insurf, surf_str = self.surf.get_kper_entry(n)
            if preserve:
                # ETSR is defined by active parameters: INETSR = count (>=1),
                # or -1 to reuse the previous period.
                recs = self.evtr_parm.get(n)
                inevtr = len(recs) if recs else -1
                evtr_str = None
            else:
                inevtr, evtr_str = self.evtr.get_kper_entry(n)
            inexdp, exdp_str = self.exdp.get_kper_entry(n)

            inievt = -1
            insgdf = -1
            ievt_str = None

            if self.netsop == 2:
                inievt, ievt_str = ievt_t2d.get_kper_entry(n)
                if inievt >= 0 and not self.parent.structured:
                    inievt = self.ievt[n].array.size

            if nseg_int > 0:
                # insgdf >= 0 if any segment data needs writing this period
                flags = [self.pxdp[i].get_kper_entry(n)[0] for i in range(nseg_int)]
                insgdf = 0 if any(v >= 0 for v in flags) else -1

            comment = f"ETS dataset 5 for stress period {n + 1}"
            if use_5a:
                # Item 5a: INETSS INETSR INETSX INIETS INSGDF
                f.write(
                    f"{insurf:10d}{inevtr:10d}{inexdp:10d}"
                    f"{inievt:10d}{insgdf:10d} #{comment}\n"
                )
            else:
                # Item 5b: INETSS INETSR INETSX
                f.write(f"{insurf:10d}{inevtr:10d}{inexdp:10d} #{comment}\n")

            if insurf >= 0:
                f.write(surf_str)
            if inevtr >= 0:
                if preserve:
                    # active-parameter records replace the ETSR array
                    write_active_array_parameters(f, self.evtr_parm[n])
                else:
                    f.write(evtr_str)
            if inexdp >= 0:
                f.write(exdp_str)
            if self.netsop == 2 and inievt >= 0:
                f.write(ievt_str)

            # (NETSEG-1) repetitions of PXDP + PETM
            if nseg_int > 0 and insgdf >= 0:
                for i in range(nseg_int):
                    _, pxdp_str = self.pxdp[i].get_kper_entry(n)
                    _, petm_str = self.petm[i].get_kper_entry(n)
                    f.write(pxdp_str)
                    f.write(petm_str)

        if close_on_exit:
            f.close()

    # ------------------------------------------------------------------
    # load
    # ------------------------------------------------------------------

    @classmethod
    def load(cls, f, model, nper=None, ext_unit_dict=None, expand_parameters=False):
        """
        Load an existing ETS package file.

        Parameters
        ----------
        f : str or file handle
            ETS file to load.
        model : MfUsg
            Model to which the package will be added.
        nper : int or None
            Number of stress periods (obtained from model if None).
        ext_unit_dict : dict or None
            External unit dictionary from parsenamefile.
        expand_parameters : bool
            When ``False`` (default) and the file is parameterized
            (``NPETS > 0``), the ETSR array parameters are *preserved*: the
            returned package keeps ``npets`` and writes the parameter syntax
            back. When ``True``, parameters are expanded to concrete ETSR
            arrays and the package writes ``NPETS=0`` (the older "Expanded
            valid write" behavior).

        Returns
        -------
        MfUsgEts
        """
        if model.verbose:
            print("loading ets package file...")

        openfile = not hasattr(f, "read")
        if openfile:
            filename = f
            f = open(filename, "r")

        # --- Skip comment lines (item 0) ---
        while True:
            line = f.readline()
            if line[0] != "#":
                break

        # --- Item 1: optional PARAMETER NPETS line ---
        npets = 0
        if "parameter" in line.lower():
            raw = line.strip().split()
            npets = int(raw[1])
            if npets > 0 and model.verbose:
                print(f"  Parameters detected. NPETS = {npets}")
            line = f.readline()

        # --- Item 2a: NETSOP IETSCB NPETS NETSEG IESFACTOR ---
        t = line.strip().split()
        netsop = int(t[0])
        ipakcb = int(t[1])
        if npets == 0 and len(t) >= 3:
            npets = int(t[2])
        netseg = int(t[3]) if len(t) >= 4 else 1
        iesfactor = int(t[4]) if len(t) >= 5 else 0

        # --- Item 2b: MXNDETS (unstructured + NETSOP==2 only) ---
        if not model.structured and netsop == 2:
            f.readline()  # read and discard; value not needed after load

        # --- Item 2c: ESFACTOR(MCOMP) ---
        esfactor = None
        if model.itrnsp and iesfactor == 1:
            line = f.readline()
            vals = line.strip().split()
            mcomp = model.mcomp if model.mcomp > 0 else 1
            esfactor = np.array([float(vals[i]) for i in range(min(mcomp, len(vals)))])

        # --- Items 3 & 4: parameter definitions ---
        pak_parms = None
        if npets > 0:
            pak_parms = mfparbc.loadarray(f, npets, model.verbose)

        # --- Stress period setup ---
        if nper is None:
            nrow, ncol, nlay, nper = model.get_nrow_ncol_nlay_nper()
        else:
            nrow, ncol, nlay, _ = model.get_nrow_ncol_nlay_nper()

        nseg_int = max(0, netseg - 1)
        use_5a = (netsop == 2) or (netseg > 1)

        # Base array shape (updated per stress period for unstructured grids)
        if not model.structured:
            u2d_shape = (1, ncol[0])
        else:
            u2d_shape = (nrow, ncol)

        surf_d, evtr_d, exdp_d, ievt_d = {}, {}, {}, {}
        evtr_parm_d = {}
        pxdp_d = [{} for _ in range(nseg_int)]
        petm_d = [{} for _ in range(nseg_int)]

        cur_surf = cur_evtr = cur_exdp = cur_ievt = []
        cur_pxdp = [[] for _ in range(nseg_int)]
        cur_petm = [[] for _ in range(nseg_int)]

        for iper in range(nper):
            line = f.readline()
            # strip inline comments before parsing
            t = line.split("#")[0].strip().split()

            insurf = int(t[0])
            inevtr = int(t[1])
            inexdp = int(t[2])
            inievt = -1
            insgdf = -1

            if use_5a:
                if len(t) >= 4:
                    inievt = int(t[3])
                if len(t) >= 5:
                    insgdf = int(t[4])

            # Adjust shape for unstructured grids
            if not model.structured:
                if netsop == 2 and inievt >= 0:
                    u2d_shape = (1, inievt)
                else:
                    u2d_shape = (1, ncol[0])

            # Item 6 / 13 – ETSS
            if insurf >= 0:
                if model.verbose:
                    print(f"   loading surf stress period {iper + 1:3d}...")
                cur_surf = Util2d.load(
                    f, model, u2d_shape, np.float32, "surf", ext_unit_dict
                )
            surf_d[iper] = cur_surf

            # Item 7 / 14 – ETSR
            if inevtr >= 0:
                if pak_parms is None or npets == 0:
                    if model.verbose:
                        print(f"   loading evtr stress period {iper + 1:3d}...")
                    cur_evtr = Util2d.load(
                        f, model, u2d_shape, np.float32, "evtr", ext_unit_dict
                    )
                else:
                    records = read_active_array_parameters(f, inevtr, pak_parms)
                    evtr_parm_d[iper] = records
                    parm_dict = {
                        name.lower(): (inst if inst else "static")
                        for name, inst in records
                    }
                    cur_evtr = mfparbc.parameter_bcfill(
                        model, u2d_shape, parm_dict, pak_parms
                    )
            evtr_d[iper] = cur_evtr

            # Item 9 / 16 – ETSX
            if inexdp >= 0:
                if model.verbose:
                    print(f"   loading exdp stress period {iper + 1:3d}...")
                cur_exdp = Util2d.load(
                    f, model, u2d_shape, np.float32, "exdp", ext_unit_dict
                )
            exdp_d[iper] = cur_exdp

            # Item 10 / 17 – IETS (NETSOP==2 only)
            if netsop == 2 and inievt >= 0:
                if model.verbose:
                    print(f"   loading ievt stress period {iper + 1:3d}...")
                t_ievt = Util2d.load(
                    f, model, u2d_shape, np.int32, "ievt", ext_unit_dict
                )
                cur_ievt = Util2d(model, u2d_shape, np.int32, t_ievt.array - 1, "ievt")
            if netsop == 2:
                ievt_d[iper] = cur_ievt

            # Items 11-12 / 18-19 – PXDP + PETM for each segment intersection
            if nseg_int > 0 and insgdf >= 0:
                for i in range(nseg_int):
                    if model.verbose:
                        print(f"   loading pxdp seg {i} stress period {iper + 1:3d}...")
                    cur_pxdp[i] = Util2d.load(
                        f, model, u2d_shape, np.float32, f"pxdp{i}", ext_unit_dict
                    )
                    if model.verbose:
                        print(f"   loading petm seg {i} stress period {iper + 1:3d}...")
                    cur_petm[i] = Util2d.load(
                        f, model, u2d_shape, np.float32, f"petm{i}", ext_unit_dict
                    )

            for i in range(nseg_int):
                pxdp_d[i][iper] = cur_pxdp[i]
                petm_d[i][iper] = cur_petm[i]

        if openfile:
            f.close()

        # Resolve unit number from ext_unit_dict
        unitnumber = None
        filenames = [None, None]
        if ext_unit_dict is not None:
            unitnumber, filenames[0] = model.get_ext_dict_attr(
                ext_unit_dict, filetype=cls._ftype()
            )
            _, filenames[1] = model.get_ext_dict_attr(ext_unit_dict, unit=ipakcb)

        # If the file is parameterized, either preserve the parameter syntax
        # (default) or expand it to concrete ETSR arrays (expand_parameters).
        if pak_parms is not None and not expand_parameters:
            npets_out = npets
            parameters = pak_parms
            evtr_parm = evtr_parm_d
        else:
            npets_out = 0 if pak_parms is not None else npets
            parameters = None
            evtr_parm = None

        return cls(
            model,
            netsop=netsop,
            ipakcb=ipakcb,
            surf=surf_d,
            evtr=evtr_d,
            exdp=exdp_d,
            ievt=ievt_d if netsop == 2 else 1,
            netseg=netseg,
            npets=npets_out,
            iesfactor=iesfactor,
            esfactor=esfactor,
            pxdp=pxdp_d if nseg_int > 0 else None,
            petm=petm_d if nseg_int > 0 else None,
            parameters=parameters,
            evtr_parm=evtr_parm,
            unitnumber=unitnumber,
            filenames=filenames,
        )

    # ------------------------------------------------------------------
    # class metadata
    # ------------------------------------------------------------------

    @staticmethod
    def _ftype():
        return "ETS"

    @staticmethod
    def _defaultunit():
        return 20
