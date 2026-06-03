"""
mfusgevt module.  Contains the MfUsgEvt class. Note that the user can access
the MfUsgEvt class as `flopy.mfusg.MfUsgEvt`.

"""

import numpy as np

from ..modflow.mfparbc import ModflowParBc as mfparbc
from ..pakbase import Package
from ..utils import Transient2d, Util2d
from ..utils.utils_def import (
    get_pak_vals_shape,
    type_from_iterable,
)
from ._usgt_parameters import (
    build_array_parameter_bc_parms,
    read_active_array_parameters,
    write_active_array_parameters,
    write_array_parameter_defs,
)


class MfUsgEvt(Package):
    """
    MODFLOW Evapotranspiration Package Class.

    Parameters
    ----------
    model : model object
        The model object (of type :class:`flopy.mfusg.MfUsgEvt`) to which
        this package will be added.
    ipakcb : int
        A flag that is used to determine if cell-by-cell budget data should be
        saved. If ipakcb is non-zero cell-by-cell budget data will be saved.
        (default is 0).
    nevtop : int
        is the recharge option code.
        1: ET is calculated only for cells in the top grid layer
        2: ET to layer defined in ievt
        3: ET to highest active cell (default is 3).
    surf : float or filename or ndarray or dict keyed on kper (zero-based)
        is the ET surface elevation. (default is 0.0, which is used for all
        stress periods).
    evtr: float or filename or ndarray or dict keyed on kper (zero-based)
        is the maximum ET flux (default is 1e-3, which is used for all
        stress periods).
    exdp : float or filename or ndarray or dict keyed on kper (zero-based)
        is the ET extinction depth (default is 1.0, which is used for all
        stress periods).
    ievt : int or filename or ndarray or dict keyed on kper (zero-based)
        is the layer indicator variable (default is 1, which is used for all
        stress periods).
    etfactor : float of array (mcomp) (default is 1.0)
        fraction of mass of the component that leaves with water
        0 = chemical component left behind in groundwater
        1 = chemical component leaves with water
        between 0 and 1 = fraction of mass of the component leaves
    iznevt : float of array (mcomp) (default is 1.0)
        array of zonal indices for applying a PET time series to zones.
        This PET input is independent of the stress period input, which
        is ignored when the zonal time series are provided.
    extension : string
        Filename extension (default is 'evt')
    unitnumber : int
        File unit number (default is None).
    filenames : str or list of str
        Filenames to use for the package and the output files. If
        filenames=None the package name will be created using the model name
        and package extension and the cbc output name will be created using
        the model name and .cbc extension (for example, mfusgtest.cbc),
        if ipakcbc is a number greater than zero. If a single string is passed
        the package will be set to the string and cbc output names will be
        created using the model name and .cbc extension, if ipakcbc is a
        number greater than zero. To define the names for all package files
        (input and output) the length of the list of strings should be 2.
        Default is None.

    Attributes
    ----------

    Methods
    -------

    See Also
    --------

    Notes
    -----
    Parameters are not supported in FloPy.

    Examples
    --------

    >>> import flopy
    >>> m = flopy.mfusg.MfUsg()
    >>> evt = flopy.mfusg.MfUsgEvt(m, nevtop=3, evtr=1.2e-4)

    """

    def __init__(
        self,
        model,
        nevtop=3,
        ipakcb=None,
        surf=0.0,
        evtr=1e-3,
        exdp=1.0,
        ievt=1,
        mxetzones=0,
        ietfactor=0,
        etfactor=0.0,
        inznevt=0,
        iznevt=0,
        npevt=0,
        parameters=None,
        evtr_parm=None,
        extension="evt",
        unitnumber=None,
        filenames=None,
        external=True,
    ):
        # set default unit number of one is not specified
        if unitnumber is None:
            unitnumber = MfUsgEvt._defaultunit()

        # set filenames
        filenames = self._prepare_filenames(filenames, 2)

        # cbc output file
        self.set_cbc_output_file(ipakcb, model, filenames[1])

        # call base package constructor
        super().__init__(
            model,
            extension=extension,
            name=self._ftype(),
            unit_number=unitnumber,
            filenames=filenames[0],
        )

        nrow, ncol, nlay, nper = self.parent.nrow_ncol_nlay_nper
        self._generate_heading()
        self.url = "evt.html"
        self.nevtop = nevtop
        if nevtop not in (1, 2, 3):
            raise ValueError(
                f"EVT NEVTOP must be 1, 2, or 3 (gwf2evt8u1.f); got {nevtop}."
            )

        self.mxetzones = int(mxetzones) if mxetzones else 0
        if self.mxetzones > 0:
            raise NotImplementedError(
                "EVT ETS zonal time-series ('ETS MXZNEVT') is not supported by "
                "MfUsgEvt. It is a dynamic, ATS-coupled execution mode, not "
                "static array I/O (gwf2evt8u1.f): (1) it requires adaptive "
                "time-stepping -- the Fortran STOPs when IATS==0; (2) the ET "
                "rates come from a separate external ETS time-series file "
                "(unit IUETS) read progressively during the run, one "
                "'Tstart Tend Factor Ets(1..MXZNEVT)' record at a time as the "
                "simulation time reaches Tend; (3) IETSOPT=1 means the "
                "time-series SUPERSEDES the array EVTR -- EVTR is recomputed each "
                "step as etsevt(IZNEVT(n))*AREA*Factor, so it is not a "
                "preservable array; and (4) per stress period an 'INEVTZONES' "
                "flag (re)reads the IZNEVT zone-index array. FloPy models static "
                "EVT I/O only; remove the ETS option (deferred -- Stage 4.6F-B)."
            )

        self.ietfactor = int(ietfactor or 0)
        # Normalize ETFACTOR to an indexable 1-D array so a scalar (e.g.
        # etfactor=2.5 with MCOMP=1) does not crash write_file's per-component
        # indexing. ETFACTOR is only written when ietfactor>0.
        self.etfactor = np.atleast_1d(etfactor) if etfactor is not None else None
        # Transport ET factor: USG-T reads MCOMP ETFACTOR values when
        # ietfactor>0 (ietfactor==0 -> implicit 0.0, ietfactor<0 -> implicit 1.0).
        if getattr(model, "itrnsp", 0) and self.ietfactor > 0:
            netf = 0 if self.etfactor is None else self.etfactor.size
            if netf != model.mcomp:
                raise ValueError(
                    f"EVT etfactor must have MCOMP={model.mcomp} value(s) when "
                    f"ietfactor>0; got {netf}."
                )

        self.external = external
        if self.external is False:
            load = True
        else:
            load = model.load

        surf_u2d_shape = get_pak_vals_shape(model, surf)
        evtr_u2d_shape = get_pak_vals_shape(model, evtr)
        exdp_u2d_shape = get_pak_vals_shape(model, exdp)
        ievt_u2d_shape = get_pak_vals_shape(model, ievt)

        self.surf = Transient2d(model, surf_u2d_shape, np.float32, surf, name="surf")
        self.evtr = Transient2d(model, evtr_u2d_shape, np.float32, evtr, name="evtr")
        self.exdp = Transient2d(model, exdp_u2d_shape, np.float32, exdp, name="exdp")
        self.ievt = Transient2d(model, ievt_u2d_shape, np.int32, ievt, name="ievt")

        # self.inznevt = [inznevt]*nper
        # self.iznevt = iznevt
        # iznevt_u2d_shape = get_pak_vals_shape(model, iznevt)
        # self.iznevt = Transient2d(model, iznevt_u2d_shape,
        # np.int32, iznevt, name="iznevt")

        # NPEVT array parameters (preserved from load, or authored from scratch).
        # Like ETS, EVT parameterizes only the EVTR (max ET-rate) array via the
        # MODFLOW array-parameter machinery (UPARARRAL/UPARARRRP/UPARARRSUB2,
        # PARTYP='EVT'); evtr_parm holds the per-SP active-parameter records.
        self.npevt = npevt
        self.parameters = parameters
        self.evtr_parm = evtr_parm if evtr_parm is not None else {}
        self.np = 0
        self.parent.add_package(self)

    def _ncells(self):
        """Maximum number of cells that have evapotranspiration (developed for
        MT3DMS SSM package).

        Returns
        -------
        ncells: int
            maximum number of evt cells

        """
        nrow, ncol, nlay, nper = self.parent.nrow_ncol_nlay_nper
        return nrow * ncol

    def _check_evtr_parm(self):
        """Validate the ``evtr_parm`` container/record structure before any use,
        so a malformed input fails with an actionable ``ValueError`` instead of an
        ``AttributeError`` (a list/string where a dict is expected) or a
        tuple-unpack error (a record that is not a ``(name, instance)`` pair).
        Semantic checks (range, duplicates, defined names, instances) are done in
        :meth:`_validate_active_params` once the structure is known to be sound.
        """
        if not isinstance(self.evtr_parm, dict):
            raise ValueError(
                "MfUsgEvt.write_file: evtr_parm must be a dict keyed by 0-based "
                f"stress period; got {type(self.evtr_parm).__name__}."
            )
        for kper, recs in self.evtr_parm.items():
            if isinstance(kper, bool) or not isinstance(kper, (int, np.integer)):
                raise ValueError(
                    f"MfUsgEvt.write_file: evtr_parm key {kper!r} must be an "
                    "integer stress period (0-based)."
                )
            if isinstance(recs, str) or not isinstance(recs, (list, tuple)):
                raise ValueError(
                    f"MfUsgEvt.write_file: evtr_parm[{kper}] must be a list of "
                    f"(name, instance) pairs; got {type(recs).__name__}."
                )
            for rec in recs:
                if (
                    isinstance(rec, str)
                    or not isinstance(rec, (list, tuple))
                    or len(rec) != 2
                ):
                    raise ValueError(
                        f"MfUsgEvt.write_file: evtr_parm[{kper}] record {rec!r} "
                        "must be a (name, instance_or_None) pair."
                    )
                name, instance = rec
                if not isinstance(name, str) or not name.strip():
                    raise ValueError(
                        f"MfUsgEvt.write_file: evtr_parm[{kper}] parameter name "
                        f"{name!r} must be a non-empty string."
                    )
                if instance is not None and not isinstance(instance, str):
                    raise ValueError(
                        f"MfUsgEvt.write_file: evtr_parm[{kper}] instance "
                        f"{instance!r} must be a string or None."
                    )

    def _resolve_parameters(self):
        """Resolve + validate EVT array parameters (EVTR) for a write.

        Accepts ``self.parameters`` as an ergonomic dict (built into a
        :class:`ModflowParBc` via the shared array-parameter machinery) or as an
        already-canonical :class:`ModflowParBc` (preserved from :meth:`load`).
        Returns ``(ModflowParBc or None, npevt)`` with ``npevt`` auto-computed
        from the definitions when omitted/0. Mirrors ``MfUsgEts`` (EVT uses the
        same UPARARRAL/UPARARRRP/UPARARRSUB2 machinery, PARTYP='EVT').
        """
        self._check_evtr_parm()
        if self.parameters is None:
            if any(self.evtr_parm.values()):
                raise ValueError(
                    "MfUsgEvt.write_file: evtr_parm activates parameters but none "
                    "are defined. Pass parameters={name: {...}} to author EVT "
                    "(EVTR) array parameters from scratch."
                )
            if self.npevt:
                raise ValueError(
                    f"MfUsgEvt.write_file: npevt={self.npevt} but no parameters "
                    "are defined; pass parameters={name: {...}} or npevt=0."
                )
            return None, 0

        # An empty parameters dict carries no definitions: fail directly rather
        # than building an empty ModflowParBc and tripping the first-period check.
        if isinstance(self.parameters, dict) and not self.parameters:
            raise ValueError(
                "MfUsgEvt.write_file: parameters is empty (no parameter "
                "definitions). Pass parameters={name: {...}} or parameters=None."
            )

        if isinstance(self.parameters, dict):
            bc_parms = build_array_parameter_bc_parms(
                self.parameters, "evt", prefix="MfUsgEvt.write_file"
            )
            params = mfparbc(bc_parms)
        else:
            params = self.parameters  # preserved ModflowParBc

        npevt = self.npevt if self.npevt else len(params.bc_parms)
        if npevt != len(params.bc_parms):
            raise ValueError(
                f"MfUsgEvt.write_file: NPEVT ({npevt}) must equal the number of "
                f"parameter definitions ({len(params.bc_parms)})."
            )
        self._validate_active_params(params)
        return params, npevt

    def _validate_active_params(self, params):
        """Validate ``evtr_parm`` activations against the (canonical) definitions.

        The first stress period must activate at least one parameter (USG-T
        reuses the previous period's EVTR when ``INEVTR<0``, so a first parametric
        period with ``INEVTR<0`` would reuse an uninitialized EVTR). Per period:
        no duplicate activation, the name must be defined, a time-varying
        parameter must name a known instance, a static one must not.
        """
        nper = self.parent.nrow_ncol_nlay_nper[3]
        if not self.evtr_parm.get(0):
            raise ValueError(
                "MfUsgEvt.write_file: the first stress period must activate EVTR "
                "parameters (evtr_parm[0]); USG-T cannot reuse a previous EVTR on "
                "the first parametric period."
            )
        for kper, recs in self.evtr_parm.items():
            if not 0 <= kper < nper:
                raise ValueError(
                    f"MfUsgEvt.write_file: evtr_parm stress period {kper} is out "
                    f"of range 0..{nper - 1}."
                )
            lowered = [name.lower() for name, _ in recs]
            if len(set(lowered)) != len(lowered):
                raise ValueError(
                    f"MfUsgEvt.write_file: a parameter is activated more than once "
                    f"in stress period {kper}: {recs}."
                )
            for name, instance in recs:
                pdef = params.bc_parms.get(name.lower())
                if pdef is None:
                    raise ValueError(
                        f"MfUsgEvt.write_file: active parameter '{name}' (stress "
                        f"period {kper}) is not defined in parameters."
                    )
                timevarying = pdef[0]["timevarying"]
                pinst = pdef[1]
                if timevarying:
                    if instance is None:
                        raise ValueError(
                            f"MfUsgEvt.write_file: parameter '{name}' is "
                            "time-varying (INSTANCES); activation must name an "
                            f"instance (stress period {kper})."
                        )
                    if str(instance).lower() not in pinst:
                        raise ValueError(
                            f"MfUsgEvt.write_file: parameter '{name}' has no "
                            f"instance '{instance}' (stress period {kper}); "
                            f"defined instances: {sorted(pinst)}."
                        )
                elif instance is not None and str(instance).lower() != "static":
                    raise ValueError(
                        f"MfUsgEvt.write_file: parameter '{name}' is static but is "
                        f"activated with instance '{instance}' (stress period "
                        f"{kper}); use None."
                    )

    def write_file(self, f=None):
        """
        Write the package file.

        Returns
        -------
        None

        """
        # Resolve + validate EVT array parameters before opening the file, so a
        # PARAMETER NPEVT header is never written without a complete, consistent
        # body (no partial file). ``params`` is a ModflowParBc (authored or
        # preserved); ``npevt`` is auto-computed when omitted.
        params, npevt = self._resolve_parameters()
        preserve = params is not None
        nrow, ncol, nlay, nper = self.parent.nrow_ncol_nlay_nper
        if f is not None:
            f_evt = f
        else:
            f_evt = open(self.fn_path, "w")
        f_evt.write(f"{self.heading}\n")
        # Item 1: optional PARAMETER NPEVT line (UPARARRAL) before item 2.
        if preserve:
            f_evt.write(f"PARAMETER {npevt:10d}\n")
        f_evt.write(f"{self.nevtop:10d}{self.ipakcb:10d}")

        # USG-T's reader takes 3 integers (NEVTOP IEVTCB IETFACTOR) whenever
        # transport (BCT) is active, so IETFACTOR must always be present then.
        if self.parent.itrnsp:
            f_evt.write(f"{self.ietfactor:10d}")
        f_evt.write("\n")

        if self.nevtop == 2:
            # USG-T validates the IEVT index: a structured layer in 1..NLAY,
            # or an unstructured node in 1..NODES (1-based in the file). Internal
            # values are 0-based, so check [0, NLAY-1] / [0, NODES-1].
            if self.parent.structured:
                hi, kind = nlay, "layer"
            else:
                from ._tabrich import node_count

                hi, kind = node_count(self.parent), "node"
            ievt = {}
            for kper, u2d in self.ievt.transient_2ds.items():
                arr = u2d.array  # 0-based internal
                if arr.size and (arr.min() < 0 or arr.max() > hi - 1):
                    raise ValueError(
                        f"EVT ievt (0-based {kind}) must be in [0, {hi - 1}] for "
                        f"NEVTOP=2; stress period {kper} is out of range."
                    )
                ievt[kper] = arr + 1  # 1-based file
            ievt = Transient2d(
                self.parent,
                self.ievt.shape,
                self.ievt.dtype,
                ievt,
                self.ievt.name,
            )
            if not self.parent.structured:
                mxndevt = np.max(
                    [u2d.array.size for kper, u2d in self.ievt.transient_2ds.items()]
                )
                f_evt.write(f"{mxndevt:10d}\n")

        # USG-T reads the MCOMP ETFACTOR array only when ietfactor>0.
        if self.parent.itrnsp and self.ietfactor > 0:
            mcomp = self.parent.mcomp
            for icomp in range(mcomp):
                f_evt.write(f"{self.etfactor[icomp]:10.2e}")
            f_evt.write("\n")

        # EVT parameter definitions (preserved from load or authored from scratch)
        if preserve:
            write_array_parameter_defs(f_evt, params)

        for n in range(nper):
            insurf, surf = self.surf.get_kper_entry(n)
            if preserve:
                # EVTR is defined by active parameters: INEVTR = count (>=1),
                # or -1 to reuse the previous period's EVTR.
                recs = self.evtr_parm.get(n)
                inevtr = len(recs) if recs else -1
                evtr = None
            else:
                inevtr, evtr = self.evtr.get_kper_entry(n)
            inexdp, exdp = self.exdp.get_kper_entry(n)
            inievt = -1
            if self.nevtop == 2:
                inievt, file_entry_ievt = ievt.get_kper_entry(n)
                if inievt >= 0 and not self.parent.structured:
                    inievt = self.ievt[n].array.size
            comment = f"Evapotranspiration dataset 5 for stress period {n + 1}"
            f_evt.write(f"{insurf:10d}{inevtr:10d}{inexdp:10d}{inievt:10d} ")
            # if self.inznevt[n] > 0:
            #     f_evt.write(f"INEVTZONES {self.inznevt[n]:10d}\n")
            f_evt.write(f"#{comment}\n")

            if insurf >= 0:
                f_evt.write(surf)
            if inevtr >= 0:
                if preserve:
                    # active-parameter records replace the EVTR array
                    write_active_array_parameters(f_evt, self.evtr_parm[n])
                else:
                    f_evt.write(evtr)
            if inexdp >= 0:
                f_evt.write(exdp)
            if self.nevtop == 2 and inievt >= 0:
                f_evt.write(file_entry_ievt)

            # if self.inznevt[n] > 0:
            #     iznevt = self.iznevt.get_kper_entry(n)
            #     f_evt.write(iznevt)

        f_evt.close()

    @classmethod
    def load(cls, f, model, nper=None, ext_unit_dict=None, expand_parameters=False):
        """
        Load an existing package.

        Parameters
        ----------
        f : filename or file handle
            File to load.
        model : model object
            The model object (of type :class:`flopy.mfusg.mf.MfUsg`) to
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

        Returns
        -------
        evt : MfUsgEvt object
            MfUsgEvt object.

        Examples
        --------

        >>> import flopy
        >>> m = flopy.mfusg.MfUsg()
        >>> evt = flopy.mfusg.mfevt.load('test.evt', m)

        """
        if model.verbose:
            print("loading evt package file...")

        openfile = not hasattr(f, "read")
        if openfile:
            filename = f
            f = open(filename, "r")

        # Dataset 0 -- header
        while True:
            line = f.readline()
            if line[0] != "#":
                break
        npar = 0
        if line.strip().lower().split()[:1] == ["parameter"]:
            raw = line.strip().split()
            if len(raw) < 2:
                raise ValueError(
                    "MfUsgEvt.load: PARAMETER line must be 'PARAMETER <NPEVT>'; "
                    f"got {line.strip()!r}."
                )
            try:
                npar = int(raw[1])
            except ValueError:
                raise ValueError(
                    "MfUsgEvt.load: PARAMETER line must be 'PARAMETER <NPEVT>' "
                    f"with an integer count; got {raw[1]!r}."
                )
            if npar > 0:
                if model.verbose:
                    print("  Parameters detected. Number of parameters = ", npar)
            line = f.readline()
        # Dataset 2
        t = line.strip().split()
        nevtop = int(t[0])
        ipakcb = int(t[1])
        ietfactor = type_from_iterable(t, 2)

        # Options. ETS zonal time-series is a dynamic ATS-coupled execution mode
        # (external IUETS time-series superseding the EVTR array), not static I/O;
        # fail here, before reading any stress-period data, so nothing is parsed
        # partially (see __init__ and Stage 4.6F-B for the full spec).
        mxetzones = 0
        if any(tok.upper() == "ETS" for tok in t):
            if openfile:
                f.close()
            raise NotImplementedError(
                "MfUsgEvt.load: EVT ETS zonal time-series ('ETS MXZNEVT') is not "
                "supported -- it is an ATS-coupled dynamic mode reading an "
                "external time-series file (IUETS) that supersedes the EVTR "
                "array; FloPy models static EVT I/O only (deferred, Stage "
                "4.6F-B)."
            )

        # dataset 2b for mfusg
        if not model.structured and nevtop == 2:
            line = f.readline()
            t = line.strip().split()
            mxndevt = int(t[0])

        # dataset 2c for mfusg
        etfactor = None
        mcomp = model.mcomp
        if mcomp > 0:
            etfactor = np.zeros(model.mcomp)
            if ietfactor < 0:
                etfactor = np.ones(mcomp)
            if ietfactor > 0:
                if mcomp > 0:
                    line = f.readline()
                    t = line.strip().split()
                    for icomp in range(mcomp):
                        etfactor[icomp] = float(t[icomp])

        # Dataset 3 and 4 - parameters data
        pak_parms = None
        if npar > 0:
            pak_parms = mfparbc.loadarray(f, npar, model.verbose)

        if nper is None:
            nrow, ncol, nlay, nper = model.get_nrow_ncol_nlay_nper()
        else:
            nrow, ncol, nlay, _ = model.get_nrow_ncol_nlay_nper()

        u2d_shape = (nrow, ncol)

        # Read data for every stress period
        surf = {}
        evtr = {}
        exdp = {}
        ievt = {}
        evtr_parm_d = {}
        # iznevt = {}
        current_surf = []
        current_evtr = []
        current_exdp = []
        current_ievt = []
        current_iznevt = []
        # inznevt = [0] * nper
        for iper in range(nper):
            line = f.readline()
            t = line.strip().split()
            # The per-SP INEVTZONES flag belongs to the ETS zonal time-series
            # mode (it (re)reads the IZNEVT zone-index array and drives the
            # external ETS file). That mode is not modeled; fail explicitly here
            # rather than mis-parsing the header (deferred, Stage 4.6F-B).
            if any(tok.upper() == "INEVTZONES" for tok in t):
                if openfile:
                    f.close()
                raise NotImplementedError(
                    "MfUsgEvt.load: the per-stress-period 'INEVTZONES' flag is "
                    "part of the EVT ETS zonal time-series mode (IZNEVT zone "
                    "array + external ETS file), which is not supported -- FloPy "
                    "models static EVT I/O only (deferred, Stage 4.6F-B)."
                )
            insurf = int(t[0])
            inevtr = int(t[1])
            inexdp = int(t[2])

            if nevtop == 2:
                inievt = int(t[3])
                if (not model.structured) and (inievt >= 0):
                    u2d_shape = (1, inievt)
            elif not model.structured:
                u2d_shape = (1, ncol[0])

            if insurf >= 0:
                if model.verbose:
                    print(f"   loading surf stress period {iper + 1:3d}...")
                t = Util2d.load(f, model, u2d_shape, np.float32, "surf", ext_unit_dict)
                current_surf = t
            surf[iper] = current_surf

            if inevtr >= 0:
                if npar == 0:
                    if model.verbose:
                        print(f"   loading evtr stress period {iper + 1:3d}...")
                    t = Util2d.load(
                        f,
                        model,
                        u2d_shape,
                        np.float32,
                        "evtr",
                        ext_unit_dict,
                    )
                else:
                    records = read_active_array_parameters(f, inevtr, pak_parms)
                    evtr_parm_d[iper] = records
                    parm_dict = {
                        name.lower(): (inst if inst else "static")
                        for name, inst in records
                    }
                    t = mfparbc.parameter_bcfill(model, u2d_shape, parm_dict, pak_parms)

                current_evtr = t
            evtr[iper] = current_evtr
            if inexdp >= 0:
                if model.verbose:
                    print(f"   loading exdp stress period {iper + 1:3d}...")
                t = Util2d.load(f, model, u2d_shape, np.float32, "exdp", ext_unit_dict)
                current_exdp = t
            exdp[iper] = current_exdp
            if nevtop == 2:
                if inievt >= 0:
                    if model.verbose:
                        print(f"   loading ievt stress period {iper + 1:3d}...")
                    t = Util2d.load(
                        f, model, u2d_shape, np.int32, "ievt", ext_unit_dict
                    )
                    current_ievt = Util2d(
                        model, u2d_shape, np.int32, t.array - 1, "ievt"
                    )
                ievt[iper] = current_ievt
            # if inznevt[iper] > 0:
            #     if model.verbose:
            #         print(f"   loading iznevt stress period {iper + 1:3d}...")
            #     current_iznevt = Util2d.load(f, model,
            #     (inznevt[iper],), np.int32, "iznevt", ext_unit_dict)
            # iznevt[iper] = current_iznevt

        if openfile:
            f.close()

        # create evt object
        args = {}
        args["ievt"] = ievt
        args["nevtop"] = nevtop
        args["evtr"] = evtr
        args["surf"] = surf
        args["exdp"] = exdp
        args["ipakcb"] = ipakcb

        args["mxetzones"] = mxetzones
        args["ietfactor"] = ietfactor
        args["etfactor"] = etfactor
        # args["inznevt"] = inznevt
        # args["iznevt"] = iznevt

        # If the file is parameterized, either preserve the parameter syntax
        # (default) or expand it to concrete EVTR arrays (expand_parameters=True
        # -> legacy "Expanded valid write", NPEVT=0 on output).
        if pak_parms is not None and not expand_parameters:
            args["npevt"] = npar
            args["parameters"] = pak_parms
            args["evtr_parm"] = evtr_parm_d

        # determine specified unit number
        unitnumber = None
        filenames = [None, None]
        if ext_unit_dict is not None:
            unitnumber, filenames[0] = model.get_ext_dict_attr(
                ext_unit_dict, filetype=cls._ftype()
            )
            _, filenames[1] = model.get_ext_dict_attr(ext_unit_dict, unit=ipakcb)

        # return evt object
        return cls(model, unitnumber=unitnumber, filenames=filenames, **args)

    @staticmethod
    def _ftype():
        return "EVT"

    @staticmethod
    def _defaultunit():
        return 22
