"""
mfusglak module.  Contains the MfUsgLak class. Note that the user can access
the MfUsgLak class as `flopy.mfusg.MfUsgLak`.

"""

import numpy as np

from ..pakbase import Package
from ..utils import Util3d, read_fixed_var, utils_def, write_fixed_var
from ..utils.util_array import Transient3d


class MfUsgLak(Package):
    """
    MODFLOW USG Transport Lake Package Class.

    Parameters
    ----------
    model : model object
        The model object (of type :class:`flopy.mfusg.MfUsg`) to which
        this package will be added.
    nlakes : int
        NLAKES Number of separate lakes.
        Sublakes of multiple-lake systems are considered separate lakes for
        input purposes. The variable NLAKES is used, with certain internal
        assumptions and approximations, to dimension arrays for the simulation.
    ipakcb : int, optional
        Toggles whether cell-by-cell budget data should be saved. If None or zero,
        budget data will not be saved (default is None).
    lwrt : int or list of ints (one per SP)
        lwrt > 0, suppresses printout from the lake package. Default is 0 (to
        print budget information)
    theta : float
        Explicit (THETA = 0.0), semi-implicit (0.0 < THETA < 1.0), or implicit
        (THETA = 1.0) solution for lake stages. SURFDEPTH is read only if
        THETA is assigned a negative value (the negative value of THETA is
        then changed to a positive value internally by the code).
        *   A new method of solving for lake stage uses only the time-weighting
            factor THETA (Merritt and Konikow, 2000, p. 52) for transient
            simulations. THETA is automatically set to a value of 1.0 for all
            steady-state stress periods. For transient stress periods, Explicit
            (THETA = 0.0), semi-implicit (0.0 < THETA < 1.0), or implicit
            (THETA = 1.0) solutions can be used to calculate lake stages. The
            option to specify negative values for THETA is supported to allow
            specification of additional variables (NSSITER, SSCNCR, SURFDEP)
            for simulations that only include transient stress periods. If
            THETA is specified as a negative value, then it is converted to a
            positive value for calculations of lake stage.
        *   In MODFLOW-2000 and later, ISS is not part of the input. Instead
            NSSITR or SSCNCR should be included if one or more stress periods
            is a steady state stress period as defined in Ss/tr in the
            Discretization file.
        *   SSCNCR and NSSITR can be read for a transient only simulation by
            placing a negative sign immediately in front of THETA. A negative
            THETA sets a flag which assumes input values for NSSITR and SSCNCR
            will follow THETA in the format as described by Merritt and Konikow
            (p. 52). A negative THETA is automatically reset to a positive
            value after values of NSSITR and SSCNCR are read.
    nssitr : int
        Maximum number of iterations for Newton's method of solution for
        equilibrium lake stages in each MODFLOW iteration for steady-state
        aquifer head solution. Only read if ISS (option flag input to DIS
        Package of MODFLOW indicating steady-state solution) is not zero or
        if THETA is specified as a negative value.
        *   NSSITR and SSCNCR may be omitted for transient solutions (ISS = 0).
        *   In MODFLOW-2000 and later, ISS is not part of the input.
            Instead NSSITR or SSCNCR should be included if one or more stress
            periods is a steady state stress period as defined in Ss/tr in the
            Discretization file.
        *   SSCNCR and NSSITR can be read for a transient only simulation by
            placing a negative sign immediately in front of THETA. A negative
            THETA sets a flag which assumes input values for NSSITR and SSCNCR
            will follow THETA in the format as described by Merritt and Konikow
            (p. 52). A negative THETA is automatically reset to a positive
            value after values of NSSITR and SSCNCR are read.
        *   If NSSITR = 0, a value of 100 will be used instead.
    sscncr : float
        Convergence criterion for equilibrium lake stage solution by Newton's
        method. Only read if ISS is not zero or if THETA is specified as a
        negative value. See notes above for nssitr.
    surfdepth : float
        The height of small topological variations (undulations) in lake-bottom
        elevations that can affect groundwater discharge to lakes. SURFDEPTH
        decreases the lakebed conductance for vertical flow across a horizontal
        lakebed caused both by a groundwater head that is between the lakebed
        and the lakebed plus SURFDEPTH and a lake stage that is also between
        the lakebed and the lakebed plus SURFDEPTH. This method provides a
        smooth transition from a condition of no groundwater discharge to a
        lake, when groundwater head is below the lakebed, to a condition of
        increasing groundwater discharge to a lake as groundwater head becomes
        greater than the elevation of the dry lakebed. The method also allows
        for the transition of seepage from a lake to groundwater when the lake
        stage decreases to the lakebed elevation. Values of SURFDEPTH ranging
        from 0.01 to 0.5 have been used successfully in test simulations.
        SURFDEP is read only if THETA is specified as a negative value.
    stages : float or list of floats
        The initial stage of each lake at the beginning of the run.
    stage_range : list of tuples (ssmn, ssmx) of length nlakes
        Where ssmn and ssmx are the minimum and maximum stages allowed for each
        lake in steady-state solution.
        *   SSMN and SSMX are not needed for a transient run and must be
            omitted when the solution is transient.
        *   When the first stress period is a steady-state stress period,
            SSMN is defined in record 3.

        For subsequent steady-state stress periods, SSMN is defined in
        record 9a.
    lakarr : array of integers (nlay, nrow, ncol)
        LKARR A value is read in for every grid cell.
        If LKARR(I,J,K) = 0, the grid cell is not a lake volume cell.
        If LKARR(I,J,K) > 0, its value is the identification number of the lake
        occupying the grid cell. LKARR(I,J,K) must not exceed the value NLAKES.
        If it does, or if LKARR(I,J,K) < 0, LKARR(I,J,K) is set to zero.
        Lake cells cannot be overlain by non-lake cells in a higher layer.
        Lake cells must be inactive cells (IBOUND = 0) and should not be
        convertible to active cells (WETDRY = 0).

        The Lake package can be used when all or some of the model layers
        containing the lake are confined.  The authors recommend using the
        Layer-Property Flow Package (LPF) for this case, although the
        BCF and HUF Packages will work too.  However, when using the BCF6
        package to define aquifer properties, lake/aquifer conductances in the
        lateral direction are based solely on the lakebed leakance (and not on
        the lateral transmissivity of the aquifer layer).  As before, when the
        BCF6 package is used, vertical lake/aquifer conductances are based on
        lakebed conductance and on the vertical hydraulic conductivity of the
        aquifer layer underlying the lake when the wet/dry option is
        implemented, and only on the lakebed leakance when the wet/dry option
        is not implemented.
    bdlknc : array of floats (nlay, nrow, ncol)
        BDLKNC A value is read in for every grid cell. The value is the lakebed
        leakance that will be assigned to lake/aquifer interfaces that occur
        in the corresponding grid cell. If the wet-dry option flag (IWDFLG) is
        not active (cells cannot rewet if they become dry), then the BDLKNC
        values are assumed to represent the combined leakances of the lakebed
        material and the aquifer material between the lake and the centers of
        the underlying grid cells, i. e., the vertical conductance values (CV)
        will not be used in the computation of conductances across lake/aquifer
        boundary faces in the vertical direction.

        IBOUND and WETDRY should be set to zero for every cell for which LKARR
        is not equal to zero. IBOUND is defined in the input to the Basic
        Package of MODFLOW. WETDRY is defined in the input to the BCF or other
        flow package of MODFLOW if the IWDFLG option is active. When used with
        the HUF package, the Lake Package has been modified to compute
        effective lake-aquifer conductance solely on the basis of the
        user-specified value of lakebed leakance; aquifer hydraulic
        conductivities are not used in this calculation. An appropriate
        informational message is now printed after the lakebed conductances
        are written to the main output file.
    sill_data : dict
        (datasets 7/8 in the documentation) Connected-lake (sublake) systems.
        A dict keyed by 0-based stress period; each value is a list of
        ``(ds8a, sillvt)`` systems (the count is dataset 7, ``NSLMS``):

        * ``ds8a`` (dataset 8a) is ``[IC, lake1, lake2, ... lakeIC]`` -- the
          number of lakes ``IC`` in the system followed by that many **1-based**
          lake numbers; the **center lake is listed first** and the rest are its
          sublakes. ``IC >= 2``.
        * ``sillvt`` (dataset 8b) is a sequence of ``IC - 1`` sill elevations,
          one per sublake in the order the sublakes appear in ``ds8a`` (the
          sill controls when the center lake is connected to that sublake).

        Example -- one system of two lakes (center 1, sublake 2)::

            sill_data = {0: [([2, 1, 2], [95.0])]}

        Datasets 7/8 are only read/written for stress periods where ``ITMP>0``
        (i.e. where ``lakarr``/``bdlknc`` are (re)specified); USG-T reuses the
        previous period and skips them otherwise. Validated in ``__init__``.
    flux_data : dict
        (dataset 9a in documentation)
        Dict of lists keyed by stress period. The list for each stress period
        is a list of lists, with each list containing the variables
        PRCPLK EVAPLK RNF WTHDRW [SSMN] [SSMX] from the documentation.
            PRCPLK : float
                The rate of precipitation per unit area at the surface of a
                lake (L/T).
            EVAPLK : float
                The rate of evaporation per unit area from the surface of a
                lake (L/T).
            RNF : float
                Overland runoff from an adjacent watershed entering the lake.
                If RNF > 0, it is specified directly as a volumetric rate, or
                flux (L3 /T). If RNF < 0, its absolute value is used as a
                dimensionless multiplier applied to the product of the lake
                precipitation rate per unit area (PRCPLK) and the surface area
                of the lake at its full stage (occupying all layer 1 lake
                cells). When RNF is entered as a dimensionless multiplier
                (RNF < 0), it is considered to be the product of two
                proportionality factors. The first is the ratio of the area of
                the basin contributing runoff to the surface area of the lake
                when it is at full stage. The second is the fraction of the
                current rainfall rate that becomes runoff to the lake. This
                procedure provides a means for the automated computation of
                runoff rate from a watershed to a lake as a function of
                varying rainfall rate. For example, if the basin area is 10
                times greater than the surface area of the lake, and 20 percent
                of the precipitation on the basin becomes overland runoff
                directly into the lake, then set RNF = -2.0.
            WTHDRW : float
                The volumetric rate, or flux (L3 /T), of water removal from a
                lake by means other than rainfall, evaporation, surface
                outflow, or groundwater seepage. A negative value indicates
                augmentation. Normally, this would be used to specify the
                rate of artificial withdrawal from a lake for human water use,
                or if negative, artificial augmentation of a lake volume for
                aesthetic or recreational purposes.
            SSMN : float
                Minimum stage allowed for each lake in steady-state solution.
                See notes on ssmn and ssmx above.
            SSMX : float
                SSMX Maximum stage allowed for each lake in steady-state
                solution.
    conc_data : dict
        (dataset 9b in documentation)
        Dict of lists keyed by stress period. The list for each stress period
        is a list of lists, with each list containing the variables
        CPPT CRNF [CAUG] from the documentation.
            CPPT : float
                concentration of solute in precipitation onto the lake surface
            CRNF : float
                concentration of solute in overland runoff directly into the lake
            CAUG : float
                concentration of solute in water used to augment the lake volume
        if
    options : list of strings
        Package options. (default is None).
    extension : string
        Filename extension (default is 'lak')
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
    >>> lak = {}
    >>> lak[0] = [[2, 3, 4, 15.6, 1050., -4]]  #this lake boundary will be
    >>>                                        #applied to all stress periods
    >>> lak = flopy.mfusg.MfUsgLak(m, nstress_period_data=strd)

    """

    def __init__(
        self,
        model,
        nlakes=1,
        ipakcb=None,
        theta=-1.0,
        nssitr=0,
        sscncr=0.0,
        surfdep=0.0,
        stages=1.0,
        stage_range=None,
        clake=None,
        tab_files=None,
        tab_units=None,
        lakarr=None,
        bdlknc=None,
        sill_data=None,
        flux_data=None,
        transportboundary=False,
        conc_data=None,
        extension="lak",
        unitnumber=None,
        filenames=None,
        options=None,
        lwrt=0,
        **kwargs,
    ):
        # set default unit number of one is not specified
        if unitnumber is None:
            unitnumber = MfUsgLak._defaultunit()

        # set filenames
        tabdata = False
        nlen = 2
        if options is not None:
            for option in options:
                if "TABLEINPUT" in option.upper():
                    tabdata = True
                    nlen += nlakes
                    break
        filenames = self._prepare_filenames(filenames, nlen)

        # cbc output file
        self.set_cbc_output_file(ipakcb, model, filenames[1])

        # table input files
        if tabdata:
            if tab_files is None:
                tab_files = filenames[2:]

        # add tab_files as external files
        if tabdata:
            # TABLEINPUT requires exactly one tabfile (and one unit) per lake;
            # otherwise write_file() crashes with IndexError on iunit_tab[n].
            if len(tab_files) != nlakes:
                raise ValueError(
                    "TABLEINPUT requires exactly one tab_file per lake: got "
                    f"{len(tab_files)} tab_files for {nlakes} lakes."
                )
            # make sure tab_files are not None
            for idx, fname in enumerate(tab_files, 1):
                if fname is None:
                    raise ValueError(
                        f"a filename must be specified for the tabfile for lake {idx}"
                    )
            # if tab_units are supplied they must also be one per lake
            if tab_units is not None and len(tab_units) != nlakes:
                raise ValueError(
                    "TABLEINPUT requires exactly one tab_unit per lake: got "
                    f"{len(tab_units)} tab_units for {nlakes} lakes."
                )
            # set unit for tab files if not passed to __init__
            if tab_units is None:
                tab_units = []
                for idx in range(len(tab_files)):
                    tab_units.append(model.next_ext_unit())
            # add tabfiles as external files
            for iu, fname in zip(tab_units, tab_files):
                model.add_external(fname, iu)

        # call base package constructor
        super().__init__(
            model,
            extension=extension,
            name=self._ftype(),
            unit_number=unitnumber,
            filenames=filenames[0],
        )

        self._generate_heading()
        self.url = "lak.html"

        if options is None:
            options = []
        self.options = options
        self.nlakes = nlakes
        self.theta = theta
        self.nssitr = nssitr
        self.sscncr = sscncr
        self.surfdep = surfdep
        self.lwrt = lwrt

        if isinstance(stages, float):
            if self.nlakes == 1:
                stages = np.array([self.nlakes], dtype=float) * stages
            else:
                stages = np.ones(self.nlakes, dtype=float) * stages
        elif isinstance(stages, list):
            stages = np.array(stages)
        if stages.shape[0] != nlakes:
            err = f"stages shape should be ({nlakes}) but is only ({stages.shape[0]})."
            raise Exception(err)
        self.stages = stages
        if stage_range is None:
            stage_range = np.ones((nlakes, 2), dtype=float)
            stage_range[:, 0] = -10000.0
            stage_range[:, 1] = 10000.0
        else:
            if isinstance(stage_range, list):
                stage_range = np.array(stage_range)
            elif isinstance(stage_range, float):
                raise Exception(
                    f"stage_range should be a list or array of size ({nlakes}, 2)"
                )

        self.dis = utils_def.get_dis(model)
        if self.dis.steady[0]:
            if stage_range.shape != (nlakes, 2):
                raise Exception(
                    "stages shape should be ({},2) but is only {}.".format(
                        nlakes, stage_range.shape
                    )
                )
        self.stage_range = stage_range

        # tabfile data
        self.tabdata = tabdata
        self.iunit_tab = tab_units

        if lakarr is None and bdlknc is None:
            err = "lakarr and bdlknc must be specified"
            raise Exception(err)
        nrow, ncol, nlay, nper = self.parent.get_nrow_ncol_nlay_nper()
        self.lakarr = Transient3d(
            model, (nlay, nrow, ncol), np.int32, lakarr, name="lakarr_"
        )
        self.bdlknc = Transient3d(
            model, (nlay, nrow, ncol), np.float32, bdlknc, name="bdlknc_"
        )

        if sill_data is not None:
            if not isinstance(sill_data, dict):
                # a bare list of systems is taken as stress period 0
                sill_data = {0: sill_data}
            # validate + canonicalize (IC/lake numbers -> int, sills -> float) so
            # write_file always emits canonical values, never a raw 2.5 / "abc".
            sill_data = self._validate_sill_data(sill_data, nper)

        if flux_data is None:
            raise ValueError(
                "LAK requires flux_data (dataset 9: PRCPLK/EVAPLK/RNF/WTHDRW per "
                "lake, keyed by 0-based stress period); USG-T reads it for every "
                "lake, so a LAK file cannot be written without it."
            )
        if flux_data is not None:
            if not isinstance(flux_data, dict):
                # convert array to a dictionary
                try:
                    flux_data = {0: flux_data}
                except:
                    err = "flux_data must be a dictionary"
                    raise Exception(err)
            for key, value in flux_data.items():
                if isinstance(value, np.ndarray):
                    td = {}
                    for k in range(value.shape[0]):
                        td[k] = value[k, :].tolist()
                    flux_data[key] = td
                    if len(list(flux_data.keys())) != nlakes:
                        raise Exception(
                            f"flux_data dictionary must have {nlakes} entries"
                        )
                elif isinstance(value, float) or isinstance(value, int):
                    td = {}
                    for k in range(self.nlakes):
                        td[k] = (np.ones(6, dtype=float) * value).tolist()
                    flux_data[key] = td
                elif isinstance(value, dict):
                    try:
                        steady = self.dis.steady[key]
                    except:
                        steady = True
                    nlen = 4
                    if steady and key > 0:
                        nlen = 6
                    for k in range(self.nlakes):
                        if k not in value:
                            raise ValueError(
                                "flux_data for stress period {} is missing lake "
                                "{} (0-based); dataset 9a needs one entry per "
                                "lake ({} lakes).".format(key + 1, k, self.nlakes)
                            )
                        td = value[k]
                        if len(td) < nlen:
                            raise ValueError(
                                "flux_data entry for stress period {} lake {} "
                                "has {} values but dataset 9a needs at least {} "
                                "({}).".format(
                                    key + 1,
                                    k,
                                    len(td),
                                    nlen,
                                    "PRCPLK EVAPLK RNF WTHDRW SSMN SSMX"
                                    if nlen == 6
                                    else "PRCPLK EVAPLK RNF WTHDRW",
                                )
                            )

        self.flux_data = flux_data
        self.sill_data = sill_data

        # Keep the transportboundary flag and the dataset-1a options in sync:
        # the header keyword drives both the dataset-9b layout and what load
        # detects, so authoring with transportboundary=True must emit the
        # TRANSPORTBOUNDARY keyword (and vice-versa).
        has_tb_opt = any("TRANSPORTBOUNDARY" in str(o).upper() for o in self.options)
        transportboundary = bool(transportboundary) or has_tb_opt
        self.transportboundary = transportboundary
        if transportboundary and not has_tb_opt:
            self.options.append("TRANSPORTBOUNDARY")
        # TRANSPORTBOUNDARY is only valid with active transport: USG-T sets
        # ILKTRNSPT from IUNITGWT and rejects the option when transport is off.
        if transportboundary and model.mcomp <= 0:
            raise ValueError(
                "LAK TRANSPORTBOUNDARY requires active transport (mcomp>0); "
                "USG-T rejects the option when transport is inactive."
            )

        mcomp = model.mcomp
        if isinstance(clake, (int, float)):
            self.clake = [[clake] * mcomp for _ in range(self.nlakes)]
        elif clake is None:
            self.clake = [[0.0] * mcomp for _ in range(self.nlakes)]
        else:
            self.clake = clake
            if mcomp > 0 and (
                len(clake) != self.nlakes or any(len(row) != mcomp for row in clake)
            ):
                raise ValueError(
                    f"LAK clake must be nlakes x mcomp ({self.nlakes} x {mcomp})."
                )

        if conc_data is not None:
            if not isinstance(conc_data, dict):
                try:
                    conc_data = {0: conc_data}
                except:
                    err = "conc_data must be a dictionary"
                    raise Exception(err)
        if model.mcomp > 0 and conc_data is None:
            raise ValueError(
                "LAK with active transport (mcomp>0) requires conc_data "
                "(dataset 9b: lake concentrations per stress period)."
            )
        self.conc_data = conc_data

        # Validate dataset 9b wherever dataset 9 (flux_data) is written, so
        # authoring fails with a clear ValueError instead of a raw KeyError or
        # TypeError inside write_file().
        if mcomp > 0:
            for kper in self.flux_data:
                if kper not in self.conc_data:
                    raise ValueError(
                        f"conc_data is missing stress period {kper} (0-based); "
                        "dataset 9b is required wherever dataset 9 (flux_data) "
                        "is written."
                    )
                cd = self.conc_data[kper]
                if not isinstance(cd, dict):
                    raise ValueError(
                        f"conc_data[{kper}] must be a dict keyed by "
                        f"(lake, component); got {type(cd).__name__}."
                    )
                for n in range(self.nlakes):
                    wthdrw = self.flux_data[kper][n][3]
                    for icomp in range(mcomp):
                        if (n, icomp) not in cd:
                            raise ValueError(
                                f"conc_data[{kper}] is missing entry "
                                f"{(n, icomp)} (lake, component; 0-based); "
                                "dataset 9b needs one entry per lake and "
                                "component."
                            )
                        val = cd[n, icomp]
                        if transportboundary:
                            if isinstance(val, (list, tuple, np.ndarray)):
                                raise ValueError(
                                    f"conc_data[{kper}][{(n, icomp)}] must be a "
                                    "single CLAKE concentration with "
                                    "TRANSPORTBOUNDARY, not a sequence."
                                )
                        else:
                            needed = 3 if wthdrw < 0 else 2
                            ok = (
                                isinstance(val, (list, tuple, np.ndarray))
                                and len(val) == needed
                            )
                            if not ok:
                                caug = ", CAUG" if needed == 3 else ""
                                sign = "<" if wthdrw < 0 else ">="
                                raise ValueError(
                                    f"conc_data[{kper}][{(n, icomp)}] must have "
                                    f"{needed} values (CPPT, CRNF{caug}) because "
                                    f"WTHDRW {sign} 0; got {val!r}."
                                )

        self.parent.add_package(self)

        return

    def _ncells(self):
        """Maximum number of cells that can have lakes (developed for
        MT3DMS SSM package).

        Returns
        -------
        ncells: int
            maximum number of lak cells

        """
        nrow, ncol, nlay, nper = self.parent.nrow_ncol_nlay_nper
        return nlay * nrow * ncol

    @staticmethod
    def _canon_int(value, what):
        """Return ``value`` as a canonical ``int`` (for IC / lake numbers), or
        raise ``ValueError``. Integer-valued floats (``2.0``) and numeric strings
        (``"2"``) are normalized; ``bool``, ``None``, non-integral floats
        (``2.5``), and non-integer strings (``"abc"``, ``"2.5"``) are rejected.
        """
        if isinstance(value, bool):
            raise ValueError(f"{what} must be an integer, not a bool ({value!r}).")
        if isinstance(value, (int, np.integer)):
            return int(value)
        if isinstance(value, (float, np.floating)):
            if float(value).is_integer():
                return int(value)
            raise ValueError(f"{what} must be a whole number; got {value!r}.")
        if isinstance(value, str):
            try:
                return int(value.strip())
            except ValueError:
                raise ValueError(f"{what} must be an integer; got {value!r}.")
        raise ValueError(
            f"{what} must be an integer; got {type(value).__name__} {value!r}."
        )

    @staticmethod
    def _canon_float(value, what):
        """Return ``value`` as a canonical ``float`` (for sill elevations), or
        raise ``ValueError``. Numeric strings are normalized; ``bool``, ``None``,
        and non-numeric strings are rejected.
        """
        if isinstance(value, bool):
            raise ValueError(f"{what} must be numeric, not a bool ({value!r}).")
        if isinstance(value, (int, np.integer, float, np.floating)):
            return float(value)
        if isinstance(value, str):
            try:
                return float(value.strip())
            except ValueError:
                raise ValueError(f"{what} must be numeric; got {value!r}.")
        raise ValueError(
            f"{what} must be numeric; got {type(value).__name__} {value!r}."
        )

    def _validate_sill_data(self, sill_data, nper):
        """Validate **and canonicalize** datasets 7/8 (connected-lake / sill
        systems) before any file is opened, so authoring fails with an actionable
        ``ValueError`` instead of a raw error -- or a silently truncated /
        non-numeric value -- inside :meth:`write_file`. Returns a canonical copy
        ``{kper: [([IC:int, lake1:int, ...], [sill1:float, ...]), ...]}`` so the
        writer always emits integer ``IC``/lake numbers and float sills.

        USG-T (``gwf2lak7u1.f``) reads ``NSLMS`` (dataset 7) and, per connected
        system, dataset 8a ``IC ISUB(1..IC)`` (the center lake first, then its
        sublakes) followed by dataset 8b ``SILLVT(1..IC-1)`` -- but only when
        ``ITMP>0`` for that stress period (the same condition under which
        ``LKARR``/``BDLKNC`` are (re)read; ``ITMP<=0`` reuses the previous
        period and skips datasets 5--8). Lake numbers are 1-based, as written in
        the file. ``IC<=0`` is the Fortran end-of-list sentinel, so a system
        needs the center lake plus at least one sublake (``IC>=2``).
        """
        canon = {}
        for kper, systems in sill_data.items():
            if not isinstance(kper, (int, np.integer)) or not 0 <= kper < nper:
                raise ValueError(
                    f"MfUsgLak: sill_data stress period {kper} is out of range "
                    f"[0, {nper}); keys are 0-based stress periods."
                )
            # datasets 7/8 are only read when ITMP>0 (LKARR/BDLKNC (re)read);
            # otherwise write_file would silently drop the sill systems.
            if self.lakarr.get_kper_entry(kper)[0] <= 0:
                raise ValueError(
                    f"MfUsgLak: sill_data given for stress period {kper} but "
                    "lakarr is not (re)specified there (ITMP<=0); USG-T only "
                    "reads datasets 7/8 when ITMP>0. Provide lakarr for that "
                    "period or move the sill systems."
                )
            if not isinstance(systems, (list, tuple)):
                raise ValueError(
                    f"MfUsgLak: sill_data[{kper}] must be a list of "
                    "(ds8a, sillvt) systems; got "
                    f"{type(systems).__name__}."
                )
            canon_systems = []
            for isys, system in enumerate(systems):
                try:
                    ds8a, sillvt = system
                    ds8a = list(ds8a)
                    sillvt = list(sillvt)
                except (TypeError, ValueError):
                    raise ValueError(
                        f"MfUsgLak: sill_data[{kper}] system {isys} must be a "
                        "(ds8a, sillvt) pair where ds8a=[IC, lake1, ... lakeIC] "
                        "and sillvt=[sill1, ... sill_(IC-1)]."
                    )
                if len(ds8a) < 1:
                    raise ValueError(
                        f"MfUsgLak: sill_data[{kper}] system {isys} ds8a is "
                        "empty; it must be [IC, lake1, ... lakeIC]."
                    )
                where = f"sill_data[{kper}] system {isys}"
                ic = self._canon_int(ds8a[0], f"MfUsgLak: {where} IC")
                lakes = [
                    self._canon_int(x, f"MfUsgLak: {where} lake number")
                    for x in ds8a[1:]
                ]
                if ic < 2:
                    raise ValueError(
                        f"MfUsgLak: {where} has IC={ic}; a connected-lake system "
                        "needs the center lake plus at least one sublake (IC>=2)."
                    )
                if ic != len(lakes):
                    raise ValueError(
                        f"MfUsgLak: {where} declares IC={ic} but lists "
                        f"{len(lakes)} lake numbers (dataset 8a is IC followed "
                        "by IC lake numbers)."
                    )
                for lake in lakes:
                    if not 1 <= lake <= self.nlakes:
                        raise ValueError(
                            f"MfUsgLak: {where} lake number {lake} is out of "
                            f"range [1, {self.nlakes}] (lake numbers are "
                            "1-based)."
                        )
                if len(set(lakes)) != len(lakes):
                    raise ValueError(
                        f"MfUsgLak: {where} repeats a lake number ({lakes}); "
                        "each lake appears once per system."
                    )
                if len(sillvt) != ic - 1:
                    raise ValueError(
                        f"MfUsgLak: {where} has {len(sillvt)} sill elevations "
                        f"but dataset 8b needs IC-1={ic - 1} (one per sublake)."
                    )
                sills = [
                    self._canon_float(x, f"MfUsgLak: {where} sill elevation")
                    for x in sillvt
                ]
                canon_systems.append(([ic, *lakes], sills))
            canon[kper] = canon_systems
        return canon

    def write_file(self):
        """
        Write the package file.

        Returns
        -------
        None

        """
        f = open(self.fn_path, "w")
        # dataset 0
        f.write(f"{self.heading}\n")

        # dataset 1a
        if len(self.options) > 0:
            for option in self.options:
                f.write(f"{option} ")
        f.write("\n")

        # dataset 1b
        f.write(
            write_fixed_var(
                [self.nlakes, self.ipakcb], free=self.parent.free_format_input
            )
        )
        # dataset 2
        steady = np.any(self.dis.steady.array)
        t = [self.theta]
        if self.theta < 0.0 or steady:
            t.append(self.nssitr)
            t.append(self.sscncr)
        if self.theta < 0.0:
            t.append(self.surfdep)
        f.write(write_fixed_var(t, free=self.parent.free_format_input))

        # dataset 3
        mcomp = self.parent.mcomp

        steady = self.dis.steady[0]
        for n in range(self.nlakes):
            ipos = [10]
            t = [self.stages[n]]
            if steady:
                ipos.append(10)
                t.append(self.stage_range[n, 0])
                ipos.append(10)
                t.append(self.stage_range[n, 1])
            if mcomp > 0:
                for icomp in range(mcomp):
                    ipos.append(10)
                    t.append(self.clake[n][icomp])
            if self.tabdata:
                ipos.append(5)
                t.append(self.iunit_tab[n])
            f.write(write_fixed_var(t, ipos=ipos, free=self.parent.free_format_input))

        ds8_keys = list(self.sill_data.keys()) if self.sill_data is not None else []
        ds9_keys = list(self.flux_data.keys()) if self.flux_data is not None else []
        nper = self.dis.steady.shape[0]
        for kper in range(nper):
            itmp, file_entry_lakarr = self.lakarr.get_kper_entry(kper)
            ibd, file_entry_bdlknc = self.bdlknc.get_kper_entry(kper)

            itmp2 = 0
            if kper in ds9_keys:
                itmp2 = 1
            elif len(ds9_keys) > 0:
                itmp2 = -1
            if isinstance(self.lwrt, list):
                tmplwrt = self.lwrt[kper]
            else:
                tmplwrt = self.lwrt
            t = [itmp, itmp2, tmplwrt]
            comment = f"Stress period {kper + 1}"
            f.write(
                write_fixed_var(t, free=self.parent.free_format_input, comment=comment)
            )

            if itmp > 0:
                f.write(file_entry_lakarr)
                f.write(file_entry_bdlknc)

                nslms = 0
                if kper in ds8_keys:
                    ds8 = self.sill_data[kper]
                    nslms = len(ds8)

                f.write(
                    write_fixed_var(
                        [nslms],
                        length=5,
                        free=self.parent.free_format_input,
                        comment="Data set 7",
                    )
                )
                if nslms > 0:
                    for n in range(nslms):
                        d1, d2 = ds8[n]
                        s = write_fixed_var(
                            d1,
                            length=5,
                            free=self.parent.free_format_input,
                            comment="Data set 8a",
                        )
                        f.write(s)
                        s = write_fixed_var(
                            d2,
                            free=self.parent.free_format_input,
                            comment="Data set 8b",
                        )
                        f.write(s)

            if itmp2 > 0:
                ds9 = self.flux_data[kper]
                # dataset 9b (concentrations) only exists with active transport
                ds9b = self.conc_data[kper] if mcomp > 0 else None
                for n in range(self.nlakes):
                    try:
                        steady = self.dis.steady[kper]
                    except:
                        steady = True
                    if kper > 0 and steady:
                        t = ds9[n]
                    else:
                        t = ds9[n][0:4]
                    s = write_fixed_var(
                        t,
                        free=self.parent.free_format_input,
                        comment="Data set 9a",
                    )
                    f.write(s)
                    if mcomp > 0:
                        if self.transportboundary:
                            # USG-T reads one line per lake with CLAKE(1:NSOL).
                            t = [ds9b[n, icomp] for icomp in range(mcomp)]
                            f.write(
                                write_fixed_var(
                                    t,
                                    free=self.parent.free_format_input,
                                    comment="Data set 9b",
                                )
                            )
                        else:
                            # classic transport: one line per (lake, component)
                            # with CPPT, CRNF [, CAUG].
                            for icomp in range(mcomp):
                                f.write(
                                    write_fixed_var(
                                        ds9b[n, icomp],
                                        free=self.parent.free_format_input,
                                        comment="Data set 9b",
                                    )
                                )

        # close the lak file
        f.close()

    @classmethod
    def load(cls, f, model, nper=None, ext_unit_dict=None):
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

        Returns
        -------
        lak : ModflowLak object
            ModflowLak object.

        Examples
        --------

        >>> import flopy
        >>> m = flopy.mfusg.MfUsg()
        >>> lak = flopy.mfusg.MfUsgLak.load('test.lak', m)

        """

        cls.dis = utils_def.get_dis(model)

        if model.verbose:
            print("loading lak package file...")

        openfile = not hasattr(f, "read")
        if openfile:
            filename = f
            f = open(filename, "r", errors="replace")

        # dataset 0 -- header
        while True:
            line = f.readline()
            if line[0] != "#":
                break

        options = []
        tabdata = False
        newline = False
        transportboundary = False
        if "TABLEINPUT" in line.upper():
            if model.verbose:
                print("   reading lak dataset 1a")
            options.append("TABLEINPUT")
            tabdata = True
            newline = True
        if "TRANSPORTBOUNDARY" in line.upper():
            options.append("TRANSPORTBOUNDARY")
            transportboundary = True
            newline = True
        if newline:
            line = f.readline()

        # read dataset 1b
        if model.verbose:
            print("   reading lak dataset 1b")
        t = line.strip().split()
        # avoid error when there is a blank line
        if t == []:
            line = f.readline()
            t = line.strip().split()

        nlakes = int(t[0])
        ipakcb = 0
        try:
            ipakcb = int(t[1])
        except:
            pass

        # read dataset 2
        line = f.readline().rstrip()
        if model.array_free_format:
            t = line.split()
        else:
            t = read_fixed_var(line, ncol=4)
        theta = float(t[0])
        nssitr, sscncr = 0, 0.0
        if theta < 0:
            try:
                nssitr = int(t[1])
            except:
                if model.verbose:
                    print("  implicit nssitr defined in file")
            try:
                sscncr = float(t[2])
            except:
                if model.verbose:
                    print("  implicit sscncr defined in file")

        surfdep = 0.0
        if theta < 0.0:
            surfdep = float(t[3])

        if nper is None:
            nrow, ncol, nlay, nper = model.get_nrow_ncol_nlay_nper()

        mcomp = model.mcomp

        if model.verbose:
            print("   reading lak dataset 3")
        stages = []
        stage_range = []
        clake = []
        if tabdata:
            tab_units = []
        else:
            tab_units = None
        for lake in range(nlakes):
            line = f.readline().rstrip()
            if model.array_free_format:
                t = line.split()
            else:
                t = read_fixed_var(line, ipos=[10, 10, 10, 5])
            stages.append(t[0])
            ipos = 1
            if cls.dis.steady[0]:
                stage_range.append((float(t[ipos]), float(t[ipos + 1])))
                ipos += 2
            if mcomp > 0:
                conc = []
                for icomp in range(mcomp):
                    conc.append(float(t[ipos]))
                    ipos += 1
                clake.append(conc)
            if tabdata:
                iu = int(t[ipos])
                tab_units.append(iu)

        lake_loc = {}
        lake_lknc = {}
        sill_data = {}
        flux_data = {}
        conc_data = {}
        lwrt = []
        for iper in range(nper):
            if model.verbose:
                print(f"   reading lak dataset 4 - for stress period {iper + 1}")
            line = f.readline().rstrip()
            if model.array_free_format:
                t = line.split()
            else:
                t = read_fixed_var(line, ncol=3)
            itmp, itmp1, tmplwrt = int(t[0]), int(t[1]), int(t[2])
            lwrt.append(tmplwrt)

            if itmp > 0:
                if model.verbose:
                    print(f"   reading lak dataset 5 - for stress period {iper + 1}")
                name = f"LKARR_StressPeriod_{iper}"
                lakarr = Util3d.load(
                    f, model, (nlay, nrow, ncol), np.int32, name, ext_unit_dict
                )
                if model.verbose:
                    print(f"   reading lak dataset 6 - for stress period {iper + 1}")
                name = f"BDLKNC_StressPeriod_{iper}"
                bdlknc = Util3d.load(
                    f,
                    model,
                    (nlay, nrow, ncol),
                    np.float32,
                    name,
                    ext_unit_dict,
                )

                lake_loc[iper] = lakarr
                lake_lknc[iper] = bdlknc

                if model.verbose:
                    print(f"   reading lak dataset 7 - for stress period {iper + 1}")
                line = f.readline().rstrip()
                t = line.split()
                nslms = int(t[0])
                ds8 = []
                if nslms > 0:
                    if model.verbose:
                        print(
                            f"   reading lak dataset 8 - for stress period {iper + 1}"
                        )
                    for i in range(nslms):
                        line = f.readline().rstrip()
                        if model.array_free_format:
                            t = line.split()
                        else:
                            ic = int(line[0:5])
                            t = read_fixed_var(line, ncol=ic + 1, length=5)
                        ic = int(t[0])
                        ds8a = [ic]
                        for j in range(1, ic + 1):
                            ds8a.append(int(t[j]))
                        line = f.readline().rstrip()
                        if model.array_free_format:
                            t = line.split()
                        else:
                            t = read_fixed_var(line, ncol=ic - 1)
                        silvt = []
                        for j in range(ic - 1):
                            silvt.append(float(t[j]))
                        ds8.append((ds8a, silvt))
                    sill_data[iper] = ds8
            if itmp1 >= 0:
                if model.verbose:
                    print(f"   reading lak dataset 9 - for stress period {iper + 1}")
                ds9 = {}
                ds9b = {}
                for n in range(nlakes):
                    line = f.readline().rstrip()
                    if model.array_free_format:
                        t = line.split()
                    else:
                        t = read_fixed_var(line, ncol=6)
                    tds = []
                    tds.append(float(t[0]))  # PRCPLK
                    tds.append(float(t[1]))  # EVAPLK
                    tds.append(float(t[2]))  # RNF
                    tds.append(float(t[3]))  # WTHDRW

                    wthdrw = float(t[3])

                    if cls.dis.steady[iper]:
                        if iper == 0:
                            tds.append(stage_range[n][0])
                            tds.append(stage_range[n][1])
                        else:
                            tds.append(float(t[4]))
                            tds.append(float(t[5]))
                    else:
                        tds.append(0.0)
                        tds.append(0.0)
                    ds9[n] = tds

                    if mcomp > 0 and not transportboundary:
                        for icomp in range(mcomp):
                            line = f.readline().rstrip()
                            if model.array_free_format:
                                t = line.split()
                            else:
                                t = read_fixed_var(line, ncol=3)
                            tds = []
                            tds.append(float(t[0]))  # CPPT
                            tds.append(float(t[1]))  # CRNF
                            if wthdrw < 0:
                                tds.append(float(t[2]))  # CAUG
                            ds9b[n, icomp] = tds

                    if mcomp > 0 and transportboundary:
                        # CLAKBC: one line per lake with CLAKE(1:NSOL).
                        line = f.readline().rstrip()
                        t = line.split()
                        for icomp in range(mcomp):
                            ds9b[n, icomp] = float(t[icomp])

                flux_data[iper] = ds9
                conc_data[iper] = ds9b

        if openfile:
            f.close()

        # convert lake data to Transient3d objects
        lake_loc = Transient3d(
            model, (nlay, nrow, ncol), np.int32, lake_loc, name="lakarr_"
        )
        lake_lknc = Transient3d(
            model, (nlay, nrow, ncol), np.float32, lake_lknc, name="bdlknc_"
        )

        # determine specified unit number
        n = 2
        if tab_units is not None:
            n += nlakes
        unitnumber = None
        filenames = [None for x in range(n)]
        if ext_unit_dict is not None:
            unitnumber, filenames[0] = model.get_ext_dict_attr(
                ext_unit_dict, filetype=cls._ftype()
            )
            if ipakcb > 0:
                iu, filenames[1] = model.get_ext_dict_attr(ext_unit_dict, unit=ipakcb)
                model.add_pop_key_list(ipakcb)

            ipos = 2
            if tab_units is not None:
                for i in range(len(tab_units)):
                    iu, filenames[ipos] = model.get_ext_dict_attr(
                        ext_unit_dict, unit=tab_units[i]
                    )
                    ipos += 1

        return cls(
            model,
            options=options,
            nlakes=nlakes,
            ipakcb=ipakcb,
            theta=theta,
            nssitr=nssitr,
            surfdep=surfdep,
            sscncr=sscncr,
            lwrt=lwrt,
            stages=stages,
            stage_range=stage_range,
            clake=clake,
            tab_units=tab_units,
            lakarr=lake_loc,
            bdlknc=lake_lknc,
            sill_data=sill_data,
            flux_data=flux_data,
            transportboundary=transportboundary,
            conc_data=conc_data,
            unitnumber=unitnumber,
            filenames=filenames,
        )

    @staticmethod
    def _ftype():
        return "LAK"

    @staticmethod
    def _defaultunit():
        return 119
