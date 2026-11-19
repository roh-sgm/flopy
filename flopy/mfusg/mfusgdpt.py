"""
Mfusgdpt module.

Contains the MfUsgDpt class. Note that the user can
access the MfUsgDpt class as `flopy.mfusg.MfUsgDpt`.
"""

import warnings

import numpy as np

from ..pakbase import Package
from ..utils import Util2d, Util3d
from ..utils.utils_def import (
    get_unitnumber_from_ext_unit_dict,
    get_util2d_shape_for_layer,
)
from .mfusg import MfUsg


class MfUsgDpt(Package):
    """Dual Porosity Transport (dpt) Package Class for MODFLOW-USG Transport.

    Parameters
    ----------
    model : model object
        The model object (of type :class:`flopy.modflow.Modflow`) to which
        this package will be added.
    ipakcb : int (0,1,-1), (default is 0)
        a flag and a unit number >0 for cell-by-cell mass flux terms.
    idptcon : int (0,1), (default is 0)
        a flag and a unit number >0 for immobile domain concentrations
    icbndimflg : int (0,1), (default is 1)
        a flag active domain for the immobile (matrix) domain the same as
        that for the mobile (fracture) domain
    iadsorbim : int (0,1,2,3), (default is 0)
        a flag for adsorption in the immobile domain(0: no adsorption,
        1: linear isotherm, 2: Freundlich isotherm, 3: Langmuir isotherm)
    idispim : int (0,1), (default is 0)
        a flag for dispersion in the immobile domain (0: no dispersion,
        1: dispersion)
    izodim : int (0,1,2), (default is 0)
        a flag for zero-order decay in the immobile domain (0: no zero-order
        decay, 1: in water, 2: on soil, read only with adsorption). The values
        3 (water and soil) and 4 (air-water interface) are valid USG-T input
        that this class does not support yet (NotImplementedError, see Notes).
    ifodim : int (0,1,2), (default is 0)
        a flag for first-order decay in the immobile domain, with the same
        values and limits as izodim
    frahk : bool, (default is False)
        a flag for fractional hydraulic conductivity in the immobile domain.
        USG-T reads it only when the flow is not dual porosity (no DPF).
    mobilesat : bool, (default is False)
        immobile domain saturation equal to initial mobile domain saturation.
        Only without DPF; mobilesat and inputsat select one value of IMSAT and
        cannot both be True.
    inputsat : bool, (default is False)
        a flag of immobile domain saturation input (SIM array). Only without
        DPF.
    icbundim : int or array of ints (nlay, nrow, ncol)
        is  cell-by-cell flag for transport simulation in immobile domain
    phif : float or array of floats (nlay, nrow, ncol)
        fraction of the total space that is occupied by the mobile domain.
        Written and read only without DPF (with DPF, PHIF comes from DPF).
    prsityim : float or array of floats (nlay, nrow, ncol)
        effective transport porosity in the immobile domain
    bulkdim : float or array of floats (nlay, nrow, ncol)
        bulk density in the immobile domain
    dlim : float or array of floats (nlay, nrow, ncol)
        longitudinal dispersivity coefficient between mobile and immobile
        domains (read with DPF and idispim != 0)
    ddtr : float or array of floats (nlay, nrow, ncol)
        mass transfer coefficient between mobile and immobile domains
    sim : float or array of floats (nlay, nrow, ncol)
        saturation in the immobile domain
    htcapsim : float or array of floats (nlay, nrow, ncol)
        soil heat capacity in the immobile domain. Not used yet: heat
        transport (IHEAT=1) is not supported by this class.
    htcondsim : float or array of floats (nlay, nrow, ncol)
        deprecated and never written: USG-T reads a single heat array in DPT
        (SOIL HEAT CAPACITY). A value other than 0.0 emits a
        DeprecationWarning.
    adsorbim : float or array of floats (nlay, nrow, ncol)
        adsorption coefficient in the immobile domain
    flichim : float or array of floats (nlay, nrow, ncol)
        Freundlich coefficient in the immobile domain
    zodrwim : float or array of floats (nlay, nrow, ncol)
        zero-order decay rate in water in the immobile domain
    zodrsim : float or array of floats (nlay, nrow, ncol)
        zero-order decay rate on soil in the immobile domain
    fodrwim : float or array of floats (nlay, nrow, ncol)
        first-order decay rate in water in the immobile domain
    fodrsim : float or array of floats (nlay, nrow, ncol)
        first-order decay rate on soil in the immobile domain
    concim : float or sequence of mcomp + imcomp values
        initial concentration in the immobile domain of each component:
        the mcomp mobile components, then the imcomp immobile components
        (USG-T order). A scalar is used for every component.
    extension : string,  (default is 'dpt').
    unitnumber : int, default is 158.
        File unit number.
    filenames : str or list of str
        Filenames to use for the package and the output files.
    add_package : bool, default is True
        Flag to add the initialised package object to the parent model object.

    Methods
    -------

    See Also
    --------

    Notes
    -----
    DPT needs the component counts set by a BCT package (``mcomp``,
    ``iheat``, ``imcomp``; or counts assigned by hand with ``itrnsp`` not 0).
    They are taken when the package is built, and ``write_file`` raises
    ``ValueError`` before anything is written if they changed or the BCT
    package was removed or replaced by one with other counts.

    Whether the flow is dual porosity is the presence of a registered DPF
    package when writing (it decides whether the name file lists DPF), and
    ``model.idpf`` (set from the name file) when loading. ``write_file``
    raises ``ValueError`` if ``model.idpf`` disagrees with the registered
    packages, for example after the DPF package was removed.

    Not supported yet (``NotImplementedError``, before any array is written
    or read): heat transport (IHEAT=1), the SD_MASS_TR option, immobile
    components (imcomp > 0) with adsorption, decay or A-W_ADSORBIM, the
    SPATIALREACT arrays of BCT CHAINDECAY, and IZODIM/IFODIM 3 without
    adsorption or 4. These are limits of this class, not invalid USG-T
    input.

    ``write_file`` builds the whole file text before opening the file, so a
    rejected package leaves an existing file or a given handle untouched.
    ``write_file(f)`` and ``load(f)`` do not close a file object passed by
    the caller.

    Examples
    --------

    >>> import flopy
    >>> ml = flopy.mfusg.MfUsg()
    >>> disu = flopy.mfusg.MfUsgDisU(model=ml, nlay=1, nodes=1,
                 iac=[1], njag=1,ja=np.array([0]), fahl=[1.0], cl12=[1.0])
    >>> bct = flopy.mfusg.MfUsgBct(ml, mcomp=1)
    >>> dpt = flopy.mfusg.MfUsgDpt(ml)"""

    def __init__(
        self,
        model,
        ipakcb=0,
        idptcon=0,
        icbndimflg=1,
        iadsorbim=0,
        idispim=0,
        izodim=0,
        ifodim=0,
        frahk=False,
        mobilesat=False,
        inputsat=False,
        icbundim=1,
        phif=0.4,
        prsityim=0.4,
        bulkdim=1.6,
        dlim=0.5,
        ddtr=0.5,
        sim=0.5,
        htcapsim=0.0,
        htcondsim=0.0,
        adsorbim=0.0,
        flichim=0.0,
        zodrwim=0.0,
        zodrsim=0.0,
        fodrwim=0.0,
        fodrsim=0.0,
        concim=0.0,
        aw_adsorbim=False,
        iarea_fnim=1,
        ikawi_fnim=1,
        awamaxim=0.0,
        awarea_x2im=0.0,
        awarea_x1im=0.0,
        awarea_x0im=0.0,
        alangawim=0.0,
        blangawim=0.0,
        extension="dpt",
        unitnumber=None,
        filenames=None,
        add_package=True,
    ):
        """Constructs the MfUsgdpt object."""
        msg = (
            "Model object must be of type flopy.mfusg.MfUsg\n"
            f"but received type: {type(model)}."
        )
        assert isinstance(model, MfUsg), msg

        # set default unit number of one is not specified
        if unitnumber is None:
            unitnumber = self._defaultunit()

        super().__init__(
            model,
            extension,
            self._ftype(),
            unitnumber,
            self._prepare_filenames(filenames),
        )

        self._generate_heading()
        self.ipakcb = ipakcb
        self.idptcon = idptcon
        self.icbndimflg = icbndimflg
        self.iadsorbim = iadsorbim
        self.idispim = idispim
        self.izodim = izodim
        self.ifodim = ifodim
        # options
        self.frahk = frahk
        self.mobilesat = mobilesat
        self.inputsat = inputsat
        self.aw_adsorbim = bool(aw_adsorbim)

        # component counts of the model (BCT), checked again by write_file
        mcomp, iheat, imcomp = self._built_counts = self._context_counts(model)
        self._check_supported(
            model,
            iheat,
            imcomp,
            self.iadsorbim,
            self.izodim,
            self.ifodim,
            self.aw_adsorbim,
        )
        self._check_imsat(self.mobilesat, self.inputsat)
        if not (np.isscalar(htcondsim) and float(htcondsim) == 0.0):
            warnings.warn(
                "MfUsgDpt: htcondsim is deprecated and is not written; USG-T "
                "reads a single heat array in DPT (SOIL HEAT CAPACITY).",
                DeprecationWarning,
                stacklevel=2,
            )
        concim = self._normalize_concim(concim, mcomp + imcomp)

        nrow, ncol, nlay, nper = self.parent.nrow_ncol_nlay_nper

        if self.icbndimflg == 0:
            self.icbundim = Util3d(
                model, (nlay, nrow, ncol), np.int32, icbundim, name="icbundim"
            )

        self.phif = Util3d(model, (nlay, nrow, ncol), np.float32, phif, name="phif")

        self.prsityim = Util3d(
            model, (nlay, nrow, ncol), np.float32, prsityim, name="prsityim"
        )

        if self.iadsorbim:
            self.bulkdim = Util3d(
                model, (nlay, nrow, ncol), np.float32, bulkdim, name="bulkdim"
            )

        self.dlim = Util3d(model, (nlay, nrow, ncol), np.float32, dlim, name="dlim")

        self.ddtr = Util3d(model, (nlay, nrow, ncol), np.float32, ddtr, name="ddtr")

        if self.inputsat:
            self.sim = Util3d(model, (nlay, nrow, ncol), np.float32, sim, name="sim")

        if isinstance(adsorbim, (int, float)):
            adsorbim = [adsorbim] * mcomp
        if isinstance(flichim, (int, float)):
            flichim = [flichim] * mcomp
        if isinstance(zodrwim, (int, float)):
            zodrwim = [zodrwim] * mcomp
        if isinstance(zodrsim, (int, float)):
            zodrsim = [zodrsim] * mcomp
        if isinstance(fodrwim, (int, float)):
            fodrwim = [fodrwim] * mcomp
        if isinstance(fodrsim, (int, float)):
            fodrsim = [fodrsim] * mcomp

        self.adsorbim = [0] * mcomp
        self.flichim = [0] * mcomp
        self.zodrwim = [0] * mcomp
        self.zodrsim = [0] * mcomp
        self.fodrwim = [0] * mcomp
        self.fodrsim = [0] * mcomp

        for icomp in range(mcomp):
            if self.iadsorbim:
                self.adsorbim[icomp] = Util3d(
                    model,
                    (nlay, nrow, ncol),
                    np.float32,
                    adsorbim[icomp],
                    name="adsorbim",
                )

            if self.iadsorbim == 2 or self.iadsorbim == 3:
                self.flichim[icomp] = Util3d(
                    model,
                    (nlay, nrow, ncol),
                    np.float32,
                    flichim[icomp],
                    name="flichim",
                )

            if self.izodim == 1 or self.izodim == 3:
                self.zodrwim[icomp] = Util3d(
                    model,
                    (nlay, nrow, ncol),
                    np.float32,
                    zodrwim[icomp],
                    name="zodrwim",
                )

            if self.iadsorbim and (self.izodim == 2 or self.izodim == 3):
                self.zodrsim[icomp] = Util3d(
                    model,
                    (nlay, nrow, ncol),
                    np.float32,
                    zodrsim[icomp],
                    name="zodrsim",
                )

            if self.ifodim == 1 or self.ifodim == 3:
                self.fodrwim[icomp] = Util3d(
                    model,
                    (nlay, nrow, ncol),
                    np.float32,
                    fodrwim[icomp],
                    name="fodrwim",
                )

            if self.iadsorbim and (self.ifodim == 2 or self.ifodim == 3):
                self.fodrsim[icomp] = Util3d(
                    model,
                    (nlay, nrow, ncol),
                    np.float32,
                    fodrsim[icomp],
                    name="fodrsim",
                )

        # mobile components, then immobile components (GWT2DPTU1AR loop
        # ICOMP = 1..NTCOMP)
        self.concim = [
            Util3d(model, (nlay, nrow, ncol), np.float32, value, name="concim")
            for value in concim
        ]

        # Immobile-domain air-water interface adsorption (A-W_ADSORBIM).
        self.iarea_fnim = iarea_fnim
        self.ikawi_fnim = ikawi_fnim
        if self.aw_adsorbim:
            # canonical integer function indices (so write_file emits "1", never
            # "1.0" or a string) and supported-branch check
            self.iarea_fnim = iarea_fnim = self._canon_int(iarea_fnim, "IAREA_FNIM")
            self.ikawi_fnim = ikawi_fnim = self._canon_int(ikawi_fnim, "IKAWI_FNIM")
            self._check_aw_adsorbim_supported(iarea_fnim, ikawi_fnim)
            if mcomp <= 0:
                raise ValueError(
                    "MfUsgDpt: A-W_ADSORBIM is a per-species (Langmuir) option "
                    "and requires active transport (mcomp>0)."
                )
            # RP1 area arrays: IAREA_FNIM==1 reads AMAX; ==4 reads X2/X1/X0.
            if iarea_fnim == 1:
                self.awamaxim = Util3d(
                    model, (nlay, nrow, ncol), np.float32, awamaxim, name="awamaxim"
                )
            elif iarea_fnim == 4:
                self.awarea_x2im = Util3d(
                    model,
                    (nlay, nrow, ncol),
                    np.float32,
                    awarea_x2im,
                    name="awarea_x2im",
                )
                self.awarea_x1im = Util3d(
                    model,
                    (nlay, nrow, ncol),
                    np.float32,
                    awarea_x1im,
                    name="awarea_x1im",
                )
                self.awarea_x0im = Util3d(
                    model,
                    (nlay, nrow, ncol),
                    np.float32,
                    awarea_x0im,
                    name="awarea_x0im",
                )
            # RP2 Langmuir A/B arrays, one pair per mobile species.
            alangawim = self._normalize_species_arrays(alangawim, mcomp, "alangawim")
            blangawim = self._normalize_species_arrays(blangawim, mcomp, "blangawim")
            self.alangawim = [0] * mcomp
            self.blangawim = [0] * mcomp
            for icomp in range(mcomp):
                self.alangawim[icomp] = Util3d(
                    model,
                    (nlay, nrow, ncol),
                    np.float32,
                    alangawim[icomp],
                    name="alangawim",
                )
                self.blangawim[icomp] = Util3d(
                    model,
                    (nlay, nrow, ncol),
                    np.float32,
                    blangawim[icomp],
                    name="blangawim",
                )

        if add_package:
            self.parent.add_package(self)

    @staticmethod
    def _context_counts(model, who="MfUsgDpt"):
        """(mcomp, iheat, imcomp) of the model, or ValueError.

        GWT2DPTU1AR takes MCOMP, MCOMPT, NTCOMP and IHEAT from the BCT module
        (gwt2dptu1.f 2.8.1:12-14), so DPT needs a BCT package that set the
        model counts, or counts assigned by hand (itrnsp not 0, no BCT)."""
        state, counts = model._transport_context()
        if state == "none":
            raise ValueError(
                f"{who}: DPT needs the component counts of a BCT package; add "
                "a BCT package first (or set itrnsp and the counts by hand)."
            )
        if state == "removed":
            raise ValueError(
                f"{who}: the BCT package was removed; its component counts are "
                "stale. Add a BCT package again."
            )
        if state == "unsynced":
            bct = model.get_package("BCT")
            reason = bct._counts_mismatch() or (
                "it was added without setting the model counts"
            )
            raise ValueError(
                f"{who}: the registered BCT package does not match the model "
                f"component counts ({reason}); build or load the BCT package "
                "with its constructor."
            )
        for name, value in zip(("mcomp", "iheat", "imcomp"), counts):
            if isinstance(value, bool) or not isinstance(value, (int, np.integer)):
                raise ValueError(f"{who}: model {name}={value!r} is not an integer")
            if value < 0:
                raise ValueError(f"{who}: model {name}={value} is negative")
        return tuple(int(v) for v in counts)

    @staticmethod
    def _check_supported(model, iheat, imcomp, iadsorbim, izodim, ifodim, aw):
        """Raise NotImplementedError for valid USG-T input that this class
        does not support yet (before any array is created, written or read).
        Line numbers are gwt2dptu1.f of USG-T 2.8.1."""
        if iheat == 1:
            raise NotImplementedError(
                "MfUsgDpt: heat transport (IHEAT=1) is not supported yet. USG-T "
                "reads BULKDIM, a SOIL HEAT CAPACITY array and a TEMPERATURE "
                "array for DPT (gwt2dptu1.f 282, 340-353, 604-606)."
            )
        if imcomp > 0 and (iadsorbim or izodim or ifodim or aw):
            raise NotImplementedError(
                "MfUsgDpt: immobile components (IMCOMP>0) are supported only "
                "without adsorption, decay and A-W_ADSORBIM: USG-T 2.8.1 reads "
                "these arrays for immobile components too (gwt2dptu1.f 388-391) "
                "and USG-T 2.7.0 does not."
            )
        bct = model.get_package("BCT")
        if (
            bct is not None
            and getattr(bct, "chaindecay", None)
            and getattr(bct, "spatialreact", None)
        ):
            raise NotImplementedError(
                "MfUsgDpt: the SPTLRCTIM arrays read with BCT CHAINDECAY and "
                "SPATIALREACT (gwt2dptu1.f 369-385) are not supported yet."
            )
        for flag, value, kind in (
            ("IZODIM", izodim, "ZODRAWIM"),
            ("IFODIM", ifodim, "FODRAWIM"),
        ):
            if value == 3 and not iadsorbim:
                raise NotImplementedError(
                    f"MfUsgDpt: {flag}=3 without adsorption is not supported: "
                    "USG-T reads the decay arrays of option 3 only when "
                    "IADSORBIM>0 (gwt2dptu1.f 439-471, 520-552). Use 1 for "
                    "decay in water."
                )
            if value == 4:
                raise NotImplementedError(
                    f"MfUsgDpt: {flag}=4 is not supported: USG-T reads the decay "
                    f"arrays of option 4 (with {kind}) only with A-W_ADSORBIM "
                    "(gwt2dptu1.f 473-504, 554-585), and "
                    f"{kind} is not modeled."
                )

    @staticmethod
    def _check_imsat(mobilesat, inputsat):
        """MOBILESAT and INPUTSAT both set IMSAT; only one may be built."""
        if mobilesat and inputsat:
            raise ValueError(
                "MfUsgDpt: mobilesat and inputsat both select the immobile "
                "domain saturation (IMSAT); set only one of MOBILESAT or "
                "INPUTSAT."
            )

    def _check_write_arrays(self, mcomp, imcomp, idpf):
        """Check the current A-W_ADSORBIM indices and that every array read
        for the current flags and counts is present; return the canonical
        (IAREA_FNIM, IKAWI_FNIM), or None without A-W_ADSORBIM."""
        aw = None
        if self.aw_adsorbim:
            aw = (
                self._canon_int(self.iarea_fnim, "IAREA_FNIM"),
                self._canon_int(self.ikawi_fnim, "IKAWI_FNIM"),
            )
            self._check_aw_adsorbim_supported(*aw)
            if mcomp <= 0:
                raise ValueError(
                    "MfUsgDpt: A-W_ADSORBIM is a per-species (Langmuir) option "
                    "and requires active transport (mcomp>0)."
                )
        need = []  # (attribute, None for one array or the list length)
        if self.icbndimflg == 0:
            need.append(("icbundim", None))
        if not idpf:
            need.append(("phif", None))
        need.append(("prsityim", None))
        if self.iadsorbim:
            need.append(("bulkdim", None))
        if idpf and self.idispim:
            need.append(("dlim", None))
        need.append(("ddtr", None))
        if self.inputsat:
            need.append(("sim", None))
        if aw is not None:
            if aw[0] == 1:
                need.append(("awamaxim", None))
            else:
                need += [
                    (n, None) for n in ("awarea_x2im", "awarea_x1im", "awarea_x0im")
                ]
            need += [("alangawim", mcomp), ("blangawim", mcomp)]
        if self.iadsorbim:
            need.append(("adsorbim", mcomp))
        if self.iadsorbim == 2 or self.iadsorbim == 3:
            need.append(("flichim", mcomp))
        if self.izodim == 1 or self.izodim == 3:
            need.append(("zodrwim", mcomp))
        if self.iadsorbim and (self.izodim == 2 or self.izodim == 3):
            need.append(("zodrsim", mcomp))
        if self.ifodim == 1 or self.ifodim == 3:
            need.append(("fodrwim", mcomp))
        if self.iadsorbim and (self.ifodim == 2 or self.ifodim == 3):
            need.append(("fodrsim", mcomp))
        need.append(("concim", mcomp + imcomp))
        for name, n in need:
            value = getattr(self, name, None)
            if n is None:
                if isinstance(value, Util3d):
                    continue
                what = "a Util3d array"
            else:
                if (
                    isinstance(value, list)
                    and len(value) == n
                    and all(isinstance(v, Util3d) for v in value)
                ):
                    continue
                what = f"a list of {n} Util3d arrays (one per component)"
            got = type(value).__name__
            if isinstance(value, (list, tuple)):
                got += f" of length {len(value)}"
            raise ValueError(
                f"MfUsgDpt: {name} must be {what} for the current flags and "
                f"component counts, got {got}; build the DPT package with "
                "these flags (or assign the arrays) before writing."
            )
        return aw

    @staticmethod
    def _normalize_concim(value, ncomp):
        """CONCIM values for the ncomp = mcomp + imcomp components: a scalar
        is used for every component; a sequence needs ncomp entries."""
        if isinstance(value, (int, float, np.integer, np.floating)) and not (
            isinstance(value, bool)
        ):
            return [value] * ncomp
        if isinstance(value, np.ndarray) and value.ndim == 0:
            return [value.item()] * ncomp
        if isinstance(value, (list, tuple, np.ndarray)):
            seq = list(value)
            if len(seq) != ncomp:
                raise ValueError(
                    f"MfUsgDpt: concim has {len(seq)} entries but needs one per "
                    f"component (mcomp + imcomp = {ncomp}): mobile components "
                    "first, then immobile components."
                )
            return seq
        raise ValueError(
            "MfUsgDpt: concim must be a number or a sequence of "
            f"{ncomp} values; got {type(value).__name__}."
        )

    @staticmethod
    def _canon_int(value, what):
        """Return ``value`` as a canonical ``int`` for an A-W_ADSORBIM function
        index, or raise ``ValueError``. Integer-valued floats (``1.0``) and
        integer strings (``"1"``) are normalized; ``bool``, ``None``,
        non-integral floats (``1.5``), and non-integer strings are rejected -- so
        ``write_file`` emits e.g. ``"1"``, never ``"1.0"`` or a string.
        """
        if isinstance(value, bool):
            raise ValueError(
                f"MfUsgDpt A-W_ADSORBIM {what} must be an integer, not a bool "
                f"({value!r})."
            )
        if isinstance(value, (int, np.integer)):
            return int(value)
        if isinstance(value, (float, np.floating)):
            if float(value).is_integer():
                return int(value)
            raise ValueError(
                f"MfUsgDpt A-W_ADSORBIM {what} must be a whole number; got {value!r}."
            )
        if isinstance(value, str):
            try:
                return int(value.strip())
            except ValueError:
                raise ValueError(
                    f"MfUsgDpt A-W_ADSORBIM {what} must be an integer; got {value!r}."
                )
        raise ValueError(
            f"MfUsgDpt A-W_ADSORBIM {what} must be an integer; got "
            f"{type(value).__name__} {value!r}."
        )

    @staticmethod
    def _normalize_species_arrays(value, mcomp, name):
        """Normalize a per-species A-W_ADSORBIM array input to a list of length
        ``mcomp``: a scalar (Python/NumPy) is broadcast to every species; a
        list/tuple/ndarray must already have length ``mcomp``; ``None`` or a
        wrong length raises a clear ``ValueError`` (never a raw
        ``IndexError``/``TypeError`` in the per-species loop).
        """
        if isinstance(value, bool):
            raise ValueError(
                f"MfUsgDpt A-W_ADSORBIM {name} must be a number or a sequence of "
                f"{mcomp} per-species values, not a bool."
            )
        if isinstance(value, (int, float, np.integer, np.floating)):
            return [value] * mcomp
        # a 0-d ndarray is a scalar; broadcast it (avoid len() on an unsized array)
        if isinstance(value, np.ndarray) and value.ndim == 0:
            return [value.item()] * mcomp
        if value is None:
            raise ValueError(
                f"MfUsgDpt A-W_ADSORBIM {name} is required (a scalar or a "
                f"sequence of {mcomp} per-species values); got None."
            )
        if isinstance(value, (list, tuple, np.ndarray)):
            seq = list(value)
            if len(seq) != mcomp:
                raise ValueError(
                    f"MfUsgDpt A-W_ADSORBIM {name} has {len(seq)} entries but "
                    f"needs one per species (mcomp={mcomp})."
                )
            return seq
        raise ValueError(
            f"MfUsgDpt A-W_ADSORBIM {name} must be a number or a sequence of "
            f"{mcomp} per-species values; got {type(value).__name__}."
        )

    @staticmethod
    def _check_aw_adsorbim_supported(iarea_fnim, ikawi_fnim):
        """Reject A-W_ADSORBIM branches that FloPy does not model, with an
        actionable ``NotImplementedError`` (before any array is written or read,
        so later items are never silently shifted).

        Implemented (array-only, no zone map / scalar constants / tables):
        ``IAREA_FNIM`` in {1 (AMAX), 4 (X2/X1/X0)} and ``IKAWI_FNIM`` in {1, 2}
        (Langmuir A/B arrays per species). The deferred branches read extra
        scalars (``ROG_SIGMA``/``SIGMA_RT``) or a zone map plus tabular functions
        (``IAREA_FNIM==5``, ``IKAWI_FNIM==4``); see ``dpt2aw_adsorb.f``.
        """
        if iarea_fnim not in (1, 4):
            raise NotImplementedError(
                f"MfUsgDpt A-W_ADSORBIM IAREA_FNIM={iarea_fnim} is not supported. "
                "Supported: 1 (AMAX array) and 4 (X2/X1/X0 arrays). IAREA_FNIM=2 "
                "(grain diameter), 3 (ROG_SIGMA*porosity), and 5 (tabular "
                "area-vs-saturation with a zone map) read extra scalars/tables "
                "and are not modeled."
            )
        if ikawi_fnim not in (1, 2):
            raise NotImplementedError(
                f"MfUsgDpt A-W_ADSORBIM IKAWI_FNIM={ikawi_fnim} is not supported. "
                "Supported: 1 and 2 (Langmuir A/B arrays per species). "
                "IKAWI_FNIM=3 (Brusseau, needs the SIGMA_RT scalar) and 4 "
                "(tabular K_AWI-vs-concentration with a zone map) are not "
                "modeled."
            )

    def _write_idpf(self):
        """IDPF of the written model: 1 when a DPF package is registered (the
        name file then lists DPF). ``model.idpf`` must agree with it."""
        model = self.parent
        has_dpf = model.get_package("DPF") is not None
        if bool(model.idpf) != has_dpf:
            if has_dpf:
                reason = "a DPF package is registered but model.idpf is 0"
            else:
                reason = (
                    f"model.idpf={model.idpf} but no DPF package is registered "
                    "(the DPF package was removed, not loaded, or idpf was set "
                    "by hand)"
                )
            raise ValueError(
                f"MfUsgDpt: {reason}. USG-T reads PHIF and the FRAHK, "
                "MOBILESAT and INPUTSAT options only without DPF (gwt2dptu1.f "
                "2.8.1:99, 263); add or remove the DPF package, or set "
                "model.idpf to match."
            )
        if has_dpf:
            used = [
                name
                for name, value in (
                    ("FRAHK", self.frahk),
                    ("MOBILESAT", self.mobilesat),
                    ("INPUTSAT", self.inputsat),
                )
                if value
            ]
            if used:
                raise ValueError(
                    f"MfUsgDpt: option(s) {', '.join(used)} apply only when the "
                    "flow is not dual porosity; USG-T does not read them with "
                    "DPF (gwt2dptu1.f 2.8.1:99)."
                )
        return int(has_dpf)

    def _file_text(self):
        """Validate and build the whole DPT file text (nothing is opened)."""
        counts = self._context_counts(self.parent)
        if counts != self._built_counts:
            raise ValueError(
                "MfUsgDpt: the model component counts (mcomp, iheat, imcomp) "
                f"changed from {self._built_counts} to {counts} after this "
                "package was built (for example a BCT package with other counts "
                "replaced it); build the DPT package again."
            )
        mcomp, iheat, imcomp = counts
        # the construction checks again, on the current flags, options and BCT
        self._check_supported(
            self.parent,
            iheat,
            imcomp,
            self.iadsorbim,
            self.izodim,
            self.ifodim,
            self.aw_adsorbim,
        )
        self._check_imsat(self.mobilesat, self.inputsat)
        idpf = self._write_idpf()
        aw = self._check_write_arrays(mcomp, imcomp, idpf)

        # Item 1: IDPTCB IDPTCON IC_BNDIM_FLG IADSORBIM IDISPIM IZODIM IFODIM,
        # 10 columns each (fixed format reads options from column 71)
        line = (
            f" {self.ipakcb:9d} {self.idptcon:9d} {self.icbndimflg:9d}"
            f" {self.iadsorbim:9d} {self.idispim:9d} {self.izodim:9d}"
            f" {self.ifodim:9d}"
        )
        if not idpf:
            if self.frahk:
                line += " FRAHK"
            if self.mobilesat:
                line += " MOBILESAT"
            if self.inputsat:
                line += " INPUTSAT"
        if aw is not None:
            line += f" A-W_ADSORBIM {aw[0]} {aw[1]}"
        entries = [line + "\n"]

        # ICBUNDIM
        if self.icbndimflg == 0:
            entries.append(self.icbundim.get_file_entry())
        # PHIF: IF(IDPF.EQ.0)
        if not idpf:
            entries.append(self.phif.get_file_entry())
        # PRSITYIM
        entries.append(self.prsityim.get_file_entry())
        # BULKDIM: IF(IADSORBIM.NE.0.OR.IHEAT.EQ.1) (heat is not supported)
        if self.iadsorbim:
            entries.append(self.bulkdim.get_file_entry())
        # DLIM: IF(IDPF.NE.0.AND.IDISPIM.NE.0)
        if idpf and self.idispim:
            entries.append(self.dlim.get_file_entry())
        # DDTTR
        entries.append(self.ddtr.get_file_entry())
        # SIM: IDPF=0 and INPUTSAT
        if self.inputsat:
            entries.append(self.sim.get_file_entry())

        # A-W_ADSORBIM RP1: area arrays (read after heat, before the species
        # loop). IAREA_FNIM==1 -> AMAX; ==4 -> X2, X1, X0.
        if aw is not None:
            if aw[0] == 1:
                entries.append(self.awamaxim.get_file_entry())
            else:  # IAREA_FNIM 4 (checked)
                entries.append(self.awarea_x2im.get_file_entry())
                entries.append(self.awarea_x1im.get_file_entry())
                entries.append(self.awarea_x0im.get_file_entry())

        # mobile components: ADSORBIM, FLICHIM, ZODRWIM, ZODRSIM, FODRWIM,
        # FODRSIM, CONCIM
        for icomp in range(mcomp):
            # A-W_ADSORBIM RP2: Langmuir A/B per species, before ADSORBIM.
            if aw is not None:
                entries.append(self.alangawim[icomp].get_file_entry())
                entries.append(self.blangawim[icomp].get_file_entry())
            if self.iadsorbim:
                entries.append(self.adsorbim[icomp].get_file_entry())
            if self.iadsorbim == 2 or self.iadsorbim == 3:
                entries.append(self.flichim[icomp].get_file_entry())
            if self.izodim == 1 or self.izodim == 3:
                entries.append(self.zodrwim[icomp].get_file_entry())
            if self.iadsorbim and (self.izodim == 2 or self.izodim == 3):
                entries.append(self.zodrsim[icomp].get_file_entry())
            if self.ifodim == 1 or self.ifodim == 3:
                entries.append(self.fodrwim[icomp].get_file_entry())
            if self.iadsorbim and (self.ifodim == 2 or self.ifodim == 3):
                entries.append(self.fodrsim[icomp].get_file_entry())
            entries.append(self.concim[icomp].get_file_entry())

        # immobile components: CONCIM only (no reactions, checked when built)
        for icomp in range(mcomp, mcomp + imcomp):
            entries.append(self.concim[icomp].get_file_entry())
        return "".join(entries)

    def write_file(self, f=None):
        """
        Write the dpt package file.

        Parameters
        ----------
        f : open file object.
            Default is None, which will result in MfUsg.fn_path being
            opened for writing. A file object is written from its current
            position and is not closed.

        The package is validated on its current flags, options, BCT package
        and arrays, and the whole text is built before the file is opened, so
        a ``ValueError`` or ``NotImplementedError`` leaves an existing file
        (or ``f``) untouched.
        """
        text = self._file_text()
        if f is None:
            with open(self.fn_path, "w") as f_obj:
                f_obj.write(text)
        else:
            f.write(text)

    @classmethod
    def load(cls, f, model, ext_unit_dict=None):
        """
        Load an existing package.

        Parameters
        ----------
        f : filename or file handle
            File to load. A file handle is not closed.
        model : model object
            The model object (of type :class:`flopy.modflow.mf.Modflow`) to
            which this package will be added. Its BCT component counts and
            ``model.idpf`` (set from the name file by ``MfUsg.load``) decide
            which arrays are read.
        ext_unit_dict : dictionary, optional
            If the arrays in the file are specified using EXTERNAL,
            or older style array control records, then `f` should be a file
            handle.  In this case ext_unit_dict is required, which can be
            constructed using the function
            :class:`flopy.utils.mfreadnam.parsenamefile`.

        Returns
        -------
        dpt : MfUsgdpt object

        Examples
        --------

        >>> import flopy
        >>> ml = flopy.mfusg.MfUsg()
        >>> dis = flopy.mfusg.MfUsgDisU.load('SeqDegEg.dis', ml)
        >>> dpt = flopy.mfusg.MfUsgdpt.load('SeqDegEg.dpt', ml)
        """
        msg = (
            "Model object must be of type flopy.mfusg.MfUsg\n"
            f"but received type: {type(model)}."
        )
        assert isinstance(model, MfUsg), msg

        if model.verbose:
            print("loading dpt package file...")

        openfile = not hasattr(f, "read")
        f_obj = open(f) if openfile else f
        try:
            kwargs = cls._load_kwargs(f_obj, model, ext_unit_dict)
        finally:
            if openfile:
                f_obj.close()

        # set package unit number
        unitnumber, filenames = get_unitnumber_from_ext_unit_dict(
            model, cls, ext_unit_dict, kwargs["ipakcb"]
        )

        return cls(model, unitnumber=unitnumber, filenames=filenames, **kwargs)

    @classmethod
    def _load_kwargs(cls, f_obj, model, ext_unit_dict):
        """Read the DPT records in the GWT2DPTU1AR order (gwt2dptu1.f)."""
        mcomp, iheat, imcomp = cls._context_counts(model, "MfUsgDpt.load")
        idpf = model.idpf
        nlay = model.nlay

        # item 1
        line = f_obj.readline()
        while line.startswith("#"):
            line = f_obj.readline()
        # URWORD reads the words in upper case
        t = line.upper().split()
        kwargs = {}

        # item 1a
        vars = {
            "ipakcb": int,
            "idptcon": int,
            "icbndimflg": int,
            "iadsorbim": int,
            "idispim": int,
            "izodim": int,
            "ifodim": int,
        }

        for i, (v, c) in enumerate(vars.items()):
            kwargs[v] = c(t[i].strip())

        # item 1a - FRAHK, MOBILESAT, INPUTSAT are read only without DPF; the
        # last of MOBILESAT and INPUTSAT sets IMSAT (gwt2dptu1.f 2.8.1:99-115)
        kwargs["frahk"] = 0
        kwargs["mobilesat"] = 0
        kwargs["inputsat"] = 0
        if not idpf:
            words = t[len(vars) :]
            kwargs["frahk"] = int("FRAHK" in words)
            imsat = [w for w in words if w in ("MOBILESAT", "INPUTSAT")]
            if imsat:
                kwargs["mobilesat"] = int(imsat[-1] == "MOBILESAT")
                kwargs["inputsat"] = int(imsat[-1] == "INPUTSAT")

        if "SD_MASS_TR" in t:
            raise NotImplementedError(
                "MfUsgDpt.load: the SD_MASS_TR option (species-dependent mass "
                "transfer rates DDTTRSD, gwt2dptu1.f 2.8.1:153-154, 590-600) is "
                "not supported yet."
            )

        # Immobile-domain air-water interface adsorption (A-W_ADSORBIM): the
        # keyword is followed by IAREA_FNIM and IKAWI_FNIM on the option line
        # (gwt2dptu1.f). Supported array-only branches are read below; the others
        # fail explicitly (before any array) rather than silently shifting reads.
        kwargs["aw_adsorbim"] = False
        if "A-W_ADSORBIM" in t:
            idx = t.index("A-W_ADSORBIM")
            try:
                a_tok, k_tok = t[idx + 1], t[idx + 2]
            except IndexError:
                raise ValueError(
                    "MfUsgDpt: A-W_ADSORBIM must be followed by two integers "
                    "IAREA_FNIM IKAWI_FNIM on the option line."
                )
            iarea_fnim = cls._canon_int(a_tok, "IAREA_FNIM")
            ikawi_fnim = cls._canon_int(k_tok, "IKAWI_FNIM")
            cls._check_aw_adsorbim_supported(iarea_fnim, ikawi_fnim)
            # A-W_ADSORBIM is a per-species option: it needs active transport.
            # Fail here (option line), before any RP1/RP2 array read, so a
            # mcomp<=0 model never hits an EOF/IndexError mid-array.
            if model.mcomp <= 0:
                raise ValueError(
                    "MfUsgDpt: A-W_ADSORBIM requires active transport (mcomp>0) "
                    "to read the per-species Langmuir arrays; model.mcomp="
                    f"{model.mcomp}."
                )
            kwargs["aw_adsorbim"] = True
            kwargs["iarea_fnim"] = iarea_fnim
            kwargs["ikawi_fnim"] = ikawi_fnim

        # unsupported input is rejected before any array is read
        cls._check_supported(
            model,
            iheat,
            imcomp,
            kwargs["iadsorbim"],
            kwargs["izodim"],
            kwargs["ifodim"],
            kwargs["aw_adsorbim"],
        )

        def read(name, dtype=np.float32):
            return cls._load_prop_arrays(f_obj, model, nlay, dtype, name, ext_unit_dict)

        # ICBUNDIM: IF(IC_BNDIM_FLG.EQ.0)
        if kwargs["icbndimflg"] == 0:
            kwargs["icbundim"] = read("icbundim", np.int32)
        # PHIF: IF(IDPF.EQ.0)
        if not idpf:
            kwargs["phif"] = read("phif")
        kwargs["prsityim"] = read("prsityim")
        # BULKDIM: IF(IADSORBIM.NE.0.OR.IHEAT.EQ.1) (heat rejected above)
        if kwargs["iadsorbim"]:
            kwargs["bulkdim"] = read("bulkdim")
        # DLIM: IF(IDPF.NE.0.AND.IDISPIM.NE.0)
        if idpf and kwargs["idispim"]:
            kwargs["dlim"] = read("dlim")
        kwargs["ddtr"] = read("ddtr")
        # SIM: IF(IDPF.EQ.0) and IMSAT=3
        if kwargs["inputsat"]:
            kwargs["sim"] = read("sim")

        # A-W_ADSORBIM RP1: area arrays (after heat, before the species loop).
        if kwargs["aw_adsorbim"]:
            if kwargs["iarea_fnim"] == 1:
                kwargs["awamaxim"] = read("awamaxim")
            elif kwargs["iarea_fnim"] == 4:
                kwargs["awarea_x2im"] = read("awarea_x2im")
                kwargs["awarea_x1im"] = read("awarea_x1im")
                kwargs["awarea_x0im"] = read("awarea_x0im")

        # species loop, mobile components
        names = ("adsorbim", "flichim", "zodrwim", "zodrsim", "fodrwim", "fodrsim")
        species = {name: [0] * mcomp for name in names}
        alangawim = [0] * mcomp
        blangawim = [0] * mcomp
        concim = []
        ads, izod, ifod = kwargs["iadsorbim"], kwargs["izodim"], kwargs["ifodim"]
        for icomp in range(mcomp):
            # A-W_ADSORBIM RP2: Langmuir A/B per species, before ADSORBIM.
            if kwargs["aw_adsorbim"]:
                alangawim[icomp] = read("alangawim")
                blangawim[icomp] = read("blangawim")
            if ads:
                species["adsorbim"][icomp] = read("adsorbim")
            if ads == 2 or ads == 3:
                species["flichim"][icomp] = read("flichim")
            if izod == 1 or izod == 3:
                species["zodrwim"][icomp] = read("zodrwim")
            if ads and (izod == 2 or izod == 3):
                species["zodrsim"][icomp] = read("zodrsim")
            if ifod == 1 or ifod == 3:
                species["fodrwim"][icomp] = read("fodrwim")
            if ads and (ifod == 2 or ifod == 3):
                species["fodrsim"][icomp] = read("fodrsim")
            concim.append(read("concim"))
        # immobile components: CONCIM only (reactions rejected above)
        for _ in range(imcomp):
            concim.append(read("concim"))

        kwargs.update(species)
        kwargs["concim"] = concim
        if kwargs["aw_adsorbim"]:
            kwargs["alangawim"] = alangawim
            kwargs["blangawim"] = blangawim
        return kwargs

    @staticmethod
    def _load_prop_arrays(f_obj, model, nlay, dtype, name, ext_unit_dict):
        if model.verbose:
            print(f"   loading {name} ...")
        prop_array = [0] * nlay
        for layer in range(nlay):
            util2d_shape = get_util2d_shape_for_layer(model, layer=layer)
            prop_array[layer] = Util2d.load(
                f_obj, model, util2d_shape, dtype, name, ext_unit_dict
            )
        return prop_array

    @staticmethod
    def _ftype():
        return "DPT"

    @staticmethod
    def _defaultunit():
        return 158
