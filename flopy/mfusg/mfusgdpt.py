"""
Mfusgdpt module.

Contains the MfUsgDpt class. Note that the user can
access the MfUsgDpt class as `flopy.mfusg.MfUsgDpt`.
"""

import numpy as np

from ..pakbase import Package
from ..utils import Util2d, Util3d
from ..utils.utils_def import (
    get_open_file_object,
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
    izodim : int (0,1), (default is 0)
        a flag for zero-order decay in the immobile domain (0: no zero-order decay,
        1: in water, 2: on soil, 3: on water and soil,4: on air-water interface)
    ifodim : int (0,1), (default is 0)
        a flag for first-order decay in the immobile domain (0: no first-order decay
        1: in water, 2: on soil, 3: on water and soil,4: on air-water interface)
    frahk : bool, (default is False)
        a flag for fractional hydraulic conductivity in the immobile domain
    mobilesat : bool, (default is False)
        immobile domain saturation equal to initial mobile domain saturation
    inputsat : bool, (default is False)
        a flag of immobile domain saturation input
    icbundim : int or array of ints (nlay, nrow, ncol)
        is  cell-by-cell flag for transport simulation in immobile domain
    phif : float or array of floats (nlay, nrow, ncol)
        fraction of the total space that is occupied by the mobile domain
    prsityim : float or array of floats (nlay, nrow, ncol)
        effective transport porosity in the immobile domain
    bulkdim : float or array of floats (nlay, nrow, ncol)
        bulk density in the immobile domain
    dlim : float or array of floats (nlay, nrow, ncol)
        longitudinal dispersivity coefficient between mobile and immobile domains
    ddtr : float or array of floats (nlay, nrow, ncol)
        mass transfer coefficient between mobile and immobile domains
    sim : float or array of floats (nlay, nrow, ncol)
        saturation in the immobile domain
    htcapsim : float or array of floats (nlay, nrow, ncol)
        heat capacity in the immobile domain
    htcondsim : float or array of floats (nlay, nrow, ncol)
        heat conductivity in the immobile domain
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
    concim : float or array of floats (nlay, nrow, ncol)
        initial concentration in the immobile domain
    extension : string,  (default is 'dpt').
    unitnumber : int, default is 58.
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

    Examples
    --------

    >>> import flopy
    >>> ml = flopy.mfusg.MfUsg()
    >>> disu = flopy.mfusg.MfUsgDisU(model=ml, nlay=1, nodes=1,
                 iac=[1], njag=1,ja=np.array([0]), fahl=[1.0], cl12=[1.0])
    >>> dpt = flopy.mfusg.MfUsgdpt(ml)"""

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
            self.unitnumber = self._defaultunit()

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

        if model.iheat:
            self.htcapsim = Util3d(
                model, (nlay, nrow, ncol), np.float32, htcapsim, name="htcapsim"
            )

            self.htcondsim = Util3d(
                model, (nlay, nrow, ncol), np.float32, htcondsim, name="htcondsim"
            )

        mcomp = model.mcomp

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

        if isinstance(concim, (int, float)):
            concim = [concim] * mcomp

        self.adsorbim = [0] * mcomp
        self.flichim = [0] * mcomp
        self.zodrwim = [0] * mcomp
        self.zodrsim = [0] * mcomp
        self.fodrwim = [0] * mcomp
        self.fodrsim = [0] * mcomp
        self.concim = [0] * mcomp

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

            if self.izodim == 1 or self.izodim == 3 or self.izodim == 4:
                self.zodrwim[icomp] = Util3d(
                    model,
                    (nlay, nrow, ncol),
                    np.float32,
                    zodrwim[icomp],
                    name="zodrwim",
                )

            if self.iadsorbim and (
                self.izodim == 2 or self.izodim == 3 or self.izodim == 4
            ):
                self.zodrsim[icomp] = Util3d(
                    model,
                    (nlay, nrow, ncol),
                    np.float32,
                    zodrsim[icomp],
                    name="zodrsim",
                )

            if self.ifodim == 1 or self.ifodim == 3 or self.ifodim == 4:
                self.fodrwim[icomp] = Util3d(
                    model,
                    (nlay, nrow, ncol),
                    np.float32,
                    fodrwim[icomp],
                    name="fodrwim",
                )

            if self.iadsorbim and (
                self.ifodim == 2 or self.ifodim == 3 or self.ifodim == 4
            ):
                self.fodrsim[icomp] = Util3d(
                    model,
                    (nlay, nrow, ncol),
                    np.float32,
                    fodrsim[icomp],
                    name="fodrsim",
                )

            self.concim[icomp] = Util3d(
                model, (nlay, nrow, ncol), np.float32, concim[icomp], name="concim"
            )

        # Immobile-domain air-water interface adsorption (A-W_ADSORBIM).
        self.aw_adsorbim = bool(aw_adsorbim)
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

    def write_file(self, f=None):
        """
        Write the dpt package file.

        Parameters
        ----------
        f : open file object.
            Default is None, which will result in MfUsg.fn_path being
            opened for writing.

        Examples
        --------
        """
        # Open file for writing
        if f is None:
            f_obj = open(self.fn_path, "w")

        #        f_obj.write(f"{self.heading}\n")

        # Item 0: IPAKCB, IDPTCON
        f_obj.write(
            f" {self.ipakcb:9d} {self.idptcon:9d} {self.icbndimflg:9d}"
            f" {self.iadsorbim:9d} {self.idispim:9d} {self.izodim:9d} {self.ifodim:9d}"
        )

        # Options
        if self.frahk:
            f_obj.write(" FRAHK")

        if self.mobilesat:
            f_obj.write(" MOBILESAT")

        if self.inputsat:
            f_obj.write(" INPUTSAT")

        if self.aw_adsorbim:
            f_obj.write(f" A-W_ADSORBIM {self.iarea_fnim} {self.ikawi_fnim}")

        f_obj.write("\n")

        # Item 1: ICBUNDIM
        if self.icbndimflg == 0:
            f_obj.write(self.icbundim.get_file_entry())

        # Item 2: PHIF
        f_obj.write(self.phif.get_file_entry())

        # Item 3: PRSITYIM
        f_obj.write(self.prsityim.get_file_entry())

        # Item 4: BULKDIM
        if self.iadsorbim:
            f_obj.write(self.bulkdim.get_file_entry())

        # Item 5: DLIM — Fortran: IF(IDPF.NE.0.AND.IDISPIM.NE.0)
        if self.parent.idpf and self.idispim:
            f_obj.write(self.dlim.get_file_entry())

        # Item 6: DDTR
        f_obj.write(self.ddtr.get_file_entry())

        # Item 7: SIM
        if self.inputsat:
            f_obj.write(self.sim.get_file_entry())

        # Item 8: HTCAPSIM, HTCONDSIM
        if self.parent.iheat == 1:
            f_obj.write(self.htcapsim.get_file_entry())
            f_obj.write(self.htcondsim.get_file_entry())

        # A-W_ADSORBIM RP1: area arrays (read after heat, before the species
        # loop). IAREA_FNIM==1 -> AMAX; ==4 -> X2, X1, X0.
        if self.aw_adsorbim:
            if self.iarea_fnim == 1:
                f_obj.write(self.awamaxim.get_file_entry())
            elif self.iarea_fnim == 4:
                f_obj.write(self.awarea_x2im.get_file_entry())
                f_obj.write(self.awarea_x1im.get_file_entry())
                f_obj.write(self.awarea_x0im.get_file_entry())

        # Item 9: ADSORBIM, FLICHIM, ZODRWIM, ZODRSIM, FODRWIM, FODRSIM, CONCIM
        mcomp = self.parent.mcomp
        for icomp in range(mcomp):
            # A-W_ADSORBIM RP2: Langmuir A/B per species, before ADSORBIM.
            if self.aw_adsorbim:
                f_obj.write(self.alangawim[icomp].get_file_entry())
                f_obj.write(self.blangawim[icomp].get_file_entry())

            if self.iadsorbim:
                f_obj.write(self.adsorbim[icomp].get_file_entry())

            if self.iadsorbim == 2 or self.iadsorbim == 3:
                f_obj.write(self.flichim[icomp].get_file_entry())

            if self.izodim == 1 or self.izodim == 3 or self.izodim == 4:
                f_obj.write(self.zodrwim[icomp].get_file_entry())

            if self.iadsorbim and (
                self.izodim == 2 or self.izodim == 3 or self.izodim == 4
            ):
                f_obj.write(self.zodrsim[icomp].get_file_entry())

            if self.ifodim == 1 or self.ifodim == 3 or self.ifodim == 4:
                f_obj.write(self.fodrwim[icomp].get_file_entry())

            if self.iadsorbim and (
                self.ifodim == 2 or self.ifodim == 3 or self.ifodim == 4
            ):
                f_obj.write(self.fodrsim[icomp].get_file_entry())

            f_obj.write(self.concim[icomp].get_file_entry())

        # close the file
        f_obj.close()

    @classmethod
    def load(cls, f, model, ext_unit_dict=None):
        """
        Load an existing package.

        Parameters
        ----------
        f : filename or file handle
            File to load.
        model : model object
            The model object (of type :class:`flopy.modflow.mf.Modflow`) to
            which this package will be added.
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

        f_obj = get_open_file_object(f, "r")

        # determine problem dimensions
        nlay = model.nlay

        # item 0
        line = f_obj.readline().upper()
        while line.startswith("#"):
            line = f_obj.readline().upper()

        t = line.split()
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
            # print(f"{v}={kwargs[v]}")

        # item 1a - options
        if "frahk" in t:
            kwargs["frahk"] = 1
        else:
            kwargs["frahk"] = 0

        if "mobilesat" in t:
            kwargs["mobilesat"] = 1
        else:
            kwargs["mobilesat"] = 0

        if "inputsat" in t:
            kwargs["inputsat"] = 1
        else:
            kwargs["inputsat"] = 0

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
                    f"to read the per-species Langmuir arrays; model.mcomp="
                    f"{model.mcomp}."
                )
            kwargs["aw_adsorbim"] = True
            kwargs["iarea_fnim"] = iarea_fnim
            kwargs["ikawi_fnim"] = ikawi_fnim

        # item 1b
        if kwargs["icbndimflg"] == 0:
            kwargs["icbundim"] = cls._load_prop_arrays(
                f_obj, model, nlay, np.int32, "icbundim", ext_unit_dict
            )

        # item 2
        kwargs["phif"] = cls._load_prop_arrays(
            f_obj, model, nlay, np.float32, "phif", ext_unit_dict
        )

        # item 3
        kwargs["prsityim"] = cls._load_prop_arrays(
            f_obj, model, nlay, np.float32, "prsityim", ext_unit_dict
        )

        # item 4
        if kwargs["iadsorbim"]:
            kwargs["bulkdim"] = cls._load_prop_arrays(
                f_obj, model, nlay, np.float32, "bulkdim", ext_unit_dict
            )

        # item 5 — Fortran: IF(IDPF.NE.0.AND.IDISPIM.NE.0)
        if model.idpf and kwargs.get("idispim", 0):
            kwargs["dlim"] = cls._load_prop_arrays(
                f_obj, model, nlay, np.float32, "dlim", ext_unit_dict
            )

        # item 6
        kwargs["ddtr"] = cls._load_prop_arrays(
            f_obj, model, nlay, np.float32, "ddtr", ext_unit_dict
        )

        # item 7
        if kwargs["inputsat"]:
            kwargs["sim"] = cls._load_prop_arrays(
                f_obj, model, nlay, np.float32, "sim", ext_unit_dict
            )

        # item 8
        if model.iheat == 1:
            kwargs["htcapsim"] = cls._load_prop_arrays(
                f_obj, model, nlay, np.float32, "htcapsim", ext_unit_dict
            )

            kwargs["htcondsim"] = cls._load_prop_arrays(
                f_obj, model, nlay, np.float32, "htcondsim", ext_unit_dict
            )

        # A-W_ADSORBIM RP1: area arrays (after heat, before the species loop).
        if kwargs["aw_adsorbim"]:
            if kwargs["iarea_fnim"] == 1:
                kwargs["awamaxim"] = cls._load_prop_arrays(
                    f_obj, model, nlay, np.float32, "awamaxim", ext_unit_dict
                )
            elif kwargs["iarea_fnim"] == 4:
                kwargs["awarea_x2im"] = cls._load_prop_arrays(
                    f_obj, model, nlay, np.float32, "awarea_x2im", ext_unit_dict
                )
                kwargs["awarea_x1im"] = cls._load_prop_arrays(
                    f_obj, model, nlay, np.float32, "awarea_x1im", ext_unit_dict
                )
                kwargs["awarea_x0im"] = cls._load_prop_arrays(
                    f_obj, model, nlay, np.float32, "awarea_x0im", ext_unit_dict
                )

        # item 9
        mcomp = model.mcomp
        adsorbim = [0] * mcomp
        flichim = [0] * mcomp
        zodrwim = [0] * mcomp
        zodrsim = [0] * mcomp
        fodrwim = [0] * mcomp
        fodrsim = [0] * mcomp
        concim = [0] * mcomp
        alangawim = [0] * mcomp
        blangawim = [0] * mcomp

        for icomp in range(mcomp):
            # A-W_ADSORBIM RP2: Langmuir A/B per species, before ADSORBIM.
            if kwargs["aw_adsorbim"]:
                alangawim[icomp] = cls._load_prop_arrays(
                    f_obj, model, nlay, np.float32, "alangawim", ext_unit_dict
                )
                blangawim[icomp] = cls._load_prop_arrays(
                    f_obj, model, nlay, np.float32, "blangawim", ext_unit_dict
                )

            if kwargs["iadsorbim"]:
                adsorbim[icomp] = cls._load_prop_arrays(
                    f_obj, model, nlay, np.float32, "adsorbim", ext_unit_dict
                )

            if kwargs["iadsorbim"] == 2 or kwargs["iadsorbim"] == 3:
                flichim[icomp] = cls._load_prop_arrays(
                    f_obj, model, nlay, np.float32, "flichim", ext_unit_dict
                )

            if kwargs["izodim"] == 1 or kwargs["izodim"] == 3 or kwargs["izodim"] == 4:
                zodrwim[icomp] = cls._load_prop_arrays(
                    f_obj, model, nlay, np.float32, "zodrwim", ext_unit_dict
                )

            if kwargs["iadsorbim"] and (
                kwargs["izodim"] == 2 or kwargs["izodim"] == 3 or kwargs["izodim"] == 4
            ):
                zodrsim[icomp] = cls._load_prop_arrays(
                    f_obj, model, nlay, np.float32, "zodrsim", ext_unit_dict
                )

            if kwargs["ifodim"] == 1 or kwargs["ifodim"] == 3 or kwargs["ifodim"] == 4:
                fodrwim[icomp] = cls._load_prop_arrays(
                    f_obj, model, nlay, np.float32, "fodrwim", ext_unit_dict
                )

            if kwargs["iadsorbim"] and (
                kwargs["ifodim"] == 2 or kwargs["ifodim"] == 3 or kwargs["ifodim"] == 4
            ):
                fodrsim[icomp] = cls._load_prop_arrays(
                    f_obj, model, nlay, np.float32, "fodrsim", ext_unit_dict
                )

            concim[icomp] = cls._load_prop_arrays(
                f_obj, model, nlay, np.float32, "concim", ext_unit_dict
            )

        kwargs["adsorbim"] = adsorbim
        kwargs["flichim"] = flichim
        kwargs["zodrwim"] = zodrwim
        kwargs["zodrsim"] = zodrsim
        kwargs["fodrwim"] = fodrwim
        kwargs["fodrsim"] = fodrsim
        kwargs["concim"] = concim
        if kwargs["aw_adsorbim"]:
            kwargs["alangawim"] = alangawim
            kwargs["blangawim"] = blangawim

        f_obj.close()
        # set package unit number
        unitnumber, filenames = get_unitnumber_from_ext_unit_dict(
            model, cls, ext_unit_dict, kwargs["ipakcb"]
        )

        return cls(model, unitnumber=unitnumber, filenames=filenames, **kwargs)

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
