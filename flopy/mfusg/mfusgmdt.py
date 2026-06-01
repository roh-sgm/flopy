"""
Mfusgmdt module.

Contains the MfUsgMdt class. Note that the user can
access the MfUsgMdt class as `flopy.mfusg.MfUsgMdt`.
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

# USG-T (gwt2mdtu1.for) reads/uses the AIOLD1MD/AIOLD2MD arrays and treats the
# time shift as active only when TSHIFTMD exceeds this threshold. Use the same
# criterion everywhere (constructor, writer, options check) so written files
# stay aligned with what the solver will read.
MDT_TSHIFT_THRESHOLD = 1.0e-10


class MfUsgMdt(Package):
    """Matrix Diffusion Transport (mdt) Package Class for MODFLOW-USG Transport.

    Parameters
    ----------
    model : model object
        The model object (of type :class:`flopy.modflow.Modflow`) to which
        this package will be added.
    ipakcb : int (0,1,-1), (default is 0)
        a flag and a unit number >0 for cell-by-cell mass flux terms.
    imdtcf : int (0,1,2,3,4,5,6,7), (default is 0)
        a flag and a unit number >0 for immobile domain concentrations
    mdflag : int, (default is 0)
        a flag for the matrix diffusion type
    volfracmd : float, (default is 0.0)
        volume fraction of fracture/high permeability domain.
        SAME AS PHIF IN DPT PACKAGE
    pormd : float, (default is 0.0)
        effective transport porosity of the mobile domain
    rhobmd : float, (default is 0.0)
        bulk density of the mobile domain
    difflenmd : float, (default is 0.0)
        diffusion length for diffusive transport within the matrix domain
    tortmd : float, (default is 0.0)
        tortuosity factor for the matrix domain
    kdmd : float, (default is 0.0)
        adsorption coefficient (Kd) of a contaminant species in immobile domain
    decaymd : float, (default is 0.0)
        first-order decay rate of a contaminant species in the immobile domain
    yieldmd : float, (default is 0.0)
        yield coefficient for chain decay
    diffmd : float, (default is 0.0)
        diffusion coefficient for the immobile domain
    aiold1md : float, (default is 0.0) (MDFLAG = 1, 3, 4, 5, 6, or 7)
        concentration integral (equation 8) of each species
    aiold2md : float, (default is 0.0) (MDFLAG = 2, 5, 6, or 7)
        concentration integral (equation 8) of each species
    frahk : bool, (default is False)
        True = hydraulic conductivity and storage terms only for mobile domain
        False = hydraulic conductivity and storage terms for matrix domains
    fradarcy : bool, (default is False)
        True = Darcy flux computed only for fracture domain
    tshiftmd : float, (default is 0.0)
        time shift for the immobile domain
    iunitAI2 :float, (default is 0.0)
        separate binary file of ai2 for the matrix domain
    crootname_md : str, (default is None)
        rootname for the concentration output binary file
    extension : str, (default is 'mdt').
    unitnumber : int, (default is 57).
        File unit number and the output files.
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
    >>> ml = flopy.mfusg.MfUsg(exe_name='USGs_1.exe')
    >>> disu = flopy.mfusg.MfUsgDisU(model=ml, nlay=1, nodes=1,
                 iac=[1], njag=1,ja=np.array([0]), fahl=[1.0], cl12=[1.0])
    >>> mdt = flopy.mfusg.MfUsgmdt(ml)"""

    def __init__(
        self,
        model,
        ipakcb=0,
        imdtcf=0,
        mdflag=0,
        volfracmd=0.0,
        pormd=0.0,
        rhobmd=0.0,
        difflenmd=0.0,
        tortmd=0.0,
        kdmd=0.0,
        decaymd=0.0,
        yieldmd=0.0,
        diffmd=0.0,
        aiold1md=0.0,
        aiold2md=0.0,
        frahk=False,
        fradarcy=False,
        tshiftmd=0.0,
        iunitAI2=0,
        crootname=None,
        extension="mdt",
        unitnumber=None,
        filenames=None,
        add_package=True,
    ):
        """Constructs the MfUsgmdt object."""
        msg = (
            "Model object must be of type flopy.mfusg.MfUsg\n"
            f"but received type: {type(model)}."
        )
        assert isinstance(model, MfUsg), msg

        if unitnumber is None:
            unitnumber = self._defaultunit()

        # call base package constructor
        super().__init__(
            model,
            extension=extension,
            name=self._ftype(),
            unit_number=unitnumber,
            filenames=self._prepare_filenames(filenames),
        )

        self._generate_heading()
        self.ipakcb = ipakcb
        self.imdtcf = imdtcf

        # options
        self.frahk = frahk
        self.fradarcy = fradarcy
        self.tshiftmd = tshiftmd
        self.iunitAI2 = iunitAI2
        self.crootname = crootname

        # FRAHK (IFRAHK=1) and FRADARCY (IFRAHK=2) map to the same flag with
        # conflicting values -- only one may be set.
        if frahk and fradarcy:
            raise ValueError(
                "MDT FRAHK and FRADARCY are mutually exclusive (both set IFRAHK)."
            )
        # gwt2mdtu1.for reads the header options only when IDPF==0 (single-domain
        # flow); with dual-porosity flow they are skipped and VOLFRACMD comes
        # from DPF/PHIF, so setting them with IDPF!=0 is invalid.
        if model.idpf and (
            frahk
            or fradarcy
            or tshiftmd > MDT_TSHIFT_THRESHOLD
            or iunitAI2 > 0
            or crootname is not None
        ):
            raise ValueError(
                "MDT header options (FRAHK/FRADARCY/TSHIFTMD/SEPARATE_AI2/"
                "MULTIFILE_MD) are only read by USG-T when IDPF==0; this model "
                "has dual-porosity flow (IDPF!=0)."
            )
        # MULTIFILE_MD output needs a positive IMDTCF unit (Fortran STOPs at 0).
        if crootname is not None and imdtcf <= 0:
            raise ValueError(
                "MDT MULTIFILE_MD (crootname) requires imdtcf>0 for the AI1/AI2 "
                "output files."
            )

        nrow, ncol, nlay, nper = self.parent.nrow_ncol_nlay_nper

        self.mdflag = Util3d(model, (nlay, nrow, ncol), np.int32, mdflag, name="mdflag")
        self.volfracmd = Util3d(
            model, (nlay, nrow, ncol), np.float32, volfracmd, name="volfracmd"
        )
        self.pormd = Util3d(model, (nlay, nrow, ncol), np.float32, pormd, name="pormd")
        self.rhobmd = Util3d(
            model, (nlay, nrow, ncol), np.float32, rhobmd, name="rhobmd"
        )
        self.difflenmd = Util3d(
            model, (nlay, nrow, ncol), np.float32, difflenmd, name="difflenmd"
        )
        self.tortmd = Util3d(
            model, (nlay, nrow, ncol), np.float32, tortmd, name="tortmd"
        )

        mcomp = model.mcomp

        def _per_comp(val, label, required):
            """Broadcast a scalar to MCOMP, or validate a per-species sequence."""
            if np.isscalar(val):
                return [val] * mcomp
            seq = list(val)
            if required and len(seq) != mcomp:
                raise ValueError(
                    f"MDT {label} must have MCOMP={mcomp} entries; got {len(seq)}."
                )
            return seq

        kdmd = _per_comp(kdmd, "kdmd", True)
        decaymd = _per_comp(decaymd, "decaymd", True)
        yieldmd = _per_comp(yieldmd, "yieldmd", True)
        diffmd = _per_comp(diffmd, "diffmd", True)
        aiold1md = _per_comp(aiold1md, "aiold1md", self.tshiftmd > MDT_TSHIFT_THRESHOLD)
        aiold2md = _per_comp(aiold2md, "aiold2md", self.tshiftmd > MDT_TSHIFT_THRESHOLD)

        self.kdmd = [0] * mcomp
        self.decaymd = [0] * mcomp
        self.yieldmd = [0] * mcomp
        self.diffmd = [0] * mcomp
        self.aiold1md = [0] * mcomp
        self.aiold2md = [0] * mcomp

        for icomp in range(mcomp):
            self.kdmd[icomp] = Util3d(
                model, (nlay, nrow, ncol), np.float32, kdmd[icomp], name="kdmd"
            )
            self.decaymd[icomp] = Util3d(
                model, (nlay, nrow, ncol), np.float32, decaymd[icomp], name="decaymd"
            )
            self.yieldmd[icomp] = Util3d(
                model, (nlay, nrow, ncol), np.float32, yieldmd[icomp], name="yieldmd"
            )
            self.diffmd[icomp] = Util3d(
                model, (nlay, nrow, ncol), np.float32, diffmd[icomp], name="diffmd"
            )
            if self.tshiftmd > MDT_TSHIFT_THRESHOLD:
                self.aiold1md[icomp] = Util3d(
                    model,
                    (nlay, nrow, ncol),
                    np.float32,
                    aiold1md[icomp],
                    name="aiold1md",
                )
                self.aiold2md[icomp] = Util3d(
                    model,
                    (nlay, nrow, ncol),
                    np.float32,
                    aiold2md[icomp],
                    name="aiold2md",
                )

        if add_package:
            self.parent.add_package(self)

    def write_file(self, f=None):
        """
        Write the mdt package file.

        Parameters
        ----------
        f : open file object.
            Default is None, which will result in MfUsg.fn_path being
            opened for writing.

        Examples
        --------
        """
        # Open file for writing (use a caller-supplied handle if given).
        opened = f is None
        f_obj = open(self.fn_path, "w") if opened else f

        f_obj.write(f"{self.heading}\n")
        f_obj.write(f"{self.ipakcb:9d} {self.imdtcf:9d}")

        # Options
        if self.frahk:
            f_obj.write(" FRAHK")

        if self.fradarcy:
            f_obj.write(" FRADARCY")

        if self.tshiftmd > MDT_TSHIFT_THRESHOLD:
            # General format so small but valid values (e.g. 1e-6) are not
            # rounded to 0.0, which would make USG-T skip the AIOLD arrays.
            f_obj.write(f" TSHIFTMD {self.tshiftmd:15.7g}")

        if self.iunitAI2 > 0:
            f_obj.write(f" SEPARATE_AI2 {self.iunitAI2:9d}")

        if self.crootname is not None:
            f_obj.write(f" MULTIFILE_MD {self.crootname}")

        f_obj.write("\n")

        # Item 1: MDFLAG, VOLFRACMD, PORMD, RHOBMD, DIFFLENMD, TORTMD
        f_obj.write(self.mdflag.get_file_entry())
        if not self.parent.idpf:
            f_obj.write(self.volfracmd.get_file_entry())
        f_obj.write(self.pormd.get_file_entry())
        f_obj.write(self.rhobmd.get_file_entry())
        f_obj.write(self.difflenmd.get_file_entry())
        f_obj.write(self.tortmd.get_file_entry())

        # Item 2 - 7: KDM, DECAYMD, YIELDMD, DIFFMD, AIOLD1MD, AIOLD2MD
        mcomp = self.parent.mcomp

        for icomp in range(mcomp):
            f_obj.write(self.kdmd[icomp].get_file_entry())
            f_obj.write(self.decaymd[icomp].get_file_entry())
            f_obj.write(self.yieldmd[icomp].get_file_entry())
            f_obj.write(self.diffmd[icomp].get_file_entry())
            if self.tshiftmd > MDT_TSHIFT_THRESHOLD:
                f_obj.write(self.aiold1md[icomp].get_file_entry())
                f_obj.write(self.aiold2md[icomp].get_file_entry())

        # close the file only if we opened it
        if opened:
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
        mdt : MfUsgmdt object

        Examples
        --------

        >>> import flopy
        >>> ml = flopy.mfusg.MfUsg()
        >>> dis = flopy.mfusg.MfUsgDisU.load('SeqDegEg.dis', ml)
        >>> mdt = flopy.mfusg.MfUsgmdt.load('SeqDegEg.mdt', ml)
        """
        msg = (
            "Model object must be of type flopy.mfusg.MfUsg\n"
            f"but received type: {type(model)}."
        )
        assert isinstance(model, MfUsg), msg

        if model.verbose:
            print("loading mdt package file...")

        f_obj = get_open_file_object(f, "r")

        nlay = model.nlay

        # item 0 + item 1a header. Keep the original case for the rootname;
        # compare keywords case-insensitively.
        line = f_obj.readline()
        while line.startswith("#"):
            line = f_obj.readline()

        traw = line.split()
        t = [tok.upper() for tok in traw]
        kwargs = {}
        kwargs["ipakcb"] = int(t[0])
        kwargs["imdtcf"] = int(t[1])

        # item 1a options -- USG-T (gwt2mdtu1.for) reads these only when IDPF==0.
        kwargs["frahk"] = 0
        kwargs["fradarcy"] = 0
        kwargs["tshiftmd"] = 0.0
        kwargs["iunitAI2"] = 0
        kwargs["crootname"] = None
        if not model.idpf:
            if "FRAHK" in t:
                kwargs["frahk"] = 1
            if "FRADARCY" in t:
                kwargs["fradarcy"] = 1
            if "TSHIFTMD" in t:
                kwargs["tshiftmd"] = float(t[t.index("TSHIFTMD") + 1])
            if "SEPARATE_AI2" in t:
                kwargs["iunitAI2"] = int(t[t.index("SEPARATE_AI2") + 1])
            if "MULTIFILE_MD" in t:
                kwargs["crootname"] = traw[t.index("MULTIFILE_MD") + 1]

        kwargs["mdflag"] = cls._load_prop_arrays(
            f_obj, model, nlay, np.int32, "mdflag", ext_unit_dict
        )

        if not model.idpf:
            kwargs["volfracmd"] = cls._load_prop_arrays(
                f_obj, model, nlay, np.float32, "volfracmd", ext_unit_dict
            )

        kwargs["pormd"] = cls._load_prop_arrays(
            f_obj, model, nlay, np.float32, "pormd", ext_unit_dict
        )

        kwargs["rhobmd"] = cls._load_prop_arrays(
            f_obj, model, nlay, np.float32, "rhobmd", ext_unit_dict
        )

        kwargs["difflenmd"] = cls._load_prop_arrays(
            f_obj, model, nlay, np.float32, "difflenmd", ext_unit_dict
        )

        kwargs["tortmd"] = cls._load_prop_arrays(
            f_obj, model, nlay, np.float32, "tortmd", ext_unit_dict
        )

        mcomp = model.mcomp

        kdmd = [0] * mcomp
        decaymd = [0] * mcomp
        yieldmd = [0] * mcomp
        diffmd = [0] * mcomp
        aiold1md = [0] * mcomp
        aiold2md = [0] * mcomp

        for icomp in range(mcomp):
            kdmd[icomp] = cls._load_prop_arrays(
                f_obj, model, nlay, np.float32, "kdmd", ext_unit_dict
            )

            decaymd[icomp] = cls._load_prop_arrays(
                f_obj, model, nlay, np.float32, "decaymd", ext_unit_dict
            )

            yieldmd[icomp] = cls._load_prop_arrays(
                f_obj, model, nlay, np.float32, "yieldmd", ext_unit_dict
            )

            diffmd[icomp] = cls._load_prop_arrays(
                f_obj, model, nlay, np.float32, "diffmd", ext_unit_dict
            )

            if kwargs["tshiftmd"] > MDT_TSHIFT_THRESHOLD:
                aiold1md[icomp] = cls._load_prop_arrays(
                    f_obj, model, nlay, np.float32, "aiold1md", ext_unit_dict
                )

                aiold2md[icomp] = cls._load_prop_arrays(
                    f_obj, model, nlay, np.float32, "aiold2md", ext_unit_dict
                )

        kwargs["kdmd"] = kdmd
        kwargs["decaymd"] = decaymd
        kwargs["yieldmd"] = yieldmd
        kwargs["diffmd"] = diffmd
        kwargs["aiold1md"] = aiold1md
        kwargs["aiold2md"] = aiold2md

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
        return "MDT"

    @staticmethod
    def _defaultunit():
        return 156
