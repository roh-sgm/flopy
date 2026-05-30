"""
Mfusgdpf module.

Contains the MfUsgDpf class. Note that the user can
access the MfUsgDpf class as `flopy.mfusg.MfUsgDpf`.
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


class MfUsgDpf(Package):
    """Dual Porosity Flow (dpf) Package Class for MODFLOW-USG Transport.

    Parameters
    ----------
    model : model object
        The model object (of type :class:`flopy.modflow.Modflow`) to which
        this package will be added.
    ipakcb : int (0,1,-1), (default is 0)
        a flag and a unit number >0 for cell-by-cell mass flux terms.
    idpfhd : int, (default is 0)
        a flag and a unit number >0 for immobile domain heads.
    idpfdd : int, (default is 0)
        a flag and a unit number >0 for immobile domain drawdown.
    iuzontabim : int or array of ints, (default is 0)
        soil type index for all groundwater flow nodes in the immobile domain
    iboundim : int, (default is 0)
        boundary variable for the immobile domain (<0 constant head, =0 no flow)
    hnewim : float or array of floats
        initial (starting) head in the immobile domain
    phif : float or array of floats
        porosity of the mobile domain.
    ddftr : float or array of floats
        dual domain flow transfer rate
    sc1im : float or array of floats
        specific storage of the immobile domain
    sc2im : float or array of floats
        specific yield or porosity of the immobile domain
    alphaim : float or array of floats
        van Genuchten alpha coefficient of the immobile domain
    betaim : float or array of floats
        van Genuchten beta coefficient of the immobile domain
    srim : float or array of floats
        van Genuchten sr coefficient of the immobile domain
    brookim : float or array of floats
        Brooks-Corey exponent for the relative permeability of the immobile domain
    bpim : float or array of floats
        Bubble point or air entry pressure head of the immobile domain
    extension : string,  (default is 'dpf').
    unitnumber : int, default is 57.
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
    >>> dpf = flopy.mfusg.MfUsgdpf(ml)"""

    def __init__(
        self,
        model,
        ipakcb=0,
        idpfhd=0,
        idpfdd=0,
        frahk=False,
        iuzontabim=0,
        iboundim=0,
        hnewim=0.0,
        phif=0.0,
        ddftr=0.0,
        sc1im=0.0,
        sc2im=0.0,
        alphaim=1.0,
        betaim=7.0,
        srim=0.05,
        brookim=6.0,
        bpim=0.0,
        extension="dpf",
        unitnumber=None,
        filenames=None,
        add_package=True,
    ):
        """Constructs the MfUsgdpf object."""
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
        self.idpfhd = idpfhd
        self.idpfdd = idpfdd
        self.frahk = frahk

        model.idpf = 1

        nrow, ncol, nlay, nper = self.parent.nrow_ncol_nlay_nper
        shape = (nlay, nrow, ncol)
        self._laycon = self._get_laycon(model)
        self._tabrich = self._has_tabrich(model)
        self._richards_layers = np.asarray(self._laycon) == 5

        if self._tabrich:
            self.iuzontabim = Util2d(
                model,
                (self._get_nodes(model),),
                np.int32,
                iuzontabim,
                name="iuzontabim",
                locat=self.unit_number[0],
            )

        self.iboundim = Util3d(
            model, shape, np.int32, iboundim, name="iboundim"
        )

        self.hnewim = Util3d(
            model, shape, np.float32, hnewim, name="hnewim"
        )

        self.phif = Util3d(model, shape, np.float32, phif, name="phif")

        self.ddftr = Util3d(model, shape, np.float32, ddftr, name="ddftr")

        self.sc1im = Util3d(model, shape, np.float32, sc1im, name="sc1im")

        self.sc2im = Util3d(model, shape, np.float32, sc2im, name="sc2im")

        if np.any(self._richards_layers) and not self._tabrich:
            self.alphaim = Util3d(
                model, shape, np.float32, alphaim, name="alphaim"
            )
            self.betaim = Util3d(model, shape, np.float32, betaim, name="betaim")
            self.srim = Util3d(model, shape, np.float32, srim, name="srim")
            self.brookim = Util3d(model, shape, np.float32, brookim, name="brookim")
            if self._has_bubblept(model):
                self.bpim = Util3d(model, shape, np.float32, bpim, name="bpim")

        if add_package:
            self.parent.add_package(self)

    def write_file(self, f=None):
        """
        Write the dpf package file.

        Parameters
        ----------
        f : open file object.
            Default is None, which will result in MfUsg.fn_path being
            opened for writing.

        """
        # Open file for writing
        if f is None:
            f_obj = open(self.fn_path, "w")
        else:
            f_obj = f

        #        f_obj.write(f"{self.heading}\n")

        # Item 0: IPAKCB, IdpfCON
        f_obj.write(f" {self.ipakcb:9d} {self.idpfhd:9d} {self.idpfdd:9d}")
        if self.frahk:
            f_obj.write(" FRAHK")
        f_obj.write(" \n")

        if self._tabrich:
            f_obj.write(self.iuzontabim.get_file_entry())
        f_obj.write(self.iboundim.get_file_entry())
        f_obj.write(self.hnewim.get_file_entry())
        f_obj.write(self.phif.get_file_entry())
        f_obj.write(self.ddftr.get_file_entry())
        f_obj.write(self.sc1im.get_file_entry())
        for layer in self._iter_convertible_layers():
            f_obj.write(self.sc2im[layer].get_file_entry())
        if np.any(self._richards_layers) and not self._tabrich:
            for layer in self._iter_richards_layers():
                f_obj.write(self.alphaim[layer].get_file_entry())
                f_obj.write(self.betaim[layer].get_file_entry())
                f_obj.write(self.srim[layer].get_file_entry())
                f_obj.write(self.brookim[layer].get_file_entry())
                if hasattr(self, "bpim"):
                    f_obj.write(self.bpim[layer].get_file_entry())

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
        dpf : MfUsgdpf object

        Examples
        --------
        >>> import flopy
        >>> ml = flopy.mfusg.MfUsg()
        >>> dis = flopy.modflow.ModflowDis.load('Test1.dis', ml)
        >>> dpf = flopy.mfusg.MfUsgdpf.load('Test1.BTN', ml)
        """
        msg = (
            "Model object must be of type flopy.mfusg.MfUsg\n"
            f"but received type: {type(model)}."
        )
        assert isinstance(model, MfUsg), msg

        if model.verbose:
            print("loading dpf package file...")

        nlay = model.nlay

        f_obj = get_open_file_object(f, "r")

        # item 0
        line = f_obj.readline().upper()
        while line.startswith("#"):
            line = f_obj.readline().upper()

        t = line.split()
        kwargs = {}

        # item 1a
        vars = {
            "ipakcb": int,
            "idpfhd": int,
            "idpfdd": int,
        }

        for i, (v, c) in enumerate(vars.items()):
            kwargs[v] = c(t[i].strip())
            # print(f"{v}={kwargs[v]}\n")
        kwargs["frahk"] = "FRAHK" in t[3:]

        laycon = cls._get_laycon(model)
        tabrich = cls._has_tabrich(model)
        richards_layers = np.asarray(laycon) == 5

        # item 1b
        if tabrich:
            kwargs["iuzontabim"] = Util2d.load(
                f_obj,
                model,
                (cls._get_nodes(model),),
                np.int32,
                "iuzontabim",
                ext_unit_dict,
            )

        kwargs["iboundim"] = cls._load_prop_arrays(
            f_obj, model, nlay, np.int32, "iboundim", ext_unit_dict
        )

        kwargs["hnewim"] = cls._load_prop_arrays(
            f_obj, model, nlay, np.float32, "hnewim", ext_unit_dict
        )

        kwargs["phif"] = cls._load_prop_arrays(
            f_obj, model, nlay, np.float32, "phif", ext_unit_dict
        )

        kwargs["ddftr"] = cls._load_prop_arrays(
            f_obj, model, nlay, np.float32, "ddftr", ext_unit_dict
        )

        kwargs["sc1im"] = cls._load_prop_arrays(
            f_obj, model, nlay, np.float32, "sc1im", ext_unit_dict
        )

        kwargs["sc2im"] = cls._load_prop_arrays(
            f_obj,
            model,
            nlay,
            np.float32,
            "sc2im",
            ext_unit_dict,
            active_layers=np.asarray(laycon) != 0,
        )

        if np.any(richards_layers) and not tabrich:
            kwargs["alphaim"] = cls._load_prop_arrays(
                f_obj,
                model,
                nlay,
                np.float32,
                "alphaim",
                ext_unit_dict,
                active_layers=richards_layers,
            )
            kwargs["betaim"] = cls._load_prop_arrays(
                f_obj,
                model,
                nlay,
                np.float32,
                "betaim",
                ext_unit_dict,
                active_layers=richards_layers,
            )
            kwargs["srim"] = cls._load_prop_arrays(
                f_obj,
                model,
                nlay,
                np.float32,
                "srim",
                ext_unit_dict,
                active_layers=richards_layers,
            )
            kwargs["brookim"] = cls._load_prop_arrays(
                f_obj,
                model,
                nlay,
                np.float32,
                "brookim",
                ext_unit_dict,
                active_layers=richards_layers,
            )
            if cls._has_bubblept(model):
                kwargs["bpim"] = cls._load_prop_arrays(
                    f_obj,
                    model,
                    nlay,
                    np.float32,
                    "bpim",
                    ext_unit_dict,
                    active_layers=richards_layers,
                )

        f_obj.close()
        # set package unit number
        unitnumber, filenames = get_unitnumber_from_ext_unit_dict(
            model, cls, ext_unit_dict, kwargs["ipakcb"]
        )

        return cls(model, unitnumber=unitnumber, filenames=filenames, **kwargs)

    @staticmethod
    def _load_prop_arrays(
        f_obj, model, nlay, dtype, name, ext_unit_dict, active_layers=None
    ):
        if model.verbose:
            print(f"   loading {name} ...")
        if active_layers is None:
            active_layers = np.ones(nlay, dtype=bool)
        prop_array = [0] * nlay
        for layer in range(nlay):
            util2d_shape = get_util2d_shape_for_layer(model, layer=layer)
            if active_layers[layer]:
                prop_array[layer] = Util2d.load(
                    f_obj, model, util2d_shape, dtype, name, ext_unit_dict
                )
            else:
                prop_array[layer] = Util2d(
                    model,
                    util2d_shape,
                    dtype,
                    0,
                    name=name,
                    locat=MfUsgDpf._defaultunit(),
                )
        return prop_array

    def _iter_convertible_layers(self):
        for layer, laycon in enumerate(self._laycon):
            if laycon != 0:
                yield layer

    def _iter_richards_layers(self):
        for layer, is_richards in enumerate(self._richards_layers):
            if is_richards:
                yield layer

    @staticmethod
    def _get_laycon(model):
        nlay = model.nlay
        flow = model.get_package("BCF6") or model.get_package("BCF")
        if flow is not None and hasattr(flow, "laycon"):
            laycon = getattr(flow.laycon, "array", flow.laycon)
            return np.asarray(laycon, dtype=np.int32)

        flow = model.get_package("LPF")
        if flow is not None and hasattr(flow, "laytyp"):
            laytyp = getattr(flow.laytyp, "array", flow.laytyp)
            return np.asarray(laytyp, dtype=np.int32)

        return np.ones(nlay, dtype=np.int32)

    @staticmethod
    def _has_tabrich(model):
        flow = model.get_package("BCF6") or model.get_package("BCF")
        return flow is not None and bool(getattr(flow, "tabrich", False))

    @staticmethod
    def _has_bubblept(model):
        flow = model.get_package("BCF6") or model.get_package("BCF")
        if flow is None:
            flow = model.get_package("LPF")
        return flow is not None and bool(getattr(flow, "bubblept", False))

    @staticmethod
    def _get_nodes(model):
        disu = model.get_package("DISU")
        if disu is not None:
            return int(disu.nodes)
        nrow, ncol, nlay, _ = model.nrow_ncol_nlay_nper
        if nrow is None:
            return int(np.asarray(ncol).sum())
        return int(nlay * nrow * ncol)

    @staticmethod
    def _ftype():
        return "DPF"

    @staticmethod
    def _defaultunit():
        return 157
