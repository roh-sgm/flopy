"""
Mfusgdpf module.

Contains the MfUsgDpf class. Note that the user can
access the MfUsgDpf class as `flopy.mfusg.MfUsgDpf`.
"""

import os
import re
import stat
from warnings import warn

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
    Items PHIF to BPIM are read layer by layer (GWF2DPFU1AR, DO 200 loop):
    for each layer PHIF, DDFTR, SC1IM, then SC2IM when the effective LAYCON
    is not 0, then ALPHAIM, BETAIM, SRIM, BROOKIM (and BPIM with BUBBLEPT)
    when the effective LAYCON is 5 and TABRICH is off. IUZONTABIM, IBOUNDIM
    and HNEWIM precede that loop, each one complete.

    Multi-layer DPF files written by earlier versions of this class stored
    these arrays property by property, which USG-Transport does not read
    correctly. Such files are not detected or reinterpreted on load; they
    must be regenerated from the original data.

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
        self.idpfhd = idpfhd
        self.idpfdd = idpfdd
        self.frahk = frahk

        model.idpf = 1

        nrow, ncol, nlay, nper = self.parent.nrow_ncol_nlay_nper
        shape = (nlay, nrow, ncol)
        self._laycon = self._get_laycon(model)
        self._tabrich = self._has_tabrich(model)
        self._bubblept = self._has_bubblept(model)
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

        # HNEWIM is read with U1DREL8 (REAL*8) when BAS DPIN/DPIO is active
        # (GWF2DPFU1AR); the precision is fixed here and checked on write
        hnewim_dtype = self._hnewim_dtype(model)
        if isinstance(hnewim, Util3d):
            layers = list(hnewim.util_2ds)
        elif isinstance(hnewim, (list, tuple)):
            layers = [u for u in hnewim if isinstance(u, Util2d)]
        else:
            layers = []
        if layers and model.get_package("BAS6") is None:
            # no BAS6 yet (e.g. MfUsgDpf.load with ext_unit_dict): keep the
            # precision of the loaded layers; write_file checks it against BAS
            hnewim_dtype = layers[0].dtype
        if any(u2d.dtype != hnewim_dtype for u2d in layers):
            # converting existing Util2d layers would either round them or
            # present already rounded values as double precision
            raise ValueError(
                f"hnewim layers must be {np.dtype(hnewim_dtype).name} Util2d "
                "objects for the current BAS precision options"
            )
        self.hnewim = Util3d(model, shape, hnewim_dtype, hnewim, name="hnewim")

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
            if self._bubblept:
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

        Notes
        -----
        HNEWIM layers with ``format.binary = True`` are written, in layer
        order, to a binary file of this package (REAL*8 with BAS DPIN/DPIO,
        REAL otherwise) that U1DREL/U1DREL8 read through a negative unit;
        the file is registered as DATA(BINARY) in the model. write_file()
        does not rewrite the name file (write_input() does). Files and
        name-file entries of other components are never replaced: when the
        default name is taken, a free ``<model>.hnewim_<i>.bin`` is used.

        """
        # checks and the binary destination are resolved before any file is
        # opened, truncated or written
        binary_layers = self._check_hnewim()
        plan = self._plan_hnewim_binary(binary_layers)
        records = self._write_hnewim_binary(binary_layers, plan)

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
        for k, u2d in enumerate(self.hnewim.util_2ds):
            f_obj.write(records[k] if k in records else u2d.get_file_entry())
        # material properties are interleaved by layer (DO 200 K=1,NLAY)
        for layer, laycon in enumerate(self._laycon):
            for name in self._layer_array_names(laycon, self._tabrich, self._bubblept):
                f_obj.write(getattr(self, name)[layer].get_file_entry())

        # close the file only if it was opened here
        if f is None:
            f_obj.close()

    def _check_hnewim(self):
        """Precision and record checks for HNEWIM; writes nothing.

        Every layer is checked, not only the Util3d container: a layer can
        be replaced after construction (e.g. ``dpf.hnewim[k] = Util2d(...)``
        keeps the dtype of the Util2d assigned).
        """
        layer_dtypes = [np.dtype(u2d.dtype) for u2d in self.hnewim.util_2ds]
        if self._hnewim_dtype(self.parent) == np.float64:
            if self.hnewim.dtype != np.float64:
                raise ValueError(
                    "BAS DPIN/DPIO is active but hnewim is stored in single "
                    "precision (DPF was created before BAS or before DPIN/DPIO "
                    "was set); create MfUsgDpf after MfUsgBas with DPIN/DPIO"
                )
            single = [k + 1 for k, dt in enumerate(layer_dtypes) if dt != np.float64]
            if single:
                raise ValueError(
                    "BAS DPIN/DPIO is active but hnewim layer(s) "
                    f"{', '.join(map(str, single))} are stored in single "
                    "precision; assign float64 values (ndarray or float64 Util2d)"
                )
        elif self.hnewim.dtype == np.float64 or np.dtype(np.float64) in layer_dtypes:
            warn(
                "hnewim is stored in double precision but BAS DPIN/DPIO is not "
                "active: the model reads it in single precision"
            )
        binary = [k for k, u2d in enumerate(self.hnewim.util_2ds) if u2d.format.binary]
        if binary and self._hnewim_dtype(self.parent) == np.float32:
            # a REAL record would round float64 values
            double = [k + 1 for k in binary if layer_dtypes[k] == np.float64]
            if double:
                raise ValueError(
                    f"binary hnewim layer(s) {', '.join(map(str, double))} hold "
                    "float64 values but BAS DPIN/DPIO is not active: the REAL "
                    "binary record would round them; assign float32 values or "
                    "activate DPIN/DPIO"
                )
        for k, u2d in enumerate(self.hnewim.util_2ds):
            expected = int(np.prod(get_util2d_shape_for_layer(self.parent, k)))
            if u2d.array.size != expected:
                raise ValueError(
                    f"hnewim layer {k + 1} has {u2d.array.size} values, "
                    f"expected {expected}"
                )
            if k not in binary:
                u2d.check_real_fields()
        return binary

    @staticmethod
    def _path_key(path, case_insensitive):
        """Normalized destination: absolute, symbolic links resolved and, on
        a case-insensitive file system, case folded."""
        key = os.path.realpath(os.path.abspath(path))
        return key.casefold() if case_insensitive else key

    def _plan_hnewim_binary(self, binary_layers):
        """Choose the binary HNEWIM file and unit without writing anything.

        Ownership is tracked by this instance (the identity, device and
        inode, of the file it created, and the name-file entry it created);
        a file or entry is never adopted because its name or path matches.
        An existing destination is rewritten only if it is still that same
        regular file (not a symbolic link, a single link) and no destination
        of another component refers to it; otherwise a free name is used.
        Destinations of packages, outputs, externals (except this instance's
        own entry) and existing files are not used.
        """
        model = self.parent
        own_fname = getattr(self, "_hnewim_bin_fname", None)
        own_unit = getattr(self, "_hnewim_bin_unit", None)
        own_entry = None
        if own_unit in model.external_units:
            idx = model.external_units.index(own_unit)
            if model.external_fnames[idx] == own_fname:
                own_entry = idx
        if not binary_layers:
            return {"stale": own_entry}
        ws = model.model_ws
        folder = os.path.join(ws, model.external_path or "")
        if not os.path.isdir(folder):
            raise FileNotFoundError(f"folder for binary hnewim not found: {folder}")
        swapped = os.path.realpath(ws).swapcase()
        case_insensitive = swapped != os.path.realpath(ws) and os.path.exists(swapped)

        def key(fname):
            fname = str(fname)
            return self._path_key(
                fname if os.path.isabs(fname) else os.path.join(ws, fname),
                case_insensitive,
            )

        destinations = [model.namefile]
        for pak in model.packagelist:
            destinations += list(pak.file_name)
            if getattr(pak, "fn_path", None):
                destinations.append(pak.fn_path)
        destinations += list(model.output_fnames)
        destinations += [
            fn for i, fn in enumerate(model.external_fnames) if i != own_entry
        ]
        taken = {key(fn) for fn in destinations}
        # files of other destinations, by identity (hard links, aliases)
        taken_ids = set()
        for fn in destinations:
            fn = str(fn)
            try:
                st = os.stat(fn if os.path.isabs(fn) else os.path.join(ws, fn))
            except OSError:
                continue
            taken_ids.add((st.st_dev, st.st_ino))
        own_id = getattr(self, "_hnewim_bin_id", None)
        for i in range(1000):
            name = f"{model.name}.hnewim{'' if i == 0 else f'_{i}'}.bin"
            fname = (
                os.path.join(model.external_path, name) if model.external_path else name
            )
            path = os.path.join(ws, fname)
            if key(fname) in taken:
                continue
            try:
                st = os.lstat(path)
            except FileNotFoundError:
                overwrite = False  # free: created exclusively
                break
            if self._is_own_file(st, own_id) and own_id not in taken_ids:
                overwrite = True  # still the file this instance created
                break
            # a file of another component (registered or not), a symbolic
            # link, or an own file now shared through another link
            continue
        else:
            raise FileExistsError("no free name found for the binary hnewim file")
        used = set(model.package_units) | set(model.output_units)
        others = {u for i, u in enumerate(model.external_units) if i != own_entry}
        if own_entry is not None and own_unit not in used and own_unit not in others:
            unit = own_unit  # repeated writes keep the same unit
        else:
            unit = None
        return {
            "fname": fname,
            "path": path,
            "unit": unit,
            "overwrite": overwrite,
            "stale": own_entry,
            "used": used | others,
            "own_id": own_id,
        }

    @staticmethod
    def _is_own_file(st, own_id):
        """True for the regular, singly linked file with this identity."""
        return (
            own_id is not None
            and stat.S_ISREG(st.st_mode)
            and (st.st_dev, st.st_ino) == own_id
            and st.st_nlink == 1
        )

    def _write_hnewim_binary(self, binary_layers, plan):
        """Write binary HNEWIM layers (U1DREL/U1DREL8 stream layout) and
        update this instance's DATA(BINARY) entry; returns the control
        records {layer: record}.

        Each layer: KSTP, KPER (int32), PERTIM, TOTIM (REAL or REAL*8),
        TEXT (16 characters), NSTRT, NNDLAY, ILAY (int32), then the layer
        values. The binary branch returns before applying CNSTNT, so the
        values written are the effective ones (multiplier applied once) and
        the record multiplier is 1.0. A previous file is never deleted.
        """
        model = self.parent
        if not binary_layers:
            if plan["stale"] is not None:
                model.remove_external(unit=self._hnewim_bin_unit)
            self._hnewim_bin_unit = self._hnewim_bin_fname = None
            return {}
        unit = plan["unit"]
        if unit is None:
            unit = model.next_ext_unit()
            while unit in plan["used"]:
                unit = model.next_ext_unit()
        real = np.float64 if self._hnewim_dtype(model) == np.float64 else np.float32
        header = np.dtype(
            [
                ("kstp", "i4"),
                ("kper", "i4"),
                ("pertim", real),
                ("totim", real),
                ("text", "S16"),
                ("nstrt", "i4"),
                ("nndlay", "i4"),
                ("ilay", "i4"),
            ]
        )
        starts = np.cumsum([0] + [u2d.array.size for u2d in self.hnewim.util_2ds])
        if plan["overwrite"]:
            # reopen without following a final symbolic link and check the
            # identity again before truncating
            fd = os.open(plan["path"], os.O_WRONLY | getattr(os, "O_NOFOLLOW", 0))
            if not self._is_own_file(os.fstat(fd), plan["own_id"]):
                os.close(fd)
                raise FileExistsError(
                    f"{plan['path']} changed since it was written; not truncated"
                )
            os.ftruncate(fd, 0)
            fbin = os.fdopen(fd, "wb")
        else:
            # a new destination is created exclusively: never over another file
            fbin = open(plan["path"], "xb")
        with fbin:
            for k in binary_layers:
                rec = np.zeros(1, dtype=header)
                rec["kstp"], rec["kper"], rec["ilay"] = 1, 1, k + 1
                rec["text"] = f"{'HNEWIM':>16s}".encode()
                rec["nstrt"], rec["nndlay"] = starts[k] + 1, starts[k + 1]
                rec.tofile(fbin)
                self.hnewim.util_2ds[k].array.astype(real).ravel().tofile(fbin)
            st = os.fstat(fbin.fileno())
        if plan["stale"] is not None:
            model.remove_external(unit=self._hnewim_bin_unit)
        model.add_external(plan["fname"], unit, binflag=True)
        self._hnewim_bin_fname, self._hnewim_bin_unit = plan["fname"], unit
        self._hnewim_bin_id = (st.st_dev, st.st_ino)  # identity, not the path
        records = {}
        for k in binary_layers:
            u2d = self.hnewim.util_2ds[k]
            if u2d.format.array_free_format:
                records[k] = f"EXTERNAL {-unit:>10d} 1.0 (BINARY) -1 #{u2d.name}\n"
            else:
                records[k] = f"{-unit:10d}{'1.0':>10s}{'':>20s}{-1:10d} #{u2d.name}\n"
        return records

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

        openfile = not hasattr(f, "read")
        f_obj = get_open_file_object(f, "r")
        try:
            kwargs = cls._load_records(f_obj, model, nlay, ext_unit_dict)
        except BaseException:
            if openfile:
                f_obj.close()
            raise
        f_obj.close()
        # set package unit number
        unitnumber, filenames = get_unitnumber_from_ext_unit_dict(
            model, cls, ext_unit_dict, kwargs["ipakcb"]
        )

        return cls(model, unitnumber=unitnumber, filenames=filenames, **kwargs)

    @classmethod
    def _load_records(cls, f_obj, model, nlay, ext_unit_dict):
        """Read items 1-6 (GWF2DPFU1AR order); returns constructor kwargs."""
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
        bubblept = cls._has_bubblept(model)
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
            f_obj,
            model,
            nlay,
            cls._hnewim_dtype(model, ext_unit_dict),
            "hnewim",
            ext_unit_dict,
            binary_multiplier=False,
            binary_negative_unit_only=True,
        )

        # material properties are interleaved by layer (DO 200 K=1,NLAY)
        names = ["phif", "ddftr", "sc1im", "sc2im"]
        if np.any(richards_layers) and not tabrich:
            names += ["alphaim", "betaim", "srim", "brookim"]
            if bubblept:
                names.append("bpim")
        layer_arrays = {name: [None] * nlay for name in names}
        for layer in range(nlay):
            util2d_shape = get_util2d_shape_for_layer(model, layer=layer)
            for name in cls._layer_array_names(laycon[layer], tabrich, bubblept):
                if model.verbose:
                    print(f"   loading {name} layer {layer + 1} ...")
                layer_arrays[name][layer] = Util2d.load(
                    f_obj, model, util2d_shape, np.float32, name, ext_unit_dict
                )
            for name in names:
                if layer_arrays[name][layer] is None:
                    layer_arrays[name][layer] = Util2d(
                        model,
                        util2d_shape,
                        np.float32,
                        0,
                        name=name,
                        locat=cls._defaultunit(),
                    )
        kwargs.update(layer_arrays)
        return kwargs

    @staticmethod
    def _load_prop_arrays(
        f_obj,
        model,
        nlay,
        dtype,
        name,
        ext_unit_dict,
        binary_multiplier=True,
        binary_negative_unit_only=False,
    ):
        if model.verbose:
            print(f"   loading {name} ...")
        prop_array = [0] * nlay
        for layer in range(nlay):
            util2d_shape = get_util2d_shape_for_layer(model, layer=layer)
            prop_array[layer] = Util2d.load(
                f_obj,
                model,
                util2d_shape,
                dtype,
                name,
                ext_unit_dict,
                binary_multiplier=binary_multiplier,
                binary_negative_unit_only=binary_negative_unit_only,
            )
        return prop_array

    @staticmethod
    def _hnewim_dtype(model, ext_unit_dict=None):
        """np.float64 when BAS DPIN or DPIO is active (IDPIN=1), else np.float32.

        Uses the BAS6 package when it is loaded; while loading, the BAS6
        options line is read from ext_unit_dict when BAS6 comes later in
        the name file.
        """
        bas = model.get_package("BAS6")
        if bas is not None:
            double = bas.double_prec or bas.double_io
        else:
            opts = []
            for item in (ext_unit_dict or {}).values():
                if item.filetype.upper() == "BAS6" and item.filehandle is not None:
                    fh = item.filehandle
                    start = fh.tell()
                    try:
                        line = fh.readline()
                        while line.startswith("#"):
                            line = fh.readline()
                    finally:
                        fh.seek(start)  # leave the BAS6 handle where it was
                    opts = re.sub(r"[^A-Z0-9\.\-\+_]", " ", line.upper()).split()
                    break
            double = "DPIN" in opts or "DPIO" in opts
        return np.float64 if double else np.float32

    @staticmethod
    def _layer_array_names(laycon, tabrich, bubblept):
        """Arrays read for one layer inside GWF2DPFU1AR's DO 200 loop."""
        names = ["phif", "ddftr", "sc1im"]
        if laycon != 0:
            names.append("sc2im")
        if laycon == 5 and not tabrich:
            names += ["alphaim", "betaim", "srim", "brookim"]
            if bubblept:
                names.append("bpim")
        return names

    @staticmethod
    def _get_flow_package(model):
        flow = model.get_package("BCF6") or model.get_package("BCF")
        if flow is None:
            flow = model.get_package("LPF")
        return flow

    @staticmethod
    def _get_laycon(model):
        """Effective LAYCON, as derived by READBCF / READLPF (gwf2bcf-lpf-u1.f).

        BCF: FloPy already stores LAYCON separately from INTERCELLT.
        LPF: LAYTYP 0 -> 0, 1..3 -> 3, 4 -> 4, 5 -> 5, <0 -> 0 with
        THICKSTRT and 3 without it. Then, with the BAS RICHARDS (or
        RICHARDS_HP) option, every LAYCON != 0 becomes 5.
        """
        nlay = model.nlay
        flow = MfUsgDpf._get_flow_package(model)
        if flow is not None and hasattr(flow, "laycon"):
            laycon = getattr(flow.laycon, "array", flow.laycon)
            laycon = np.array(laycon, dtype=np.int32)
        elif flow is not None and hasattr(flow, "laytyp"):
            laytyp = getattr(flow.laytyp, "array", flow.laytyp)
            laytyp = np.asarray(laytyp, dtype=np.int32)
            options = str(getattr(flow, "options", "")).upper().split()
            laycon = np.where(laytyp == 0, 0, 3).astype(np.int32)
            laycon[laytyp == 4] = 4
            laycon[laytyp == 5] = 5
            if "THICKSTRT" in options:
                laycon[laytyp < 0] = 0
        else:
            return np.ones(nlay, dtype=np.int32)

        bas = model.get_package("BAS6")
        if bas is not None and getattr(bas, "richards", False):
            laycon[laycon != 0] = 5
        return laycon

    @staticmethod
    def _has_tabrich(model):
        flow = MfUsgDpf._get_flow_package(model)
        return flow is not None and bool(getattr(flow, "tabrich", False))

    @staticmethod
    def _has_bubblept(model):
        flow = MfUsgDpf._get_flow_package(model)
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
