"""
mfusgghb module. Contains the MfUsgGhb class.
"""

import numpy as np

from ..modflow.mfghb import ModflowGhb
from ..pakbase import Package
from ..utils import MfList
from ..utils.recarray_utils import create_empty_recarray


class MfUsgGhb(ModflowGhb):
    """MODFLOW-USG General-Head Boundary Package (GHB) with transport support.

    Supports unstructured grids (node-based) and AUX concentration variables.

    Parameters
    ----------
    model : MfUsg
        The model object to which this package will be added.
    stress_period_data : dict, optional
        Dictionary keyed by zero-based stress period index. Each value is a
        recarray with fields matching ``dtype``. For unstructured grids the
        required base fields are ``node``, ``bhead``, and ``cond``.
    dtype : np.dtype, optional
        Custom dtype. If None, defaults to the unstructured base dtype plus
        any AUX fields declared in ``options``.
    ipakcb : int
        Unit number for cell-by-cell budget output (default 0 = not written).
    options : list of str, optional
        Package options (e.g. ``["AUX C01"]``).
    extension : str
        Filename extension (default ``"ghb"``).
    unitnumber : int, optional
        File unit number.
    filenames : str or list of str, optional
        Package filename(s).
    add_package : bool
        Add package to model on construction (default True).

    Examples
    --------
    >>> import flopy
    >>> m = flopy.mfusg.MfUsg()
    >>> ghb = flopy.mfusg.MfUsgGhb.load("MDV.ghb", m, nper=12)
    """

    def __init__(
        self,
        model,
        ipakcb=None,
        stress_period_data=None,
        dtype=None,
        options=None,
        extension="ghb",
        unitnumber=None,
        filenames=None,
        add_package=True,
    ):
        if unitnumber is None:
            unitnumber = ModflowGhb._defaultunit()
        if options is None:
            options = []

        # Prepare two filename slots (package input + CBC output) and register
        # the CBC output file — mirrors ModflowGhb.__init__ exactly.
        filenames = self._prepare_filenames(filenames, 2)
        self.set_cbc_output_file(ipakcb, model, filenames[1])

        Package.__init__(
            self,
            model,
            extension=extension,
            name=self._ftype(),
            unit_number=unitnumber,
            filenames=filenames[0],
        )
        self._generate_heading()
        self.url = "ghb.html"
        self.np = 0
        self.options = options

        if dtype is not None:
            self.dtype = dtype
        else:
            self.dtype = self.get_default_dtype(structured=self.parent.structured)

        self._update_aux_options()
        self.stress_period_data = MfList(self, stress_period_data)

        if add_package:
            self.parent.add_package(self)

    # ------------------------------------------------------------------
    # Static helpers
    # ------------------------------------------------------------------

    @staticmethod
    def get_default_dtype(structured=True):
        if structured:
            return ModflowGhb.get_default_dtype(structured=True)
        return np.dtype([
            ("node", int),
            ("bhead", np.float32),
            ("cond", np.float32),
        ])

    @staticmethod
    def get_empty(ncells=0, aux_names=None, structured=True):
        dtype = MfUsgGhb.get_default_dtype(structured=structured)
        if aux_names is not None:
            dtype = Package.add_to_dtype(dtype, aux_names, np.float32)
        return create_empty_recarray(ncells, dtype, default_value=-1.0e10)

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _update_aux_options(self):
        """Append AUX entries to options for dtype fields beyond the base fields."""
        n_base = len(self.get_default_dtype(structured=self.parent.structured).names)
        for name in self.dtype.names[n_base:]:
            token = f"AUX {name.upper()}"
            if not any(name.lower() in o.lower() for o in self.options):
                self.options.append(token)

    # ------------------------------------------------------------------
    # write_file
    # ------------------------------------------------------------------

    def write_file(self):
        """Write the package file in MODFLOW-USG-T GHB format."""
        if self.parent.structured:
            ModflowGhb.write_file(self)
            return
        nper = self.parent.nper
        n_base = len(self.get_default_dtype(structured=self.parent.structured).names)

        with open(self.fn_path, "w") as f:
            f.write(f"{self.heading}\n")

            # Header: MXACTB IGHBCB [AUX ...]
            line = f" {self.stress_period_data.mxact:9d} {self.ipakcb}"
            for opt in self.options:
                line += f" {opt}"
            f.write(line + "\n")

            for kper in range(nper):
                if kper in self.stress_period_data.data:
                    kdata = self.stress_period_data[kper]
                    f.write(f" {len(kdata)} 0    Stress Period {kper + 1}\n")
                    for rec in kdata:
                        row = (f" {int(rec['node']) + 1}"
                               f"  {float(rec['bhead']):.6f}"
                               f"  {float(rec['cond']):.6e}")
                        for name in self.dtype.names[n_base:]:
                            row += f" {float(rec[name]):.6e}"
                        f.write(row + "\n")
                else:
                    f.write(f" -1 0    Stress Period {kper + 1}\n")

    # ------------------------------------------------------------------
    # load
    # ------------------------------------------------------------------

    @classmethod
    def load(cls, f, model, nper=None, ext_unit_dict=None, check=True):
        """Load a MODFLOW-USG-T GHB package from file.

        Parameters
        ----------
        f : str or file-like
            Path to the .ghb file, or an open file handle.
        model : MfUsg
            Model to attach the package to.
        nper : int, optional
            Number of stress periods. Defaults to ``model.nper``.

        Returns
        -------
        MfUsgGhb
        """
        if model.structured:
            return ModflowGhb.load(f, model, nper=nper, ext_unit_dict=ext_unit_dict, check=check)

        if model.verbose:
            print("loading mfusg ghb package file...")
        if nper is None:
            nper = model.nper

        openfile = not hasattr(f, "read")
        if openfile:
            f = open(f, "r")

        # Skip comment lines, find the header line
        line = f.readline()
        while line.startswith("#"):
            line = f.readline()

        # Parse header: MXACTB IGHBCB [AUX varname ...]
        tokens = line.split()
        ipakcb = int(tokens[1]) if len(tokens) > 1 else 0
        options = []
        aux_names = []
        i = 2
        while i < len(tokens):
            if tokens[i].upper() == "AUX" and i + 1 < len(tokens):
                aux_names.append(tokens[i + 1])
                options.append(f"AUX {tokens[i + 1]}")
                i += 2
            else:
                options.append(tokens[i])
                i += 1

        dtype = cls.get_default_dtype(structured=model.structured)
        if aux_names:
            dtype = Package.add_to_dtype(dtype, aux_names, np.float32)

        spd = {}
        current = None

        for kper in range(nper):
            line = f.readline()
            if not line:
                break
            nact = int(line.split()[0])

            if nact < 0:
                if current is not None:
                    spd[kper] = current.copy()
            else:
                current = create_empty_recarray(nact, dtype, default_value=0.0)
                for idx in range(nact):
                    vals = f.readline().split()
                    for j, name in enumerate(dtype.names):
                        kind = dtype[name].kind
                        current[idx][name] = (
                            int(vals[j]) if kind in ("i", "u") else float(vals[j])
                        )
                    current[idx]["node"] -= 1
                spd[kper] = current

        if openfile:
            f.close()

        unitnumber = None
        filenames = [None]
        if ext_unit_dict is not None:
            unitnumber_ext, fname_ext = model.get_ext_dict_attr(
                ext_unit_dict, filetype=cls._ftype()
            )
            if unitnumber_ext is not None:
                unitnumber = unitnumber_ext
            if fname_ext is not None:
                filenames = [fname_ext]

        return cls(
            model,
            stress_period_data=spd,
            dtype=dtype,
            ipakcb=ipakcb,
            options=options,
            extension="ghb",
            unitnumber=unitnumber,
            filenames=filenames,
        )
