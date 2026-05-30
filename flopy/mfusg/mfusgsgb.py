"""
mfusgsgb module. Contains the MfUsgSgb class.

Specified Gradient Boundary (SGB) Package -- USG-Transport.

Reference: Fortran subroutines ``GLO2SGBU1AR`` / ``GLO2SGBU1RP`` in
``glo2sgbu1.f`` of USG-Transport 2.7. SGB applies a specified hydraulic
gradient at boundary nodes (distinct from GHB, which specifies a
head-to-head conductance). It is a node-based (unstructured) list package::

    Item 1 (header):   MXACTS ISGBCB [AUX <name> ...] [NOPRINT]
    Per stress period: ITMP NP
                       ITMP rows of:  NODE  GRADIENT  [aux ...]

``ITMP < 0`` reuses the previous stress period's list. ``ULSTRDU`` reads a
single node id per row (``RLIST(1)``) followed by the gradient and any AUX
values, so each record is ``NODE GRADIENT [aux ...]``.

Internal ``node`` values are 0-based; the file is written 1-based.
Named parameters (``NPSGB > 0``) are not supported in this version.
"""

import numpy as np

from ..pakbase import Package
from ..utils import MfList
from ..utils.recarray_utils import create_empty_recarray
from .mfusg import MfUsg


class MfUsgSgb(Package):
    """MODFLOW-USG Specified Gradient Boundary (SGB) Package.

    Node-based package for unstructured USG-Transport grids, with optional
    AUX concentration variables for transport.

    Parameters
    ----------
    model : MfUsg
        The model object to which this package will be added.
    ipakcb : int, optional
        Unit number for cell-by-cell budget output (``ISGBCB``; default
        ``None`` = not written).
    stress_period_data : dict, optional
        Dictionary keyed by zero-based stress period index. Each value is a
        recarray with fields matching ``dtype``. The required base fields are
        ``node`` (0-based) and ``gradient``.
    dtype : np.dtype, optional
        Custom dtype. If None, defaults to the base dtype plus any AUX fields
        declared in ``options``.
    options : list of str, optional
        Package options (e.g. ``["AUX C01"]`` or ``["NOPRINT"]``).
    extension : str
        Filename extension (default ``"sgb"``).
    unitnumber : int, optional
        File unit number.
    filenames : str or list of str, optional
        Package filename(s).
    add_package : bool
        Add package to model on construction (default True).

    Examples
    --------
    >>> import flopy
    >>> m = flopy.mfusg.MfUsg(structured=False)
    >>> spd = {0: [[0, 0.01]]}  # node 0 (0-based), gradient 0.01
    >>> sgb = flopy.mfusg.MfUsgSgb(m, stress_period_data=spd)
    """

    def __init__(
        self,
        model,
        ipakcb=None,
        stress_period_data=None,
        dtype=None,
        options=None,
        extension="sgb",
        unitnumber=None,
        filenames=None,
        add_package=True,
    ):
        assert isinstance(model, MfUsg), (
            "Model object must be of type flopy.mfusg.MfUsg\n"
            f"but received type: {type(model)}."
        )

        if unitnumber is None:
            unitnumber = MfUsgSgb._defaultunit()
        if options is None:
            options = []

        # Two filename slots (package input + CBC output); register the CBC
        # output file -- mirrors MfUsgGhb.__init__.
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
        self.url = "sgb.html"
        self.np = 0
        self.options = options

        if dtype is not None:
            self.dtype = dtype
        else:
            self.dtype = self.get_default_dtype()

        self._update_aux_options()
        self.stress_period_data = MfList(self, stress_period_data)

        if add_package:
            self.parent.add_package(self)

    # ------------------------------------------------------------------
    # Static helpers
    # ------------------------------------------------------------------

    @staticmethod
    def get_default_dtype():
        """Node-based SGB dtype: node (0-based) and gradient (float64)."""
        return np.dtype([("node", int), ("gradient", np.float64)])

    @staticmethod
    def get_empty(ncells=0, aux_names=None):
        dtype = MfUsgSgb.get_default_dtype()
        if aux_names is not None:
            dtype = Package.add_to_dtype(dtype, aux_names, np.float64)
        return create_empty_recarray(ncells, dtype, default_value=-1.0e10)

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _update_aux_options(self):
        """Append AUX entries to options for dtype fields beyond the base fields."""
        n_base = len(self.get_default_dtype().names)
        for name in self.dtype.names[n_base:]:
            if not any(name.lower() in o.lower() for o in self.options):
                self.options.append(f"AUX {name.upper()}")

    # ------------------------------------------------------------------
    # write_file
    # ------------------------------------------------------------------

    def write_file(self):
        """Write the package file in MODFLOW-USG-T SGB format."""
        nper = self.parent.nper
        n_base = len(self.get_default_dtype().names)

        with open(self.fn_path, "w") as f:
            f.write(f"{self.heading}\n")

            # Item 1: MXACTS ISGBCB [AUX ...] [NOPRINT]
            line = f" {self.stress_period_data.mxact:9d} {self.ipakcb}"
            for opt in self.options:
                line += f" {opt}"
            f.write(line + "\n")

            for kper in range(nper):
                if kper in self.stress_period_data.data:
                    kdata = self.stress_period_data[kper]
                    f.write(f" {len(kdata)} 0    Stress Period {kper + 1}\n")
                    for rec in kdata:
                        row = (
                            f" {int(rec['node']) + 1}"
                            f"  {float(rec['gradient']):.6e}"
                        )
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
        """Load a MODFLOW-USG-T SGB package from file."""
        if model.verbose:
            print("loading mfusg sgb package file...")
        if nper is None:
            nper = model.nper

        openfile = not hasattr(f, "read")
        if openfile:
            f = open(f, "r")

        # Skip comment lines, find the header line
        line = f.readline()
        while line.startswith("#"):
            line = f.readline()

        # Parse header: MXACTS ISGBCB [AUX varname ...] [NOPRINT]
        tokens = line.split()
        if tokens and tokens[0].upper() == "PARAMETER":
            raise NotImplementedError(
                "MfUsgSgb does not support named SGB parameters (NPSGB > 0)."
            )
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

        dtype = cls.get_default_dtype()
        if aux_names:
            dtype = Package.add_to_dtype(dtype, aux_names, np.float64)

        spd = {}
        current = None

        for kper in range(nper):
            line = f.readline()
            if not line:
                break
            sp_tokens = line.split()
            itmp = int(sp_tokens[0])
            np_sp = int(sp_tokens[1]) if len(sp_tokens) > 1 else 0
            if np_sp > 0:
                raise NotImplementedError(
                    "MfUsgSgb does not support active SGB parameters (NP > 0)."
                )

            if itmp < 0:
                if current is not None:
                    spd[kper] = current.copy()
            else:
                current = create_empty_recarray(itmp, dtype, default_value=0.0)
                for idx in range(itmp):
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
            ipakcb=ipakcb,
            stress_period_data=spd,
            dtype=dtype,
            options=options,
            extension="sgb",
            unitnumber=unitnumber,
            filenames=filenames,
        )

    @staticmethod
    def _ftype():
        return "SGB"

    @staticmethod
    def _defaultunit():
        return 159
