"""
mfusgriv module. Contains the MfUsgRiv class.
"""

import numpy as np

from ..modflow.mfriv import ModflowRiv
from ..pakbase import Package
from ..utils import MfList
from ..utils.recarray_utils import create_empty_recarray


class MfUsgRiv(ModflowRiv):
    """MODFLOW-USG River Package (RIV) with transport support.

    Supports unstructured grids (node-based), AUX concentration variables,
    and an optional trailing reach-ID column (``irch``).

    Parameters
    ----------
    model : MfUsg
        The model object to which this package will be added.
    stress_period_data : dict, optional
        Dictionary keyed by zero-based stress period index. Each value is a
        recarray with base fields ``node``, ``stage``, ``cond``, ``rbot``,
        plus any AUX fields and optionally ``irch``.
    dtype : np.dtype, optional
        Custom dtype. Defaults to the unstructured base dtype plus declared
        AUX fields.
    irdflag : int
        Read-data flag written to the package header (default 50).
    options : list of str, optional
        Package options (e.g. ``["AUX C01"]``).
    extension : str
        Filename extension (default ``"riv"``).
    unitnumber : int, optional
        File unit number.
    filenames : str or list of str, optional
        Package filename(s).
    add_package : bool
        Add package to model on construction (default True).

    Notes
    -----
    The ``irch`` field (reach identifier) is NOT declared as an AUX variable
    in the file header. It is written and read as a trailing integer column.

    Examples
    --------
    >>> import flopy
    >>> m = flopy.mfusg.MfUsg()
    >>> riv = flopy.mfusg.MfUsgRiv.load("MDV.riv", m, nper=12)
    """

    # Fields beyond the base dtype that are written positionally, not as AUX.
    _NON_AUX_TRAILING = frozenset({"irch"})

    def __init__(
        self,
        model,
        ipakcb=None,
        stress_period_data=None,
        dtype=None,
        irdflag=50,
        options=None,
        extension="riv",
        unitnumber=None,
        filenames=None,
        add_package=True,
    ):
        if unitnumber is None:
            unitnumber = ModflowRiv._defaultunit()
        if options is None:
            options = []

        # Prepare two filename slots (package input + CBC output) and register
        # the CBC output file — mirrors ModflowRiv.__init__ exactly.
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
        self.url = "riv.html"
        self.np = 0
        self.irdflag = irdflag
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
            return ModflowRiv.get_default_dtype(structured=True)
        return np.dtype([
            ("node", int),
            ("stage", np.float64),
            ("cond", np.float32),
            ("rbot", np.float32),
        ])

    @staticmethod
    def get_empty(ncells=0, aux_names=None, structured=True):
        dtype = MfUsgRiv.get_default_dtype(structured=structured)
        if aux_names is not None:
            dtype = Package.add_to_dtype(dtype, aux_names, np.float32)
        return create_empty_recarray(ncells, dtype, default_value=-1.0e10)

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _update_aux_options(self):
        """Append AUX entries to options for non-trailing extra dtype fields."""
        n_base = len(self.get_default_dtype(structured=self.parent.structured).names)
        for name in self.dtype.names[n_base:]:
            if name.lower() in self._NON_AUX_TRAILING:
                continue
            token = f"AUX {name.upper()}"
            if not any(name.lower() in o.lower() for o in self.options):
                self.options.append(token)

    # ------------------------------------------------------------------
    # write_file
    # ------------------------------------------------------------------

    def write_file(self):
        """Write the package file in MODFLOW-USG-T RIV format."""
        if self.parent.structured:
            ModflowRiv.write_file(self)
            return
        nper = self.parent.nper
        n_base = len(self.get_default_dtype(structured=self.parent.structured).names)

        with open(self.fn_path, "w") as f:
            f.write(f"{self.heading}\n")

            # Header: MXACTR IRDFLAG [AUX ...]
            line = f" {self.stress_period_data.mxact} {self.irdflag}"
            for opt in self.options:
                line += f" {opt}"
            f.write(line + "\n")

            for kper in range(nper):
                if kper in self.stress_period_data.data:
                    kdata = self.stress_period_data[kper]
                    f.write(f" {len(kdata)} 0    Stress Period {kper + 1}\n")
                    for rec in kdata:
                        row = (f" {int(rec['node']) + 1}"
                               f"  {float(rec['stage']):.6f}"
                               f"  {float(rec['cond']):.6e}"
                               f"  {float(rec['rbot']):.6f}")
                        for name in self.dtype.names[n_base:]:
                            if name.lower() in self._NON_AUX_TRAILING:
                                row += f" {int(rec[name])}"
                            else:
                                row += f" {float(rec[name]):.6e}"
                        f.write(row + "\n")
                else:
                    f.write(f" -1 0    Stress Period {kper + 1}\n")

    # ------------------------------------------------------------------
    # load
    # ------------------------------------------------------------------

    @classmethod
    def load(cls, f, model, nper=None, ext_unit_dict=None, check=True):
        """Load a MODFLOW-USG-T RIV package from file.

        Parameters
        ----------
        f : str or file-like
            Path to the .riv file, or an open file handle.
        model : MfUsg
            Model to attach the package to.
        nper : int, optional
            Number of stress periods. Defaults to ``model.nper``.

        Returns
        -------
        MfUsgRiv
        """
        if model.structured:
            return ModflowRiv.load(f, model, nper=nper, ext_unit_dict=ext_unit_dict, check=check)

        if model.verbose:
            print("loading mfusg riv package file...")
        if nper is None:
            nper = model.nper

        openfile = not hasattr(f, "read")
        if openfile:
            f = open(f, "r")

        # Skip comment lines, find the header line
        line = f.readline()
        while line.startswith("#"):
            line = f.readline()

        # Parse header: MXACTR IRDFLAG [AUX varname ...]
        tokens = line.split()
        irdflag = int(tokens[1])
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

        n_base = len(cls.get_default_dtype(structured=model.structured).names)
        n_expected = n_base + len(aux_names)

        # irch detection is deferred to first non-empty data row
        has_irch = None
        spd = {}
        current = None

        for kper in range(nper):
            line = f.readline()
            if not line:
                break
            nact = int(line.split()[0])
            # second token is CLN flag (0 = none), ignored

            if nact < 0:
                # Reuse previous SP data
                if current is not None:
                    spd[kper] = current.copy()
                continue

            # On first non-empty SP, peek at the first data row to detect irch
            if has_irch is None and nact > 0:
                first_line = f.readline()
                first_vals = first_line.split()
                has_irch = len(first_vals) > n_expected
                if has_irch:
                    dtype = np.dtype(list(dtype.descr) + [("irch", np.int32)])

                current = create_empty_recarray(nact, dtype, default_value=0.0)
                _fill_row(current, 0, first_vals, dtype)
                for idx in range(1, nact):
                    _fill_row(current, idx, f.readline().split(), dtype)
            else:
                current = create_empty_recarray(nact, dtype, default_value=0.0)
                for idx in range(nact):
                    _fill_row(current, idx, f.readline().split(), dtype)

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
            irdflag=irdflag,
            options=options,
            extension="riv",
            unitnumber=unitnumber,
            filenames=filenames,
        )


def _fill_row(arr, i, vals, dtype):
    """Fill one row of a recarray from a list of string tokens."""
    for j, name in enumerate(dtype.names):
        kind = dtype[name].kind
        arr[i][name] = int(vals[j]) if kind in ("i", "u") else float(vals[j])
    arr[i]["node"] -= 1
