"""
mfusgriv module. Contains the MfUsgRiv class.
"""

from warnings import warn

import numpy as np

from ..modflow.mfriv import ModflowRiv
from ..pakbase import Package
from ..utils import MfList
from ..utils.recarray_utils import create_empty_recarray
from ..utils.utils_def import get_unitnumber_from_ext_unit_dict


class MfUsgRiv(ModflowRiv):
    """MODFLOW-USG River Package (RIV) with transport support.

    Supports unstructured grids (node-based) and AUX concentration variables.

    Parameters
    ----------
    model : MfUsg
        The model object to which this package will be added.
    stress_period_data : dict, optional
        Dictionary keyed by zero-based stress period index. Each value is a
        recarray with base fields ``node``, ``stage``, ``cond``, ``rbot``,
        plus any AUX fields.
    dtype : np.dtype, optional
        Custom dtype. Defaults to the unstructured base dtype plus declared
        AUX fields.
    ipakcb : int, optional
        IRIVCB, the second header field (GWF2RIV7U1AR): > 0 saves
        cell-by-cell flows on that unit, < 0 prints them when ICBCFL is not
        0, 0 (or None) neither.
    irdflag : int, optional
        Not part of the USG-T RIV input. Kept for backward compatibility of
        the call signature; it is ignored (it was formerly written where the
        engine reads IRIVCB). A value different from ``ipakcb`` warns.
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
    On unstructured grids the list is read by ULSTRDU: node, stage,
    conductance, bottom elevation and the declared AUX values, in the order
    of the AUX/AUXILIARY declarations. Values after them are not read by
    USG-T; on load they are ignored with a warning. An ``irch`` field that is
    not declared (the former undeclared trailing column) is rejected on
    write; a variable declared as ``AUX IRCH`` (any case) is an ordinary AUX
    variable and is kept.

    Examples
    --------
    >>> import flopy
    >>> m = flopy.mfusg.MfUsg()
    >>> riv = flopy.mfusg.MfUsgRiv.load("MDV.riv", m, nper=12)
    """

    # Extra fields that are rejected unless declared as AUX (the former
    # undeclared trailing column, which USG-T does not read).
    _UNDECLARED_REJECTED = frozenset({"irch"})

    def __init__(
        self,
        model,
        ipakcb=None,
        stress_period_data=None,
        dtype=None,
        irdflag=None,
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
        if irdflag is not None and irdflag != self.ipakcb:
            warn(
                f"MfUsgRiv: irdflag={irdflag} is ignored; the second header field "
                f"is IRIVCB, taken from ipakcb ({self.ipakcb})"
            )
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

    @staticmethod
    def _declared_aux(options):
        """AUX/AUXILIARY names declared in options, in order (whole tokens)."""
        names = []
        for opt in options:
            tokens = str(opt).split()
            i = 0
            while i < len(tokens):
                if tokens[i].upper() in ("AUX", "AUXILIARY") and i + 1 < len(tokens):
                    names.append(tokens[i + 1])
                    i += 2
                else:
                    i += 1
        return names

    def _update_aux_options(self):
        """Declare as AUX the extra dtype fields that are not declared yet
        (whole-token, case-insensitive match); an undeclared 'irch' is not
        declared automatically (it is rejected on write)."""
        n_base = len(self.get_default_dtype(structured=self.parent.structured).names)
        declared = {n.casefold() for n in self._declared_aux(self.options)}
        for name in self.dtype.names[n_base:]:
            if name.casefold() in declared:
                continue
            if name.casefold() in self._UNDECLARED_REJECTED:
                continue
            self.options.append(f"AUX {name.upper()}")
            declared.add(name.casefold())

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
        # checks done before the file is opened
        fields = {n.casefold(): n for n in self.dtype.names[n_base:]}
        declared = self._declared_aux(self.options)
        undeclared = [
            n
            for c, n in fields.items()
            if c in self._UNDECLARED_REJECTED
            and c not in {d.casefold() for d in declared}
        ]
        if undeclared:
            raise ValueError(
                f"MfUsgRiv: field(s) {undeclared} are not declared as AUX; an "
                "undeclared trailing column is not read by USG-T RIV "
                "(GWF2RIV7U1RP/ULSTRDU read node, stage, cond, rbot and the "
                "declared AUX values). Remove the field, or declare it "
                "explicitly (e.g. options=['AUX IRCH']) if it is an AUX variable"
            )
        missing = [d for d in declared if d.casefold() not in fields]
        if missing:
            raise ValueError(f"MfUsgRiv: AUX variable(s) {missing} have no data field")
        # AUX values are read in the order of the declarations
        aux_fields = [fields[d.casefold()] for d in declared]

        with open(self.fn_path, "w") as f:
            f.write(f"{self.heading}\n")

            # Header: MXACTR IRIVCB [options] (GWF2RIV7U1AR)
            line = f" {self.stress_period_data.mxact} {self.ipakcb}"
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
                        for name in aux_fields:
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

        # Parse header: MXACTR IRIVCB [options] (GWF2RIV7U1AR)
        tokens = line.split()
        ipakcb = int(tokens[1])
        options = []
        aux_names = []
        i = 2
        while i < len(tokens):
            if tokens[i].upper() in ("AUX", "AUXILIARY") and i + 1 < len(tokens):
                aux_names.append(tokens[i + 1])
                options.append(f"{tokens[i]} {tokens[i + 1]}")
                i += 2
            else:
                options.append(tokens[i])
                i += 1

        dtype = cls.get_default_dtype(structured=model.structured)
        if aux_names:
            dtype = Package.add_to_dtype(dtype, aux_names, np.float32)

        n_base = len(cls.get_default_dtype(structured=model.structured).names)
        n_expected = n_base + len(aux_names)

        spd = {}
        current = None
        extra_values = False

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

            current = create_empty_recarray(nact, dtype, default_value=0.0)
            for idx in range(nact):
                vals = f.readline().split()
                # USG-T (ULSTRDU) does not read values after the AUX values
                extra_values |= len(vals) > n_expected
                _fill_row(current, idx, vals, dtype)

            spd[kper] = current

        if openfile:
            f.close()
        if extra_values:
            warn(
                "MfUsgRiv.load: values after node, stage, cond, rbot and the "
                "declared AUX values are not read by USG-T and were ignored"
            )

        unitnumber, filenames = get_unitnumber_from_ext_unit_dict(
            model, cls, ext_unit_dict, ipakcb
        )

        return cls(
            model,
            ipakcb=ipakcb,
            stress_period_data=spd,
            dtype=dtype,
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
