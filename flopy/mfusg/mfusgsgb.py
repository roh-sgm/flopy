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

Named parameters (``NPSGB > 0``) are **definition-preserving only** as of Stage
4.4C (review follow-up): the leading ``PARAMETER NPSGB MXS`` record and the
per-parameter definitions (``UPARLSTRP`` header + ``NLST`` rows) round-trip
(load -> write -> reload), but **active SGB parameters are unsupported**. USG-T
2.7 defines SGB parameters as ``PARTYP='SGB'`` (``UPARLSTRP``, glo2sgbu1.f:97)
but activates them as ``PTYP='G'`` (``UPARLSTSUB``, glo2sgbu1.f:185); since a
parameter has a single type, any active SGB parameter trips a "Parameter type
conflict" (parutl7.f:684/800) and aborts the run. Therefore a per-stress-period
``NP>0`` raises ``NotImplementedError`` on load/write, parameter ``INSTANCES``
(``NUMINST>0``) raise ``NotImplementedError``, and inconsistent preserved state
(including ``MXS<=0`` or ``MXS`` below the total definition rows) raises
``ValueError``. Non-parametric authoring is unchanged.
"""

import numpy as np

from ..pakbase import Package
from ..utils import MfList
from ..utils.recarray_utils import create_empty_recarray
from ._usgt_list import begin_list_block
from ._usgt_parameters import (
    read_list_parameter_count,
    read_list_parameter_header,
    write_list_parameter_count,
    write_list_parameter_header,
)
from .mfusg import MfUsg

# USG-T 2.7 has a parameter-type mismatch for SGB: definitions are read with
# PARTYP='SGB' (UPARLSTRP, glo2sgbu1.f:97) but activations require PTYP='G'
# (UPARLSTSUB, glo2sgbu1.f:185). Since a parameter has a single PARTYP, any
# active SGB parameter trips "Parameter type conflict" (parutl7.f:684/800) and
# aborts the run. So active SGB parameters are not executable in USG-T 2.7;
# FloPy preserves SGB parameter *definitions* only.
_SGB_ACTIVE_PARAM_MSG = (
    "MfUsgSgb: active SGB parameters (per-stress-period NP>0) are not supported. "
    "USG-T 2.7 defines SGB parameters as PARTYP='SGB' (UPARLSTRP, "
    "glo2sgbu1.f:97) but activates them as PTYP='G' (UPARLSTSUB, "
    "glo2sgbu1.f:185), so any active SGB parameter aborts the Fortran with a "
    "'Parameter type conflict' (parutl7.f:684/800). Only SGB parameter "
    "definitions (no activations) are preserved."
)


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
        parameters=None,
        mxs=0,
        active_params=None,
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

        # Preserved SGB list-parameter state (set by load when NPSGB > 0):
        #   parameters: {name: {"partyp", "parval", "nlst", "data"}} (rows 0-based)
        #   mxs:        MXS from the leading PARAMETER NPSGB MXS record
        #   active_params: {kper: [name, ...]} -- always empty from load (active
        #     SGB parameters are unsupported; see _SGB_ACTIVE_PARAM_MSG). Kept so
        #     a manually constructed non-empty value is rejected on write.
        self.parameters = parameters
        self.mxs = mxs
        self.active_params = active_params if active_params is not None else {}

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

    def _write_sgb_rows(self, f, data, n_base):
        """Write SGB list rows: NODE (1-based) GRADIENT [aux ...]."""
        for rec in data:
            row = f" {int(rec['node']) + 1}  {float(rec['gradient']):.6e}"
            for name in self.dtype.names[n_base:]:
                row += f" {float(rec[name]):.6e}"
            f.write(row + "\n")

    def _validate_parameter_write(self):
        """Validate preserved SGB parameter state before a parameterized write.

        Active SGB parameters are unsupported (USG-T 2.7 defines them as
        PARTYP='SGB' but activates them as PTYP='G', so any activation aborts the
        run) and raise ``NotImplementedError``. Parameter *definitions* (with no
        activations) are preserved; inconsistent state raises ``ValueError`` so a
        ``PARAMETER`` header is never written without a complete, consistent body
        (including ``MXS`` large enough for the definition rows).
        """
        if any(self.active_params.values()):
            raise NotImplementedError(_SGB_ACTIVE_PARAM_MSG)
        if not self.parameters:
            return
        total = 0
        for name, pdef in self.parameters.items():
            missing = [k for k in ("partyp", "parval", "nlst", "data") if k not in pdef]
            if missing:
                raise ValueError(
                    f"MfUsgSgb.write_file: parameter '{name}' is missing keys "
                    f"{missing} (each definition needs partyp, parval, nlst, "
                    "data)."
                )
            if len(pdef["data"]) != pdef["nlst"]:
                raise ValueError(
                    f"MfUsgSgb.write_file: parameter '{name}' declares nlst="
                    f"{pdef['nlst']} but carries {len(pdef['data'])} rows."
                )
            total += pdef["nlst"]
        if self.mxs <= 0:
            raise ValueError(
                "MfUsgSgb.write_file: MXS (the PARAMETER NPSGB MXS record) must "
                f"be > 0 when SGB parameter definitions are present; got "
                f"mxs={self.mxs}. A from-scratch parameterized SGB without a "
                "valid MXS does not represent a loadable file."
            )
        if self.mxs < total:
            raise ValueError(
                f"MfUsgSgb.write_file: MXS ({self.mxs}) must be >= the total "
                f"number of parameter list entries ({total})."
            )

    def write_file(self):
        """Write the package file in MODFLOW-USG-T SGB format."""
        preserve = bool(self.parameters) or any(self.active_params.values())
        if preserve:
            # Validate before opening the file so a parameterized header is never
            # written without a complete, consistent body.
            self._validate_parameter_write()
        nper = self.parent.nper
        n_base = len(self.get_default_dtype().names)

        with open(self.fn_path, "w") as f:
            f.write(f"{self.heading}\n")

            # Item 1a (optional): PARAMETER NPSGB MXS
            if preserve:
                write_list_parameter_count(f, len(self.parameters), self.mxs)

            # Item 1: MXACTS ISGBCB [AUX ...] [NOPRINT]
            line = f" {self.stress_period_data.mxact:9d} {self.ipakcb}"
            for opt in self.options:
                line += f" {opt}"
            f.write(line + "\n")

            # Items 2-3: parameter definitions (UPARLSTRP header + NLST rows)
            if preserve:
                for name, pdef in self.parameters.items():
                    write_list_parameter_header(
                        f, name, pdef["partyp"], pdef["parval"], pdef["nlst"]
                    )
                    self._write_sgb_rows(f, pdef["data"], n_base)

            # Per stress period: ITMP NP. Active SGB parameters are unsupported
            # (see _validate_parameter_write), so NP is always 0 and only
            # non-parametric rows are written; parameter definitions above are
            # preserved but never activated.
            for kper in range(nper):
                if kper in self.stress_period_data.data:
                    kdata = self.stress_period_data[kper]
                    f.write(f" {len(kdata)} 0    Stress Period {kper + 1}\n")
                    self._write_sgb_rows(f, kdata, n_base)
                else:
                    f.write(f" -1 0    Stress Period {kper + 1}\n")

    # ------------------------------------------------------------------
    # load
    # ------------------------------------------------------------------

    @staticmethod
    def _read_sgb_rows(f, count, dtype, model, ext_unit_dict):
        """Read ``count`` SGB list rows (``NODE GRADIENT [aux ...]``) into a
        0-based recarray, honoring leading SFAC / OPEN-CLOSE / EXTERNAL controls.

        SFAC is inert on the gradient (the Fortran ISCLOC=2 scales an internal
        dummy column, not the gradient), so it is not applied to the stored data.
        """
        data = create_empty_recarray(count, dtype, default_value=0.0)
        if count == 0:
            return data
        source, _sfac, first_line, to_close = begin_list_block(
            f, model, ext_unit_dict, package="SGB"
        )
        for idx in range(count):
            row = first_line if idx == 0 else source.readline()
            vals = row.split()
            for j, name in enumerate(dtype.names):
                kind = dtype[name].kind
                data[idx][name] = int(vals[j]) if kind in ("i", "u") else float(vals[j])
            data[idx]["node"] -= 1
        if to_close is not None:
            source.close()
        return data

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

        # Optional item 1a: PARAMETER NPSGB MXS (UPARLSTAL); then the header line.
        npsgb, mxs = read_list_parameter_count(line)
        if npsgb > 0:
            line = f.readline()

        # Parse header: MXACTS ISGBCB [AUX varname ...] [NOPRINT]
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

        dtype = cls.get_default_dtype()
        if aux_names:
            dtype = Package.add_to_dtype(dtype, aux_names, np.float64)

        # Items 2-3: NPSGB list-parameter definitions (UPARLSTRP header + rows).
        parameters = None
        if npsgb > 0:
            parameters = {}
            for _ in range(npsgb):
                name, partyp, parval, nlst, numinst = read_list_parameter_header(
                    f.readline()
                )
                if numinst > 0:
                    raise NotImplementedError(
                        "MfUsgSgb.load: SGB parameter INSTANCES (NUMINST>0) are "
                        "not supported. SGB active parameters are unsupported in "
                        "USG-T 2.7 (PARTYP 'SGB' vs 'G' mismatch), so instanced "
                        "definitions could never be activated either."
                    )
                data = cls._read_sgb_rows(f, nlst, dtype, model, ext_unit_dict)
                parameters[name] = {
                    "partyp": partyp,
                    "parval": parval,
                    "nlst": nlst,
                    "data": data,
                }

        spd = {}
        active_params = {}
        current = None

        for kper in range(nper):
            line = f.readline()
            if not line:
                break
            sp_tokens = line.split()
            itmp = int(sp_tokens[0])
            np_sp = int(sp_tokens[1]) if len(sp_tokens) > 1 else 0

            # Non-parametric list (ITMP<0 reuses the previous period's rows).
            if itmp < 0:
                if current is not None:
                    spd[kper] = current.copy()
            else:
                current = cls._read_sgb_rows(f, itmp, dtype, model, ext_unit_dict)
                spd[kper] = current

            # Active SGB parameters (per-SP NP>0) are not executable in USG-T 2.7
            # (PARTYP='SGB' definition vs PTYP='G' activation -> type conflict).
            if np_sp > 0:
                raise NotImplementedError(_SGB_ACTIVE_PARAM_MSG)

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
            parameters=parameters,
            mxs=mxs,
            active_params=active_params,
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
