"""
mfusgdrt module. Contains the MfUsgDrt class.

Drain Return (DRT8) Package -- USG-Transport.

Reference: Fortran subroutines ``GWF2DRT8U1AR`` / ``GWF2DRT8U1RP`` /
``SGWF2DRT8LR`` in ``gwf2drt8u.f`` of USG-Transport 2.7. DRT8 extends the
MODFLOW drain with return flow: a proportion of the drained water can be
returned to a single recipient node, or spread over several recipient nodes,
optionally with a transport concentration-change type.

File layout (unstructured, free format)::

    Item 1:  MXADRT IDRTCB NPDRT MXL [options]
             options: [AUX <name> ...] [RETURNFLOW] [NOPRINT]
                      [SPREAD <MXSPREADNDS>] [CHANGEC]
    Per SP:  ITMP
             ITMP drain lines, each:
                NODE  EL  COND  NR  [PROP]  [IDCHNGTYP]  [aux ...]
             NR > 0 : single recipient node (inline)
             NR < 0 : -NR spreading recipient nodes follow this line as a
                      U1DINT block
             NR = 0 : no return flow

``ITMP < 0`` reuses the previous stress period's drains and recipients.

Internal ``node`` and recipient values are 0-based; the file is 1-based.
A recipient list of length 1 is written inline (``NR > 0``); length > 1 is
written as a spreading block (``NR < 0``). The two forms are physically
equivalent for a single node, so spreading-of-one normalizes to inline.

The spreading ``U1DINT`` block is read with the full set of control records
(Stage 4.5B): ``INTERNAL`` / ``CONSTANT`` (inline) and ``EXTERNAL <unit>`` /
``OPEN/CLOSE <fname>`` (resolved via ``ext_unit_dict`` / ``model.model_ws``),
matching USG-T's ``U1DINT`` (utl7u1.f). On write the recipients are always
emitted inline (``INTERNAL (FREE)``, ``Expanded valid write``) -- EXTERNAL /
OPEN-CLOSE are not preserved.

Notes
-----
The USG-T transport extensions (RETURNFLOW recipient nodes, ``CHANGEC`` /
``IDCHNGTYP``, ``SPREAD``) are implemented for **unstructured** grids only.
Constructing ``MfUsgDrt`` on a structured (DIS) model raises
``NotImplementedError`` — use :class:`flopy.modflow.ModflowDrt` for classic
structured DRT. (On load, ``MfUsgDrt.load`` delegates structured files to
``ModflowDrt.load``, which returns a base-class object.)

Named parameters (``NPDRT > 0``) are **preserved** (load -> write -> reload) as of
Stage 4.4D. DRT is internally consistent in USG-T 2.7: definitions and
activations both use ``PARTYP='DRT'`` (``UPARLSTRP`` at ``gwf2drt8u.f:135`` and
``SGWF2DRT8LS`` at ``gwf2drt8u.f:1108``), so an active DRT parameter is valid
(unlike SGB, whose ``'SGB'`` vs ``'G'`` mismatch aborts the run). The
``PARAMETER`` count lives in item 1 (``MXADRT IDRTCB NPDRT MXL``); each parameter
owns ``NLST`` drain rows (with their RETURNFLOW recipients / spreading blocks),
and its value scales ``COND`` (Fortran ``IPVL=5``). Per stress period, ``ITMP NP``
gives the non-parametric count and the number of active parameters; each
activation is a parameter name.

Not supported (explicit failure rather than partial write):

* Parameter ``INSTANCES`` (``NUMINST>0``): the Fortran supports them, but FloPy
  does not yet (instances combined with per-row recipient lists) -> load raises
  ``NotImplementedError``; the from-scratch authoring API has no instance field.
* An ``EXTERNAL`` spreading-node unit that cannot be resolved via
  ``ext_unit_dict`` -> ``NotImplementedError`` (inline the nodes instead).

From-scratch ``NPDRT>0`` parameter **authoring** is supported (Stage 4.6B): pass
``parameters={name: {...}}`` and ``active_params={kper: [name, ...]}`` (see the
``parameters`` argument). ``active_params`` that reference an undefined name (or
no ``parameters`` at all) raise ``ValueError``. Activated ``SPREAD`` (``NR<0``)
recipients still round-trip structurally but are not execution-guaranteed (USG-T
copies ``DRTF`` but not ``NodDRT`` on activation).
"""

import numpy as np

from ..modflow.mfdrt import ModflowDrt
from ..pakbase import Package
from ._usgt_list import begin_list_block
from ._usgt_parameters import (
    check_parameter_name,
    check_parval,
    read_active_list_parameters,
    read_list_parameter_header,
    resolve_list_parameter,
    write_active_list_parameters,
    write_list_parameter_header,
)
from ._usgt_returnflow import read_u1dint_list, write_u1dint_list
from .mfusg import MfUsg


class MfUsgDrt(ModflowDrt):
    """MODFLOW-USG Drain Return (DRT8) Package with transport extensions.

    Parameters
    ----------
    model : MfUsg
        The model object to which this package will be added.
    ipakcb : int, optional
        Unit number for cell-by-cell budget output (``IDRTCB``).
    stress_period_data : dict, optional
        Dictionary keyed by zero-based stress period. Each value is a recarray
        matching ``dtype``. Base fields are ``node`` (0-based), ``elev`` and
        ``cond``; with ``RETURNFLOW``, ``rfprop``; with ``CHANGEC``,
        ``idchngtyp``; plus any AUX fields.
    recipient_nodes : dict, optional
        Dictionary keyed by zero-based stress period. Each value is a list
        (one entry per drain record, in the same order as
        ``stress_period_data``) of 0-based recipient node lists.
    dtype : np.dtype, optional
        Custom dtype. If None, derived from the active options.
    options : list of str, optional
        Package options, e.g. ``["RETURNFLOW", "CHANGEC", "AUX C01"]``. A
        ``SPREAD`` token is added automatically on write when needed.
    parameters : dict, optional
        Named ``NPDRT`` list-parameter definitions, keyed by parameter name.
        Each value is a dict (Stage 4.6B from-scratch authoring; also the form
        ``load`` preserves):

        * ``parval`` (required) — the parameter value (scales ``COND``);
        * ``data`` (required) — the parameter's drain rows: a recarray matching
          the active ``dtype``, or any array-like the dtype can build (e.g. a
          list of tuples), nodes 0-based;
        * ``partyp`` (optional) — defaults to ``"DRT"`` and must be ``"DRT"``;
        * ``nlst`` (optional) — number of rows; computed from ``data`` if
          omitted, validated if given;
        * ``recipient_nodes`` (optional, ``RETURNFLOW`` only) — one 0-based
          recipient list per row; defaults to empty lists if omitted.
    mxl : int, optional
        ``MXL`` (item 1): max parameter list entries. Computed as the total of
        all definition ``nlst`` when omitted/0; validated to be >= that total if
        given.
    active_params : dict, optional
        ``{kper: [name, ...]}`` — parameters activated per zero-based stress
        period (``ITMP NP`` records). Names must be defined in ``parameters``
        (case-insensitive); duplicates within a period and out-of-range ``kper``
        raise ``ValueError``.
    extension : str
        Filename extension (default ``"drt"``).
    unitnumber : int, optional
        File unit number.
    filenames : str or list of str, optional
        Package filename(s).
    add_package : bool
        Add package to model on construction (default True).
    """

    def __init__(
        self,
        model,
        ipakcb=None,
        stress_period_data=None,
        recipient_nodes=None,
        dtype=None,
        options=None,
        parameters=None,
        mxl=0,
        active_params=None,
        extension="drt",
        unitnumber=None,
        filenames=None,
        add_package=True,
    ):
        assert isinstance(model, MfUsg), (
            "Model object must be of type flopy.mfusg.MfUsg\n"
            f"but received type: {type(model)}."
        )

        # MfUsgDrt implements the *unstructured* USG-T DRT8 format only. For
        # structured (DIS) grids the classic MODFLOW DRT has no USG-T return-flow
        # extensions, so authoring goes through the base class instead. (On
        # load, MfUsgDrt.load delegates structured files to ModflowDrt.load.)
        if model.structured:
            raise NotImplementedError(
                "MfUsgDrt implements the unstructured USG-T DRT8 format; for "
                "structured (DIS) grids use flopy.modflow.ModflowDrt instead."
            )

        if unitnumber is None:
            unitnumber = ModflowDrt._defaultunit()
        if options is None:
            options = []

        filenames = self._prepare_filenames(filenames, 2)
        self.set_cbc_output_file(ipakcb, model, filenames[1])

        # Build the package directly (mirrors MfUsgGhb: we do not call
        # ModflowDrt.__init__ because the USG-T storage model differs).
        Package.__init__(
            self,
            model,
            extension=extension,
            name=self._ftype(),
            unit_number=unitnumber,
            filenames=filenames[0],
        )
        self._generate_heading()
        self.url = "drt.html"
        # Strip any SPREAD token; the value is recomputed from the data on write.
        self.options = self._strip_spread(options)

        self.returnflow = any(
            o.split()[0].upper() == "RETURNFLOW" for o in self.options
        )
        self.changec = self.returnflow and any(
            o.split()[0].upper() == "CHANGEC" for o in self.options
        )

        aux_names = self._aux_names_from_options(self.options)
        if dtype is not None:
            self.dtype = dtype
        else:
            self.dtype = self.get_usg_dtype(
                returnflow=self.returnflow, changec=self.changec, aux_names=aux_names
            )

        self.stress_period_data = {
            kper: self._to_recarray(val)
            for kper, val in (stress_period_data or {}).items()
        }
        self.recipient_nodes = recipient_nodes or {}

        # Preserved DRT list-parameter state (set by load when NPDRT > 0):
        #   parameters: {name: {"partyp", "parval", "nlst", "data", "recipient_nodes"}}
        #     (rows 0-based; recipient_nodes is one 0-based list per definition row)
        #   mxl:        MXL from item 1 (max parameter list entries)
        #   active_params: {kper: [name, ...]} active parameters per stress period
        self.parameters = parameters
        self.mxl = mxl
        self.active_params = active_params if active_params is not None else {}

        if add_package:
            self.parent.add_package(self)

    # ------------------------------------------------------------------
    # Static / helper constructors
    # ------------------------------------------------------------------

    @staticmethod
    def get_usg_dtype(returnflow=True, changec=False, aux_names=None):
        """Build the unstructured USG-T DRT dtype for the active options."""
        cols = [("node", int), ("elev", np.float64), ("cond", np.float64)]
        if returnflow:
            cols.append(("rfprop", np.float64))
            if changec:
                cols.append(("idchngtyp", int))
        dtype = np.dtype(cols)
        if aux_names:
            dtype = Package.add_to_dtype(dtype, aux_names, np.float64)
        return dtype

    @staticmethod
    def _aux_names_from_options(options):
        names = []
        for opt in options:
            toks = opt.split()
            if toks and toks[0].upper() in ("AUX", "AUXILIARY") and len(toks) > 1:
                names.append(toks[1])
        return names

    @staticmethod
    def _strip_spread(options):
        kept = []
        skip = False
        for opt in options:
            if skip:
                skip = False
                continue
            toks = opt.split()
            if toks and toks[0].upper() == "SPREAD":
                # token may be "SPREAD" alone (value next) or "SPREAD n"
                if len(toks) == 1:
                    skip = True
                continue
            kept.append(opt)
        return kept

    def _base_field_count(self):
        n = 3  # node, elev, cond
        if self.returnflow:
            n += 1  # rfprop
            if self.changec:
                n += 1  # idchngtyp
        return n

    def _aux_field_names(self):
        return list(self.dtype.names[self._base_field_count() :])

    def _to_recarray(self, data):
        if data is None:
            return None
        if isinstance(data, np.recarray):
            return data
        return np.array(data, dtype=self.dtype).view(np.recarray)

    # ------------------------------------------------------------------
    # write_file
    # ------------------------------------------------------------------

    def _validated_recipients(self, kper, nrec):
        """Recipient lists for a stress period, validated against its records.

        With RETURNFLOW active there must be exactly one recipient list per
        drain record. Omitting ``recipient_nodes`` for a period means every
        record has zero recipients. Mismatched lengths raise so spreading /
        return-flow metadata cannot be silently dropped or shifted on write.
        """
        if not self.returnflow:
            return [[] for _ in range(nrec)]
        recips = self.recipient_nodes.get(kper)
        if recips is None:
            return [[] for _ in range(nrec)]
        if len(recips) != nrec:
            raise ValueError(
                f"MfUsgDrt: recipient_nodes[{kper}] has {len(recips)} entries "
                f"but stress_period_data[{kper}] has {nrec} records; provide "
                "exactly one recipient list per record (or omit the period for "
                "all-zero recipients)."
            )
        return recips

    @staticmethod
    def _check_param_name(name, what):
        """Validate a parameter/active name (shared helper; see
        :func:`flopy.mfusg._usgt_parameters.check_parameter_name`)."""
        check_parameter_name(name, what, prefix="MfUsgDrt.write_file")

    @staticmethod
    def _check_parval(name, parval):
        """Validate ``parval`` (shared helper; see
        :func:`flopy.mfusg._usgt_parameters.check_parval`)."""
        return check_parval(parval, name, prefix="MfUsgDrt.write_file")

    def _normalize_param(self, name, pdef):
        """Validate + canonicalize one DRT parameter definition for writing.

        Accepts an ergonomic from-scratch dict *or* an already-canonical loaded
        one and returns ``{partyp, parval, nlst, data (recarray),
        recipient_nodes}``. ``data`` may be a recarray or any array-like the
        active ``dtype`` can build; ``nlst`` and ``recipient_nodes`` are computed
        when omitted and validated when given. ``partyp`` defaults to ``DRT`` and
        must be ``DRT`` (the Fortran activates DRT params as ``PARTYP='DRT'``).
        Raises ``ValueError`` for invalid input, before any file is opened.
        """
        self._check_param_name(name, "parameter name")
        if not isinstance(pdef, dict):
            raise ValueError(
                f"MfUsgDrt.write_file: parameter '{name}' must be a dict, got "
                f"{type(pdef).__name__}."
            )
        partyp = pdef.get("partyp", "DRT")
        if str(partyp).upper() != "DRT":
            raise ValueError(
                f"MfUsgDrt.write_file: parameter '{name}' partyp must be 'DRT'; "
                f"got {partyp!r}."
            )
        parval = self._check_parval(name, pdef.get("parval"))
        data = pdef.get("data")
        if data is None or len(data) == 0:
            raise ValueError(
                f"MfUsgDrt.write_file: parameter '{name}' needs non-empty 'data' "
                "drain rows."
            )
        if not isinstance(data, np.recarray):
            data = np.array(data, dtype=self.dtype).view(np.recarray)
        if np.any(np.asarray(data["node"]) < 0):
            raise ValueError(
                f"MfUsgDrt.write_file: parameter '{name}' has a negative node; "
                "drain nodes are 0-based and must be non-negative."
            )
        nlst = pdef.get("nlst")
        if nlst is None:
            nlst = len(data)
        elif nlst != len(data):
            raise ValueError(
                f"MfUsgDrt.write_file: parameter '{name}' declares nlst={nlst} "
                f"but carries {len(data)} rows."
            )
        recips = pdef.get("recipient_nodes")
        if recips is None:
            recips = [[] for _ in range(nlst)]
        if len(recips) != nlst:
            raise ValueError(
                f"MfUsgDrt.write_file: parameter '{name}' has {len(recips)} "
                f"recipient lists but nlst={nlst} (need one per definition row)."
            )
        if not self.returnflow and any(len(r) > 0 for r in recips):
            raise ValueError(
                f"MfUsgDrt.write_file: parameter '{name}' has recipient_nodes but "
                "the RETURNFLOW option is not enabled."
            )
        for ri, row in enumerate(recips):
            for nd in row:
                if not float(nd).is_integer() or int(nd) < 0:
                    raise ValueError(
                        f"MfUsgDrt.write_file: parameter '{name}' recipient_nodes "
                        f"row {ri} has invalid node {nd!r}; recipients are 0-based "
                        "non-negative integers."
                    )
        return {
            "partyp": "DRT",
            "parval": parval,
            "nlst": nlst,
            "data": data,
            "recipient_nodes": [list(r) for r in recips],
        }

    def _validate_parameter_write(self):
        """Normalize + validate DRT parameter state for a parameterized write.

        Supports both preserved (loaded) and **from-scratch** parameter dicts
        (Stage 4.6B). Returns ``(params, mxl)`` with canonical definitions and the
        resolved ``MXL`` (``= sum(nlst)`` when ``mxl`` is omitted). All validation
        runs before the file is opened, so a parameterized item 1 is never written
        without a complete, consistent body and no partial file is produced.
        ``INSTANCES`` remains unsupported (load raises; the authoring API has no
        instance field).
        """
        if not self.parameters:
            raise ValueError(
                "MfUsgDrt.write_file: active_params reference parameters but none "
                "are defined. Pass parameters={name: {...}} to author DRT "
                "parameter definitions from scratch."
            )
        params = {}
        total = 0
        for name, pdef in self.parameters.items():
            params[name] = self._normalize_param(name, pdef)
            total += params[name]["nlst"]
        lowered_defs = [name.lower() for name in params]
        if len(set(lowered_defs)) != len(lowered_defs):
            raise ValueError(
                "MfUsgDrt.write_file: duplicate parameter definition name "
                f"(case-insensitive): {list(params)}. The Fortran upper-cases "
                "PARNAM, so definition names must be unique ignoring case."
            )
        mxl = self.mxl if self.mxl else total
        if mxl < total:
            raise ValueError(
                f"MfUsgDrt.write_file: MXL ({mxl}) must be >= the total number of "
                f"parameter list entries ({total})."
            )
        defined = {name.lower() for name in params}
        nper = self.parent.nper
        for kper, names in self.active_params.items():
            if not 0 <= kper < nper:
                raise ValueError(
                    f"MfUsgDrt.write_file: active_params stress period {kper} is "
                    f"out of range 0..{nper - 1}."
                )
            for nm in names:
                self._check_param_name(nm, f"active parameter (stress period {kper})")
            lowered = [nm.lower() for nm in names]
            if len(set(lowered)) != len(lowered):
                raise ValueError(
                    f"MfUsgDrt.write_file: a parameter is activated more than "
                    f"once in stress period {kper}: {names}. USG-T aborts when a "
                    "parameter is already active this stress period."
                )
            for nm in names:
                if nm.lower() not in defined:
                    raise ValueError(
                        f"MfUsgDrt.write_file: active parameter '{nm}' (stress "
                        f"period {kper}) is not defined in parameters."
                    )
        return params, mxl

    def write_file(self):
        """Write the package file in MODFLOW-USG-T DRT8 format.

        Unstructured only; ``__init__`` rejects structured models, so there is
        no structured-delegation path here.
        """
        preserve = bool(self.parameters) or any(self.active_params.values())
        params, mxl = ({}, 0)
        if preserve:
            # Normalizes from-scratch / preserved definitions and validates the
            # whole parameterized state before the file is opened.
            params, mxl = self._validate_parameter_write()
        nper = self.parent.nper
        aux_names = self._aux_field_names()
        has_aux = len(aux_names) > 0

        mxadrt = self._max_active_drains(params)
        mxspread = self._max_spread_nodes(params)
        npdrt = len(params)

        with open(self.fn_path, "w") as f:
            f.write(f"{self.heading}\n")

            # Item 1: MXADRT IDRTCB NPDRT MXL [options]
            line = f" {mxadrt:9d} {self.ipakcb} {npdrt} {mxl}"
            for opt in self.options:
                line += f" {opt}"
            if mxspread > 0:
                line += f" SPREAD {mxspread}"
            f.write(line + "\n")

            # Items 2-3: parameter definitions (UPARLSTRP header + NLST rows,
            # each with its RETURNFLOW recipients / spreading block).
            for name, pdef in params.items():
                write_list_parameter_header(
                    f, name, pdef["partyp"], pdef["parval"], pdef["nlst"]
                )
                recips = pdef["recipient_nodes"]
                for i, rec in enumerate(pdef["data"]):
                    self._write_drain_line(f, rec, recips[i], has_aux, aux_names)

            for kper in range(nper):
                active = self.active_params.get(kper, []) if preserve else []
                if kper not in self.stress_period_data:
                    f.write(f" -1 {len(active)}    Stress Period {kper + 1}\n")
                    write_active_list_parameters(f, active)
                    continue
                recarray = self.stress_period_data[kper]
                recips = self._validated_recipients(kper, len(recarray))
                f.write(f" {len(recarray)} {len(active)}    Stress Period {kper + 1}\n")
                for i, rec in enumerate(recarray):
                    self._write_drain_line(f, rec, recips[i], has_aux, aux_names)
                # Active-parameter records (SGWF2DRT8LS) follow the non-param rows.
                write_active_list_parameters(f, active)

    def _max_active_drains(self, params):
        """Max active drain-return cells in any stress period (the Fortran's
        NDRTCL, which must not exceed MXADRT).

        Each active parameter contributes its NLST rows on top of the
        non-parametric drains: ``NDRTCL = NDRTNP + sum(nlst of active params)``
        (gwf2drt8u.f: SGWF2DRT8LS does ``NDRTCL = NDRTCL + NLST`` and aborts if
        ``NDRTCL > MXADRT``). ``ITMP<0`` reuses the previous period's
        non-parametric count. ``params`` is the canonical (normalized) parameter
        dict from ``_validate_parameter_write`` (empty when not parameterized).
        """
        mx = 0
        prev_nonparam = 0
        for kper in range(self.parent.nper):
            if kper in self.stress_period_data:
                prev_nonparam = len(self.stress_period_data[kper])
            nonparam = prev_nonparam  # ITMP<0 reuse keeps the previous count
            active = 0
            for name in self.active_params.get(kper, []):
                pdef = resolve_list_parameter(params, name)
                if pdef is not None:
                    active += pdef["nlst"]
            mx = max(mx, nonparam + active)
        return mx

    def _max_spread_nodes(self, params):
        """Max total spreading-recipient nodes (non-parametric per SP, and the
        canonical parameter definitions, which the Fortran also stores in
        NodDRT)."""
        mx = 0
        for kper, recarray in self.stress_period_data.items():
            recips = self.recipient_nodes.get(kper, [])
            total = sum(len(r) for r in recips if len(r) > 1)
            mx = max(mx, total)
        if params:
            def_total = sum(
                len(r)
                for pdef in params.values()
                for r in pdef.get("recipient_nodes", [])
                if len(r) > 1
            )
            mx = max(mx, def_total)
        return mx

    def _write_drain_line(self, f, rec, rnodes, has_aux, aux_names):
        n = len(rnodes)
        line = (
            f" {int(rec['node']) + 1}"
            f"  {float(rec['elev']):.6e}"
            f"  {float(rec['cond']):.6e}"
        )
        if self.returnflow:
            if n == 0:
                nr = 0
            elif n == 1:
                nr = int(rnodes[0]) + 1  # inline single recipient
            else:
                nr = -n  # spreading ground
            line += f"  {nr}"
            if nr != 0 or has_aux:
                line += f"  {float(rec['rfprop']):.6e}"
            if self.changec:
                line += f"  {int(rec['idchngtyp'])}"
        for name in aux_names:
            line += f" {float(rec[name]):.6e}"
        f.write(line + "\n")
        # Spreading block immediately follows the line (interleaved).
        if self.returnflow and n > 1:
            write_u1dint_list(f, [int(x) + 1 for x in rnodes])

    @classmethod
    def _read_drain_rows(
        cls, f, count, returnflow, changec, naux, dtype, model, ext_unit_dict
    ):
        """Read ``count`` drain rows (+ interleaved spreading U1DINT blocks) into
        a recarray and a per-row recipient list. Honors leading SFAC / EXTERNAL /
        OPEN-CLOSE list controls; SFAC scales COND (Fortran ISCLOC=5)."""
        if count == 0:
            return np.array([], dtype=dtype).view(np.recarray), []
        source, sfac, first_line, to_close = begin_list_block(
            f, model, ext_unit_dict, package="DRT"
        )
        records = []
        recip_lists = []
        for idx in range(count):
            row = first_line if idx == 0 else source.readline()
            rec, recips = cls._parse_drain_tokens(
                row.split(), returnflow, changec, naux, source, model, ext_unit_dict
            )
            records.append(rec)
            recip_lists.append(recips)
        if to_close is not None:
            source.close()
        recarray = np.array(records, dtype=dtype).view(np.recarray)
        if sfac != 1.0 and len(recarray) > 0:
            recarray["cond"] = recarray["cond"] * sfac
        return recarray, recip_lists

    # ------------------------------------------------------------------
    # load
    # ------------------------------------------------------------------

    @classmethod
    def load(cls, f, model, nper=None, ext_unit_dict=None, check=True):
        """Load a MODFLOW-USG-T DRT8 package from file."""
        if model.structured:
            return ModflowDrt.load(
                f, model, nper=nper, ext_unit_dict=ext_unit_dict, check=check
            )

        if model.verbose:
            print("loading mfusg drt package file...")
        if nper is None:
            nper = model.nper

        openfile = not hasattr(f, "read")
        if openfile:
            f = open(f, "r")

        line = f.readline()
        while line.startswith("#"):
            line = f.readline()

        options, aux_names, ipakcb, returnflow, changec, npdrt, mxl = cls._parse_header(
            line
        )
        dtype = cls.get_usg_dtype(
            returnflow=returnflow, changec=changec, aux_names=aux_names
        )
        naux = len(aux_names)

        # Items 2-3: NPDRT parameter definitions (UPARLSTRP header + NLST rows).
        parameters = None
        if npdrt > 0:
            parameters = {}
            for _ in range(npdrt):
                name, partyp, parval, nlst, numinst = read_list_parameter_header(
                    f.readline()
                )
                if numinst > 0:
                    raise NotImplementedError(
                        "MfUsgDrt.load: DRT parameter INSTANCES (NUMINST>0) are "
                        "supported by the Fortran but not yet by FloPy (instances "
                        "combined with per-row recipient lists)."
                    )
                data, recips = cls._read_drain_rows(
                    f, nlst, returnflow, changec, naux, dtype, model, ext_unit_dict
                )
                parameters[name] = {
                    "partyp": partyp,
                    "parval": parval,
                    "nlst": nlst,
                    "data": data,
                    "recipient_nodes": recips,
                }

        spd = {}
        recipient_nodes = {}
        active_params = {}
        prev_recarray = None
        prev_recips = None

        for kper in range(nper):
            line = f.readline()
            if not line:
                break
            parts = line.split()
            itmp = int(parts[0])
            # The per-SP header is "ITMP NP" only when parameters exist; with
            # NPDRT==0 it is just "ITMP" (parts[1:] may be an inline comment).
            np_sp = int(parts[1]) if (npdrt > 0 and len(parts) > 1) else 0

            if itmp < 0:
                # Reuse the previous period's non-parametric drains + recipients.
                if prev_recarray is not None:
                    spd[kper] = prev_recarray.copy()
                    recipient_nodes[kper] = [list(r) for r in prev_recips]
            else:
                recarray, recip_lists = cls._read_drain_rows(
                    f, itmp, returnflow, changec, naux, dtype, model, ext_unit_dict
                )
                spd[kper] = recarray
                recipient_nodes[kper] = recip_lists
                prev_recarray = recarray
                prev_recips = recip_lists

            # Active-parameter records (SGWF2DRT8LS) for this stress period.
            if np_sp > 0:
                if not parameters:
                    raise NotImplementedError(
                        "MfUsgDrt.load: active DRT parameters (NP>0) without "
                        "NPDRT parameter definitions are not supported."
                    )
                active_params[kper] = read_active_list_parameters(f, np_sp)

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
            recipient_nodes=recipient_nodes,
            dtype=dtype,
            options=options,
            parameters=parameters,
            mxl=mxl,
            active_params=active_params,
            extension="drt",
            unitnumber=unitnumber,
            filenames=filenames,
        )

    @staticmethod
    def _parse_header(line):
        """Parse DRT item 1: MXADRT IDRTCB NPDRT MXL [options].

        Returns (options, aux_names, ipakcb, returnflow, changec, npdrt, mxl).
        NPDRT/MXL come from item 1 directly (DRT has no separate PARAMETER line).
        """
        tokens = line.split()
        ipakcb = int(tokens[1]) if len(tokens) > 1 else 0
        npdrt = int(tokens[2]) if len(tokens) > 2 else 0
        mxl = int(tokens[3]) if len(tokens) > 3 else 0

        options = []
        aux_names = []
        returnflow = False
        changec = False
        i = 4
        while i < len(tokens):
            t = tokens[i].upper()
            if t in ("AUX", "AUXILIARY") and i + 1 < len(tokens):
                aux_names.append(tokens[i + 1])
                options.append(f"AUX {tokens[i + 1]}")
                i += 2
            elif t == "RETURNFLOW":
                returnflow = True
                options.append("RETURNFLOW")
                i += 1
            elif t == "CHANGEC":
                changec = True
                options.append("CHANGEC")
                i += 1
            elif t == "NOPRINT":
                options.append("NOPRINT")
                i += 1
            elif t == "SPREAD":
                # value recomputed on write; drop the token (and its value)
                if i + 1 < len(tokens):
                    i += 2
                else:
                    i += 1
            else:
                options.append(tokens[i])
                i += 1

        changec = changec and returnflow
        return options, aux_names, ipakcb, returnflow, changec, npdrt, mxl

    @staticmethod
    def _parse_drain_tokens(
        toks, returnflow, changec, naux, f, model=None, ext_unit_dict=None
    ):
        """Parse one drain line (+ spreading block); return (record, recipients)."""
        idx = 0
        node = int(toks[idx]) - 1
        idx += 1
        elev = float(toks[idx])
        idx += 1
        cond = float(toks[idx])
        idx += 1
        rec = [node, elev, cond]
        recips = []
        if returnflow:
            nr = int(toks[idx])
            idx += 1
            if nr != 0 or naux > 0:
                rfprop = float(toks[idx])
                idx += 1
            else:
                rfprop = 0.0
            rec.append(rfprop)
            if changec:
                rec.append(int(toks[idx]))
                idx += 1
            if nr > 0:
                recips = [nr - 1]  # inline single recipient
            elif nr < 0:
                nodes_1based = read_u1dint_list(  # spreading block
                    f, -nr, model, ext_unit_dict, package="DRT"
                )
                recips = [n - 1 for n in nodes_1based]
        for j in range(naux):
            rec.append(float(toks[idx + j]))
        return tuple(rec), recips

    @staticmethod
    def _ftype():
        return "DRT"
