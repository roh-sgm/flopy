"""
mfusgqrt module. Contains the MfUsgQrt class.

Sink with Return Flow (QRT) Package -- USG-Transport.

Reference: Fortran subroutines ``GWF2QRT8U1AR`` / ``GWF2QRT8U1RP`` /
``SGWF2QRT8LR`` in ``gwf2QRT8u.f`` of USG-Transport 2.7. QRT8 extracts water
at sink nodes and optionally returns a proportion of it to one or more
recipient nodes (analogous to DRT, but for sinks).

File layout (unstructured, free format)::

    Item 1:  MXAQRT MXRTCELLS IQRTCB NPQRT MXL [options]
             options: [AUX <name> ...] [RETURNFLOW] [AUTOFLOWREDUCE]
                      [CHANGEC] [NOPRINT] [IUNIT_AFR_QRT <unit>]
    Per SP:  ITMP
             ITMP sink lines, each:
                NODE  Q  [NumRT  Rfprop]  [IQCHNGTYP]  [aux ...]
             then, for each sink with NumRT > 0 (in sink order), a U1DINT
             block of NumRT recipient node ids.

``ITMP < 0`` reuses the previous stress period's sinks and recipients.

Internal ``node`` and recipient values are 0-based; the file is 1-based.

Named parameters (``NPQRT > 0``) are **structurally preserved** (load -> write ->
reload) as of Stage 4.4E, but **not execution-guaranteed**. QRT is type-consistent
(definitions and activations both use ``PARTYP='QRT'``: ``UPARLSTRP`` at
``gwf2QRT8u.f:183`` and ``SGWF2QRT8LS`` at ``gwf2QRT8u.f:1118``), so an active QRT
parameter is type-valid. However, two USG-T 2.7 Fortran issues mean FloPy can only
promise faithful round-trip, not guaranteed execution, for activated QRT
parameters:

* **the parameter value scales the wrong field.** ``SGWF2QRT8LS`` uses
  ``IPVL1=IPVL2=5`` (``gwf2QRT8u.f:1119-1120``), i.e. it scales ``QRTF(5)=NumRT``
  (the *recipient count*), not ``QRTF(4)=Q``. ``SFAC`` correctly scales ``Q``
  (``ISCLOC=4``), so this is a Fortran bug: with ``PARVAL != 1.0`` the recipient
  count is corrupted. FloPy does **not** apply the parameter value to ``Q``.
* **recipient nodes are not copied on activation.** ``SGWF2QRT8LS`` copies
  ``QRTF`` but not ``NodQRT``, so an activated parameter's RETURNFLOW recipients
  are not resolved at run time (same limitation as DRT).

The ``NPQRT``/``MXL`` count is in item 1 (no separate ``PARAMETER`` line); each
parameter owns ``NLST`` sink rows with their recipient ``U1DINT`` blocks (read
after the rows, in sink order). Per stress period, ``ITMP NP`` gives the
non-parametric count and the number of active parameters.

Not supported (explicit failure rather than partial write):

* Parameter ``INSTANCES`` (``NUMINST>0``): the Fortran supports them, FloPy does
  not yet -> load raises ``NotImplementedError``.
* From-scratch parameter authoring (active parameters with no loaded
  definitions) -> ``NotImplementedError``.
* The ``TRANSIENTQ`` transient-flow time-series option.
* ``EXTERNAL`` / ``OPEN/CLOSE`` recipient-node lists (see
  :mod:`flopy.mfusg._usgt_returnflow`).
"""

import numpy as np

from ..pakbase import Package
from ._usgt_list import begin_list_block
from ._usgt_parameters import (
    read_active_list_parameters,
    read_list_parameter_header,
    resolve_list_parameter,
    write_active_list_parameters,
    write_list_parameter_header,
)
from ._usgt_returnflow import read_u1dint_list, write_u1dint_list
from .mfusg import MfUsg


class MfUsgQrt(Package):
    """MODFLOW-USG Sink with Return Flow (QRT) Package.

    Parameters
    ----------
    model : MfUsg
        The model object to which this package will be added.
    ipakcb : int, optional
        Unit number for cell-by-cell budget output (``IQRTCB``).
    stress_period_data : dict, optional
        Dictionary keyed by zero-based stress period. Each value is a recarray
        matching ``dtype``. Base fields are ``node`` (0-based) and ``q``; with
        the ``RETURNFLOW`` option, ``rfprop`` (return proportion); with
        ``CHANGEC``, ``iqchngtyp``; plus any AUX fields.
    recipient_nodes : dict, optional
        Dictionary keyed by zero-based stress period. Each value is a list
        (one entry per sink record, in the same order as ``stress_period_data``)
        of 0-based recipient node lists. Only used with ``RETURNFLOW``.
    dtype : np.dtype, optional
        Custom dtype. If None, derived from the active options.
    options : list of str, optional
        Package options, e.g. ``["RETURNFLOW", "CHANGEC", "AUX C01"]``.
    extension : str
        Filename extension (default ``"qrt"``).
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
        extension="qrt",
        unitnumber=None,
        filenames=None,
        add_package=True,
    ):
        assert isinstance(model, MfUsg), (
            "Model object must be of type flopy.mfusg.MfUsg\n"
            f"but received type: {type(model)}."
        )

        if unitnumber is None:
            unitnumber = MfUsgQrt._defaultunit()
        if options is None:
            options = []

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
        self.url = "qrt.html"
        self.options = options

        # Derive feature flags from the options block.
        self.returnflow = any(o.split()[0].upper() == "RETURNFLOW" for o in options)
        # CHANGEC (return-flow solute change type) is only read with RETURNFLOW.
        self.changec = self.returnflow and any(
            o.split()[0].upper() == "CHANGEC" for o in options
        )

        aux_names = self._aux_names_from_options(options)
        if dtype is not None:
            self.dtype = dtype
        else:
            self.dtype = self.get_default_dtype(
                returnflow=self.returnflow, changec=self.changec, aux_names=aux_names
            )

        self.stress_period_data = {
            kper: self._to_recarray(val)
            for kper, val in (stress_period_data or {}).items()
        }
        self.recipient_nodes = recipient_nodes or {}

        # Preserved QRT list-parameter state (set by load when NPQRT > 0):
        #   parameters: {name: {"partyp","parval","nlst","data","recipient_nodes"}}
        #     (rows 0-based; recipient_nodes is one 0-based list per definition row)
        #   mxl:        MXL from item 1 (max parameter list entries)
        #   active_params: {kper: [name, ...]} active parameters per stress period
        # See module docstring: structurally preserved, not execution-guaranteed.
        self.parameters = parameters
        self.mxl = mxl
        self.active_params = active_params if active_params is not None else {}

        if add_package:
            self.parent.add_package(self)

    # ------------------------------------------------------------------
    # Static / helper constructors
    # ------------------------------------------------------------------

    @staticmethod
    def get_default_dtype(returnflow=True, changec=False, aux_names=None):
        """Build the QRT dtype for the active options."""
        cols = [("node", int), ("q", np.float64)]
        if returnflow:
            cols.append(("rfprop", np.float64))
            if changec:
                cols.append(("iqchngtyp", int))
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

    def _base_field_count(self):
        n = 2  # node, q
        if self.returnflow:
            n += 1  # rfprop
            if self.changec:
                n += 1  # iqchngtyp
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
        """Return the recipient lists for a stress period, validated against it.

        With RETURNFLOW active there must be exactly one recipient list per
        stress record. Omitting ``recipient_nodes`` for a period is allowed and
        means every record has zero recipients. Mismatched lengths raise so
        return-flow metadata cannot be silently dropped or shifted on write.
        """
        if not self.returnflow:
            return [[] for _ in range(nrec)]
        recips = self.recipient_nodes.get(kper)
        if recips is None:
            return [[] for _ in range(nrec)]
        if len(recips) != nrec:
            raise ValueError(
                f"MfUsgQrt: recipient_nodes[{kper}] has {len(recips)} entries "
                f"but stress_period_data[{kper}] has {nrec} records; provide "
                "exactly one recipient list per record (or omit the period for "
                "all-zero recipients)."
            )
        return recips

    def _validate_parameter_write(self):
        """Validate preserved QRT parameter state before a parameterized write.

        Raises ``NotImplementedError`` for from-scratch authoring (active
        parameters with no loaded definitions) and ``ValueError`` for
        inconsistent definitions, so a parameterized item 1 is never written
        without a complete, consistent body. Validation runs before the file is
        opened, so no partial file is produced.
        """
        if not self.parameters:
            raise NotImplementedError(
                "MfUsgQrt.write_file cannot author QRT parameter definitions "
                "from scratch (active parameters without loaded definitions). "
                "Parameter preservation is supported for files read by "
                "MfUsgQrt.load; for from-scratch input use no parameters."
            )
        total = 0
        for name, pdef in self.parameters.items():
            missing = [
                k
                for k in ("partyp", "parval", "nlst", "data", "recipient_nodes")
                if k not in pdef
            ]
            if missing:
                raise ValueError(
                    f"MfUsgQrt.write_file: parameter '{name}' is missing keys "
                    f"{missing} (needs partyp, parval, nlst, data, "
                    "recipient_nodes)."
                )
            if len(pdef["data"]) != pdef["nlst"]:
                raise ValueError(
                    f"MfUsgQrt.write_file: parameter '{name}' declares nlst="
                    f"{pdef['nlst']} but carries {len(pdef['data'])} rows."
                )
            if len(pdef["recipient_nodes"]) != pdef["nlst"]:
                raise ValueError(
                    f"MfUsgQrt.write_file: parameter '{name}' has "
                    f"{len(pdef['recipient_nodes'])} recipient lists but nlst="
                    f"{pdef['nlst']} (need one per definition row)."
                )
            total += pdef["nlst"]
        if self.mxl <= 0:
            raise ValueError(
                "MfUsgQrt.write_file: MXL (item 1) must be > 0 when QRT "
                f"parameter definitions are present; got mxl={self.mxl}."
            )
        if self.mxl < total:
            raise ValueError(
                f"MfUsgQrt.write_file: MXL ({self.mxl}) must be >= the total "
                f"number of parameter list entries ({total})."
            )
        defined = {name.lower() for name in self.parameters}
        for kper, names in self.active_params.items():
            lowered = [nm.lower() for nm in names]
            if len(set(lowered)) != len(lowered):
                raise ValueError(
                    f"MfUsgQrt.write_file: a parameter is activated more than "
                    f"once in stress period {kper}: {names}. USG-T aborts when a "
                    "parameter is already active this stress period."
                )
            for nm in names:
                if nm.lower() not in defined:
                    raise ValueError(
                        f"MfUsgQrt.write_file: active parameter '{nm}' (stress "
                        f"period {kper}) is not defined in parameters."
                    )

    def _max_active_sinks(self):
        """Max active sink-return cells in any stress period (the Fortran's
        NQRTCL, which must not exceed MXAQRT).

        Each active parameter contributes its NLST rows on top of the
        non-parametric sinks (gwf2QRT8u.f: SGWF2QRT8LS does
        ``NQRTCL = NQRTCL + NLST`` and aborts if ``NQRTCL > MXAQRT``). ``ITMP<0``
        reuse carries the previous period's non-parametric count.
        """
        mx = 0
        prev_nonparam = 0
        for kper in range(self.parent.nper):
            if kper in self.stress_period_data:
                prev_nonparam = len(self.stress_period_data[kper])
            nonparam = prev_nonparam
            active = 0
            for name in self.active_params.get(kper, []):
                pdef = resolve_list_parameter(self.parameters, name)
                if pdef is not None:
                    active += pdef["nlst"]
            mx = max(mx, nonparam + active)
        return mx

    def _max_rt_cells(self):
        """Max total recipient nodes (NodQRT): the max non-parametric total per
        stress period, and the parameter definitions (read into NodQRT in AR)."""
        mx = 0
        for kper, recs in self.recipient_nodes.items():
            if kper in self.stress_period_data:
                mx = max(mx, sum(len(r) for r in recs))
        if self.parameters:
            def_total = sum(
                len(r)
                for pdef in self.parameters.values()
                for r in pdef.get("recipient_nodes", [])
            )
            mx = max(mx, def_total)
        return mx

    def write_file(self):
        """Write the package file in MODFLOW-USG-T QRT format."""
        preserve = bool(self.parameters) or any(self.active_params.values())
        if preserve:
            self._validate_parameter_write()
        nper = self.parent.nper
        aux_names = self._aux_field_names()
        has_aux = len(aux_names) > 0

        mxaqrt = self._max_active_sinks()
        mxrtcells = self._max_rt_cells()
        npqrt = len(self.parameters) if preserve else 0
        mxl = self.mxl if preserve else 0

        with open(self.fn_path, "w") as f:
            f.write(f"{self.heading}\n")

            # Item 1: MXAQRT MXRTCELLS IQRTCB NPQRT MXL [options]
            line = f" {mxaqrt:9d} {mxrtcells:9d} {self.ipakcb} {npqrt} {mxl}"
            for opt in self.options:
                line += f" {opt}"
            f.write(line + "\n")

            # Items 2-3: parameter definitions (UPARLSTRP header + NLST sink rows
            # + their recipient U1DINT blocks).
            if preserve:
                for name, pdef in self.parameters.items():
                    write_list_parameter_header(
                        f, name, pdef["partyp"], pdef["parval"], pdef["nlst"]
                    )
                    self._write_sink_block(
                        f, pdef["data"], pdef["recipient_nodes"], has_aux, aux_names
                    )

            for kper in range(nper):
                active = self.active_params.get(kper, []) if preserve else []
                if kper not in self.stress_period_data:
                    f.write(f" -1 {len(active)}    Stress Period {kper + 1}\n")
                    write_active_list_parameters(f, active)
                    continue
                recarray = self.stress_period_data[kper]
                recips = self._validated_recipients(kper, len(recarray))
                f.write(f" {len(recarray)} {len(active)}    Stress Period {kper + 1}\n")
                self._write_sink_block(f, recarray, recips, has_aux, aux_names)
                # Active-parameter records (SGWF2QRT8LS) follow the non-param rows.
                write_active_list_parameters(f, active)

    def _write_sink_block(self, f, recarray, recips, has_aux, aux_names):
        """Write the sink rows then their recipient U1DINT blocks (sink order)."""
        for i, rec in enumerate(recarray):
            self._write_sink_line(f, rec, len(recips[i]), has_aux, aux_names)
        if self.returnflow:
            for i in range(len(recarray)):
                if len(recips[i]) > 0:
                    write_u1dint_list(f, [int(n) + 1 for n in recips[i]])

    def _write_sink_line(self, f, rec, numrt, has_aux, aux_names):
        line = f" {int(rec['node']) + 1}  {float(rec['q']):.6e}"
        if self.returnflow:
            line += f"  {numrt}"
            if numrt > 0 or has_aux:
                line += f"  {float(rec['rfprop']):.6e}"
            if self.changec:
                line += f"  {int(rec['iqchngtyp'])}"
        for name in aux_names:
            line += f" {float(rec[name]):.6e}"
        f.write(line + "\n")

    # ------------------------------------------------------------------
    # load
    # ------------------------------------------------------------------

    @classmethod
    def load(cls, f, model, nper=None, ext_unit_dict=None, check=True):
        """Load a MODFLOW-USG-T QRT package from file."""
        if model.verbose:
            print("loading mfusg qrt package file...")
        if nper is None:
            nper = model.nper

        openfile = not hasattr(f, "read")
        if openfile:
            f = open(f, "r")

        line = f.readline()
        while line.startswith("#"):
            line = f.readline()

        options, aux_names, ipakcb, returnflow, changec, npqrt, mxl = cls._parse_header(
            line
        )
        dtype = cls.get_default_dtype(
            returnflow=returnflow, changec=changec, aux_names=aux_names
        )
        naux = len(aux_names)

        # Items 2-3: NPQRT parameter definitions (UPARLSTRP header + NLST rows +
        # recipient blocks). Structurally preserved -- see module docstring.
        parameters = None
        if npqrt > 0:
            parameters = {}
            for _ in range(npqrt):
                name, partyp, parval, nlst, numinst = read_list_parameter_header(
                    f.readline()
                )
                if numinst > 0:
                    raise NotImplementedError(
                        "MfUsgQrt.load: QRT parameter INSTANCES (NUMINST>0) are "
                        "supported by the Fortran but not yet by FloPy (instances "
                        "combined with per-row recipient lists)."
                    )
                data, recips = cls._read_sink_rows(
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
            # NPQRT==0 it is just "ITMP" (parts[1:] may be an inline comment).
            np_sp = int(parts[1]) if (npqrt > 0 and len(parts) > 1) else 0

            if itmp < 0:
                if prev_recarray is not None:
                    spd[kper] = prev_recarray.copy()
                    recipient_nodes[kper] = [list(r) for r in prev_recips]
            else:
                recarray, recip_lists = cls._read_sink_rows(
                    f, itmp, returnflow, changec, naux, dtype, model, ext_unit_dict
                )
                spd[kper] = recarray
                recipient_nodes[kper] = recip_lists
                prev_recarray = recarray
                prev_recips = recip_lists

            # Active-parameter records (SGWF2QRT8LS) for this stress period.
            if np_sp > 0:
                if not parameters:
                    raise NotImplementedError(
                        "MfUsgQrt.load: active QRT parameters (NP>0) without "
                        "NPQRT parameter definitions are not supported."
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
            extension="qrt",
            unitnumber=unitnumber,
            filenames=filenames,
        )

    @classmethod
    def _read_sink_rows(
        cls, f, count, returnflow, changec, naux, dtype, model, ext_unit_dict
    ):
        """Read ``count`` sink rows + their trailing recipient U1DINT blocks into
        a recarray and a per-row recipient list. Honors leading SFAC / EXTERNAL /
        OPEN-CLOSE list controls; SFAC scales Q (Fortran ISCLOC=4)."""
        if count == 0:
            return np.array([], dtype=dtype).view(np.recarray), []
        source, sfac, first_line, to_close = begin_list_block(
            f, model, ext_unit_dict, package="QRT"
        )
        records = []
        numrt_list = []
        for idx in range(count):
            row = first_line if idx == 0 else source.readline()
            rec, numrt = cls._parse_sink_tokens(row.split(), returnflow, changec, naux)
            records.append(rec)
            numrt_list.append(numrt)
        recip_lists = []
        for i in range(count):
            if returnflow and numrt_list[i] > 0:
                nodes_1based = read_u1dint_list(source, numrt_list[i])
                recip_lists.append([n - 1 for n in nodes_1based])
            else:
                recip_lists.append([])
        if to_close is not None:
            source.close()
        recarray = np.array(records, dtype=dtype).view(np.recarray)
        if sfac != 1.0 and len(recarray) > 0:
            recarray["q"] = recarray["q"] * sfac
        return recarray, recip_lists

    @staticmethod
    def _parse_header(line):
        """Parse QRT item 1: MXAQRT MXRTCELLS IQRTCB NPQRT MXL [options].

        Returns (options, aux_names, ipakcb, returnflow, changec, npqrt, mxl).
        NPQRT/MXL come from item 1 directly (QRT has no separate PARAMETER line).
        """
        tokens = line.split()
        ipakcb = int(tokens[2]) if len(tokens) > 2 else 0
        npqrt = int(tokens[3]) if len(tokens) > 3 else 0
        mxl = int(tokens[4]) if len(tokens) > 4 else 0

        options = []
        aux_names = []
        returnflow = False
        changec = False
        i = 5
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
            elif t == "AUTOFLOWREDUCE":
                options.append("AUTOFLOWREDUCE")
                i += 1
            elif t == "NOPRINT":
                options.append("NOPRINT")
                i += 1
            elif t == "IUNIT_AFR_QRT" and i + 1 < len(tokens):
                options.append(f"IUNIT_AFR_QRT {tokens[i + 1]}")
                i += 2
            elif t == "TRANSIENTQ":
                raise NotImplementedError(
                    "MfUsgQrt does not support the TRANSIENTQ transient-flow "
                    "time-series option."
                )
            else:
                options.append(tokens[i])
                i += 1

        changec = changec and returnflow
        return options, aux_names, ipakcb, returnflow, changec, npqrt, mxl

    @staticmethod
    def _parse_sink_tokens(toks, returnflow, changec, naux):
        """Parse one sink data line; return (record tuple, NumRT)."""
        idx = 0
        node = int(toks[idx]) - 1
        idx += 1
        q = float(toks[idx])
        idx += 1
        rec = [node, q]
        numrt = 0
        if returnflow:
            numrt = int(toks[idx])
            idx += 1
            if numrt > 0 or naux > 0:
                rfprop = float(toks[idx])
                idx += 1
            else:
                rfprop = 0.0
            rec.append(rfprop)
            if changec:
                rec.append(int(toks[idx]))
                idx += 1
        for j in range(naux):
            rec.append(float(toks[idx + j]))
        return tuple(rec), numrt

    @staticmethod
    def _ftype():
        return "QRT"

    @staticmethod
    def _defaultunit():
        return 160
