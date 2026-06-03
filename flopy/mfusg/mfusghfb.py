"""
mfusghfb module.  Contains the MfUsgHfb class.
"""

import numpy as np

from ..modflow.mfhfb import ModflowHfb
from ..pakbase import Package
from ..utils.recarray_utils import create_empty_recarray
from ._usgt_list import begin_list_block
from ._usgt_parameters import (
    check_parameter_name,
    check_parval,
    read_active_list_parameters,
    read_list_parameter_header,
    write_active_list_parameters,
    write_list_parameter_header,
)


class MfUsgHfb(ModflowHfb):
    """MODFLOW-USG Hydraulic Flow Barrier (HFB6) Package.

    Extends :class:`ModflowHfb` with TRANSIENT_HFB support and proper
    registration for USG-T models.  The HFB applies only to the
    Groundwater Flow Process; it has no transport AUX variables.

    Parameters
    ----------
    model : MfUsg
        Model object.
    nhfbnp : int
        Number of non-parametric HFB barriers (dataset 4).
    hfb_data : np.recarray, optional
        Static barrier definitions.  Fields depend on grid type:
        - Unstructured: ``node1``, ``node2``, ``hydchr``
        - Structured:   ``k``, ``irow1``, ``icol1``, ``irow2``, ``icol2``, ``hydchr``
        Indices are **zero-based** (converted to 1-based on write).
    nphfb : int
        Number of HFB barriers defined by named parameters (dataset 1).
        Auto-computed as ``len(parameters)`` when omitted/0; validated if given.
    mxfb : int
        Maximum number of barrier list entries defined by parameters (``MXFBP``
        in the Fortran). Auto-computed as the total of all definition ``nlst``
        when omitted/0; validated to be ``>=`` that total if given.
    nacthfb : int
        Number of active HFB parameters (dataset 5). Auto-computed as
        ``len(acthfb_names)`` when omitted/0; validated if given.
    parameters : dict or None
        HFB list-parameter definitions, keyed by parameter name (Stage 4.6C-A
        from-scratch authoring; also the form :meth:`load` preserves). Each value
        is a dict:

        * ``parval`` (required) — the value (scales the barrier ``HYDCHR``);
        * ``data`` (required) — the parameter's barrier rows: a recarray matching
          the active dtype, or any array-like the dtype can build (zero-based
          indices, written 1-based); unstructured ``(node1, node2, hydchr)`` or
          structured ``(k, irow1, icol1, irow2, icol2, hydchr)``;
        * ``partyp`` (optional) — defaults to / must be ``"HFB"``;
        * ``nlst`` (optional) — number of rows; computed from ``data`` if
          omitted, validated if given.
    acthfb_names : list of str or None
        Names of the active parameters (dataset 6), in file order. Must be
        defined in *parameters* (case-insensitive), with no duplicates.
    options : list of str, optional
        Extra options written to the header line (e.g. ``['NOPRINT']``).
        Do **not** include ``TRANSIENT_HFB`` here; use *transient=True*.
    transient : bool
        If True, write ``TRANSIENT_HFB`` in the header and include
        per-stress-period IHFBRD + optional updated barrier data.
    stress_period_data : dict, optional
        Only used when *transient=True*.  Keys are zero-based stress
        period indices.  Values are either:
        - ``None`` / missing key → IHFBRD = -1 (reuse previous SP data)
        - np.recarray with same dtype as *hfb_data* → IHFBRD = 1 and NHFBNP
          records are written.  Indices are zero-based.
    extension : str
        File extension (default ``'hfb'``).
    unitnumber : int, optional
    filenames : str or list of str, optional

    Notes
    -----
    HFB uses MODFLOW **list** parameters (``UPARLSTRP`` / ``UPARLSTSUB`` in
    ``parutl7.f``): each parameter owns NLST barrier rows and its value scales
    the barrier ``HYDCHR`` factor. Parameterized HFB files are **preserved**
    (Stage 4.4B: load -> write -> reload keeps ``NPHFB``, the definition blocks,
    and the active-parameter records) and, as of Stage 4.6C-A, **authorable from
    scratch**: pass ``parameters={name: {parval, data, ...}}`` + ``acthfb_names``
    (the counts ``nphfb``/``mxfb``/``nacthfb`` are auto-computed). Definitions are
    normalized and validated before the file is opened (``ValueError``, no
    partial file): parameter/active names must be single whitespace-free tokens
    of <=10 chars and unique case-insensitively; ``parval`` a number or single
    token; barrier indices non-negative.

    Still not supported, failing explicitly: parameter ``INSTANCES`` (the Fortran
    aborts) -> ``NotImplementedError`` on load; and ``TRANSIENT_HFB`` combined
    with parameters (``NPHFB > 0``), which would re-read/redefine parameters each
    stress period -> ``NotImplementedError``. Non-parametric HFB (static and
    ``TRANSIENT_HFB``) is unchanged.

    Examples
    --------
    >>> import flopy, numpy as np
    >>> m = flopy.mfusg.MfUsg()
    >>> hfb = flopy.mfusg.MfUsgHfb.load("model.hfb", m, nper=12)
    """

    def __init__(
        self,
        model,
        nhfbnp=0,
        hfb_data=None,
        nacthfb=0,
        nphfb=0,
        mxfb=0,
        parameters=None,
        acthfb_names=None,
        options=None,
        transient=False,
        stress_period_data=None,
        extension="hfb",
        unitnumber=None,
        filenames=None,
        add_package=True,
    ):
        if unitnumber is None:
            unitnumber = ModflowHfb._defaultunit()
        if options is None:
            options = []

        # Build hfb_data recarray
        dtype = MfUsgHfb.get_default_dtype(structured=model.structured)
        if hfb_data is None and stress_period_data:
            first = next(
                (v for _, v in sorted(stress_period_data.items()) if v is not None),
                None,
            )
            if first is not None:
                hfb_data = first
                nhfbnp = len(first)
        if hfb_data is None:
            hfb_data = create_empty_recarray(0, dtype)
            nhfbnp = 0
        elif not isinstance(hfb_data, np.recarray):
            hfb_data = np.array(hfb_data, dtype=dtype).view(np.recarray)
        if nhfbnp == 0:
            nhfbnp = len(hfb_data)

        # Parent __init__ handles Package registration, self.hfb_data, etc.
        # add_package is intercepted here so we can control it independently.
        super().__init__(
            model,
            nphfb=nphfb,
            mxfb=mxfb,
            nhfbnp=nhfbnp,
            hfb_data=hfb_data,
            nacthfb=nacthfb,
            options=options,
            extension=extension,
            unitnumber=unitnumber,
            filenames=self._prepare_filenames(filenames),
        )

        self.transient = transient
        self.stress_period_data = stress_period_data or {}

        # Preserved HFB list-parameter definitions and active-parameter records
        # (set by load when NPHFB > 0).
        self.parameters = parameters
        self.acthfb_names = acthfb_names if acthfb_names is not None else []

        # The parent always adds itself; remove then re-add only if requested.
        if not add_package:
            model.pop_package(self._ftype())

    # ------------------------------------------------------------------
    # write_file
    # ------------------------------------------------------------------

    def write_file(self):
        """Write HFB6 package file honouring TRANSIENT_HFB when set."""
        # Parameter intent: definitions, an active list, or a declared NPHFB.
        preserve = bool(self.parameters) or bool(self.acthfb_names) or self.nphfb > 0
        params, nphfb, mxfb, nacthfb = ({}, 0, 0, 0)
        if preserve:
            if self.transient:
                raise NotImplementedError(
                    "MfUsgHfb.write_file does not support TRANSIENT_HFB combined "
                    "with HFB parameters (NPHFB > 0): under TRANSIENT_HFB the "
                    "Fortran re-reads the parameter definitions every stress "
                    "period via UPARLSTRP with ITERP=1, so the second period "
                    "would re-define the same parameter names and USG-T aborts "
                    "with 'Duplicate parameter name' (parutl7.f:604-609). Use "
                    "non-transient parameterized HFB, or NPHFB=0 for transient "
                    "barriers."
                )
            # Normalize from-scratch / preserved definitions, auto-compute the
            # counts, and validate the whole parameterized state before opening
            # the file, so a parameterized header is never written without a
            # complete, consistent body and no partial file is produced.
            params, nphfb, mxfb, nacthfb = self._validate_parameter_write()
        structured = self.parent.structured
        nper = self.parent.nper

        with open(self.fn_path, "w") as f:
            f.write(f"{self.heading}\n")

            # Dataset 1: NPHFB MXFBP NHFBNP [OPTIONS]
            header = f"{nphfb:10d}{mxfb:10d}{self.nhfbnp:10d}"
            for opt in self.options:
                header += f"  {opt}"
            if self.transient:
                header += "  TRANSIENT_HFB"
            f.write(header + "\n")

            # Parameterized HFB (non-transient): dataset 2-3 parameter
            # definitions, then dataset 4 non-parametric barriers, then
            # dataset 5-6 active-parameter records.
            if preserve:
                for name, pdef in params.items():
                    write_list_parameter_header(
                        f, name, pdef["partyp"], pdef["parval"], pdef["nlst"]
                    )
                    self._write_hfb_rows(f, pdef["data"], structured)
                self._write_hfb_rows(f, self.hfb_data, structured)
                f.write(f"{nacthfb:10d}\n")
                write_active_list_parameters(f, self.acthfb_names)
                return

            if not self.transient:
                self._write_hfb_rows(f, self.hfb_data, structured)
                return

            # Fortran reads IHFBRD at the start of each stress period.  When
            # IHFBRD > 0 it then reads NHFBNP HFB rows; IHFBRD is a flag, not a
            # row count.
            for kper in range(nper):
                sp_data = self.stress_period_data.get(kper, None)
                if sp_data is None and kper == 0 and self.nhfbnp > 0:
                    sp_data = self.hfb_data

                if sp_data is None:
                    f.write("-1\n")
                    continue

                if not isinstance(sp_data, np.recarray):
                    sp_data = np.array(sp_data, dtype=self.hfb_data.dtype).view(
                        np.recarray
                    )
                if len(sp_data) != self.nhfbnp:
                    raise ValueError(
                        "Transient HFB stress_period_data entries must contain "
                        f"NHFBNP={self.nhfbnp} rows; got {len(sp_data)} for kper {kper}."
                    )
                f.write("1\n")
                self._write_hfb_rows(f, sp_data, structured)

    def _normalize_param(self, name, pdef, structured):
        """Validate + canonicalize one HFB parameter definition for writing.

        Accepts an ergonomic from-scratch dict *or* an already-canonical loaded
        one and returns ``{partyp, parval, nlst, data (recarray)}``. ``data`` may
        be a recarray or any array-like the active barrier ``dtype`` can build;
        ``nlst`` is computed from ``len(data)`` when omitted, validated when
        given. ``partyp`` defaults to / must be ``HFB`` (the Fortran activates HFB
        params as ``PARTYP='HFB'``). Barrier indices must be non-negative 0-based
        integers. Raises ``ValueError`` for invalid input, before any file is
        opened.
        """
        check_parameter_name(name, "parameter name", prefix="MfUsgHfb.write_file")
        if not isinstance(pdef, dict):
            raise ValueError(
                f"MfUsgHfb.write_file: parameter '{name}' must be a dict, got "
                f"{type(pdef).__name__}."
            )
        partyp = pdef.get("partyp", "HFB")
        if str(partyp).upper().strip() != "HFB":
            raise ValueError(
                f"MfUsgHfb.write_file: parameter '{name}' partyp must be 'HFB'; "
                f"got {partyp!r}."
            )
        parval = check_parval(pdef.get("parval"), name, prefix="MfUsgHfb.write_file")
        data = pdef.get("data")
        if data is None or len(data) == 0:
            raise ValueError(
                f"MfUsgHfb.write_file: parameter '{name}' needs non-empty 'data' "
                "barrier rows."
            )
        dtype = MfUsgHfb.get_default_dtype(structured=structured)
        if not isinstance(data, np.recarray):
            data = np.array(data, dtype=dtype).view(np.recarray)
        idx_fields = (
            ("k", "irow1", "icol1", "irow2", "icol2")
            if structured
            else ("node1", "node2")
        )
        for fld in idx_fields:
            if np.any(np.asarray(data[fld]) < 0):
                raise ValueError(
                    f"MfUsgHfb.write_file: parameter '{name}' has a negative "
                    f"'{fld}'; barrier indices are 0-based and must be "
                    "non-negative."
                )
        nlst = pdef.get("nlst")
        if nlst is None:
            nlst = len(data)
        elif nlst != len(data):
            raise ValueError(
                f"MfUsgHfb.write_file: parameter '{name}' declares nlst={nlst} "
                f"but carries {len(data)} barrier rows."
            )
        # Preserve the original partyp text (the Fortran upper-cases it anyway),
        # so a loaded file's casing round-trips unchanged.
        return {"partyp": partyp, "parval": parval, "nlst": nlst, "data": data}

    def _validate_parameter_write(self):
        """Normalize + validate HFB parameter state for a parameterized write.

        Supports both preserved (loaded) and **from-scratch** parameter dicts
        (Stage 4.6C-A). Returns ``(params, nphfb, mxfb, nacthfb)`` with canonical
        definitions and the resolved counts (``NPHFB = len(params)``,
        ``MXFBP = sum(nlst)``, ``NACTHFB = len(acthfb_names)`` when omitted/0).
        All validation runs before the file is opened, so a ``NPHFB>0`` header is
        never written without a complete, consistent body and no partial file is
        produced. ``INSTANCES`` and ``TRANSIENT_HFB`` + parameters remain
        unsupported (handled by ``load`` / ``write_file``).
        """
        if not self.parameters:
            raise ValueError(
                "MfUsgHfb.write_file: NPHFB>0 / acthfb_names reference parameters "
                "but none are defined. Pass parameters={name: {...}} to author "
                "HFB parameter definitions from scratch."
            )
        structured = self.parent.structured
        params = {}
        total = 0
        for name, pdef in self.parameters.items():
            params[name] = self._normalize_param(name, pdef, structured)
            total += params[name]["nlst"]
        lowered_defs = [name.lower() for name in params]
        if len(set(lowered_defs)) != len(lowered_defs):
            raise ValueError(
                "MfUsgHfb.write_file: duplicate parameter definition name "
                f"(case-insensitive): {list(params)}. The Fortran upper-cases "
                "PARNAM, so definition names must be unique ignoring case."
            )
        nphfb = self.nphfb if self.nphfb else len(params)
        if nphfb != len(params):
            raise ValueError(
                f"MfUsgHfb.write_file: NPHFB ({nphfb}) must equal the number of "
                f"parameter definitions ({len(params)})."
            )
        mxfb = self.mxfb if self.mxfb else total
        if mxfb < total:
            raise ValueError(
                f"MfUsgHfb.write_file: MXFBP ({mxfb}) must be >= the total number "
                f"of parameter barrier rows ({total})."
            )
        for nm in self.acthfb_names:
            check_parameter_name(
                nm, "active parameter name", prefix="MfUsgHfb.write_file"
            )
        nacthfb = self.nacthfb if self.nacthfb else len(self.acthfb_names)
        if nacthfb != len(self.acthfb_names):
            raise ValueError(
                f"MfUsgHfb.write_file: nacthfb ({nacthfb}) must equal the number "
                f"of active parameter names ({len(self.acthfb_names)})."
            )
        lowered = [nm.lower() for nm in self.acthfb_names]
        if len(set(lowered)) != len(lowered):
            raise ValueError(
                "MfUsgHfb.write_file: a parameter is activated more than once in "
                f"the HFB active list: {self.acthfb_names}. USG-T aborts when a "
                "parameter is already activated (SGWF2HFB7SUB)."
            )
        defined = {name.lower() for name in params}
        for nm in self.acthfb_names:
            if nm.lower() not in defined:
                raise ValueError(
                    f"MfUsgHfb.write_file: active parameter '{nm}' is not "
                    "defined in parameters."
                )
        return params, nphfb, mxfb, nacthfb

    @staticmethod
    def _write_hfb_rows(f, data, structured):
        for row in data:
            if structured:
                vals = [
                    f"{int(row['k']) + 1:10d}",
                    f"{int(row['irow1']) + 1:10d}",
                    f"{int(row['icol1']) + 1:10d}",
                    f"{int(row['irow2']) + 1:10d}",
                    f"{int(row['icol2']) + 1:10d}",
                    f"{float(row['hydchr']):13.6g}",
                ]
            else:
                vals = [
                    f"{int(row['node1']) + 1:10d}",
                    f"{int(row['node2']) + 1:10d}",
                    f"{float(row['hydchr']):13.6g}",
                ]
            f.write("".join(vals) + "\n")

    @staticmethod
    def _read_hfb_rows(f, nrows, dtype, structured, model=None, ext_unit_dict=None):
        data = create_empty_recarray(nrows, dtype)
        if nrows == 0:
            return data
        # A barrier list may begin with MODFLOW list controls (SFAC, OPEN/CLOSE,
        # EXTERNAL) per SGWF2HFB7RL/RLU. begin_list_block consumes them and
        # returns the source to read the rows from, the SFAC scale, and the first
        # already-read row. SFAC scales the barrier FACTOR/HYDCHR (HFB(6,II) =
        # FACTOR*SFAC in the Fortran).
        source, sfac, line, to_close = begin_list_block(
            f, model=model, ext_unit_dict=ext_unit_dict, package="HFB"
        )
        for i in range(nrows):
            if i > 0:
                line = source.readline()
            t = line.split()
            if structured:
                data[i] = (
                    int(t[0]) - 1,
                    int(t[1]) - 1,
                    int(t[2]) - 1,
                    int(t[3]) - 1,
                    int(t[4]) - 1,
                    float(t[5]) * sfac,
                )
            else:
                data[i] = (int(t[0]) - 1, int(t[1]) - 1, float(t[2]) * sfac)
        if to_close is not None:
            to_close.close()
        return data

    # ------------------------------------------------------------------
    # load
    # ------------------------------------------------------------------

    @classmethod
    def load(cls, f, model, nper=None, ext_unit_dict=None):
        """Load an HFB6 package from a USG-T file.

        Parameters
        ----------
        f : str or file-like
        model : MfUsg
        nper : int, optional
        ext_unit_dict : dict, optional

        Returns
        -------
        MfUsgHfb
        """
        if model.verbose:
            print("loading mfusg hfb6 package file...")

        if nper is None:
            nper = model.nper
        structured = model.structured

        openfile = not hasattr(f, "read")
        if openfile:
            f = open(f, "r")

        # Skip comments
        line = f.readline()
        while line.startswith("#"):
            line = f.readline()

        # Dataset 1: NPHFB MXFB NHFBNP [OPTIONS]
        tokens = line.split()
        nphfb = int(tokens[0])
        mxfb = int(tokens[1])
        nhfbnp = int(tokens[2])
        options = []
        transient = False
        for tok in tokens[3:]:
            if tok.upper() == "TRANSIENT_HFB":
                transient = True
            elif tok.upper() not in ("OPTIONS",):
                options.append(tok)

        dtype = cls.get_default_dtype(structured=structured)
        stress_period_data = {}
        parameters = None
        acthfb_names = []
        nacthfb = 0

        if nphfb > 0:
            if transient:
                raise NotImplementedError(
                    "MfUsgHfb.load does not support TRANSIENT_HFB combined with "
                    "HFB parameters (NPHFB > 0): the Fortran re-reads/redefines "
                    "the parameter definitions every stress period (UPARLSTRP "
                    "ITERP=1). This combination is out of scope."
                )
            # Dataset 2-3: NPHFB list-parameter definitions (header + NLST rows).
            parameters = {}
            for _ in range(nphfb):
                name, partyp, parval, nlst, numinst = read_list_parameter_header(
                    f.readline()
                )
                if numinst > 0:
                    raise NotImplementedError(
                        "MfUsgHfb.load: HFB parameter INSTANCES are not "
                        "supported -- USG-T itself rejects them: gwf2hfb7u1.f:"
                        "135-137 writes 'INSTANCES ARE NOT SUPPORTED FOR HFB' "
                        "and calls USTOP when NUMINST>0."
                    )
                data = cls._read_hfb_rows(
                    f, nlst, dtype, structured, model, ext_unit_dict
                )
                parameters[name] = {
                    "partyp": partyp,
                    "parval": parval,
                    "nlst": nlst,
                    "data": data,
                }
            # Dataset 4: barriers not defined by parameters.
            hfb_data = cls._read_hfb_rows(
                f, nhfbnp, dtype, structured, model, ext_unit_dict
            )
            # Dataset 5-6: number of active parameters and their names.
            nacthfb = int(f.readline().split()[0])
            acthfb_names = read_active_list_parameters(f, nacthfb)
        elif transient:
            hfb_data = create_empty_recarray(0, dtype)
            for kper in range(nper):
                line = f.readline()
                if not line:
                    break
                ihfbrd = int(line.split()[0])
                if ihfbrd > 0:
                    sp_arr = cls._read_hfb_rows(
                        f, nhfbnp, dtype, structured, model, ext_unit_dict
                    )
                    stress_period_data[kper] = sp_arr
                    if len(hfb_data) == 0:
                        hfb_data = sp_arr.copy()
        else:
            hfb_data = cls._read_hfb_rows(
                f, nhfbnp, dtype, structured, model, ext_unit_dict
            )

        if openfile:
            f.close()

        unitnumber = None
        filenames = [None]
        if ext_unit_dict is not None:
            unitnumber, filenames[0] = model.get_ext_dict_attr(
                ext_unit_dict, filetype=cls._ftype()
            )

        return cls(
            model,
            nphfb=nphfb,
            mxfb=mxfb,
            nhfbnp=nhfbnp,
            hfb_data=hfb_data,
            nacthfb=nacthfb,
            parameters=parameters,
            acthfb_names=acthfb_names,
            options=options,
            transient=transient,
            stress_period_data=stress_period_data,
            extension="hfb",
            unitnumber=unitnumber,
            filenames=filenames,
        )
