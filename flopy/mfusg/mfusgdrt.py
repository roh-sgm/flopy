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

Notes
-----
The USG-T transport extensions (RETURNFLOW recipient nodes, ``CHANGEC`` /
``IDCHNGTYP``, ``SPREAD``) are implemented for **unstructured** grids only.
Constructing ``MfUsgDrt`` on a structured (DIS) model raises
``NotImplementedError`` — use :class:`flopy.modflow.ModflowDrt` for classic
structured DRT. (On load, ``MfUsgDrt.load`` delegates structured files to
``ModflowDrt.load``, which returns a base-class object.)

Not supported (explicit failure rather than partial write):

* Named parameters (``NPDRT > 0``).
* ``EXTERNAL`` / ``OPEN/CLOSE`` spreading-node lists.
"""

import numpy as np

from ..modflow.mfdrt import ModflowDrt
from ..pakbase import Package
from ._usgt_list import begin_list_block
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
        return list(self.dtype.names[self._base_field_count():])

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

    def write_file(self):
        """Write the package file in MODFLOW-USG-T DRT8 format.

        Unstructured only; ``__init__`` rejects structured models, so there is
        no structured-delegation path here.
        """
        nper = self.parent.nper
        aux_names = self._aux_field_names()
        has_aux = len(aux_names) > 0

        mxadrt = max((len(v) for v in self.stress_period_data.values()), default=0)
        mxspread = self._max_spread_nodes()

        with open(self.fn_path, "w") as f:
            f.write(f"{self.heading}\n")

            # Item 1: MXADRT IDRTCB NPDRT MXL [options]
            line = f" {mxadrt:9d} {self.ipakcb} 0 0"
            for opt in self.options:
                line += f" {opt}"
            if mxspread > 0:
                line += f" SPREAD {mxspread}"
            f.write(line + "\n")

            for kper in range(nper):
                if kper not in self.stress_period_data:
                    f.write(f" -1    Stress Period {kper + 1}\n")
                    continue
                recarray = self.stress_period_data[kper]
                recips = self._validated_recipients(kper, len(recarray))
                f.write(f" {len(recarray)} 0    Stress Period {kper + 1}\n")
                for i, rec in enumerate(recarray):
                    self._write_drain_line(f, rec, recips[i], has_aux, aux_names)

    def _max_spread_nodes(self):
        """Max total spreading-recipient nodes in any stress period."""
        mx = 0
        for kper, recarray in self.stress_period_data.items():
            recips = self.recipient_nodes.get(kper, [])
            total = sum(len(r) for r in recips if len(r) > 1)
            mx = max(mx, total)
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

        options, aux_names, ipakcb, returnflow, changec = cls._parse_header(line)
        dtype = cls.get_usg_dtype(
            returnflow=returnflow, changec=changec, aux_names=aux_names
        )
        naux = len(aux_names)

        spd = {}
        recipient_nodes = {}
        prev_recarray = None
        prev_recips = None

        for kper in range(nper):
            line = f.readline()
            if not line:
                break
            itmp = int(line.split()[0])

            if itmp < 0:
                if prev_recarray is not None:
                    spd[kper] = prev_recarray.copy()
                    recipient_nodes[kper] = [list(r) for r in prev_recips]
                continue

            # Honor leading SFAC / EXTERNAL / OPEN-CLOSE list controls; drain
            # rows and their interleaved spreading U1DINT blocks are read from
            # the (possibly redirected) source. SFAC scales COND (Fortran
            # ISCLOC=5).
            source, sfac, first_line, to_close = (f, 1.0, None, None)
            if itmp > 0:
                source, sfac, first_line, to_close = begin_list_block(
                    f, model, ext_unit_dict, package="DRT"
                )

            records = []
            recip_lists = []
            for idx in range(itmp):
                row = first_line if idx == 0 else source.readline()
                rec, recips = cls._parse_drain_tokens(
                    row.split(), returnflow, changec, naux, source
                )
                records.append(rec)
                recip_lists.append(recips)

            if to_close is not None:
                source.close()

            recarray = np.array(records, dtype=dtype).view(np.recarray)
            if sfac != 1.0 and len(recarray) > 0:
                recarray["cond"] = recarray["cond"] * sfac
            spd[kper] = recarray
            recipient_nodes[kper] = recip_lists
            prev_recarray = recarray
            prev_recips = recip_lists

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
            extension="drt",
            unitnumber=unitnumber,
            filenames=filenames,
        )

    @staticmethod
    def _parse_header(line):
        """Parse DRT item 1: MXADRT IDRTCB NPDRT MXL [options]."""
        tokens = line.split()
        npdrt = int(tokens[2]) if len(tokens) > 2 else 0
        if npdrt > 0:
            raise NotImplementedError(
                "MfUsgDrt does not support named DRT parameters (NPDRT > 0)."
            )
        ipakcb = int(tokens[1]) if len(tokens) > 1 else 0

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
        return options, aux_names, ipakcb, returnflow, changec

    @staticmethod
    def _parse_drain_tokens(toks, returnflow, changec, naux, f):
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
                nodes_1based = read_u1dint_list(f, -nr)  # spreading block
                recips = [n - 1 for n in nodes_1based]
        for j in range(naux):
            rec.append(float(toks[idx + j]))
        return tuple(rec), recips

    @staticmethod
    def _ftype():
        return "DRT"
