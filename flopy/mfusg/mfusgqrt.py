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

Not supported in this version (explicit failure rather than partial write):

* Named parameters (``NPQRT > 0``).
* The ``TRANSIENTQ`` transient-flow time-series option.
* ``EXTERNAL`` / ``OPEN/CLOSE`` recipient-node lists (see
  :mod:`flopy.mfusg._usgt_returnflow`).
"""

import numpy as np

from ..pakbase import Package
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

    def write_file(self):
        """Write the package file in MODFLOW-USG-T QRT format."""
        nper = self.parent.nper
        aux_names = self._aux_field_names()
        has_aux = len(aux_names) > 0

        mxaqrt = max((len(v) for v in self.stress_period_data.values()), default=0)
        mxrtcells = 0
        for kper, recs in self.recipient_nodes.items():
            if kper in self.stress_period_data:
                mxrtcells = max(mxrtcells, sum(len(r) for r in recs))

        with open(self.fn_path, "w") as f:
            f.write(f"{self.heading}\n")

            # Item 1: MXAQRT MXRTCELLS IQRTCB NPQRT MXL [options]
            line = f" {mxaqrt:9d} {mxrtcells:9d} {self.ipakcb} 0 0"
            for opt in self.options:
                line += f" {opt}"
            f.write(line + "\n")

            for kper in range(nper):
                if kper not in self.stress_period_data:
                    f.write(f" -1    Stress Period {kper + 1}\n")
                    continue
                recarray = self.stress_period_data[kper]
                recips = self.recipient_nodes.get(kper, [[]] * len(recarray))
                f.write(f" {len(recarray)} 0    Stress Period {kper + 1}\n")
                for i, rec in enumerate(recarray):
                    rnodes = recips[i] if i < len(recips) else []
                    self._write_sink_line(f, rec, len(rnodes), has_aux, aux_names)
                # Recipient-node U1DINT blocks, in sink order, skipping NumRT==0
                if self.returnflow:
                    for i in range(len(recarray)):
                        rnodes = recips[i] if i < len(recips) else []
                        if len(rnodes) > 0:
                            write_u1dint_list(f, [int(n) + 1 for n in rnodes])

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

        options, aux_names, ipakcb, returnflow, changec = cls._parse_header(line)
        dtype = cls.get_default_dtype(
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

            records = []
            numrt_list = []
            for _ in range(itmp):
                toks = f.readline().split()
                rec, numrt = cls._parse_sink_tokens(
                    toks, returnflow, changec, naux
                )
                records.append(rec)
                numrt_list.append(numrt)

            recip_lists = []
            for i in range(itmp):
                if returnflow and numrt_list[i] > 0:
                    nodes_1based = read_u1dint_list(f, numrt_list[i])
                    recip_lists.append([n - 1 for n in nodes_1based])
                else:
                    recip_lists.append([])

            recarray = np.array(records, dtype=dtype).view(np.recarray)
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
            extension="qrt",
            unitnumber=unitnumber,
            filenames=filenames,
        )

    @staticmethod
    def _parse_header(line):
        """Parse QRT item 1: MXAQRT MXRTCELLS IQRTCB NPQRT MXL [options]."""
        tokens = line.split()
        npqrt = int(tokens[3]) if len(tokens) > 3 else 0
        if npqrt > 0:
            raise NotImplementedError(
                "MfUsgQrt does not support named QRT parameters (NPQRT > 0)."
            )
        ipakcb = int(tokens[2]) if len(tokens) > 2 else 0

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
        return options, aux_names, ipakcb, returnflow, changec

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
