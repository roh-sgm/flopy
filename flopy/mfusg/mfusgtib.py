"""mfusgtib module.  Contains the MfUsgTib class. Note that the user can access
the MfUsgTib class as `flopy.mfusg.MfUsgTib`.

Transient Ibound (TIB) Package — USG-Transport.

Reference: Fortran subroutine ``GWF2TIB1RP`` in ``glo2basu1.f`` of USG-Transport.
Per-stress-period free-format layout:

    NIB0 NIB1 NIBM1 [NICB0 NICB1 NICBM1]   <- 3 ints, or 6 with transport (BCT)
    [NIB0  node ids via U1DINT       -> IBOUND := 0   (inactivate flow)]
    [NIB1  lines: ICELL [HEAD v|AVHEAD]  -> IBOUND := 1   (activate flow)]
    [NIBM1 lines: ICELL [HEAD v|AVHEAD]  -> IBOUND := -1  (prescribed head)]
    [NICB0 node ids via U1DINT       -> ICBUND := 0   (inactivate transport)]
    [NICB1 lines: ICELL [CONC c..|AVCONC] -> ICBUND := 1   (activate transport)]
    [NICBM1 lines: ICELL [CONC c..|AVCONC] -> ICBUND := -1 (prescribed conc)]

A count ``<= 0`` skips that block for the stress period (the corresponding
IBOUND/ICBUND values simply persist from the previous period; TIB does not
re-zero them). The header is 6 integers whenever the model has an active BCT
(transport) package and 3 integers otherwise — this mirrors ``INBCT`` in the
Fortran and is independent of whether any transport cells are listed.

This class supports two modes:

* **Semantic** — build from Python/numpy with ``stress_period_data`` (the
  primary authoring workflow; no existing ``.tib`` required) or recover it from
  a file with ``MfUsgTib.load(..., parse=True)``. Internal node numbers are
  0-based; the writer emits 1-based file ids.
* **Raw round-trip** — ``MfUsgTib.load`` defaults to ``parse=False`` and keeps
  the file body verbatim (``raw_body``) for byte-exact rewrite. This is the safe
  fallback for files using constructs the semantic parser does not model (e.g.
  ``EXTERNAL``/``OPEN-CLOSE`` U1DINT node lists or inline annotations) and is
  what ``MfUsg.load`` uses so arbitrary real models round-trip without loss.
"""

from __future__ import annotations

import io

import numpy as np

from ..pakbase import Package
from ._usgt_returnflow import read_u1dint_list, write_u1dint_list
from .mfusg import MfUsg

_SPD_KEYS = ("ib0", "ib1", "ibm1", "icb0", "icb1", "icbm1")


def _fmt_num(v):
    """Format a float so it reloads exactly (repr round-trips Python floats)."""
    return repr(float(v))


class MfUsgTib(Package):
    """MODFLOW-USG Transient Ibound (TIB) Package Class.

    Parameters
    ----------
    model : flopy.mfusg.MfUsg
        Parent model.
    stress_period_data : dict, optional
        Semantic per-stress-period data keyed by zero-based ``kper``. Each value
        is a dict with any of the keys below (missing keys mean "no cells of
        that kind this period"). All node numbers are **0-based**; the writer
        converts to 1-based file ids.

        ``"ib0"``  : array-like of node ids to inactivate (``IBOUND := 0``).
        ``"icb0"`` : array-like of node ids to inactivate transport
            (``ICBUND := 0``); requires a BCT package on the model.
        ``"ib1"``  : list of ``(node, head)`` to activate (``IBOUND := 1``).
        ``"ibm1"`` : list of ``(node, head)`` for prescribed head
            (``IBOUND := -1``).
            ``head`` is a float (written as ``HEAD <value>``), the string
            ``"AVHEAD"`` (head averaged from connected active cells), or ``None``
            (no keyword — reuse the cell's existing head; USG-T requires the cell
            to already be active in that case).
        ``"icb1"``  : list of ``(node, conc)`` to activate transport
            (``ICBUND := 1``); requires BCT.
        ``"icbm1"`` : list of ``(node, conc)`` for prescribed concentration
            (``ICBUND := -1``); requires BCT.
            ``conc`` is a sequence of ``MCOMP`` floats (written as
            ``CONC c1 c2 ...``), the string ``"AVCONC"``, or ``None``.
    blocks : dict[int, str], optional
        Legacy hand-built per-SP raw text keyed by zero-based ``kper`` (each
        value begins at that period's header line). Retained for backward
        compatibility; prefer ``stress_period_data``.
    raw_body : str, optional
        Verbatim file body (set by ``load`` in raw round-trip mode).
    extension : str, optional
        File extension (default "tib").
    unitnumber : int, optional
        File unit number (default 102, matching common GMS/GV8 convention).
    filenames : str or list of str, optional
        Package filename override.

    Notes
    -----
    Provide exactly one of ``stress_period_data``, ``blocks``, or ``raw_body``.
    With none of them the package writes an all-zero header for every stress
    period (a valid no-op TIB).
    """

    def __init__(
        self,
        model,
        stress_period_data=None,
        blocks=None,
        raw_body=None,
        extension="tib",
        unitnumber=None,
        filenames=None,
    ):
        msg = (
            "Model object must be of type flopy.mfusg.MfUsg\n"
            f"but received type: {type(model)}."
        )
        assert isinstance(model, MfUsg), msg

        if unitnumber is None:
            unitnumber = MfUsgTib._defaultunit()

        filenames = self._prepare_filenames(filenames, num=1)

        super().__init__(
            model,
            extension=extension,
            name=self._ftype(),
            unit_number=unitnumber,
            filenames=filenames,
        )
        self._generate_heading()

        if stress_period_data is not None:
            self.stress_period_data = self._normalize_spd(stress_period_data, model)
        else:
            self.stress_period_data = None
        self.blocks = dict(blocks) if blocks else {}
        self.raw_body = raw_body
        self.parent.add_package(self)

    # ------------------------------------------------------------------ #
    # Authoring / validation
    # ------------------------------------------------------------------ #
    @staticmethod
    def _normalize_spd(stress_period_data, model):
        """Validate and normalize semantic stress-period data (0-based nodes)."""
        transport = getattr(model, "bct", None) is not None
        ncomp = int(getattr(model.bct, "mcomp", 0)) if transport else 0

        out = {}
        for kper, sp in stress_period_data.items():
            unknown = set(sp) - set(_SPD_KEYS)
            if unknown:
                raise ValueError(
                    f"Unknown TIB stress_period_data keys: {sorted(unknown)}; "
                    f"valid keys are {list(_SPD_KEYS)}."
                )
            d = {}
            for key in ("ib0", "icb0"):
                vals = sp.get(key)
                if vals is not None and len(vals) > 0:
                    arr = np.asarray(vals, dtype=int).ravel()
                    if (arr < 0).any():
                        raise ValueError(
                            f"TIB '{key}' nodes must be 0-based (>= 0); "
                            f"got minimum {int(arr.min())}."
                        )
                    d[key] = arr
            for key in ("ib1", "ibm1"):
                recs = [MfUsgTib._norm_head_rec(r, key) for r in (sp.get(key) or [])]
                if recs:
                    d[key] = recs
            for key in ("icb1", "icbm1"):
                recs = [
                    MfUsgTib._norm_conc_rec(r, key, ncomp, transport)
                    for r in (sp.get(key) or [])
                ]
                if recs:
                    d[key] = recs

            if not transport and any(k in d for k in ("icb0", "icb1", "icbm1")):
                raise ValueError(
                    "TIB transport data (icb0/icb1/icbm1) requires an active "
                    "BCT package on the model."
                )
            if d:
                out[int(kper)] = d
        return out

    @staticmethod
    def _norm_head_rec(rec, key):
        node, head = rec
        node = int(node)
        if node < 0:
            raise ValueError(f"TIB '{key}' node must be 0-based (>= 0); got {node}.")
        if head is None:
            return (node, None)
        if isinstance(head, str):
            if head.upper() != "AVHEAD":
                raise ValueError(
                    f"TIB '{key}' head option must be a number, 'AVHEAD', or "
                    f"None; got {head!r}."
                )
            return (node, "AVHEAD")
        return (node, float(head))

    @staticmethod
    def _norm_conc_rec(rec, key, ncomp, transport):
        node, conc = rec
        node = int(node)
        if node < 0:
            raise ValueError(f"TIB '{key}' node must be 0-based (>= 0); got {node}.")
        if conc is None:
            return (node, None)
        if isinstance(conc, str):
            if conc.upper() != "AVCONC":
                raise ValueError(
                    f"TIB '{key}' conc option must be a sequence of values, "
                    f"'AVCONC', or None; got {conc!r}."
                )
            return (node, "AVCONC")
        arr = np.asarray(conc, dtype=float).ravel()
        if transport and ncomp and len(arr) != ncomp:
            raise ValueError(
                f"TIB '{key}' CONC record for node {node} has {len(arr)} "
                f"value(s); expected MCOMP = {ncomp}."
            )
        return (node, arr)

    # ------------------------------------------------------------------ #
    # Writing
    # ------------------------------------------------------------------ #
    def write_file(self, check=False):
        """Write the TIB package file."""
        nper = self.parent.nper
        with open(self.fn_path, "w") as f:
            f.write(f"{self.heading}\n")
            if self.raw_body is not None:
                f.write(self.raw_body)
                if self.raw_body and not self.raw_body.endswith("\n"):
                    f.write("\n")
                return
            if self.stress_period_data is not None:
                self._write_semantic(f)
                return
            # Legacy raw-blocks path; also emits all-zero headers for empty SPs.
            for kper in range(nper):
                block = self.blocks.get(kper, "")
                if not block:
                    inbct = getattr(self.parent, "bct", None)
                    n = 6 if inbct is not None else 3
                    f.write((" 0" * n) + "\n")
                else:
                    if not block.endswith("\n"):
                        block = block + "\n"
                    f.write(block)

    def _write_semantic(self, f):
        """Write semantic ``stress_period_data`` as TIB stress-period records."""
        transport = getattr(self.parent, "bct", None) is not None
        ncomp = int(getattr(self.parent.bct, "mcomp", 0)) if transport else 0
        spd = self.stress_period_data or {}

        for kper in range(self.parent.nper):
            sp = spd.get(kper, {})
            ib0 = sp.get("ib0", [])
            ib1 = sp.get("ib1", [])
            ibm1 = sp.get("ibm1", [])
            icb0 = sp.get("icb0", [])
            icb1 = sp.get("icb1", [])
            icbm1 = sp.get("icbm1", [])

            if transport:
                f.write(
                    f" {len(ib0)} {len(ib1)} {len(ibm1)}"
                    f" {len(icb0)} {len(icb1)} {len(icbm1)}\n"
                )
            else:
                f.write(f" {len(ib0)} {len(ib1)} {len(ibm1)}\n")

            if len(ib0) > 0:
                write_u1dint_list(f, [int(n) + 1 for n in ib0])
            for node, head in ib1:
                f.write(self._fmt_head_record(node, head))
            for node, head in ibm1:
                f.write(self._fmt_head_record(node, head))

            if transport:
                if len(icb0) > 0:
                    write_u1dint_list(f, [int(n) + 1 for n in icb0])
                for node, conc in icb1:
                    f.write(self._fmt_conc_record(node, conc))
                for node, conc in icbm1:
                    f.write(self._fmt_conc_record(node, conc))

    @staticmethod
    def _fmt_head_record(node, head):
        n1 = int(node) + 1
        if head is None:
            return f" {n1}\n"
        if isinstance(head, str):
            return f" {n1} {head}\n"
        return f" {n1} HEAD {_fmt_num(head)}\n"

    @staticmethod
    def _fmt_conc_record(node, conc):
        n1 = int(node) + 1
        if conc is None:
            return f" {n1}\n"
        if isinstance(conc, str):
            return f" {n1} {conc}\n"
        vals = " ".join(_fmt_num(c) for c in conc)
        return f" {n1} CONC {vals}\n"

    # ------------------------------------------------------------------ #
    # Loading
    # ------------------------------------------------------------------ #
    @classmethod
    def load(cls, f, model, nper=None, ext_unit_dict=None, check=False, parse=False):
        """Load an existing TIB file.

        By default (``parse=False``) the file body is preserved verbatim for a
        byte-exact round-trip (``raw_body``). With ``parse=True`` the body is
        parsed into semantic ``stress_period_data``; if it uses a construct the
        parser does not model the loader falls back to the raw round-trip rather
        than writing partial data.
        """
        msg = (
            "Model object must be of type flopy.mfusg.MfUsg\n"
            f"but received type: {type(model)}."
        )
        assert isinstance(model, MfUsg), msg

        if model.verbose:
            print("loading tib package file...")

        if hasattr(f, "read"):
            fh = f
        else:
            fh = open(f, "r")

        try:
            lines = fh.readlines()
        finally:
            if not hasattr(f, "read"):
                fh.close()

        # Strip the leading comment header (lines starting with "#").
        i = 0
        while i < len(lines) and lines[i].lstrip().startswith("#"):
            i += 1
        raw_body = "".join(lines[i:])

        stress_period_data = None
        if parse:
            try:
                stress_period_data = cls._parse_semantic(lines, model)
            except (NotImplementedError, ValueError, IndexError) as e:
                if model.verbose:
                    print(f"  TIB semantic parse failed ({e}); keeping raw round-trip.")
                stress_period_data = None

        # Preserve original unit + filename from NAM.
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

        if stress_period_data is not None:
            return cls(
                model,
                stress_period_data=stress_period_data,
                unitnumber=unitnumber,
                filenames=filenames,
            )
        return cls(
            model,
            raw_body=raw_body,
            unitnumber=unitnumber,
            filenames=filenames,
        )

    @classmethod
    def _parse_semantic(cls, lines, model):
        """Parse raw TIB lines into semantic 0-based ``stress_period_data``.

        Raises on any unsupported construct or inconsistency so ``load`` can
        fall back to the raw round-trip instead of writing partial data.
        """
        transport = getattr(model, "bct", None) is not None
        ncomp = int(getattr(model.bct, "mcomp", 0)) if transport else 0
        nper = model.nper

        filtered = [
            ln for ln in lines if ln.strip() and not ln.lstrip().startswith("#")
        ]
        fh = io.StringIO("".join(filtered))

        out = {}
        for kper in range(nper):
            header = fh.readline()
            if not header:
                break
            toks = header.split()
            if transport:
                if len(toks) < 6:
                    raise ValueError(
                        f"TIB SP{kper + 1} header needs 6 integers (transport "
                        f"active); got {header.strip()!r}."
                    )
                nib0, nib1, nibm1, nicb0, nicb1, nicbm1 = (int(t) for t in toks[:6])
            else:
                if len(toks) < 3:
                    raise ValueError(
                        f"TIB SP{kper + 1} header needs 3 integers; got "
                        f"{header.strip()!r}."
                    )
                nib0, nib1, nibm1 = (int(t) for t in toks[:3])
                nicb0 = nicb1 = nicbm1 = 0

            sp = {}
            if nib0 > 0:
                sp["ib0"] = np.array(
                    [v - 1 for v in read_u1dint_list(fh, nib0)], dtype=int
                )
            if nib1 > 0:
                sp["ib1"] = [cls._read_head_record(fh) for _ in range(nib1)]
            if nibm1 > 0:
                sp["ibm1"] = [cls._read_head_record(fh) for _ in range(nibm1)]
            if transport:
                if nicb0 > 0:
                    sp["icb0"] = np.array(
                        [v - 1 for v in read_u1dint_list(fh, nicb0)], dtype=int
                    )
                if nicb1 > 0:
                    sp["icb1"] = [
                        cls._read_conc_record(fh, ncomp) for _ in range(nicb1)
                    ]
                if nicbm1 > 0:
                    sp["icbm1"] = [
                        cls._read_conc_record(fh, ncomp) for _ in range(nicbm1)
                    ]
            if sp:
                out[kper] = sp

        rest = fh.read()
        if rest.strip():
            raise ValueError(
                "Trailing TIB content not consumed by the semantic parser."
            )
        return out

    @staticmethod
    def _read_head_record(fh):
        line = fh.readline()
        if not line:
            raise ValueError("Unexpected EOF reading a TIB head record.")
        toks = line.split()
        node = int(toks[0]) - 1
        if len(toks) == 1:
            return (node, None)
        kw = toks[1].upper()
        if kw == "HEAD":
            return (node, float(toks[2]))
        if kw == "AVHEAD":
            return (node, "AVHEAD")
        raise NotImplementedError(f"Unsupported TIB head option: {toks[1]!r}.")

    @staticmethod
    def _read_conc_record(fh, ncomp):
        line = fh.readline()
        if not line:
            raise ValueError("Unexpected EOF reading a TIB conc record.")
        toks = line.split()
        node = int(toks[0]) - 1
        if len(toks) == 1:
            return (node, None)
        kw = toks[1].upper()
        if kw == "CONC":
            vals = [float(t) for t in toks[2 : 2 + ncomp]]
            if len(vals) != ncomp:
                raise ValueError(
                    f"TIB CONC record needs {ncomp} value(s); got {len(vals)}."
                )
            return (node, np.array(vals))
        if kw == "AVCONC":
            return (node, "AVCONC")
        raise NotImplementedError(f"Unsupported TIB conc option: {toks[1]!r}.")

    @staticmethod
    def _ftype():
        return "TIB"

    @staticmethod
    def _defaultunit():
        return 102
