"""MfUsgDpt (USG-T DPT): read sequence, components, DPF coupling, handles.

Expected read sequences are written by hand from GWT2DPTU1AR (USG-T 2.8.1
gwt2dptu1.f; line numbers below are 2.8.1 unless noted, 2.7.0 in
parentheses). They never come from the FloPy writer. Every array is one
CONSTANT record with a distinct value, so the sequence of values identifies
the arrays. Engine tests (requires_usgt) compare a FloPy file with a
hand-written file on the same model.
"""

import builtins
import io
import re
import struct
import warnings
from pathlib import Path

import numpy as np
import pytest

from autotest.usgt_exe import USGT_EXE, requires_usgt
from flopy.mfusg import (
    MfUsg,
    MfUsgBas,
    MfUsgBcf,
    MfUsgBct,
    MfUsgDisU,
    MfUsgDpf,
    MfUsgDpt,
    MfUsgLpf,
    MfUsgOc,
    MfUsgSms,
)
from flopy.modflow import ModflowDis

# distinct values per array
PHIF, PRS, BULK, DLIM, DDTR, SIM = 0.41, 0.21, 1.61, 0.51, 0.111, 0.91
ADS, ZW, FW, C1, CIM2 = 1.1, 3.1, 5.1, 11.0, 31.0
SENTINEL = "SENTINEL DPT FILE\n"
BCT_KW = {"mcomp": 1, "icbndflg": 1, "idisp": 0, "prsity": 0.2, "conc": 21.0}
IMCOMP_BCT = {"imcomp": 1, "imconc": 2.0}
BASE = {
    "phif": PHIF,
    "prsityim": PRS,
    "bulkdim": BULK,
    "dlim": DLIM,
    "ddtr": DDTR,
    "sim": SIM,
    "concim": C1,
}


def _model(ws, name, structured=False, dpf=False, bct=None, manual=None):
    m = MfUsg(modelname=name, model_ws=str(ws), structured=structured)
    if structured:
        ModflowDis(m, nlay=1, nrow=1, ncol=3, nper=1, steady=True)
    else:
        MfUsgDisU(
            m,
            nodes=3,
            nlay=1,
            njag=3,
            iac=np.ones(3, dtype=int),
            ja=np.arange(1, 4),
            cl12=np.ones(3),
            fahl=np.ones(3),
            area=1.0,
            top=[0.0] * 3,
            bot=[-1.0] * 3,
        )
    MfUsgBas(m)
    if dpf:
        MfUsgBcf(m)
        MfUsgDpf(m)
    if bct is not None:
        MfUsgBct(m, **{**BCT_KW, **bct})
    if manual is not None:  # manual counts mode (no BCT), as the A-W tests
        m.itrnsp = 1
        m.mcomp, m.iheat, m.imcomp = manual
    return m


def _values(text):
    """(header tokens, CONSTANT values) of a DPT file text."""
    lines = [ln for ln in text.splitlines() if ln.strip() and not ln.startswith("#")]
    head = lines[0].split()
    vals = []
    for ln in lines[1:]:
        t = ln.split()
        assert t[0].upper() == "CONSTANT", ln
        vals.append(round(float(t[1]), 6))
    return head, vals


def _written(dpt):
    dpt.write_file()
    return _values(Path(dpt.fn_path).read_text())


def _hand(ws, name, header, values):
    f = Path(ws) / f"{name}_hand.dpt"
    f.write_text(header + "\n" + "".join(f"CONSTANT {v}\n" for v in values))
    return f


def _assert_rejected_untouched(dpt, exc, match):
    """write_file raises before opening: the file and a handle stay intact."""
    p = Path(dpt.fn_path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(SENTINEL)
    with pytest.raises(exc, match=match):
        dpt.write_file()
    assert p.read_text() == SENTINEL
    buf = io.StringIO("PREFIX\n")
    buf.seek(0, io.SEEK_END)
    pos = buf.tell()
    with pytest.raises(exc, match=match):
        dpt.write_file(f=buf)
    assert not buf.closed and buf.getvalue() == "PREFIX\n" and buf.tell() == pos


def _first(util):
    return round(float(util.array.ravel()[0]), 6)


# ---------------------------------------------------------------- controls


@pytest.mark.parametrize("structured", [False, True], ids=["disu", "structured"])
def test_dpt_base_sequence(function_tmpdir, structured):
    m = _model(function_tmpdir, "cb", structured=structured, bct={})
    dpt = MfUsgDpt(m, **BASE)
    head, vals = _written(dpt)
    assert head == ["0", "0", "1", "0", "0", "0", "0"]
    assert vals == [PHIF, PRS, DDTR, C1]  # 263, 274, 301, 604


def test_dpt_adsorption_decay_sequence(function_tmpdir):
    m = _model(function_tmpdir, "cad", bct={})
    dpt = MfUsgDpt(
        m,
        iadsorbim=1,
        izodim=1,
        ifodim=3,
        adsorbim=ADS,
        zodrwim=ZW,
        fodrwim=FW,
        fodrsim=6.1,
        **BASE,
    )
    assert _written(dpt)[1] == [PHIF, PRS, BULK, DDTR, ADS, ZW, FW, 6.1, C1]


@pytest.mark.parametrize("flag", ["izodim", "ifodim"])
def test_dpt_decay2_without_adsorption_reads_nothing(function_tmpdir, flag):
    # 439-440 / 520-521: option 2 reads only with IADSORBIM>0
    m = _model(function_tmpdir, "c2", bct={})
    dpt = MfUsgDpt(m, **{flag: 2}, zodrsim=4.1, fodrsim=6.1, **BASE)
    assert _written(dpt)[1] == [PHIF, PRS, DDTR, C1]


def test_dpt_inputsat_without_dpf_serialized(function_tmpdir):
    # serialization only (111-114, 324-335); not promoted as executable
    m = _model(function_tmpdir, "cis", bct={})
    dpt = MfUsgDpt(m, inputsat=True, **BASE)
    head, vals = _written(dpt)
    assert head[7:] == ["INPUTSAT"] and vals == [PHIF, PRS, DDTR, SIM, C1]


def test_dpt_load_basic(function_tmpdir):
    m = _model(function_tmpdir, "clb", bct={})
    f = _hand(function_tmpdir, "clb", "0 0 1 1 0 0 0", [PHIF, PRS, BULK, DDTR, ADS, C1])
    d = MfUsgDpt.load(str(f), m)
    got = [d.phif, d.prsityim, d.bulkdim, d.ddtr, d.adsorbim[0], d.concim[0]]
    assert [_first(u) for u in got] == [PHIF, PRS, BULK, DDTR, ADS, C1]


EX4_DPT = (
    "         0        33         0         1         0         0         0\n"
    + "CONSTANT          1                                #icbundim"
    + " " * 22
    + "\n"
    + "CONSTANT    4.000000E-01                           #phif"
    + " " * 26
    + "\n"
    + "CONSTANT    2.000000E-01                           #prsityim"
    + " " * 22
    + "\n"
    + "CONSTANT    1.600000E+00                           #bulkdim"
    + " " * 23
    + "\n"
    + "CONSTANT    1.000000E-01                           #ddtr"
    + " " * 26
    + "\n"
    + "CONSTANT    1.000000E-01                           #adsorbim"
    + " " * 22
    + "\n"
    + "CONSTANT    0.000000E+00                           #concim"
    + " " * 24
    + "\n"
)


def test_dpt_ex4_bytes(function_tmpdir, example_data_path):
    # golden from the develop writer (USGT-DPT-A2 04_tests/golden_ex4)
    ws = example_data_path / "mfusg_transport" / "Ex4_Dual_Domain"
    m = MfUsg.load(
        "DualDomain.nam",
        model_ws=str(ws),
        check=False,
        load_only=["DIS", "BAS6", "BCF6", "BCT", "DPT"],
    )
    m.change_model_ws(str(function_tmpdir))
    m.dpt.write_file()
    assert Path(m.dpt.fn_path).read_text() == EX4_DPT


# ---------------------------------------------------------------- unit number


def test_dpt_default_unit_number_in_name_file(function_tmpdir):
    m = _model(function_tmpdir, "un", bct={})
    dpt = MfUsgDpt(m, **BASE)
    assert dpt.unit_number[0] == 158
    entries = m.get_name_file_entries()
    assert re.search(r"^DPT\s+158\s+un\.dpt$", entries, re.M), entries


# ---------------------------------------------------------------- handles


def test_dpt_write_handle_appends_at_cursor_and_stays_open(function_tmpdir):
    m = _model(function_tmpdir, "wh", bct={})
    dpt = MfUsgDpt(m, **BASE)
    buf = io.StringIO("PREFIX\n")
    buf.seek(0, io.SEEK_END)
    dpt.write_file(f=buf)
    assert not buf.closed
    text = buf.getvalue()
    assert text.startswith("PREFIX\n")
    assert _values(text[len("PREFIX\n") :])[1] == [PHIF, PRS, DDTR, C1]
    assert buf.tell() == len(text)
    assert not Path(dpt.fn_path).exists()  # the handle is used, not fn_path


def test_dpt_load_handle_not_closed_cursor_after_data(function_tmpdir):
    m = _model(function_tmpdir, "lh", bct={})
    body = (
        f"0 0 1 0 0 0 0\nCONSTANT {PHIF}\nCONSTANT {PRS}\nCONSTANT {DDTR}\n"
        f"CONSTANT {C1}\n"
    )
    buf = io.StringIO(body + "TRAILER\n")
    d = MfUsgDpt.load(buf, m)
    assert not buf.closed
    assert buf.tell() == len(body)
    assert _first(d.concim[0]) == C1


def test_dpt_load_file_closed_after_error(function_tmpdir, monkeypatch):
    m = _model(function_tmpdir, "lfe", bct={})
    f = _hand(function_tmpdir, "lfe", "0 0 1 0 0 0 0 SD_MASS_TR", [PHIF])
    opened = []
    real_open = builtins.open

    def spy(*args, **kwargs):
        fh = real_open(*args, **kwargs)
        opened.append(fh)
        return fh

    monkeypatch.setattr(builtins, "open", spy)
    with pytest.raises(NotImplementedError, match="SD_MASS_TR"):
        MfUsgDpt.load(str(f), m)
    monkeypatch.undo()
    assert opened and all(fh.closed for fh in opened)


# ---------------------------------------------------------------- counts and lifecycle


def test_dpt_counts_needed_at_construction(function_tmpdir):
    m = _model(function_tmpdir, "nc")  # no BCT, itrnsp=0
    with pytest.raises(ValueError, match="BCT"):
        MfUsgDpt(m, **BASE)


@pytest.mark.parametrize(
    "case", ["replaced_other_counts", "removed", "mutated", "added_without_ctor"]
)
def test_dpt_lifecycle_rejected_before_open(function_tmpdir, case):
    m = _model(function_tmpdir, f"lc_{case}", bct={})
    dpt = MfUsgDpt(m, **BASE)
    if case == "replaced_other_counts":
        MfUsgBct(m, **{**BCT_KW, "mcomp": 2})
    elif case == "removed":
        m.remove_package("BCT")
    elif case == "mutated":
        m.mcomp = 2
    else:  # registered BCT that did not set the counts (unsynced)
        bct = MfUsgBct(m, **{**BCT_KW, "mcomp": 2}, add_package=False)
        m.add_package(bct)
    _assert_rejected_untouched(dpt, ValueError, "BCT|count")


def test_dpt_lifecycle_replaced_same_counts_writes(function_tmpdir):
    m = _model(function_tmpdir, "lcs", bct={})
    dpt = MfUsgDpt(m, **BASE)
    MfUsgBct(m, **BCT_KW)
    assert _written(dpt)[1] == [PHIF, PRS, DDTR, C1]


def test_dpt_manual_counts_mutated_rejected(function_tmpdir):
    m = _model(function_tmpdir, "lcm", manual=(1, 0, 0))
    dpt = MfUsgDpt(m, **BASE)
    m.mcomp = 2
    _assert_rejected_untouched(dpt, ValueError, "count")


# ---------------------------------------------------------------- DPF coupling


@pytest.mark.parametrize(
    "idispim, expected",
    [(0, [PRS, DDTR, C1]), (1, [PRS, DLIM, DDTR, C1])],
    ids=["idispim0", "idispim1"],
)
def test_dpt_dpf_layout_without_phif(function_tmpdir, idispim, expected):
    # 263 (PHIF only IDPF=0), 292 (DLIM IDPF!=0 and IDISPIM!=0)
    m = _model(function_tmpdir, f"dpf{idispim}", dpf=True, bct={})
    dpt = MfUsgDpt(m, idispim=idispim, **BASE)
    head, vals = _written(dpt)
    assert head == ["0", "0", "1", "0", str(idispim), "0", "0"]
    assert vals == expected


def test_dpt_dpf_added_after_dpt(function_tmpdir):
    m = _model(function_tmpdir, "dpa", bct={})
    dpt = MfUsgDpt(m, **BASE)
    assert _written(dpt)[1] == [PHIF, PRS, DDTR, C1]
    MfUsgBcf(m)
    MfUsgDpf(m)
    assert _written(dpt)[1] == [PRS, DDTR, C1]


def test_dpt_dpf_replaced_writes(function_tmpdir):
    m = _model(function_tmpdir, "dpx", dpf=True, bct={})
    dpt = MfUsgDpt(m, **BASE)
    MfUsgDpf(m)  # replaces the registered DPF package
    assert _written(dpt)[1] == [PRS, DDTR, C1]


@pytest.mark.parametrize("opt", ["frahk", "mobilesat", "inputsat"])
def test_dpt_dpf_single_domain_options_rejected(function_tmpdir, opt):
    # 99: options parsed only with IDPF=0
    m = _model(function_tmpdir, f"dpo_{opt}", dpf=True, bct={})
    dpt = MfUsgDpt(m, **{opt: True}, **BASE)
    _assert_rejected_untouched(dpt, ValueError, "DPF")


def test_dpt_dpf_removed_stale_idpf_rejected(function_tmpdir):
    m = _model(function_tmpdir, "dpr", dpf=True, bct={})
    dpt = MfUsgDpt(m, idispim=1, **BASE)
    m.remove_package("DPF")
    assert m.idpf == 1  # set by MfUsgDpf, not reset when it is removed
    _assert_rejected_untouched(dpt, ValueError, "idpf")


def test_dpt_manual_idpf_without_dpf_rejected(function_tmpdir):
    m = _model(function_tmpdir, "dpm", bct={})
    m.idpf = 1
    dpt = MfUsgDpt(m, **BASE)
    _assert_rejected_untouched(dpt, ValueError, "idpf")


def test_dpt_dpf_registered_idpf_zero_rejected(function_tmpdir):
    m = _model(function_tmpdir, "dpz", dpf=True, bct={})
    m.idpf = 0
    dpt = MfUsgDpt(m, **BASE)
    _assert_rejected_untouched(dpt, ValueError, "idpf")


def test_dpt_dpf_load_without_phif(function_tmpdir):
    m = _model(function_tmpdir, "dpl", dpf=True, bct={})
    f = _hand(function_tmpdir, "dpl", "0 0 1 0 1 0 0", [PRS, DLIM, DDTR, C1])
    d = MfUsgDpt.load(str(f), m)
    got = [d.prsityim, d.dlim, d.ddtr, d.concim[0]]
    assert [_first(u) for u in got] == [PRS, DLIM, DDTR, C1]


def test_dpt_dpf_load_ignores_single_domain_options(function_tmpdir):
    # with IDPF=1 the option words are not parsed (99)
    m = _model(function_tmpdir, "dpi", dpf=True, bct={})
    f = _hand(function_tmpdir, "dpi", "0 0 1 0 0 0 0 FRAHK INPUTSAT", [PRS, DDTR, C1])
    d = MfUsgDpt.load(str(f), m)
    assert (d.frahk, d.mobilesat, d.inputsat) == (0, 0, 0)
    assert [_first(u) for u in (d.prsityim, d.ddtr, d.concim[0])] == [PRS, DDTR, C1]


def _write_dpf_model(ws, swap_nam=False):
    m = _model(ws, "nm", dpf=True, bct={})
    MfUsgDpt(m, idispim=1, **BASE)
    m.write_input()
    nam = Path(ws) / "nm.nam"
    if swap_nam:  # list DPT before DPF in the name file
        lines = nam.read_text().splitlines(keepends=True)
        dpt = [ln for ln in lines if ln.startswith("DPT")]
        rest = [ln for ln in lines if not ln.startswith("DPT")]
        i = next(k for k, ln in enumerate(rest) if ln.startswith("DPF"))
        nam.write_text("".join(rest[:i] + dpt + rest[i:]))
    return nam


@pytest.mark.parametrize("swap_nam", [False, True], ids=["dpf_first", "dpt_first"])
def test_dpt_nam_load_with_dpf(function_tmpdir, swap_nam):
    nam = _write_dpf_model(function_tmpdir, swap_nam)
    m = MfUsg.load(nam.name, model_ws=str(function_tmpdir), check=False)
    d = m.get_package("DPT")
    got = [d.prsityim, d.dlim, d.ddtr, d.concim[0]]
    assert [_first(u) for u in got] == [PRS, DLIM, DDTR, C1]


def test_dpt_nam_load_without_loading_dpf(function_tmpdir):
    # model.idpf comes from the name file; the DPF object is not required
    nam = _write_dpf_model(function_tmpdir)
    m = MfUsg.load(
        nam.name,
        model_ws=str(function_tmpdir),
        check=False,
        load_only=["DISU", "BAS6", "BCT", "DPT"],
    )
    assert m.idpf == 1 and m.get_package("DPF") is None
    d = m.get_package("DPT")
    assert _first(d.dlim) == DLIM and _first(d.concim[0]) == C1
    # writing it alone is ambiguous (the name file would not list DPF)
    _assert_rejected_untouched(d, ValueError, "idpf")


# ---------------------------------------------------------------- options


@pytest.mark.parametrize(
    "opts, flags",
    [
        ("FRAHK", (1, 0, 0)),
        ("frahk", (1, 0, 0)),
        ("MOBILESAT", (0, 1, 0)),
        ("mobilesat", (0, 1, 0)),
        ("Inputsat", (0, 0, 1)),
        ("INPUTSAT MOBILESAT", (0, 1, 0)),  # 108-114: the last IMSAT wins
        ("MOBILESAT INPUTSAT", (0, 0, 1)),
        ("FRAHK inputsat mobilesat", (1, 1, 0)),
    ],
    ids=[
        "frahk_upper",
        "frahk_lower",
        "mobilesat_upper",
        "mobilesat_lower",
        "inputsat_mixed",
        "inputsat_then_mobilesat",
        "mobilesat_then_inputsat",
        "frahk_lower_pair",
    ],
)
def test_dpt_load_options(function_tmpdir, opts, flags):
    m = _model(function_tmpdir, "lo", bct={})
    inputsat = flags[2]
    arrays = [PHIF, PRS, DDTR] + ([SIM] if inputsat else []) + [C1]
    f = _hand(function_tmpdir, "lo", f"0 0 1 0 0 0 0 {opts}", arrays)
    d = MfUsgDpt.load(str(f), m)
    assert (int(d.frahk), int(d.mobilesat), int(d.inputsat)) == flags
    assert _first(d.concim[0]) == C1
    if inputsat:
        assert _first(d.sim) == SIM


def test_dpt_mobilesat_and_inputsat_ambiguous_when_built(function_tmpdir):
    # two booleans for one IMSAT value: ambiguous on construction only
    m = _model(function_tmpdir, "mx", bct={})
    with pytest.raises(ValueError, match="IMSAT"):
        MfUsgDpt(m, mobilesat=True, inputsat=True, **BASE)


# ---------------------------------------------------------------- immobile components


@pytest.mark.parametrize("mode", ["bct", "manual"])
@pytest.mark.parametrize(
    "concim, expected",
    [([C1, CIM2], [C1, CIM2]), (C1, [C1, C1])],
    ids=["list", "scalar"],
)
def test_dpt_imcomp_concim_written(function_tmpdir, mode, concim, expected):
    # 364 loop ICOMP=1..NTCOMP; (376) immobile -> CONCIM only: common format
    # with IADSORBIM=0, IZODIM=0, IFODIM=0
    if mode == "bct":
        kw = {"bct": IMCOMP_BCT}
    else:
        kw = {"manual": (1, 0, 1)}
    m = _model(function_tmpdir, f"im_{mode}", **kw)
    dpt = MfUsgDpt(m, **{**BASE, "concim": concim})
    assert _written(dpt)[1] == [PHIF, PRS, DDTR, *expected]


@pytest.mark.parametrize(
    "concim", [[C1], [C1, CIM2, 1.0], [], True], ids=["short", "long", "empty", "bool"]
)
def test_dpt_concim_length_checked(function_tmpdir, concim):
    m = _model(function_tmpdir, "iml", bct=IMCOMP_BCT)
    with pytest.raises(ValueError, match="concim"):
        MfUsgDpt(m, **{**BASE, "concim": concim})


def test_dpt_imcomp_load(function_tmpdir):
    m = _model(function_tmpdir, "imld", bct=IMCOMP_BCT)
    body = (
        f"0 0 1 0 0 0 0\nCONSTANT {PHIF}\nCONSTANT {PRS}\nCONSTANT {DDTR}\n"
        f"CONSTANT {C1}\nCONSTANT {CIM2}\n"
    )
    buf = io.StringIO(body)
    d = MfUsgDpt.load(buf, m)
    assert [_first(c) for c in d.concim] == [C1, CIM2]
    assert buf.tell() == len(body)


@pytest.mark.parametrize(
    "kw",
    [{"iadsorbim": 1}, {"izodim": 1}, {"ifodim": 1}, {"aw_adsorbim": True}],
    ids=["ads", "zod", "fod", "aw"],
)
def test_dpt_imcomp_reactions_not_supported(function_tmpdir, kw):
    # 388-391 vs (376): 2.8.1 reads reaction arrays for immobile components
    m = _model(function_tmpdir, "imr", bct=IMCOMP_BCT)
    with pytest.raises(NotImplementedError, match="IMCOMP"):
        MfUsgDpt(m, **kw, **BASE)


def test_dpt_imcomp_reactions_load_rejected_after_header(function_tmpdir):
    m = _model(function_tmpdir, "imrl", bct=IMCOMP_BCT)
    head = "0 0 1 1 0 0 0\n"
    buf = io.StringIO(head + f"CONSTANT {PHIF}\n")
    with pytest.raises(NotImplementedError, match="IMCOMP"):
        MfUsgDpt.load(buf, m)
    assert not buf.closed and buf.tell() == len(head)


# ---------------------------------------------------------------- heat, htcondsim

HEAT = {
    "iheat": 1,
    "bulkd": 1.6,
    "htcondw": 0.58,
    "rhow": 1000.0,
    "htcapw": 4186.0,
    "htcaps": 51.0,
    "htconds": 52.0,
    "heat": 31.0,
}


def test_dpt_heat_not_supported_construction(function_tmpdir):
    m = _model(function_tmpdir, "ht", bct=HEAT)
    with pytest.raises(NotImplementedError, match="IHEAT"):
        MfUsgDpt(m, **BASE)


def test_dpt_heat_not_supported_load_after_header(function_tmpdir):
    m = _model(function_tmpdir, "htl", bct=HEAT)
    head = "0 0 1 0 0 0 0\n"
    buf = io.StringIO(head + f"CONSTANT {PHIF}\n")
    with pytest.raises(NotImplementedError, match="IHEAT"):
        MfUsgDpt.load(buf, m)
    assert not buf.closed and buf.tell() == len(head)


def test_dpt_htcondsim_deprecated_not_written(function_tmpdir):
    # USG-T reads one heat array (340-353)
    m = _model(function_tmpdir, "hc", bct={})
    with pytest.warns(DeprecationWarning, match="htcondsim"):
        dpt = MfUsgDpt(m, htcondsim=0.82, **BASE)
    assert 0.82 not in _written(dpt)[1]


def test_dpt_htcondsim_default_no_warning(function_tmpdir):
    m = _model(function_tmpdir, "hd", bct={})
    with warnings.catch_warnings():
        warnings.simplefilter("error", DeprecationWarning)
        MfUsgDpt(m, **BASE)


# ---------------------------------------------------------------- SD_MASS_TR, chain


@pytest.mark.parametrize("opt", ["SD_MASS_TR", "sd_mass_tr"], ids=["upper", "lower"])
def test_dpt_sd_mass_tr_load_rejected_after_header(function_tmpdir, opt):
    # 153-154 (DDTTR omitted, DDTTRSD per mobile species 590-600)
    m = _model(function_tmpdir, "sd", bct={})
    head = f"0 0 1 0 0 0 0 {opt}\n"
    buf = io.StringIO(head + f"CONSTANT {PHIF}\n")
    with pytest.raises(NotImplementedError, match="SD_MASS_TR"):
        MfUsgDpt.load(buf, m)
    assert not buf.closed and buf.tell() == len(head)


def test_dpt_chain_spatialreact_not_supported(function_tmpdir):
    # 175-177, 369-385: SPTLRCTIM read when ICHAIN!=0 and ISPRCT!=0
    bct = {
        "chaindecay": True,
        "spatialreact": True,
        "nparent": [0],
        "jparent": [0],
        "stotio": [0],
    }
    m = _model(function_tmpdir, "ch", bct=bct)
    with pytest.raises(NotImplementedError, match="SPTLRCTIM"):
        MfUsgDpt(m, **BASE)


# ---------------------------------------------------------------- decay options 3/4


@pytest.mark.parametrize(
    "kw",
    [
        {"izodim": 3},
        {"ifodim": 3},
        {"izodim": 4},
        {"ifodim": 4},
        {"izodim": 4, "iadsorbim": 1},
        {"ifodim": 4, "iadsorbim": 1},
    ],
    ids=["izod3", "ifod3", "izod4", "ifod4", "izod4_ads", "ifod4_ads"],
)
def test_dpt_decay_options_not_supported(function_tmpdir, kw):
    # option 3 reads arrays only with IADSORBIM>0 (439-471, 520-552); option 4
    # only with A-W_ADSORBIM (473-504, 554-585)
    m = _model(function_tmpdir, "dz", bct={})
    with pytest.raises(NotImplementedError, match=r"IZODIM|IFODIM"):
        MfUsgDpt(m, **kw, **BASE)


@pytest.mark.parametrize("flag", ["izodim", "ifodim"])
def test_dpt_aw_decay4_not_supported(function_tmpdir, flag):
    # ZODRAWIM/FODRAWIM (494-502, 575-583) not modeled
    m = _model(function_tmpdir, "aw4", bct={})
    with pytest.raises(NotImplementedError, match=r"ZODRAWIM|FODRAWIM"):
        MfUsgDpt(m, aw_adsorbim=True, iarea_fnim=1, ikawi_fnim=1, **{flag: 4}, **BASE)


@pytest.mark.parametrize(
    "header", ["0 0 1 0 0 3 0", "0 0 1 0 0 0 4"], ids=["izod3", "ifod4"]
)
def test_dpt_decay_options_load_rejected_after_header(function_tmpdir, header):
    m = _model(function_tmpdir, "dzl", bct={})
    buf = io.StringIO(header + f"\nCONSTANT {PHIF}\n")
    with pytest.raises(NotImplementedError, match=r"IZODIM|IFODIM"):
        MfUsgDpt.load(buf, m)
    assert buf.tell() == len(header) + 1


# ------------------------------------------------- write-time contract (v2)

# write_file checks the current flags, options, BCT package and arrays, not
# only the state at construction: a mutation outside the supported contract
# is rejected before the file is opened; a supported mutation is written.
AW = {
    "aw_adsorbim": True,
    "iarea_fnim": 1,
    "ikawi_fnim": 1,
    "awamaxim": 2.5,
    "alangawim": 3.5,
    "blangawim": 4.5,
}
CHAIN_BCT = {
    "chaindecay": True,
    "spatialreact": True,
    "nparent": [0],
    "jparent": [0],
    "stotio": [0],
}


@pytest.mark.parametrize(
    "built, mutation, exc, match",
    [
        ({}, {"izodim": 4}, NotImplementedError, "IZODIM=4"),
        ({}, {"ifodim": 4}, NotImplementedError, "IFODIM=4"),
        ({}, {"izodim": 3}, NotImplementedError, "IZODIM=3"),
        ({}, {"ifodim": 3}, NotImplementedError, "IFODIM=3"),
        ({"iadsorbim": 1}, {"ifodim": 4}, NotImplementedError, "IFODIM=4"),
        ({}, {"mobilesat": True, "inputsat": True}, ValueError, "IMSAT"),
        ({"mobilesat": True}, {"inputsat": True}, ValueError, "IMSAT"),
        ({}, {"iadsorbim": 1}, ValueError, "bulkdim"),
        ({}, {"izodim": 1}, ValueError, "zodrwim"),
        ({}, {"icbndimflg": 0}, ValueError, "icbundim"),
        ({}, {"inputsat": True}, ValueError, "sim must"),
        ({}, {"concim": []}, ValueError, "concim"),
    ],
    ids=[
        "izod4",
        "ifod4",
        "izod3",
        "ifod3",
        "ads_ifod4",
        "both_imsat",
        "inputsat_added_to_mobilesat",
        "adsorption_without_arrays",
        "izod1_without_arrays",
        "icbundim_missing",
        "sim_missing",
        "concim_emptied",
    ],
)
def test_dpt_write_rechecks_mutations(function_tmpdir, built, mutation, exc, match):
    m = _model(function_tmpdir, "wm", bct={})
    dpt = MfUsgDpt(m, **built, **BASE)
    for key, value in mutation.items():
        setattr(dpt, key, value)
    _assert_rejected_untouched(dpt, exc, match)


@pytest.mark.parametrize(
    "mutation",
    [{"iadsorbim": 1}, {"izodim": 1}, AW],
    ids=["adsorption", "decay", "aw"],
)
def test_dpt_write_rechecks_imcomp_reactions(function_tmpdir, mutation):
    # 388-391 (2.8.1) against 2.7.0:376, as at construction
    m = _model(function_tmpdir, "wi", bct=IMCOMP_BCT)
    dpt = MfUsgDpt(m, **BASE)
    for key, value in mutation.items():
        setattr(dpt, key, value)
    _assert_rejected_untouched(dpt, NotImplementedError, "IMCOMP")


@pytest.mark.parametrize("how", ["replaced", "modified"])
def test_dpt_write_rechecks_bct_chain_spatialreact(function_tmpdir, how):
    # 369-385: the current BCT package is checked, with the same counts
    m = _model(function_tmpdir, f"wc_{how}", bct={})
    dpt = MfUsgDpt(m, **BASE)
    if how == "replaced":
        MfUsgBct(m, **{**BCT_KW, **CHAIN_BCT})
    else:
        bct = m.get_package("BCT")
        bct.chaindecay = 1
        bct.spatialreact = 1
    assert m._transport_context() == ("bct", (1, 0, 0))
    _assert_rejected_untouched(dpt, NotImplementedError, "SPTLRCTIM")


@pytest.mark.parametrize(
    "mutation, exc, match",
    [
        ({"iarea_fnim": 2}, NotImplementedError, "IAREA_FNIM=2"),
        ({"ikawi_fnim": 3}, NotImplementedError, "IKAWI_FNIM=3"),
        ({"izodim": 4}, NotImplementedError, "ZODRAWIM"),
        ({"iarea_fnim": 4}, ValueError, "awarea_x2im"),
        ({"ikawi_fnim": 1.5}, ValueError, "IKAWI_FNIM"),
        ({"iarea_fnim": True}, ValueError, "IAREA_FNIM"),
        ({"alangawim": []}, ValueError, "alangawim"),
    ],
    ids=[
        "iarea2",
        "ikawi3",
        "izod4",
        "iarea4_without_arrays",
        "ikawi_fraction",
        "iarea_bool",
        "langmuir_emptied",
    ],
)
def test_dpt_write_rechecks_aw_mutations(function_tmpdir, mutation, exc, match):
    m = _model(function_tmpdir, "wa", bct={})
    dpt = MfUsgDpt(m, **AW, **BASE)
    for key, value in mutation.items():
        setattr(dpt, key, value)
    _assert_rejected_untouched(dpt, exc, match)


@pytest.mark.parametrize(
    "counts, match",
    [((1, 0, 0), "awamaxim"), ((0, 0, 0), r"mcomp>0")],
    ids=["mcomp1", "mcomp0"],
)
def test_dpt_write_rechecks_aw_enabled_after_build(function_tmpdir, counts, match):
    m = _model(function_tmpdir, "wae", manual=counts)
    dpt = MfUsgDpt(m, **BASE)
    dpt.aw_adsorbim = True
    _assert_rejected_untouched(dpt, ValueError, match)


@pytest.mark.parametrize(
    "built, mutation, head, vals",
    [
        (
            {},
            {"ipakcb": 50, "mobilesat": True},
            ["50", "0", "1", "0", "0", "0", "0", "MOBILESAT"],
            [PHIF, PRS, DDTR, C1],
        ),
        (
            {},
            {"prsityim": 0.33, "concim": [5.0]},
            ["0", "0", "1", "0", "0", "0", "0"],
            [PHIF, 0.33, DDTR, 5.0],
        ),
        (
            {},
            {"izodim": 2, "ifodim": 2},
            ["0", "0", "1", "0", "0", "2", "2"],
            [PHIF, PRS, DDTR, C1],
        ),
        (
            {"iadsorbim": 1, "izodim": 1, "adsorbim": ADS, "zodrwim": ZW},
            {"iadsorbim": 0},
            ["0", "0", "1", "0", "0", "1", "0"],
            [PHIF, PRS, DDTR, ZW, C1],
        ),
        (
            AW,
            {"ikawi_fnim": 2.0, "iarea_fnim": "1"},
            ["0", "0", "1", "0", "0", "0", "0", "A-W_ADSORBIM", "1", "2"],
            [PHIF, PRS, DDTR, 2.5, 3.5, 4.5, C1],
        ),
    ],
    ids=[
        "header_and_option",
        "arrays_assigned",
        "decay2_no_arrays",
        "adsorption_off",
        "aw_ikawi2",
    ],
)
def test_dpt_write_supported_mutations(function_tmpdir, built, mutation, head, vals):
    m = _model(function_tmpdir, "wv", bct={})
    dpt = MfUsgDpt(m, **{**BASE, **built})
    for key, value in mutation.items():
        setattr(dpt, key, value)
    assert _written(dpt) == (head, vals)


def test_dpt_write_supported_mutation_dpf_idispim(function_tmpdir):
    # 292: DLIM read with DPF once IDISPIM is set
    m = _model(function_tmpdir, "wvd", dpf=True, bct={})
    dpt = MfUsgDpt(m, **BASE)
    dpt.idispim = 1
    assert _written(dpt) == (["0", "0", "1", "0", "1", "0", "0"], [PRS, DLIM, DDTR, C1])


# ------------------------------------------------------------ engine (requires_usgt)


def _eng_model(ws, dpf=False, bct_kw=None):
    ml = MfUsg(modelname="m", model_ws=str(ws), structured=False, exe_name=USGT_EXE)
    MfUsgDisU(
        ml,
        nodes=5,
        nlay=1,
        njag=13,
        iac=np.array([2, 3, 3, 3, 2]),
        ja=np.array([0, 1, 1, 0, 2, 2, 1, 3, 3, 2, 4, 4, 3]),
        ivsd=0,
        cl12=5.0,
        fahl=100.0,
        top=10.0,
        bot=0.0,
        area=100.0,
        nper=1,
        perlen=10.0,
        nstp=2,
        steady=False,
    )
    MfUsgBas(ml, ibound=1, strt=9.5)
    MfUsgLpf(ml, laytyp=0, hk=1.0, ss=1.0e-3)
    MfUsgSms(ml, linmeth=1)
    if dpf:
        MfUsgDpf(ml, iboundim=1, phif=0.3, ddftr=0.0, sc1im=1.0e-4)
    kw = {
        "itrnsp": 1,
        "mcomp": 1,
        "prsity": 0.2,
        "idisp": 1,
        "dl": 0.0,
        "dt": 0.0,
        "diffnc": 0.0,
        "conc": 10.0,
    }
    MfUsgBct(ml, **{**kw, **(bct_kw or {})})
    MfUsgOc(ml, stress_period_data={(0, 1): ["save head", "save concentration"]})
    return ml


def _run(ml):
    success, buff = ml.run_model(silent=True, report=True)
    assert success and any("normal termination" in b.lower() for b in buff), buff[-5:]
    ws = Path(ml.model_ws)
    return {f: (ws / f).read_bytes() for f in ("m.con", "m.cim")}


def _records(data):
    """(text, values) of unstructured single-precision CON/CIM records."""
    out, off = [], 0
    while off < len(data):
        _, _, _, _, text, n1, n2, _ = struct.unpack_from("<iiff16siii", data, off)
        off += 44
        n = n2 - n1 + 1
        out.append((text.decode().strip(), struct.unpack_from(f"<{n}f", data, off)))
        off += 4 * n
    return out


def _flopy_run(ws, dpf=False, bct_kw=None, unit=None, **dpt_kw):
    ml = _eng_model(ws, dpf, bct_kw)
    MfUsgDpt(ml, idptcon=77, unitnumber=unit, **dpt_kw)
    ml.add_output_file(77, fname="m.cim", binflag=True, package="DPT")
    ml.write_input()
    return _run(ml)


def _hand_run(ws, hand_text, dpf=False, bct_kw=None):
    ml = _eng_model(ws, dpf, bct_kw)
    MfUsgDpt(ml, filenames="manual.dpt")
    ml.add_output_file(77, fname="m.cim", binflag=True, package="DPT")
    ml.write_input(
        SelPackList=[p.name[0] for p in ml.packagelist if p.name[0] != "DPT"]
    )
    (Path(ml.model_ws) / "manual.dpt").write_text(hand_text)
    return _run(ml)


BASE_RUN = {"phif": 0.4, "prsityim": 0.2, "ddtr": 0.05, "concim": 0.0}


@requires_usgt
def test_dpt_engine_base_flopy_equals_hand(function_tmpdir):
    # default unit number (158)
    fa = _flopy_run(function_tmpdir / "flopy", **BASE_RUN)
    fb = _hand_run(
        function_tmpdir / "hand",
        "0 77 1 0 0 0 0\nCONSTANT 0.4\nCONSTANT 0.2\nCONSTANT 0.05\nCONSTANT 0.0\n",
    )
    assert fa == fb


@requires_usgt
def test_dpt_engine_dpf_flopy_equals_hand(function_tmpdir):
    fa = _flopy_run(function_tmpdir / "flopy", dpf=True, unit=158, **BASE_RUN)
    fb = _hand_run(
        function_tmpdir / "hand",
        "0 77 1 0 0 0 0\nCONSTANT 0.2\nCONSTANT 0.05\nCONSTANT 0.0\n",
        dpf=True,
    )
    assert fa == fb


@requires_usgt
@pytest.mark.parametrize("cim2", [7.0, 3.0], ids=["cim7", "cim3"])
def test_dpt_engine_imcomp_flopy_equals_hand(function_tmpdir, cim2):
    fa = _flopy_run(
        function_tmpdir / "flopy",
        bct_kw=IMCOMP_BCT,
        unit=158,
        **{**BASE_RUN, "concim": [0.0, cim2]},
    )
    fb = _hand_run(
        function_tmpdir / "hand",
        "0 77 1 0 0 0 0\nCONSTANT 0.4\nCONSTANT 0.2\nCONSTANT 0.05\n"
        f"CONSTANT 0.0\nCONSTANT {cim2}\n",
        bct_kw=IMCOMP_BCT,
    )
    assert fa == fb
    recs = _records(fa["m.cim"])
    # 2.8.1 saves one record per component (IMMOBILE CONC01/02); 2.7.0 saves
    # the mobile component only: no equivalence is forced
    if len(recs) == 2:
        assert recs[1][0] == "IMMOBILE CONC02" and np.allclose(recs[1][1], cim2)


@requires_usgt
def test_dpt_engine_sensitivity_ddtr(function_tmpdir):
    out = [
        _flopy_run(function_tmpdir / f"s{i}", unit=158, **{**BASE_RUN, "ddtr": v})
        for i, v in enumerate((0.05, 0.0))
    ]
    assert out[0]["m.con"] != out[1]["m.con"] and out[0]["m.cim"] != out[1]["m.cim"]
    assert np.allclose(_records(out[1]["m.cim"])[0][1], 0.0)  # no exchange


@requires_usgt
def test_dpt_engine_sensitivity_concim_mobile(function_tmpdir):
    out = [
        _flopy_run(function_tmpdir / f"c{i}", unit=158, **{**BASE_RUN, "concim": v})
        for i, v in enumerate((0.0, 5.0))
    ]
    assert out[0]["m.con"] != out[1]["m.con"] and out[0]["m.cim"] != out[1]["m.cim"]


@requires_usgt
def test_dpt_engine_sensitivity_concim_immobile(function_tmpdir):
    out = [
        _flopy_run(
            function_tmpdir / f"i{i}",
            bct_kw=IMCOMP_BCT,
            unit=158,
            **{**BASE_RUN, "concim": [0.0, v]},
        )
        for i, v in enumerate((7.0, 3.0))
    ]
    # the immobile component does not exchange with the mobile domain
    assert out[0]["m.con"] == out[1]["m.con"]
    recs = [_records(o["m.cim"]) for o in out]
    if len(recs[0]) == 2:  # 2.8.1 writes the immobile component
        assert recs[0][0] == recs[1][0] and recs[0][1] != recs[1][1]
