"""Transport data of the CLN nodes (MfUsgCln) for USG-T 2.8.1.

USG-T 2.8.1 reads them in SGWT2BCT1CLN (glo2btnu1.f :2034-2165), after the
CLN IBOUND and STRT: ICBUND (IC_IBOUND_FLG=0), DLX and ATXY (IDISP not 0),
then for each of the NTCOMP components (:2101) MCOMP SPTLRCT arrays
(SPATIALREACT, :2107-2118), ZODRW (IZOD 1 or 3, :2124) and FODRW (IFOD 1 or 3,
:2134) for the mobile species only (:2122-2123), and CONC (TEMP for the
temperature, :2146). Components: species, temperature (IHEAT=1), immobile.

The hand-written CLN files below use the FloPy flow part (written while no
BCT is registered) followed by a transport part written by hand from those
lines. Engine tests (requires_usgt) compare FloPy files with hand-written
files and check that the results depend on the data; they use a one-layer
DISU grid.
"""

import builtins
import io
import struct
from pathlib import Path

import numpy as np
import pytest

from autotest.usgt_exe import USGT_EXE, requires_usgt
from flopy.mfusg import (
    MfUsg,
    MfUsgBas,
    MfUsgBct,
    MfUsgCln,
    MfUsgDisU,
    MfUsgLpf,
    MfUsgOc,
    MfUsgSms,
)
from flopy.modflow import ModflowChd
from flopy.utils import Util2d

HEAT = {
    "iheat": 1,
    "htcondw": 0.58,
    "rhow": 1000.0,
    "htcapw": 4186.0,
    "htcaps": 800.0,
    "htconds": 2.0,
    "heat": 20.0,
    "bulkd": 1600.0,
}
IMCOMP = {"imcomp": 1, "imconc": 2.0}
ADS = {"iadsorb": 1, "adsorb": 0.0, "bulkd": 1.6}
AW = {
    "aw_adsorb": True,
    "iarea_fn": 1,
    "ikawi_fn": 1,
    "awamax": 100.0,
    "alangaw": 0.0,
    "blangaw": 0.0,
}
BCT_KW = {
    "itrnsp": 1,
    "mcomp": 1,
    "idisp": 1,
    "dl": 0.0,
    "dt": 0.0,
    "diffnc": 0.0,
    "prsity": 0.2,
    "conc": 1.0,
}
NCLN = 2  # CLN nodes


def _model(ws, bct=None, exe=None, chain=False):
    """5-node DISU line (one layer), CHD at both ends; BCT when bct is a dict.

    chain=True marks the BCT with CHAINDECAY after construction: the chain
    data of MfUsgBct are not used by CLN (they are being changed in another
    lane), only the ``chaindecay`` flag is."""
    ml = MfUsg(
        modelname="m", model_ws=str(ws), structured=False, exe_name=exe or "mfusgt"
    )
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
        nstp=5,
        steady=False,
    )
    MfUsgBas(ml, ibound=1, strt=9.5)
    if bct is not None:
        b = MfUsgBct(ml, **{**BCT_KW, **bct})
        if chain:
            b.chaindecay = 1
    return ml


def _cln(ml, **kw):
    kw = {"dll": 0.0, "dlm": 0.0, **kw}  # DLX and ATXY, read with IDISP not 0
    return MfUsgCln(
        ml,
        ncln=1,
        iclnnds=-1,
        nndcln=NCLN,
        nclngwc=2,
        node_prop=[
            [1, 1, 0, 10.0, 0.5, 1.57, 0, 0],
            [2, 1, 0, 10.0, 0.5, 1.57, 0, 0],
        ],
        cln_gwc=[[1, 1, 0, 0.0, 10.0, 1.0, 0], [2, 5, 0, 0.0, 10.0, 1.0, 0]],
        cln_circ=[[1, 0.5, 3.23e10]],
        strt=9.5,
        unitnumber=[71, 0, 0, 0, 0, 37, 0, 0],
        **kw,
    )


def _flow_text(ws):
    """CLN file written by FloPy without BCT: flow part only."""
    ws = Path(ws) / "flow_part"
    cln = _cln(_model(ws))
    cln.write_file()
    return Path(cln.fn_path).read_text()


def _hand(arrays):
    """Transport part written by hand: CONSTANT for a scalar, else INTERNAL."""
    out = []
    for label, value in arrays:
        if np.isscalar(value):
            out.append(f"CONSTANT {value!r}   #{label}\n")
        else:
            out.append(f"INTERNAL 1.0 (FREE) -1   #{label}\n")
            out.append(" ".join(repr(float(v)) for v in value) + "\n")
    return "".join(out)


def _transport(text):
    """(label, values) of each array after the CLN STRT record."""
    lines = text.splitlines()
    i = next(k for k, ln in enumerate(lines) if "#strt" in ln) + 1
    out = []
    while i < len(lines):
        tok = lines[i].split()
        if not tok:
            i += 1
            continue
        label = lines[i].split("#", 1)[1].split()[0] if "#" in lines[i] else ""
        if tok[0].upper() == "CONSTANT":
            out.append((label, [float(tok[1])] * NCLN))
            i += 1
        elif tok[0].upper() == "INTERNAL":
            mult, vals = float(tok[1]), []
            i += 1
            while len(vals) < NCLN:
                vals += [mult * float(v) for v in lines[i].split()]
                i += 1
            out.append((label, vals))
        else:
            raise AssertionError(f"unexpected line {lines[i]!r}")
    return out


def _expect(arrays):
    return [
        (lab, [float(v)] * NCLN if np.isscalar(v) else list(v)) for lab, v in arrays
    ]


def _written(cln):
    cln.write_file()
    return _transport(Path(cln.fn_path).read_text())


DISP = [("dll", 0.0), ("dlm", 0.0)]
ZOD = {"zodrw": 1.0e-3, "zodrs": 0.0, "zodraw": 0.0}
FOD = {"fodrw": 2.0e-3, "fodrs": 0.0, "fodraw": 0.0}


def _decay_bct(which, value):
    """BCT keywords for IZOD or IFOD equal to value, with the adsorption that
    MfUsgBct needs for the arrays of that value."""
    kw = {which: value, **(ZOD if which == "izod" else FOD)}
    if value in (2, 3):
        kw.update(ADS)
    if value == 4:
        kw.update(AW)
    return kw


# (bct keywords, chain flag, cln keywords, expected transport arrays)
CASES = {
    "mcomp1": ({}, False, {"conc": [5.0]}, DISP + [("conc", 5.0)]),
    "mcomp1_icbund": (
        {"icbndflg": 0},
        False,
        {"icbund": 1, "conc": [5.0]},
        [("icbund", 1)] + DISP + [("conc", 5.0)],
    ),
    "mcomp2": (
        {"mcomp": 2, "conc": [1.0, 0.0]},
        False,
        {"conc": [5.0, 6.0]},
        DISP + [("conc", 5.0), ("conc", 6.0)],
    ),
    "heat": (
        HEAT,
        False,
        {"conc": [5.0, 25.0]},
        DISP + [("conc", 5.0), ("conc", 25.0)],
    ),
    "imcomp": (
        IMCOMP,
        False,
        {"conc": [5.0, 3.0]},
        DISP + [("conc", 5.0), ("conc", 3.0)],
    ),
    "heat_imcomp": (
        {**HEAT, **IMCOMP},
        False,
        {"conc": [5.0, 25.0, 3.0]},
        DISP + [("conc", 5.0), ("conc", 25.0), ("conc", 3.0)],
    ),
    "decay_heat": (
        {**HEAT, "izod": 1, "ifod": 1, "zodrw": 1.0e-3, "fodrw": 2.0e-3},
        False,
        {"zodrw": 0.01, "fodrw": 0.02, "conc": [5.0, 25.0]},
        DISP + [("zodrw", 0.01), ("fodrw", 0.02), ("conc", 5.0), ("conc", 25.0)],
    ),
    "sptlrct_mcomp2": (
        {"mcomp": 2, "conc": [1.0, 0.0], "spatialreact": True},
        True,
        {"sptlrct": [[0.0, 0.0], [[0.5, 0.5], [0.25, 0.25]]], "conc": [5.0, 6.0]},
        DISP
        + [("sptlrct", 0.0), ("sptlrct", 0.0), ("conc", 5.0)]
        + [("sptlrct", 0.5), ("sptlrct", 0.25), ("conc", 6.0)],
    ),
    "sptlrct_heat": (
        {"mcomp": 2, "conc": [1.0, 0.0], "spatialreact": True, **HEAT},
        True,
        {
            "sptlrct": [[0.0, 0.0], [[0.5, 0.5], [0.25, 0.25]], 0.0],
            "conc": [5.0, 6.0, 25.0],
        },
        DISP
        + [("sptlrct", 0.0), ("sptlrct", 0.0), ("conc", 5.0)]
        + [("sptlrct", 0.5), ("sptlrct", 0.25), ("conc", 6.0)]
        + [("sptlrct", 0.0), ("sptlrct", 0.0), ("conc", 25.0)],
    ),
}
for _w, _lab in (("izod", "zodrw"), ("ifod", "fodrw")):
    for _v in range(5):
        _arr = [(_lab, 0.01)] if _v in (1, 3) else []
        CASES[f"{_w}{_v}"] = (
            _decay_bct(_w, _v),
            False,
            {_lab: 0.01, "conc": [5.0]},
            DISP + _arr + [("conc", 5.0)],
        )


# ------------------------------------------------------------------ write


@pytest.mark.parametrize("case", list(CASES))
def test_cln_transport_written_as_engine_reads(function_tmpdir, case):
    bct, chain, kw, arrays = CASES[case]
    ml = _model(function_tmpdir, bct=bct, chain=chain)
    assert _written(_cln(ml, **kw)) == _expect(arrays)


@pytest.mark.parametrize("case", [c for c in CASES if c != "mcomp1"])
def test_cln_transport_hand_file_loads_and_rewrites(function_tmpdir, case):
    """Control independent of the writer: a hand-written file is loaded
    component by component and written back with the same arrays."""
    bct, chain, _, arrays = CASES[case]
    path = Path(function_tmpdir) / "hand.cln"
    path.write_text(_flow_text(function_tmpdir) + _hand(arrays))
    ml = _model(function_tmpdir / "load", bct=bct, chain=chain)
    cln = MfUsgCln.load(str(path), ml)
    ncomp = ml.ntcomp
    assert len(cln.conc) == ncomp
    concs = [a for a in _expect(arrays) if a[0] == "conc"]
    for icomp in range(ncomp):
        assert np.allclose(cln.conc[icomp].array, concs[icomp][1])
    assert _written(cln) == _expect(arrays)


def test_cln_transport_scalar_conc_for_all_components(function_tmpdir):
    ml = _model(function_tmpdir, bct={**HEAT, **IMCOMP})
    got = _written(_cln(ml, conc=4.0))
    assert got == _expect(DISP + [("conc", 4.0)] * 3)


def test_cln_transport_mcomp1_old_forms_kept(function_tmpdir):
    """Control: the MCOMP=1 forms used so far (tutorial03) give the same file."""
    ml = _model(function_tmpdir, bct={"spatialreact": True}, chain=True)
    got = _written(_cln(ml, sptlrct=[0.5], conc=[[0.0, 1.0]]))
    assert got == _expect(DISP + [("sptlrct", 0.5), ("conc", [0.0, 1.0])])


def test_cln_transport_longer_conc_extra_ignored(function_tmpdir):
    """Control (behaviour kept from upstream): extra entries are not written."""
    ml = _model(function_tmpdir, bct={})
    assert _written(_cln(ml, conc=[5.0, 99.0])) == _expect(DISP + [("conc", 5.0)])


# SPTLRCT entries of one component: one array (scalar, array or flat list,
# as upstream reads it) used for all parent positions, or a list of arrays,
# one per parent position (exactly MCOMP).
SPR2 = {"mcomp": 2, "conc": [1.0, 0.0], "spatialreact": True}


def test_cln_sptlrct_per_node_mcomp1_kept(function_tmpdir):
    """Upstream form: a flat list is one array per node."""
    ml = _model(function_tmpdir, bct={"spatialreact": True}, chain=True)
    got = _written(_cln(ml, sptlrct=[[0.1, 0.2]], conc=[5.0]))
    assert got == _expect(DISP + [("sptlrct", [0.1, 0.2]), ("conc", 5.0)])


def test_cln_sptlrct_per_node_used_for_all_parents(function_tmpdir):
    ml = _model(function_tmpdir, bct=SPR2, chain=True)
    got = _written(_cln(ml, sptlrct=[[0.1, 0.2], (0.3, 0.4)], conc=[5.0, 6.0]))
    assert got == _expect(
        DISP
        + [("sptlrct", [0.1, 0.2])] * 2
        + [("conc", 5.0)]
        + [("sptlrct", [0.3, 0.4])] * 2
        + [("conc", 6.0)]
    )


def test_cln_sptlrct_one_array_per_parent(function_tmpdir):
    ml = _model(function_tmpdir, bct=SPR2, chain=True)
    sptlrct = [[[0.1, 0.2], [0.3, 0.4]], [np.array([0.5, 0.6]), 0.7]]
    got = _written(_cln(ml, sptlrct=sptlrct, conc=[5.0, 6.0]))
    assert got == _expect(
        DISP
        + [("sptlrct", [0.1, 0.2]), ("sptlrct", [0.3, 0.4]), ("conc", 5.0)]
        + [("sptlrct", [0.5, 0.6]), ("sptlrct", 0.7), ("conc", 6.0)]
    )


@pytest.mark.parametrize("nparents", [1, 3], ids=["fewer", "more"])
def test_cln_sptlrct_parent_count_must_be_mcomp(function_tmpdir, nparents):
    ml = _model(function_tmpdir, bct=SPR2, chain=True)
    with pytest.raises(ValueError, match=r"sptlrct\[1\]"):
        _cln(ml, sptlrct=[0.0, [[0.1, 0.2]] * nparents], conc=[5.0, 6.0])
    assert ml.get_package("CLN") is None


def test_cln_sptlrct_load_write_keeps_parent_arrays(function_tmpdir):
    arrays = DISP + [
        ("sptlrct", [0.1, 0.2]),
        ("sptlrct", [0.3, 0.4]),
        ("conc", 5.0),
        ("sptlrct", [0.5, 0.6]),
        ("sptlrct", [0.7, 0.8]),
        ("conc", 6.0),
    ]
    path = Path(function_tmpdir) / "hand.cln"
    path.write_text(_flow_text(function_tmpdir) + _hand(arrays))
    ml = _model(function_tmpdir / "load", bct=SPR2, chain=True)
    cln = MfUsgCln.load(str(path), ml)
    assert [len(s) for s in cln.sptlrct] == [2, 2]
    assert all(isinstance(u, Util2d) for s in cln.sptlrct for u in s)
    assert np.allclose(cln.sptlrct[1][1].array.ravel(), [0.7, 0.8])
    assert _written(cln) == _expect(arrays)


def test_cln_per_node_conc_and_decay_kept(function_tmpdir):
    """Control: a per-node list in conc, zodrw or fodrw keeps its meaning."""
    bct = {"izod": 1, "zodrw": 1.0e-3, "ifod": 1, "fodrw": 2.0e-3}
    ml = _model(function_tmpdir, bct=bct)
    kw = {"zodrw": [[0.01, 0.02]], "fodrw": [[0.03, 0.04]], "conc": [[0.1, 0.2]]}
    got = _written(_cln(ml, **kw))
    assert got == _expect(
        DISP + [("zodrw", [0.01, 0.02]), ("fodrw", [0.03, 0.04]), ("conc", [0.1, 0.2])]
    )


@pytest.mark.parametrize(
    "bct,chain,kw,match",
    [
        (HEAT, False, {"conc": [[5.0, 5.0]]}, "conc"),
        (HEAT, False, {}, "conc"),
        (IMCOMP, False, {"conc": [5.0]}, "conc"),
        (
            {"mcomp": 2, "conc": [1.0, 0.0], "izod": 1, "zodrw": 1.0e-3},
            False,
            {"zodrw": [0.01], "conc": [5.0, 6.0]},
            "zodrw",
        ),
        (
            {"mcomp": 2, "conc": [1.0, 0.0], "ifod": 3, **ADS, **FOD},
            False,
            {"conc": [5.0, 6.0]},
            "fodrw",
        ),
        (
            {"mcomp": 2, "conc": [1.0, 0.0], "spatialreact": True},
            True,
            {"sptlrct": [[[0.0, 0.0]], [[0.5, 0.5]]], "conc": [5.0, 6.0]},
            "sptlrct",
        ),
        (
            {"mcomp": 2, "conc": [1.0, 0.0], "spatialreact": True, **HEAT},
            True,
            {"sptlrct": [[0.0, 0.0], [0.5, 0.0]], "conc": [5.0, 6.0, 25.0]},
            "sptlrct",
        ),
    ],
    ids=[
        "heat_conc_short",
        "heat_conc_missing",
        "imcomp_conc_short",
        "zodrw_short",
        "fodrw_missing",
        "sptlrct_parents_short",
        "sptlrct_components_short",
    ],
)
def test_cln_transport_short_input_rejected(function_tmpdir, bct, chain, kw, match):
    ml = _model(function_tmpdir, bct=bct, chain=chain)
    with pytest.raises(ValueError, match=match):
        _cln(ml, **kw)
    assert ml.get_package("CLN") is None


# ------------------------------------------------------- write-time checks


def _unsynced_bct(ml):
    """Register a BCT that did not set the model counts."""
    ml.add_package(MfUsgBct(ml, **{**BCT_KW, "mcomp": 2}, add_package=False))


CHANGES = {
    "mcomp": lambda ml: MfUsgBct(ml, **{**BCT_KW, "mcomp": 2, "conc": [1.0, 0.0]}),
    "iheat": lambda ml: MfUsgBct(ml, **{**BCT_KW, **HEAT}),
    "imcomp": lambda ml: MfUsgBct(ml, **{**BCT_KW, **IMCOMP}),
    "izod": lambda ml: MfUsgBct(ml, **{**BCT_KW, "izod": 1, "zodrw": 1.0e-3}),
    "ifod": lambda ml: MfUsgBct(ml, **{**BCT_KW, "ifod": 3, **ADS, **FOD}),
    "idisp": lambda ml: MfUsgBct(ml, **{**BCT_KW, "idisp": 0}),
    "icbndflg": lambda ml: MfUsgBct(ml, **{**BCT_KW, "icbndflg": 0}),
    "unsynced": _unsynced_bct,
}


@pytest.mark.parametrize("change", list(CHANGES))
def test_cln_write_rechecks_bct_layout(function_tmpdir, change):
    ml = _model(function_tmpdir, bct={})
    cln = _cln(ml, conc=[5.0])
    cln.write_file()
    before = Path(cln.fn_path).read_bytes()
    CHANGES[change](ml)
    with pytest.raises(ValueError, match="BCT"):
        cln.write_file()
    assert Path(cln.fn_path).read_bytes() == before


def test_cln_write_bct_added_after_cln_rejected(function_tmpdir):
    ml = _model(function_tmpdir)
    cln = _cln(ml)
    MfUsgBct(ml, **BCT_KW)
    with pytest.raises(ValueError, match="BCT"):
        cln.write_file()
    assert not Path(cln.fn_path).exists()


def test_cln_write_bct_replaced_same_layout(function_tmpdir):
    """Control: an equivalent BCT keeps the CLN package valid."""
    ml = _model(function_tmpdir, bct={})
    cln = _cln(ml, conc=[5.0])
    MfUsgBct(ml, **BCT_KW)
    assert _written(cln) == _expect(DISP + [("conc", 5.0)])


def test_cln_write_bct_removed_no_transport(function_tmpdir):
    """Control: without BCT in the name file USG-T reads no CLN transport."""
    ml = _model(function_tmpdir, bct={})
    cln = _cln(ml, conc=[5.0])
    ml.remove_package("BCT")
    assert _written(cln) == []


@pytest.mark.parametrize("how", ["construction", "load"])
def test_cln_unsynced_bct_rejected(function_tmpdir, how):
    ml = _model(function_tmpdir)
    _unsynced_bct(ml)
    hand = _flow_text(function_tmpdir) + _hand(DISP + [("conc", 5.0)] * 2)
    with pytest.raises(ValueError, match="BCT"):
        if how == "construction":
            _cln(ml, conc=[5.0, 6.0])
        else:
            MfUsgCln.load(io.StringIO(hand), ml)
    assert ml.get_package("CLN") is None


# ----------------------------------------------------------------- guards
# Probes P1 and P2 of lane O3-C2 (USG-T 2.8.1 and 2.7.0, gfortran builds):
# SIGSEGV when CLN reads SPTLRCT without CHAINDECAY (glo2btnu1.f :412-419,
# :2107) or ZODRW/FODRW with IZOD/IFOD=3 and IADSORB=0 (:305-355, :2124,
# :2134), because those arrays are not allocated.

GUARDS = {
    "spatialreact_no_chain": (
        {"spatialreact": True},
        {"sptlrct": 0.5},
        DISP + [("sptlrct", 0.5), ("conc", 5.0)],
        "SPATIALREACT",
    ),
    "izod3_no_adsorption": (
        {"izod": 3, **ZOD},
        {"zodrw": 0.01},
        DISP + [("zodrw", 0.01), ("conc", 5.0)],
        "IZOD",
    ),
    "ifod3_no_adsorption": (
        {"ifod": 3, **FOD},
        {"fodrw": 0.01},
        DISP + [("fodrw", 0.01), ("conc", 5.0)],
        "IFOD",
    ),
}


@pytest.mark.parametrize("case", list(GUARDS))
def test_cln_guard_construction(function_tmpdir, case):
    bct, kw, _, match = GUARDS[case]
    ml = _model(function_tmpdir, bct=bct)
    with pytest.raises(NotImplementedError, match=match):
        _cln(ml, conc=[5.0], **kw)
    assert ml.get_package("CLN") is None


@pytest.mark.parametrize("case", list(GUARDS))
def test_cln_guard_write_keeps_file(function_tmpdir, case):
    bct, kw, _, match = GUARDS[case]
    ml = _model(function_tmpdir, bct={})
    cln = _cln(ml, conc=[5.0])
    cln.write_file()
    before = Path(cln.fn_path).read_bytes()
    MfUsgBct(ml, **{**BCT_KW, **bct})
    with pytest.raises(NotImplementedError, match=match):
        cln.write_file()
    assert Path(cln.fn_path).read_bytes() == before


@pytest.mark.parametrize("case", list(GUARDS))
def test_cln_guard_load_before_opening(function_tmpdir, case, monkeypatch):
    bct, _, arrays, match = GUARDS[case]
    path = Path(function_tmpdir) / "hand.cln"
    path.write_text(_flow_text(function_tmpdir) + _hand(arrays))
    ml = _model(function_tmpdir / "load", bct=bct)
    opened = []
    real_open = builtins.open

    def tracking_open(file, *args, **kwargs):
        opened.append(Path(file).name if isinstance(file, (str, Path)) else file)
        return real_open(file, *args, **kwargs)

    monkeypatch.setattr(builtins, "open", tracking_open)
    with pytest.raises(NotImplementedError, match=match):
        MfUsgCln.load(str(path), ml)
    monkeypatch.undo()
    assert path.name not in opened
    # a handle of the caller is left open and unread
    with open(path) as f:
        with pytest.raises(NotImplementedError, match=match):
            MfUsgCln.load(f, ml)
        assert not f.closed
        assert f.tell() == 0
    assert ml.get_package("CLN") is None


# ------------------------------------------------------ engine (requires_usgt)


def _eng_model(ws, bct, cln_kw, chain=False):
    ml = _model(ws, exe=USGT_EXE)
    MfUsgLpf(ml, laytyp=0, hk=1.0, ss=1.0e-3)
    MfUsgSms(ml, linmeth=1)
    ModflowChd(ml, stress_period_data={0: [[0, 10.0, 10.0], [4, 9.0, 9.0]]})
    b = MfUsgBct(ml, **{**BCT_KW, **bct})
    if chain:
        b.chaindecay = 1
    cln = _cln(ml, **cln_kw)
    MfUsgOc(ml, stress_period_data={(0, 4): ["save head", "save concentration"]})
    return ml, cln


def _records(data, cln):
    """Labels and values of the concentration records: GW (unstructured
    header: first node, last node, layer) or CLN (NCLNNDS, 1, 1)."""
    out, off = [], 0
    while off < len(data):
        a, b = struct.unpack_from("<ii", data, off + 32)
        n = a * b if cln else b - a + 1
        label = data[off + 16 : off + 32].decode().strip()
        out.append((label, np.array(struct.unpack_from(f"<{n}f", data, off + 44))))
        off += 44 + 4 * n
    return out


def _run(ml, skip=(), texts=None):
    ml.write_input(
        SelPackList=[p.name[0] for p in ml.packagelist if p.name[0] not in skip]
    )
    ws = Path(ml.model_ws)
    for name, text in (texts or {}).items():
        (ws / name).write_text(text)
    success, buff = ml.run_model(silent=True, report=True)
    assert success and any("normal termination" in str(b).lower() for b in buff)
    return (
        _records((ws / "m.clncn").read_bytes(), cln=True),
        _records((ws / "m.con").read_bytes(), cln=False),
    )


def _flopy_run(ws, bct, cln_kw, chain=False, bct_text=None):
    """Returns the CLN and GW records and the transport arrays FloPy wrote."""
    ml, cln = _eng_model(ws, bct, cln_kw, chain=chain)
    if bct_text is None:
        out = _run(ml)
    else:
        out = _run(ml, skip=("BCT",), texts={"m.bct": bct_text})
    return out + (_transport(Path(cln.fn_path).read_text()),)


def _hand_run(ws, bct, cln_kw, arrays, chain=False, bct_text=None):
    """Same model; CLN file = FloPy flow part + transport part by hand."""
    ml, _ = _eng_model(ws, bct, cln_kw, chain=chain)
    texts = {"m.cln": _flow_text(ws) + _hand(arrays)}
    skip = ["CLN"]
    if bct_text is not None:
        texts["m.bct"] = bct_text
        skip.append("BCT")
    return _run(ml, skip=skip, texts=texts)


def _same(a, b):
    return [x[0] for x in a] == [x[0] for x in b] and all(
        np.array_equal(x[1], y[1]) for x, y in zip(a, b)
    )


# (bct keywords, cln keywords, hand-written arrays, changed cln keywords)
ENGINE = {
    "heat": (
        HEAT,
        {"conc": [5.0, 25.0]},
        DISP + [("conc", 5.0), ("conc", 25.0)],
        {"conc": [5.0, 40.0]},
    ),
    "imcomp": (
        IMCOMP,
        {"conc": [5.0, 3.0]},
        DISP + [("conc", 5.0), ("conc", 3.0)],
        {"conc": [5.0, 4.0]},
    ),
    "ifod1": (
        {"ifod": 1, "fodrw": 0.01},
        {"fodrw": 0.2, "conc": [5.0]},
        DISP + [("fodrw", 0.2), ("conc", 5.0)],
        {"fodrw": 0.0, "conc": [5.0]},
    ),
    "ifod1_heat": (
        {**HEAT, "ifod": 1, "fodrw": 0.01},
        {"fodrw": 0.2, "conc": [5.0, 25.0]},
        DISP + [("fodrw", 0.2), ("conc", 5.0), ("conc", 25.0)],
        {"fodrw": 0.0, "conc": [5.0, 25.0]},
    ),
    "izod3_adsorption": (
        {"izod": 3, **ADS, **ZOD},
        {"zodrw": 0.05, "conc": [5.0]},
        DISP + [("zodrw", 0.05), ("conc", 5.0)],
        {"zodrw": 0.0, "conc": [5.0]},
    ),
}


@requires_usgt
@pytest.mark.parametrize("case", list(ENGINE))
def test_cln_engine_flopy_equals_hand_and_sensitivity(function_tmpdir, case):
    bct, kw, arrays, changed = ENGINE[case]
    flopy_out = _flopy_run(function_tmpdir / "flopy", bct, kw)
    assert flopy_out[2] == _expect(arrays)
    hand_out = _hand_run(function_tmpdir / "hand", bct, kw, arrays)
    assert _same(flopy_out[0], hand_out[0]) and _same(flopy_out[1], hand_out[1])
    ntcomp = len([a for a in arrays if a[0] == "conc"])
    assert len(flopy_out[0]) == ntcomp  # one CLN record per component
    other = _flopy_run(function_tmpdir / "changed", bct, changed)
    assert not _same(flopy_out[0], other[0])


# Hand-written BCT for SPTLRCT (FloPy cannot write CHAINDECAY with SPATIALREACT
# yet): glo2btnu1.f item 1 (:146-170) with the options on the same line
# (:64-141), PRSITY (:1562), ANGLEX (:1594), DL and DT (:1626, :1637), then
# per species (:1732): NPARENT (:1744), JPARENT (:1751) and one SPTLRCT per
# parent and layer (:1765-1778), FODRW (IFOD=1, :1918-1923) and CONC.
# Species 2 is the daughter of species 1.
CHAIN_BCT = (
    "# BCT written by hand\n"
    " 1 0 2 1 1 0 0 -999.0 1.0e-6 1 0 0.0 0 1 0 0 0 0 0 CHAINDECAY SPATIALREACT\n"
    "CONSTANT 0.2   PRSITY\n"
    "CONSTANT 0.0   ANGLEX\n"
    "CONSTANT 0.0   DL\n"
    "CONSTANT 0.0   DT\n"
    "0   NPARENT species 1\n"
    "CONSTANT 0.1   FODRW species 1\n"
    "CONSTANT 1.0   CONC species 1\n"
    "1   NPARENT species 2\n"
    "CONSTANT 1   JPARENT species 2\n"
    "CONSTANT 0.5   SPTLRCT species 2 parent 1\n"
    "CONSTANT 0.0   FODRW species 2\n"
    "CONSTANT 0.0   CONC species 2\n"
)
CHAIN_BCT_KW = {
    "mcomp": 2,
    "ifod": 1,
    "fodrw": [0.1, 0.0],
    "conc": [1.0, 0.0],
    "spatialreact": True,
}


@requires_usgt
def test_cln_engine_sptlrct_end_to_end(function_tmpdir):
    kw = {
        "sptlrct": [0.0, [[0.8, 0.8], [0.0, 0.0]]],
        "fodrw": [0.1, 0.0],
        "conc": [5.0, 0.0],
    }
    arrays = DISP + [
        ("sptlrct", 0.0),
        ("sptlrct", 0.0),
        ("fodrw", 0.1),
        ("conc", 5.0),
        ("sptlrct", 0.8),
        ("sptlrct", 0.0),
        ("fodrw", 0.0),
        ("conc", 0.0),
    ]
    args = (CHAIN_BCT_KW, kw)
    flopy_out = _flopy_run(function_tmpdir / "flopy", *args, True, CHAIN_BCT)
    assert flopy_out[2] == _expect(arrays)
    hand_out = _hand_run(function_tmpdir / "hand", *args, arrays, True, CHAIN_BCT)
    assert _same(flopy_out[0], hand_out[0]) and _same(flopy_out[1], hand_out[1])
    assert [x[0] for x in flopy_out[0]] == ["CONC01", "CONC02"]
    # the daughter in the CLN depends on the CLN generation factor
    kw0 = {**kw, "sptlrct": [[0.0, 0.0], [0.0, 0.0]]}
    other = _flopy_run(function_tmpdir / "changed", CHAIN_BCT_KW, kw0, True, CHAIN_BCT)
    assert _same(flopy_out[0][:1], other[0][:1])
    assert not np.array_equal(flopy_out[0][1][1], other[0][1][1])
