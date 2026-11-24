"""MfUsgBct CHAINDECAY, SPATIALREACT, MULTIFILE and REACTION1 (USG-T 2.8.1).

Contract, from glo2btnu1.f (GWT2BCT1AR and the two readers SGWT2BCT1S,
structured, and SGWT2BCT1G, unstructured); line numbers are 2.8.1
(structured / unstructured):

- CHAINDECAY (82-85): first item of every species, NPARENT in a list-directed
  record (1158 / 1744); with NPARENT > 0, JPARENT with NPARENT integers
  (1165 / 1751) and STOTIO with NPARENT reals (1174 / 1760).
- SPATIALREACT (86-89): with CHAINDECAY, instead of STOTIO, for each parent
  one array per layer (1179-1193 / 1765-1778), used per node by GWT2BCT1GEN
  (3936-3939).
- MULTIFILE (68-72): the root name is read in upper case (URWORD,
  utl7u1.f 89-95); the engine writes ROOT//' CONC01'//'.CON' (214, GET_TEXTC).
- REACTION1 (133-137): GWT2RXN1AR reads its block from the BCT file before
  the arrays (431, 459-463); not implemented, guard on load and construction.

The structured reader of USG-T 2.8.1 and 2.7.0 (local gfortran builds) stops
with a segmentation fault on CHAINDECAY with NPARENT > 0 (JPARENT is read into
ITEMP, allocated only while ICBUND is read, 912-923, and used at 1165), so the
structured cases are checked for input/output only, and writing them warns.

Expected files are written by hand from this contract, never by the FloPy
writer, and are read back with the reader of test_usgt_bct_b1.py. Engine
tests (requires_usgt) use the 5-node DISU model of test_usgt_bct_b1.py.
"""

import io
import struct
import warnings
from pathlib import Path

import numpy as np
import pytest

from autotest.test_usgt_bct_b1 import (
    PRS,
    SENTINEL,
    _consts,
    _hand,
    _head,
    _model,
    _parse,
    _state,
    _written,
)
from autotest.usgt_exe import USGT_EXE, requires_usgt
from flopy.mfusg import (
    MfUsg,
    MfUsgBas,
    MfUsgBct,
    MfUsgDisU,
    MfUsgLpf,
    MfUsgOc,
    MfUsgSms,
)
from flopy.modflow import ModflowDis

GRIDS = pytest.mark.parametrize("structured", [True, False], ids=["dis", "disu"])
CHAIN = " CHAINDECAY"
SPR = " CHAINDECAY SPATIALREACT"
INT_I = "INTERNAL 1 (FREE) -1\n"
INT_R = "INTERNAL 1.0 (FREE) -1\n"


def _bct(ml, **kw):
    base = {
        "mcomp": 3,
        "itvd": 0,
        "idisp": 0,
        "prsity": PRS,
        "conc": [21.0, 22.0, 23.0],
    }
    return MfUsgBct(ml, **{**base, **kw})


def _quiet_write(bct, ws, name="written.bct"):
    """Write, ignoring the structured CHAINDECAY warning (tested apart)."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        return _written(bct, ws, name)


# ------------------------------------------------------------------ CHAINDECAY

# species 1 without parents, 2 with one parent, 3 with parents 1 and 2
CHAIN_KW = {
    "chaindecay": True,
    "nparent": [0, 1, 2],
    "jparent": [0, [1], [1, 2]],
    "stotio": [0, [0.4], [0.25, 0.75]],
}
HAND_CHAIN = (
    _head(mcomp=3, opts=CHAIN)
    + _consts(PRS)
    + "0\n"
    + _consts(21.0)
    + "1\n"
    + INT_I
    + " 1\n"
    + INT_R
    + " 0.4\n"
    + _consts(22.0)
    + "2\n"
    + INT_I
    + " 1 2\n"
    + INT_R
    + " 0.25 0.75\n"
    + _consts(23.0)
)


@GRIDS
def test_bct_chaindecay_written(function_tmpdir, structured):
    bct = _bct(_model(function_tmpdir, structured=structured), **CHAIN_KW)
    assert _parse(_quiet_write(bct, function_tmpdir)) == _parse(HAND_CHAIN)


@GRIDS
def test_bct_chaindecay_loaded(function_tmpdir, structured):
    ml = _model(function_tmpdir, structured=structured)
    bct = MfUsgBct.load(_hand(function_tmpdir, HAND_CHAIN), ml)
    assert list(bct.nparent) == [0, 1, 2]
    assert bct.jparent[0] == 0 and bct.stotio[0] == 0
    assert bct.jparent[1].array.tolist() == [1]
    assert bct.jparent[2].array.tolist() == [1, 2]
    assert np.allclose(bct.stotio[1].array, [0.4])
    assert np.allclose(bct.stotio[2].array, [0.25, 0.75])
    assert _parse(_quiet_write(bct, function_tmpdir)) == _parse(HAND_CHAIN)


def test_bct_chaindecay_scalar_per_species(function_tmpdir):
    """One parent per species given as scalars (07c notebook, Ex7 file form)."""
    kw = {
        "mcomp": 4,
        "conc": 0.0,
        "chaindecay": True,
        "nparent": [0, 1, 1, 1],
        "jparent": [0, 1, 2, 3],
        "stotio": [0, 0.795, 0.737, 0.640],
    }
    bct = _bct(_model(function_tmpdir, structured=False), **kw)
    hand = _head(mcomp=4, opts=CHAIN) + _consts(PRS) + "0\n" + _consts(0.0)
    for parent, s in ((1, 0.795), (2, 0.737), (3, 0.640)):
        hand += "1\n" + _consts(parent, s, 0.0)
    assert _parse(_written(bct, function_tmpdir)) == _parse(hand)


def _without(kw, name):
    return {k: v for k, v in kw.items() if k != name}


@pytest.mark.parametrize(
    "kw, match",
    [
        (_without(CHAIN_KW, "nparent"), "nparent"),
        ({**CHAIN_KW, "nparent": [0, 1]}, "nparent"),
        ({**CHAIN_KW, "nparent": [0, 1, -1]}, "nparent"),
        ({**CHAIN_KW, "nparent": [0, 1, 4]}, "nparent"),
        ({**CHAIN_KW, "nparent": [0, 1, 2.0]}, "nparent"),
        ({**CHAIN_KW, "jparent": [0, [1], [1]]}, "jparent"),
        ({**CHAIN_KW, "jparent": [0, [1], [1, 4]]}, "jparent"),
        ({**CHAIN_KW, "jparent": [0, [1], [0, 2]]}, "jparent"),
        ({**CHAIN_KW, "stotio": [0, [0.4], [0.25]]}, "stotio"),
        (_without(CHAIN_KW, "stotio"), "stotio"),
        ({**_without(CHAIN_KW, "stotio"), "spatialreact": True}, "sptlrct"),
        (
            {**CHAIN_KW, "spatialreact": True, "sptlrct": [0, [0.1], [0.1]]},
            "sptlrct",
        ),
    ],
    ids=[
        "nparent_missing",
        "nparent_length",
        "nparent_negative",
        "nparent_above_mcomp",
        "nparent_float",
        "jparent_count",
        "jparent_above_mcomp",
        "jparent_zero",
        "stotio_count",
        "stotio_missing",
        "sptlrct_missing",
        "sptlrct_count",
    ],
)
def test_bct_chaindecay_invalid(function_tmpdir, kw, match):
    """Rejected before the model is changed."""
    ml = _model(function_tmpdir)
    before = _state(ml)
    with pytest.raises(ValueError, match=match):
        _bct(ml, ipakcb=55, **kw)
    assert _state(ml) == before


def _set(**values):
    def change(bct):
        for name, value in values.items():
            setattr(bct, name, value)

    return change


@pytest.mark.parametrize(
    "kw, change, match",
    [
        ({}, _set(chaindecay=1), "nparent"),
        (CHAIN_KW, _set(nparent=[0, 1]), "nparent"),
        (CHAIN_KW, _set(nparent=[0, 2, 2]), "jparent"),
        (CHAIN_KW, _set(spatialreact=1), "sptlrct"),
    ],
    ids=["chaindecay_added", "nparent_length", "nparent_changed", "spatialreact"],
)
@pytest.mark.parametrize("target", ["file", "handle"])
def test_bct_chaindecay_write_recheck(function_tmpdir, target, kw, change, match):
    """Changed after construction: rejected before anything is written."""
    bct = _bct(_model(function_tmpdir, structured=False), **kw)
    bct.fn_path = str(function_tmpdir / "keep.bct")
    Path(bct.fn_path).write_text(SENTINEL)
    change(bct)
    handle = io.StringIO() if target == "handle" else None
    with pytest.raises(ValueError, match=match):
        bct.write_file(f=handle)
    assert Path(bct.fn_path).read_text() == SENTINEL
    if handle is not None:
        assert handle.getvalue() == ""


@pytest.mark.parametrize(
    "structured, nparent, warned",
    [(True, [0, 1, 2], True), (True, [0, 0, 0], False), (False, [0, 1, 2], False)],
    ids=["dis_parents", "dis_no_parents", "disu_parents"],
)
def test_bct_chaindecay_structured_warning(
    function_tmpdir, structured, nparent, warned
):
    """Structured CHAINDECAY with parents: written as the format says, with a
    warning (the local 2.8.1 and 2.7.0 builds stop while reading it)."""
    kw = {**CHAIN_KW, "nparent": nparent}
    if not any(nparent):
        kw = {**kw, "jparent": [0, 0, 0], "stotio": [0, 0, 0]}
    bct = _bct(_model(function_tmpdir, structured=structured), **kw)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        _written(bct, function_tmpdir)
    found = [w for w in caught if "1165" in str(w.message)]
    assert bool(found) == warned
    if warned:
        assert issubclass(found[0].category, UserWarning)
        assert "912-923" in str(found[0].message)


@GRIDS
def test_bct_chaindecay_without_parents_or_jparent(function_tmpdir, structured):
    """NPARENT=0 for every species needs neither jparent nor stotio (USG-T
    reads only NPARENT): built, written as the format says and reloaded."""
    ml = _model(function_tmpdir, structured=structured)
    bct = _bct(ml, chaindecay=True, nparent=[0, 0, 0])
    hand = (
        _head(mcomp=3, opts=CHAIN)
        + _consts(PRS)
        + "0\n"
        + _consts(21.0)
        + "0\n"
        + _consts(22.0)
        + "0\n"
        + _consts(23.0)
    )
    assert _parse(_written(bct, function_tmpdir)) == _parse(hand)
    loaded = MfUsgBct.load(_hand(function_tmpdir, hand), ml)
    assert list(loaded.nparent) == [0, 0, 0]
    assert loaded.jparent == [0, 0, 0] and loaded.stotio == [0, 0, 0]
    assert _parse(_written(loaded, function_tmpdir, "rewritten.bct")) == _parse(hand)


# ---------------------------------------------------------------- SPATIALREACT


def _model2(ws, structured):
    """MfUsg with 2 layers: DIS 2x1x3 or DISU with 3 and 2 nodes, free format."""
    ml = MfUsg(modelname="m2", model_ws=str(ws), structured=structured)
    if structured:
        ModflowDis(ml, nlay=2, nrow=1, ncol=3, top=0.0, botm=[-1.0, -2.0])
    else:
        MfUsgDisU(
            ml,
            nodes=5,
            nlay=2,
            nodelay=[3, 2],
            njag=5,
            iac=np.ones(5, dtype=int),
            ja=np.arange(1, 6),
            cl12=np.ones(5),
            fahl=np.ones(5),
            area=1.0,
            top=[0.0] * 5,
            bot=[-1.0] * 5,
        )
    MfUsgBas(ml, ibound=1, strt=1.0)
    ml.free_format_input = True
    return ml


# species 3 with parents 1 and 2; one array per layer for each parent
SPR_KW = {
    "chaindecay": True,
    "spatialreact": True,
    "nparent": [0, 0, 2],
    "jparent": [0, 0, [1, 2]],
    "sptlrct": [0, 0, [[0.1, 0.2], [0.3, 0.4]]],
}


def _hand_spr(sptlrct_text):
    return (
        _head(mcomp=3, opts=SPR)
        + _consts(PRS, PRS)
        + "0\n"
        + _consts(21.0, 21.0)
        + "0\n"
        + _consts(22.0, 22.0)
        + "2\n"
        + INT_I
        + " 1 2\n"
        + sptlrct_text
        + _consts(23.0, 23.0)
    )


@GRIDS
def test_bct_spatialreact_written(function_tmpdir, structured):
    bct = _bct(_model2(function_tmpdir, structured), **SPR_KW)
    hand = _hand_spr(_consts(0.1, 0.2, 0.3, 0.4))
    assert _parse(_quiet_write(bct, function_tmpdir)) == _parse(hand)
    assert bct.stotio == [0, 0, 0]


@GRIDS
def test_bct_spatialreact_loaded(function_tmpdir, structured):
    # parent 2, layer 2 with one value per node (DIS 3 columns, DISU 2 nodes)
    layer2 = " 0.41 0.42 0.43\n" if structured else " 0.41 0.42\n"
    hand = _hand_spr(_consts(0.1, 0.2, 0.3) + INT_R + layer2)
    ml = _model2(function_tmpdir, structured)
    bct = MfUsgBct.load(_hand(function_tmpdir, hand), ml)
    assert list(bct.nparent) == [0, 0, 2]
    assert bct.jparent[2].array.tolist() == [1, 2]
    assert bct.stotio == [0, 0, 0]
    spr = bct.sptlrct[2]
    assert len(spr) == 2 and all(len(p) == 2 for p in spr)
    assert np.allclose(spr[0][0].array, 0.1) and np.allclose(spr[0][1].array, 0.2)
    assert np.allclose(spr[1][0].array, 0.3)
    assert np.allclose(spr[1][1].array.ravel(), [float(v) for v in layer2.split()])
    assert _parse(_quiet_write(bct, function_tmpdir)) == _parse(hand)


def _model_state(ml):
    """Packages, counts and the output files a package registers."""
    return (
        _state(ml),
        list(ml.output_fnames),
        list(ml.output_units),
        list(ml.output_binflag),
        list(ml.output_packages),
    )


@pytest.mark.parametrize(
    "parents, match",
    [
        ([[0.1], [0.3, 0.4]], "1 layer values for 2 layers"),
        ([[np.ones(4), 0.2], [0.3, 0.4]], "does not match"),
    ],
    ids=["layer_count", "node_count"],
)
@GRIDS
def test_bct_spatialreact_layer_shape_invalid(
    function_tmpdir, structured, parents, match
):
    """A per-layer SPTLRCT value the arrays cannot take is rejected before the
    model is changed (no cbc output file registered, no package added)."""
    ml = _model2(function_tmpdir, structured)
    before = _model_state(ml)
    with pytest.raises(Exception, match=match):
        _bct(ml, ipakcb=55, **{**SPR_KW, "sptlrct": [0, 0, parents]})
    assert _model_state(ml) == before


# ------------------------------------------------------------------ MULTIFILE


def test_bct_multifile_roundtrip(function_tmpdir):
    """Written as given; loaded in upper case, as the engine reads it."""
    bct = _bct(
        _model(function_tmpdir, structured=False), mcomp=1, conc=21.0, crootname="mroot"
    )
    text = _written(bct, function_tmpdir, "a.bct")
    assert text.splitlines()[1].split()[19:] == ["MULTIFILE", "mroot"]
    ml = _model(function_tmpdir, structured=False)
    loaded = MfUsgBct.load(str(function_tmpdir / "a.bct"), ml)
    again = _written(loaded, function_tmpdir, "b.bct")
    assert again.splitlines()[1].split()[19:] == ["MULTIFILE", "MROOT"]
    assert _parse(again)[2:] == _parse(text)[2:]


# ------------------------------------------------------------------ REACTION1


def test_bct_reaction1_load_guard(function_tmpdir):
    """Guard (not support): raised from item 1, before anything else is read."""
    hand = _head(opts=" REACTION1") + "NOT A BCT RECORD\n"
    ml = _model(function_tmpdir)
    before = _state(ml)
    with pytest.raises(NotImplementedError, match="REACTION1"):
        MfUsgBct.load(_hand(function_tmpdir, hand), ml)
    assert _state(ml) == before


def test_bct_reaction1_construction_guard(function_tmpdir):
    ml = _model(function_tmpdir)
    before = _state(ml)
    with pytest.raises(NotImplementedError, match="REACTION1"):
        _bct(ml, mcomp=1, conc=21.0, reaction1=True)
    assert _state(ml) == before


# ------------------------------------------------------------ engine (requires_usgt)


def _eng_run(ws, mcomp, bct_kw=None, hand=None, structured=False):
    """Run the 5-node model (DISU, or DIS 1x1x5 with IDISP=1); BCT from FloPy
    (bct_kw) or hand-written. Returns (normal termination, model)."""
    ml = MfUsg(
        modelname="m", model_ws=str(ws), structured=structured, exe_name=USGT_EXE
    )
    if structured:
        ModflowDis(
            ml,
            nlay=1,
            nrow=1,
            ncol=5,
            delr=10.0,
            delc=10.0,
            top=10.0,
            botm=0.0,
            nper=1,
            perlen=10.0,
            nstp=2,
            steady=False,
        )
    else:
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
    kw = {"mcomp": mcomp, "prsity": 0.2, "idisp": 1, "dl": 0.0, "dt": 0.0}
    bct = MfUsgBct(ml, **{**kw, **(bct_kw or {})})
    MfUsgOc(ml, stress_period_data={(0, 1): ["save head", "save concentration"]})
    if hand is None:
        ml.write_input()
    else:
        ml.write_input(
            SelPackList=[p.name[0] for p in ml.packagelist if p.name[0] != "BCT"]
        )
        Path(bct.fn_path).write_text("# hand-written BCT\n" + hand)
    success, buff = ml.run_model(silent=True, report=True)
    return success and any("normal termination" in b.lower() for b in buff), ml


def _records(path, structured=False):
    """Last record of each label of a single-precision concentration file.

    Stream files (no record markers) and sequential unformatted files (a
    4-byte length before and after the 44-byte header and the values; the
    local gfortran builds write the MULTIFILE files this way) are read.
    """
    data, out, off = Path(path).read_bytes(), {}, 0
    seq = struct.unpack_from("<i", data, 0)[0] == 44
    while off < len(data):
        h = off + 4 if seq else off
        text = data[h + 16 : h + 32].decode().strip()
        a, b, _ = struct.unpack_from("<iii", data, h + 32)
        n = a * b if structured else b - a + 1
        v = h + 44 + (8 if seq else 0)
        if seq:
            assert struct.unpack_from("<ii", data, h + 44) == (44, 4 * n)
        out[text] = np.array(struct.unpack_from(f"<{n}f", data, v))
        off = v + 4 * n + (4 if seq else 0)
    return out


def _con(ml):
    return (Path(ml.model_ws) / "m.con").read_bytes()


def _eng_head(mcomp, opts=""):
    return f" 1 0 {mcomp} 1 1 0 0 -999.0 1.0e-6 1 0 0.0 0 1 0 0 0 0 0{opts}\n"


ENG_GRID = _consts(0.2, 0.0, 0.0, 0.0)  # PRSITY, ANGLEX, DL, DT
# species 1 and 2 decay (FODRW) and feed species 3
ENG_SPECIES = {"ifod": 1, "fodrw": 0.05, "conc": [10.0, 20.0, 0.0]}


def _eng_chain(stotio=(0.2, 0.9)):
    return {
        **ENG_SPECIES,
        "chaindecay": True,
        "nparent": [0, 0, 2],
        "jparent": [0, 0, [1, 2]],
        "stotio": [0, 0, list(stotio)],
    }


def _eng_hand(parents_text, opts=CHAIN):
    return (
        _eng_head(3, opts)
        + ENG_GRID
        + "0\n"
        + _consts(0.05, 10.0)
        + "0\n"
        + _consts(0.05, 20.0)
        + "2\n"
        + INT_I
        + " 1 2\n"
        + parents_text
        + _consts(0.05, 0.0)
    )


@requires_usgt
def test_bct_engine_chaindecay_flopy_equals_hand(function_tmpdir):
    fa = _eng_run(function_tmpdir / "flopy", 3, _eng_chain())
    fb = _eng_run(function_tmpdir / "hand", 3, hand=_eng_hand(INT_R + " 0.2 0.9\n"))
    assert fb[0], "hand-written file did not run"
    assert fa[0] and _con(fa[1]) == _con(fb[1])


@requires_usgt
def test_bct_engine_chaindecay_stotio_sensitivity(function_tmpdir):
    """Each STOTIO coefficient changes species 3; the parents do not change."""
    runs = {"ref": (0.2, 0.9), "first": (0.9, 0.9), "second": (0.2, 0.2)}
    out = {k: _eng_run(function_tmpdir / k, 3, _eng_chain(s)) for k, s in runs.items()}
    assert all(ok for ok, _ in out.values())
    c = {k: _records(Path(ml.model_ws) / "m.con") for k, (_, ml) in out.items()}
    for k in ("first", "second"):
        assert np.all(np.abs(c[k]["CONC03"] - c["ref"]["CONC03"]) > 1e-3), k
        for parent in ("CONC01", "CONC02"):
            assert np.array_equal(c[k][parent], c["ref"][parent])
    # the generation term is linear in the coefficients
    x = c["second"]["CONC03"] / 0.2  # parents 1 + 2 at unit coefficient
    assert np.allclose(c["first"]["CONC03"], 0.9 * x, rtol=1e-4)


@requires_usgt
def test_bct_engine_structured_chaindecay_without_parents(function_tmpdir):
    """Structured reader with CHAINDECAY and NPARENT=0 (MCOMP=1): FloPy equals
    hand and equals the run without CHAINDECAY."""
    kw = {"ifod": 1, "fodrw": 0.05, "conc": 10.0}
    chain = {**kw, "chaindecay": True, "nparent": [0], "jparent": [0], "stotio": [0]}
    hand = _eng_head(1, CHAIN) + _consts(0.2, 0.0, 0.0) + "0\n" + _consts(0.05, 10.0)
    fa = _eng_run(function_tmpdir / "flopy", 1, chain, structured=True)
    fb = _eng_run(function_tmpdir / "hand", 1, hand=hand, structured=True)
    fc = _eng_run(function_tmpdir / "plain", 1, kw, structured=True)
    assert fb[0] and fc[0], "hand-written or plain file did not run"
    assert fa[0] and _con(fa[1]) == _con(fb[1]) == _con(fc[1])


def _eng_spr(first, second=0.9):
    return {
        **_eng_chain(),
        "spatialreact": True,
        "stotio": None,
        "sptlrct": [0, 0, [[first], [second]]],
    }


@requires_usgt
def test_bct_engine_spatialreact_flopy_equals_hand(function_tmpdir):
    first = [0.2, 0.2, 0.2, 0.2, 0.5]
    hand = _eng_hand(INT_R + " 0.2 0.2 0.2 0.2 0.5\n" + _consts(0.9), SPR)
    fa = _eng_run(function_tmpdir / "flopy", 3, _eng_spr(np.array(first)))
    fb = _eng_run(function_tmpdir / "hand", 3, hand=hand)
    assert fb[0], "hand-written file did not run"
    assert fa[0] and _con(fa[1]) == _con(fb[1])


@requires_usgt
def test_bct_engine_spatialreact_sensitivity(function_tmpdir):
    """Uniform SPTLRCT equals STOTIO; the values are used per node."""
    out = {
        "stotio": _eng_run(function_tmpdir / "stotio", 3, _eng_chain()),
        "uniform": _eng_run(function_tmpdir / "uniform", 3, _eng_spr(0.2)),
        "node5": _eng_run(
            function_tmpdir / "node5", 3, _eng_spr(np.array([0.2] * 4 + [0.5]))
        ),
    }
    assert all(ok for ok, _ in out.values())
    assert _con(out["uniform"][1]) == _con(out["stotio"][1])
    c = {
        k: _records(Path(ml.model_ws) / "m.con")["CONC03"] for k, (_, ml) in out.items()
    }
    assert np.array_equal(c["node5"][:4], c["uniform"][:4])
    assert c["node5"][4] > c["uniform"][4] + 1e-3


@requires_usgt
def test_bct_engine_multifile(function_tmpdir):
    """MULTIFILE reaches the engine: one file per species named from the root
    in upper case; FloPy, hand-written and loaded-rewritten files agree."""
    kw = {**ENG_SPECIES, "mcomp": 2, "conc": [10.0, 20.0], "crootname": "mroot"}
    names = ["MROOT CONC01.CON", "MROOT CONC02.CON"]
    hand = _eng_head(2, " MULTIFILE mroot") + ENG_GRID + _consts(0.05, 10.0, 0.05, 20.0)
    fa = _eng_run(function_tmpdir / "flopy", 2, kw)
    fb = _eng_run(function_tmpdir / "hand", 2, hand=hand)
    assert fa[0] and fb[0]
    for _, ml in (fa, fb):
        found = sorted(p.name for p in Path(ml.model_ws).iterdir() if "CONC0" in p.name)
        assert found == names
    # loaded from the FloPy file and rewritten (root now upper case)
    ws = function_tmpdir / "reloaded"
    ok, ml = _eng_run(ws, 2, kw)
    assert ok
    bct = MfUsgBct.load(str(Path(fa[1].model_ws) / "m.bct"), ml)
    assert bct.options.split()[-2:] == ["MULTIFILE", "MROOT"]
    bct.write_file()
    success, buff = ml.run_model(silent=True, report=True)
    assert success and any("normal termination" in b.lower() for b in buff)
    for name in names:
        ref = (Path(fa[1].model_ws) / name).read_bytes()
        assert (Path(fb[1].model_ws) / name).read_bytes() == ref
        assert (ws / name).read_bytes() == ref
    # the species files hold the values of the single concentration file
    plain = _eng_run(function_tmpdir / "plain", 2, _without(kw, "crootname"))
    assert plain[0]
    single = _records(Path(plain[1].model_ws) / "m.con")
    for i, name in enumerate(names, start=1):
        per = _records(Path(fa[1].model_ws) / name)
        assert np.array_equal(per[f"CONC0{i}"], single[f"CONC0{i}"])
