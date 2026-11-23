"""MfUsgMdt (USG-T MDT): per-component read sequence, guards, handles.

Expected read sequences are written by hand from GWT2MDTU1AR (USG-T 2.8.1
gwt2mdtu1.for; the 2.7.0 reader is the same). They never come from the FloPy
writer. MDT has one reader for structured and unstructured grids. The
per-species arrays KDMD, DECAYMD, YIELDMD, DIFFMD (and AIOLD1, AIOLD2 with
TSHIFTMD > 1e-10) are read for ICOMP = 1..NTCOMP (line 251): the MCOMP mobile
species, then the temperature when IHEAT is 1. Every array is one CONSTANT
record with a distinct value, so the sequence of values identifies the
arrays. Engine tests (requires_usgt) compare a FloPy file with a hand-written
file on the same model.
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
    MfUsgDisU,
    MfUsgLpf,
    MfUsgMdt,
    MfUsgOc,
    MfUsgSms,
)
from flopy.modflow import ModflowDis

# distinct values per array; per-component arrays as [species 1, second
# component (temperature or species 2)]
MDFLAG, VOLFRAC, POR, RHOB, DIFFLEN, TORT = 2, 0.41, 0.21, 1.61, 0.51, 0.71
KD, DECAY, YIELD, DIFF = [1.1, 1.2], [2.1, 2.2], [3.1, 3.2], [4.1, 4.2]
AI1, AI2 = [5.1, 5.2], [6.1, 6.2]
COMMON = [MDFLAG, VOLFRAC, POR, RHOB, DIFFLEN, TORT]
SENTINEL = "SENTINEL MDT FILE\n"
BCT_KW = {"mcomp": 1, "icbndflg": 1, "idisp": 0, "prsity": 0.2, "conc": 21.0}
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
IMCOMP_BCT = {"imcomp": 1, "imconc": 2.0}
BASE = {
    "mdflag": MDFLAG,
    "volfracmd": VOLFRAC,
    "pormd": POR,
    "rhobmd": RHOB,
    "difflenmd": DIFFLEN,
    "tortmd": TORT,
}


def _kw(ncomp, tshift=False):
    """MfUsgMdt arguments with one distinct value per component."""
    kw = dict(
        BASE,
        kdmd=KD[:ncomp],
        decaymd=DECAY[:ncomp],
        yieldmd=YIELD[:ncomp],
        diffmd=DIFF[:ncomp],
    )
    if tshift:
        kw.update(tshiftmd=2.0, aiold1md=AI1[:ncomp], aiold2md=AI2[:ncomp])
    return kw


def _sequence(ncomp, tshift=False):
    """Values in the GWT2MDTU1AR order (IDPF=0), written by hand."""
    seq = list(COMMON)
    for i in range(ncomp):
        seq += [KD[i], DECAY[i], YIELD[i], DIFF[i]]
        if tshift:
            seq += [AI1[i], AI2[i]]
    return seq


def _model(ws, name, structured=False, bct=None, manual=None):
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
    if bct is not None:
        MfUsgBct(m, **{**BCT_KW, **bct})
    if manual is not None:  # counts assigned by hand, no BCT
        m.itrnsp = 1
        m.mcomp, m.iheat, m.imcomp = manual
    return m


def _values(text):
    """(header tokens, CONSTANT values) of an MDT file text."""
    lines = [ln for ln in text.splitlines() if ln.strip() and not ln.startswith("#")]
    head = lines[0].split("#")[0].split()
    vals = []
    for ln in lines[1:]:
        t = ln.split()
        assert t[0].upper() == "CONSTANT", ln
        vals.append(float(t[1]))
    return head, vals


def _written(mdt):
    mdt.write_file()
    return _values(Path(mdt.fn_path).read_text())[1]


def _hand_text(values, header="0 0"):
    return header + "\n" + "".join(f"CONSTANT {v}\n" for v in values)


def _same(values, expected):
    return len(values) == len(expected) and np.allclose(values, expected, rtol=1e-6)


def _assert_rejected_untouched(mdt, exc, match):
    """write_file raises before opening: the file and a handle stay intact."""
    p = Path(mdt.fn_path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(SENTINEL)
    with pytest.raises(exc, match=match):
        mdt.write_file()
    assert p.read_text() == SENTINEL
    buf = io.StringIO("PREFIX\n")
    buf.seek(0, io.SEEK_END)
    pos = buf.tell()
    with pytest.raises(exc, match=match):
        mdt.write_file(f=buf)
    assert not buf.closed and buf.getvalue() == "PREFIX\n" and buf.tell() == pos


def _first(util):
    return float(util.array.ravel()[0])


# ------------------------------------------------------- heat (IHEAT=1): NTCOMP


@pytest.mark.parametrize("structured", [False, True], ids=["disu", "structured"])
def test_mdt_heat_sequence(function_tmpdir, structured):
    m = _model(function_tmpdir, "hs", structured=structured, bct=HEAT)
    mdt = MfUsgMdt(m, **_kw(2))
    assert len(mdt.kdmd) == len(mdt.diffmd) == 2
    assert _same(_written(mdt), _sequence(2))


def test_mdt_heat_hand_file_loads_and_rewrites(function_tmpdir):
    m = _model(function_tmpdir, "hl", bct=HEAT)
    hand = function_tmpdir / "hand.mdt"
    hand.write_text(_hand_text(_sequence(2)))
    mdt = MfUsgMdt.load(str(hand), m)
    for name, expected in (
        ("kdmd", KD),
        ("decaymd", DECAY),
        ("yieldmd", YIELD),
        ("diffmd", DIFF),
    ):
        arrays = getattr(mdt, name)
        assert len(arrays) == 2, name
        assert np.allclose([_first(a) for a in arrays], expected), name
    assert _same(_written(mdt), _values(hand.read_text())[1])


def test_mdt_heat_tshift_sequence(function_tmpdir):
    m = _model(function_tmpdir, "ht", bct=HEAT)
    hand = function_tmpdir / "hand.mdt"
    hand.write_text(_hand_text(_sequence(2, tshift=True), "0 0 TSHIFTMD 2.0"))
    built = MfUsgMdt(m, **_kw(2, tshift=True))
    assert _same(_written(built), _sequence(2, tshift=True))
    loaded = MfUsgMdt.load(str(hand), _model(function_tmpdir, "ht2", bct=HEAT))
    assert np.allclose([_first(a) for a in loaded.aiold1md], AI1)
    assert np.allclose([_first(a) for a in loaded.aiold2md], AI2)
    assert _same(_written(loaded), _sequence(2, tshift=True))


def test_mdt_heat_scalars_broadcast_to_temperature(function_tmpdir):
    m = _model(function_tmpdir, "hb", bct=HEAT)
    mdt = MfUsgMdt(m, **BASE, kdmd=1.1, decaymd=2.1, yieldmd=3.1, diffmd=4.1)
    assert _same(_written(mdt), COMMON + [1.1, 2.1, 3.1, 4.1] * 2)


def test_mdt_heat_manual_counts(function_tmpdir):
    m = _model(function_tmpdir, "hm", manual=(1, 1, 0))
    assert _same(_written(MfUsgMdt(m, **_kw(2))), _sequence(2))


def test_mdt_heat_species_only_length_rejected(function_tmpdir):
    m = _model(function_tmpdir, "hw", bct=HEAT)
    with pytest.raises(ValueError, match=r"kdmd.*2 entries"):
        MfUsgMdt(m, **_kw(1))


# --------------------------------------------------- no heat: MCOMP (unchanged)


@pytest.mark.parametrize("mode", ["bct", "manual", "none"])
def test_mdt_no_heat_sequence(function_tmpdir, mode):
    if mode == "bct":
        m = _model(function_tmpdir, "nb", bct={"mcomp": 2})
    elif mode == "manual":
        m = _model(function_tmpdir, "nm", manual=(2, 0, 0))
    else:  # no BCT, itrnsp=0: the counts of the model are used as they are
        m = _model(function_tmpdir, "nn")
        m.mcomp = 2
    assert _same(_written(MfUsgMdt(m, **_kw(2))), _sequence(2))


def test_mdt_no_transport_default_model(function_tmpdir):
    # a new model has mcomp=0: no per-species array, no BCT needed
    m = _model(function_tmpdir, "nd")
    assert _same(_written(MfUsgMdt(m, **BASE)), COMMON)


@pytest.mark.parametrize(
    "folder, nam, mcomp",
    [
        ("DiscreteFracture", "USG_discrete_fracture.nam", 1),
        ("SandTank", "usg_sand_tank.nam", 1),
        ("Multispecies", "USG_Multispecies.nam", 4),
    ],
)
def test_mdt_ex7_shipped_file_rewritten(
    function_tmpdir, example_data_path, folder, nam, mcomp
):
    # the Ex7 MDT files were not written by FloPy
    ws = example_data_path / "mfusg_transport" / "Ex7_Matrix_Diffusion" / folder
    m = MfUsg.load(
        nam, model_ws=str(ws), check=False, load_only=["DISU", "BAS6", "BCT", "MDT"]
    )
    assert m.mcomp == mcomp and len(m.mdt.kdmd) == mcomp
    shipped = _values(Path(m.mdt.fn_path).read_text())[1]
    m.change_model_ws(str(function_tmpdir))
    assert _same(_written(m.mdt), shipped)


# --------------------------------------------------------- IMCOMP > 0: guarded


@pytest.mark.parametrize("mode", ["bct", "manual"])
def test_mdt_imcomp_rejected_at_construction(function_tmpdir, mode):
    if mode == "bct":
        m = _model(function_tmpdir, "ic", bct=IMCOMP_BCT)
    else:
        m = _model(function_tmpdir, "ic", manual=(1, 0, 1))
    with pytest.raises(NotImplementedError, match="IMCOMP"):
        MfUsgMdt(m, **_kw(1))
    assert m.get_package("MDT") is None


def test_mdt_imcomp_rejected_at_load_before_reading(function_tmpdir):
    m = _model(function_tmpdir, "il", bct=IMCOMP_BCT)
    buf = io.StringIO(_hand_text(_sequence(2)))
    with pytest.raises(NotImplementedError, match="IMCOMP"):
        MfUsgMdt.load(buf, m)
    assert not buf.closed and buf.tell() == 0


def test_mdt_imcomp_rejected_at_write(function_tmpdir):
    m = _model(function_tmpdir, "iw", manual=(1, 0, 0))
    mdt = MfUsgMdt(m, **_kw(1))
    m.imcomp = 1
    _assert_rejected_untouched(mdt, NotImplementedError, "IMCOMP")


# ------------------------------------------------------- counts and lifecycle


def _unsynced_bct(m):
    """Register a BCT that did not set the model counts."""
    bct = MfUsgBct(m, **{**BCT_KW, "mcomp": 2}, add_package=False)
    m.add_package(bct)


@pytest.mark.parametrize(
    "case", ["mcomp", "iheat", "bct_replaced", "bct_removed", "bct_unsynced"]
)
def test_mdt_write_rechecks_counts(function_tmpdir, case):
    if case in ("mcomp", "iheat"):
        m = _model(function_tmpdir, f"wr_{case}", manual=(1, 0, 0))
    else:
        m = _model(function_tmpdir, f"wr_{case}", bct={})
    mdt = MfUsgMdt(m, **_kw(1))
    if case == "mcomp":
        m.mcomp = 2
    elif case == "iheat":
        m.iheat = 1
    elif case == "bct_replaced":
        MfUsgBct(m, **{**BCT_KW, "mcomp": 2})
    elif case == "bct_removed":
        m.remove_package("BCT")
    else:
        _unsynced_bct(m)
    _assert_rejected_untouched(mdt, ValueError, "BCT|count")


def test_mdt_write_bct_replaced_same_counts(function_tmpdir):
    m = _model(function_tmpdir, "ws", bct={})
    mdt = MfUsgMdt(m, **_kw(1))
    MfUsgBct(m, **BCT_KW)
    assert _same(_written(mdt), _sequence(1))


def test_mdt_write_rechecks_array_lists(function_tmpdir):
    m = _model(function_tmpdir, "wl", bct={})
    mdt = MfUsgMdt(m, **_kw(1))
    mdt.diffmd = []
    _assert_rejected_untouched(mdt, ValueError, "diffmd")


@pytest.mark.parametrize("case", ["bct_removed", "bct_unsynced"])
@pytest.mark.parametrize("how", ["construction", "load"])
def test_mdt_stale_bct_counts_rejected(function_tmpdir, case, how):
    m = _model(function_tmpdir, "st", bct={})
    if case == "bct_removed":
        m.remove_package("BCT")
    else:
        _unsynced_bct(m)
    with pytest.raises(ValueError, match="BCT"):
        if how == "construction":
            MfUsgMdt(m, **_kw(1))
        else:
            MfUsgMdt.load(io.StringIO(_hand_text(_sequence(1))), m)


# --------------------------------------------------------------- file handles


def test_mdt_load_handle_not_closed_cursor_after_data(function_tmpdir):
    m = _model(function_tmpdir, "lh", bct=HEAT)
    body = _hand_text(_sequence(2))
    buf = io.StringIO(body + "TRAILER\n")
    mdt = MfUsgMdt.load(buf, m)
    assert not buf.closed
    assert buf.tell() == len(body)
    assert np.isclose(_first(mdt.diffmd[1]), DIFF[1])


def test_mdt_load_file_closed_after_error(function_tmpdir, monkeypatch):
    m = _model(function_tmpdir, "le", bct={})
    f = function_tmpdir / "bad.mdt"
    f.write_text(_hand_text(_sequence(1), "0 notanumber"))
    opened = []
    real_open = builtins.open

    def spy(*args, **kwargs):
        fh = real_open(*args, **kwargs)
        opened.append(fh)
        return fh

    monkeypatch.setattr(builtins, "open", spy)
    with pytest.raises(ValueError):
        MfUsgMdt.load(str(f), m)
    monkeypatch.undo()
    assert opened and all(fh.closed for fh in opened)


# ------------------------------------------------------------ engine (requires_usgt)

ENG_COMMON = [2, 0.5, 0.3, 1.6, 0.5, 0.5]
ENG = {
    "mdflag": 2,
    "volfracmd": 0.5,
    "pormd": 0.3,
    "rhobmd": 1.6,
    "difflenmd": 0.5,
    "tortmd": 0.5,
    "kdmd": [0.1, 0.2],
    "decaymd": [1.0e-3, 2.0e-3],
    "yieldmd": [1.0, 0.5],
    "diffmd": [1.0e-2, 2.0e-2],
}


def _eng_values(kw):
    seq = list(ENG_COMMON)
    for i in range(len(kw["kdmd"])):
        seq += [kw["kdmd"][i], kw["decaymd"][i], kw["yieldmd"][i], kw["diffmd"][i]]
    return seq


def _eng_model(ws):
    """5-node DISU line, one solute and heat (IHEAT=1), no boundaries."""
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
    MfUsgBct(
        ml,
        itrnsp=1,
        mcomp=1,
        idisp=1,
        dl=0.0,
        dt=0.0,
        diffnc=0.0,
        prsity=0.2,
        conc=10.0,
        **HEAT,
    )
    MfUsgOc(ml, stress_period_data={(0, 1): ["save head", "save concentration"]})
    return ml


def _records(data):
    """Values of the unstructured single-precision concentration records."""
    out, off = [], 0
    while off < len(data):
        n1, n2 = struct.unpack_from("<ii", data, off + 32)
        n = n2 - n1 + 1
        out.append(np.array(struct.unpack_from(f"<{n}f", data, off + 44)))
        off += 44 + 4 * n
    return out


def _run(ml):
    success, buff = ml.run_model(silent=True, report=True)
    normal = any("normal termination" in str(b).lower() for b in buff)
    ws = Path(ml.model_ws)
    return success and normal, ws


def _flopy_run(ws, **mdt_kw):
    ml = _eng_model(ws)
    MfUsgMdt(ml, **mdt_kw)
    ml.write_input()
    ok, ws = _run(ml)
    assert ok
    return _records((ws / "m.con").read_bytes()), (ws / "m.hds").read_bytes()


def _hand_run(ws, values):
    """Same model, MDT file written by hand (not by the FloPy writer)."""
    ml = _eng_model(ws)
    mdt = MfUsgMdt(ml)
    ml.write_input(
        SelPackList=[p.name[0] for p in ml.packagelist if p.name[0] != "MDT"]
    )
    Path(mdt.fn_path).write_text(_hand_text(values))
    return _run(ml)


@requires_usgt
def test_mdt_engine_heat_hand_file_runs(function_tmpdir):
    # control of the hand-written file itself: NTCOMP = 2 blocks
    ok, ws = _hand_run(function_tmpdir, _eng_values(ENG))
    assert ok
    recs = _records((ws / "m.con").read_bytes())
    assert len(recs) == 2  # solute, temperature
    # mass and heat diffuse into the matrix
    assert np.all(recs[0] < 10.0) and np.all(recs[1] < 20.0)
    # the engine echoes the arrays of both components
    lst = (ws / "m.list").read_text(errors="replace")
    assert lst.count("mdT ADSORPTION COEF") == lst.count("mdT MATRIX DIFFSN COEFF") == 2


@requires_usgt
def test_mdt_engine_heat_species_only_file_fails(function_tmpdir):
    # negative control: a file with the MCOMP species blocks only
    ok, ws = _hand_run(function_tmpdir, _eng_values(ENG)[:10])
    assert not ok
    # the listing, which can be cut short or missing when the engine stops,
    # has no complete echo of the arrays of component 2 (gwt2mdtu1.for:251-303),
    # no transport solution and no end of run
    listing = ws / "m.list"
    lst = listing.read_text(errors="replace").upper() if listing.is_file() else ""
    assert lst.count("MDT MATRIX DIFFSN COEFF") < 2
    assert "TRANSPORT SOLUTION COMPLETE" not in lst
    assert "ELAPSED RUN TIME" not in lst


@requires_usgt
def test_mdt_engine_heat_flopy_equals_hand(function_tmpdir):
    con, hds = _flopy_run(function_tmpdir / "flopy", **ENG)
    text = (function_tmpdir / "flopy" / "m.mdt").read_text()
    assert _same(_values(text)[1], _eng_values(ENG))
    ok, ws = _hand_run(function_tmpdir / "hand", _eng_values(ENG))
    assert ok
    hand = _records((ws / "m.con").read_bytes())
    assert len(con) == len(hand) == 2
    assert all(np.array_equal(a, b) for a, b in zip(con, hand))
    assert hds == (ws / "m.hds").read_bytes()


@requires_usgt
@pytest.mark.parametrize("block", [0, 1], ids=["solute", "temperature"])
def test_mdt_engine_heat_sensitivity(function_tmpdir, block):
    base = _flopy_run(function_tmpdir / "base", **ENG)[0]
    diffmd = list(ENG["diffmd"])
    diffmd[block] *= 5.0
    changed = _flopy_run(function_tmpdir / "changed", **{**ENG, "diffmd": diffmd})[0]
    assert np.all(changed[block] < base[block])  # more diffusion into the matrix
    if block == 1:  # the temperature block does not change the solute
        assert np.array_equal(changed[0], base[0])
