"""MfUsgBct header and read conditions (USG-T 2.8.1 glo2btnu1.f).

Contract, from GWT2BCT1AR and the two readers SGWT2BCT1S (structured,
lines 849-1469) and SGWT2BCT1G (unstructured, lines 1471-2032); line numbers
are 2.8.1 (structured / unstructured). The 2.7.0 readers are the same for
these items.

- IFMBC /= 0: one record MBEGWUNF MBEGWUNT MBECLNUNF MBECLNUNT (357-379).
- Zero-order decay, per species (1262-1344 / 1840-1915): IZOD=1 ZODRW;
  IZOD=2 ZODRS only with IADSORB > 0; IZOD=3 ZODRW ZODRS only with
  IADSORB > 0; IZOD=4 ZODRW ZODRS ZODRAW only with A-W_ADSORB; nothing
  otherwise. IFOD is the same with FODRW FODRS FODRAW (1347-1429 / 1918-1993).
- SOLUBILITY: SOLLIM(MCOMP) and SOLSLOPE(MCOMP), once, after the A-W arrays
  and before the species loop (1139-1142 / 1723-1726).
- BULKD: structured only with IADSORB /= 0 (964); unstructured with
  IADSORB /= 0 or IHEAT = 1 (1578). The structured reader still uses BULKD
  for the heat capacity term (1114), so IHEAT = 1 with IADSORB = 0 on a
  structured grid is rejected by MfUsgBct (guard, not support).
- VARIABLE_SOLUBILITY (94-98, 1237-1259 / 1817-1837) is not implemented:
  guard on load and on construction.

Expected files are written by hand from this contract, never by the FloPy
writer, and are read back with the small reader below. Engine tests
(requires_usgt) compare a FloPy file with a hand-written file on the same model.
"""

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
    MfUsgOc,
    MfUsgSms,
)
from flopy.modflow import ModflowDis
from flopy.utils.mfreadnam import NamData

# distinct values per array
PRS, BULK, AWA, ALA, BLA, ADS, CONC = 0.21, 1.61, 7.1, 8.1, 9.1, 1.1, 21.0
HTCAPS, HTCONDS, TEMP = 51.0, 52.0, 31.0
HEAT = {
    "iheat": 1,
    "htcondw": 0.58,
    "rhow": 1000.0,
    "htcapw": 4200.0,
    "htcaps": HTCAPS,
    "htconds": HTCONDS,
    "heat": TEMP,
}
HEAT_LINE = " 0.58 1000.0 4200.0\n"
SENTINEL = "SENTINEL BCT FILE\n"
GRIDS = pytest.mark.parametrize("structured", [True, False], ids=["dis", "disu"])


def _model(ws, name="m", structured=True):
    """MfUsg with a 1-layer, 3-cell DIS (or DISU) grid and BAS, free format."""
    ml = MfUsg(modelname=name, model_ws=str(ws), structured=structured)
    if structured:
        ModflowDis(ml, nlay=1, nrow=1, ncol=3, nper=1, steady=True)
    else:
        MfUsgDisU(
            ml,
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
    MfUsgBas(ml, ibound=1, strt=1.0)
    ml.free_format_input = True
    return ml


def _bct(ml, **kw):
    base = {"mcomp": 1, "itvd": 0, "idisp": 0, "prsity": PRS, "conc": CONC}
    return MfUsgBct(ml, **{**base, **kw})


def _head(mcomp=1, iadsorb=0, izod=0, ifod=0, ifmbc=0, iheat=0, opts=""):
    """Hand-written item 1 (options on the same record)."""
    return (
        f" 1 0 {mcomp} 1 0 {iadsorb} 0 -999.0 1.0e-6 0 0 0.0"
        f" {izod} {ifod} {ifmbc} {iheat} 0 0 0{opts}\n"
    )


def _consts(*values):
    return "".join(f"CONSTANT {v}\n" for v in values)


def _parse(text):
    """Item 1 and the records that follow, in file order (no FloPy reader).

    CONSTANT records give one value, INTERNAL records the values of the next
    line (times CNSTNT), any other record its numbers.
    """
    lines = [
        ln for ln in text.splitlines() if ln.strip() and not ln.lstrip().startswith("#")
    ]
    t = lines[0].split()
    out = [[float(v) for v in t[:19]], [v.upper() for v in t[19:]]]
    i = 1
    while i < len(lines):
        t = lines[i].split()
        key = t[0].upper()
        if key == "CONSTANT":
            out.append(round(float(t[1]), 5))
        elif key == "INTERNAL":
            i += 1
            c = float(t[1])
            out.append(("internal", [round(c * float(v), 5) for v in lines[i].split()]))
        else:
            out.append(("record", [round(float(v), 5) for v in t]))
        i += 1
    return out


def _written(bct, ws, name="written.bct"):
    bct.fn_path = str(Path(ws) / name)
    bct.write_file()
    return Path(bct.fn_path).read_text()


def _hand(ws, text, name="hand.bct"):
    p = Path(ws) / name
    p.write_text("# hand-written BCT\n" + text)
    return str(p)


def _value(arr):
    """Value of a per-layer list of constant Util2d."""
    return round(float(arr[0].array.mean()), 5)


def _state(ml):
    return (
        [type(p).__name__ for p in ml.packagelist],
        (ml.itrnsp, ml.mcomp, ml.iheat, ml.imcomp),
    )


# ---------------------------------------------------------------- IFMBC /= 0


def _hand_ifmbc(ifmbc):
    return _head(ifmbc=ifmbc) + " 71 72 0 0\n" + _consts(PRS, CONC)


@pytest.mark.parametrize("ifmbc", [2, -1, 1])
def test_bct_ifmbc_units_record_loaded(function_tmpdir, ifmbc):
    """The MBE units record is read with any IFMBC /= 0 (ifmbc=1: control)."""
    ml = _model(function_tmpdir)
    path = _hand(function_tmpdir, _hand_ifmbc(ifmbc))
    pk = ml.mfnam_packages
    eud = {
        150: NamData("BCT", path, None, pk),
        71: NamData("DATA", str(function_tmpdir / "m.mbf"), None, pk),
        72: NamData("DATA", str(function_tmpdir / "m.mbt"), None, pk),
    }
    bct = MfUsgBct.load(path, ml, ext_unit_dict=eud)
    assert bct.ifmbc == ifmbc
    assert list(bct.unit_number[2:6]) == [71, 72, 0, 0]
    assert (_value(bct.prsity), _value(bct.conc[0])) == (PRS, CONC)
    assert _parse(_written(bct, function_tmpdir)) == _parse(_hand_ifmbc(ifmbc))


# ------------------------------------------------------------ IZOD / IFOD 1-4

DECAY = {"izod": ("zodrw", "zodrs", "zodraw"), "ifod": ("fodrw", "fodrs", "fodraw")}
DVAL = {
    "zodrw": 3.1,
    "zodrs": 4.1,
    "zodraw": 5.1,
    "fodrw": 3.2,
    "fodrs": 4.2,
    "fodraw": 5.2,
}
# (IZOD or IFOD, IADSORB, A-W_ADSORB) -> arrays read (W water, S soil, A air-water)
DECAY_CASES = {
    "1": ((1, 0, False), "W"),
    "1_ads_aw": ((1, 1, True), "W"),
    "2_noads": ((2, 0, False), ""),
    "2_ads": ((2, 1, False), "S"),
    "3_ads": ((3, 1, False), "WS"),
    "4_aw_ads": ((4, 1, True), "WSA"),
    # the cases below differ from the reader before this change
    "3_noads": ((3, 0, False), ""),
    "3_noads_aw": ((3, 0, True), ""),
    "4_noaw_noads": ((4, 0, False), ""),
    "4_noaw_ads": ((4, 1, False), ""),
    "4_aw_noads": ((4, 0, True), "WSA"),
}
DECAY_PARAMS = pytest.mark.parametrize("case", list(DECAY_CASES))
DECAY_FLAGS = pytest.mark.parametrize("flag", ["izod", "ifod"])


def _decay_read(flag, case):
    return [n for n, k in zip(DECAY[flag], "WSA") if k in DECAY_CASES[case][1]]


def _hand_decay(flag, case, mcomp=1):
    """Hand-written file: species-independent arrays, then per species the
    A-W partitioning arrays, ADSORB, the decay arrays read and CONC."""
    (order, iadsorb, aw), _ = DECAY_CASES[case]
    vals = [PRS]
    if iadsorb:
        vals.append(BULK)
    if aw:
        vals.append(AWA)
    for ic in range(mcomp):
        if aw:
            vals += [ALA + ic, BLA + ic]
        if iadsorb:
            vals.append(ADS + ic)
        vals += [DVAL[n] + ic for n in _decay_read(flag, case)]
        vals.append(CONC + ic)
    head = _head(
        mcomp=mcomp,
        iadsorb=iadsorb,
        opts=" A-W_ADSORB 1 1" if aw else "",
        **{flag: order},
    )
    return head + _consts(*vals)


def _decay_kw(flag, case, mcomp=1):
    """FloPy arguments: every decay array is passed, read or not."""
    (order, iadsorb, aw), _ = DECAY_CASES[case]

    def species(v):
        return [v + ic for ic in range(mcomp)]

    kw = {"mcomp": mcomp, flag: order, "iadsorb": iadsorb, "conc": species(CONC)}
    kw.update(bulkd=BULK, adsorb=species(ADS))
    kw.update({n: species(DVAL[n]) for n in DECAY[flag]})
    if aw:
        kw.update(aw_adsorb=True, iarea_fn=1, ikawi_fn=1, awamax=AWA)
        kw.update(alangaw=species(ALA), blangaw=species(BLA))
    return kw


@GRIDS
@DECAY_FLAGS
@DECAY_PARAMS
def test_bct_decay_arrays_written(function_tmpdir, case, flag, structured):
    ml = _model(function_tmpdir, structured=structured)
    bct = _bct(ml, **_decay_kw(flag, case))
    assert _parse(_written(bct, function_tmpdir)) == _parse(_hand_decay(flag, case))
    assert [n for n in DECAY[flag] if hasattr(bct, n)] == _decay_read(flag, case)


@DECAY_FLAGS
@DECAY_PARAMS
def test_bct_decay_arrays_written_two_species(function_tmpdir, case, flag):
    ml = _model(function_tmpdir, structured=False)
    bct = _bct(ml, **_decay_kw(flag, case, mcomp=2))
    hand = _hand_decay(flag, case, mcomp=2)
    assert _parse(_written(bct, function_tmpdir)) == _parse(hand)


@GRIDS
@DECAY_FLAGS
@DECAY_PARAMS
def test_bct_decay_arrays_loaded(function_tmpdir, case, flag, structured):
    """Hand-written file -> load -> every array at its place -> same file."""
    ml = _model(function_tmpdir, structured=structured)
    hand = _hand_decay(flag, case, mcomp=2)
    bct = MfUsgBct.load(_hand(function_tmpdir, hand), ml)
    read = _decay_read(flag, case)
    assert [n for n in DECAY[flag] if hasattr(bct, n)] == read
    for ic in range(2):
        assert [_value(getattr(bct, n)[ic]) for n in read] == [
            round(DVAL[n] + ic, 5) for n in read
        ]
        assert _value(bct.conc[ic]) == CONC + ic
    assert _parse(_written(bct, function_tmpdir)) == _parse(hand)


@DECAY_FLAGS
@pytest.mark.parametrize("case", ["3_noads", "4_noaw_noads", "4_noaw_ads"])
def test_bct_decay_arrays_not_read_not_required(function_tmpdir, case, flag):
    """Arrays the engine does not read are not required and not kept."""
    ml = _model(function_tmpdir)
    kw = {k: v for k, v in _decay_kw(flag, case).items() if k not in DECAY[flag]}
    bct = _bct(ml, **kw)
    assert not [n for n in DECAY[flag] if hasattr(bct, n)]
    assert _parse(_written(bct, function_tmpdir)) == _parse(_hand_decay(flag, case))


def _set(**values):
    def change(bct):
        for name, value in values.items():
            setattr(bct, name, value)

    return change


@pytest.mark.parametrize(
    "kw, change, match",
    [
        ({}, _set(izod=1), "zodrw"),
        ({}, _set(ifod=1), "fodrw"),
        ({"izod": 1, "zodrw": 3.1}, _set(izod=3, iadsorb=1), "zodrs"),
        ({"ifod": 4, "fodrw": 3.2}, _set(aw_adsorb=1), "fodrw, fodrs, fodraw"),
        ({}, _set(solubility=1), "sollim, solslope"),
    ],
    ids=["izod", "ifod", "izod_iadsorb", "ifod_aw", "solubility"],
)
@pytest.mark.parametrize("target", ["file", "handle"])
def test_bct_write_rejects_flags_without_arrays(
    function_tmpdir, target, kw, change, match
):
    """Flags changed after construction: rejected before anything is written."""
    ml = _model(function_tmpdir)
    bct = _bct(ml, **kw)
    bct.fn_path = str(function_tmpdir / "keep.bct")
    Path(bct.fn_path).write_text(SENTINEL)
    change(bct)
    handle = io.StringIO() if target == "handle" else None
    with pytest.raises(ValueError, match=match):
        bct.write_file(f=handle)
    assert Path(bct.fn_path).read_text() == SENTINEL
    if handle is not None:
        assert handle.getvalue() == ""


# ----------------------------------------------------------------- SOLUBILITY

SOL = " SOLUBILITY"
INTERNAL = "INTERNAL 1.0 (FREE) -1\n"
# name: (structured, FloPy arguments, hand-written file)
SOLUBILITY_CASES = {
    "one_species": (
        True,
        {"sollim": 5.0, "solslope": 0.5},
        _head(opts=SOL) + _consts(PRS, 5.0, 0.5, CONC),
    ),
    "two_species": (
        False,
        {"mcomp": 2, "sollim": [5.0, 7.0], "solslope": 0.5, "conc": [CONC, 22.0]},
        _head(mcomp=2, opts=SOL)
        + _consts(PRS)
        + INTERNAL
        + " 5.0 7.0\n"
        + _consts(0.5, CONC, 22.0),
    ),
    # per-species arrays (ADSORB, FODRW) come after SOLLIM and SOLSLOPE
    "adsorb_ifod1": (
        True,
        {
            "iadsorb": 1,
            "ifod": 1,
            "bulkd": BULK,
            "adsorb": ADS,
            "fodrw": 3.2,
            "sollim": 5.0,
            "solslope": 0.5,
        },
        _head(iadsorb=1, ifod=1, opts=SOL)
        + _consts(PRS, BULK, 5.0, 0.5, ADS, 3.2, CONC),
    ),
    # after the heat arrays, before CONC and temperature
    "heat_disu": (
        False,
        {**HEAT, "bulkd": BULK, "sollim": 5.0, "solslope": 0.5},
        _head(iheat=1, opts=SOL)
        + HEAT_LINE
        + _consts(PRS, BULK, HTCAPS, HTCONDS, 5.0, 0.5, CONC, TEMP),
    ),
    # after AWAMAX (AW_ADSORB1RP1), before ALANGAW and BLANGAW (AW_ADSORB1RP2)
    "aw_adsorb": (
        True,
        {
            "aw_adsorb": True,
            "iarea_fn": 1,
            "ikawi_fn": 1,
            "awamax": AWA,
            "alangaw": ALA,
            "blangaw": BLA,
            "sollim": 5.0,
            "solslope": 0.5,
        },
        _head(opts=" A-W_ADSORB 1 1" + SOL)
        + _consts(PRS, AWA, 5.0, 0.5, ALA, BLA, CONC),
    ),
}
SOLUBILITY_PARAMS = pytest.mark.parametrize("case", list(SOLUBILITY_CASES))


@SOLUBILITY_PARAMS
def test_bct_solubility_written(function_tmpdir, case):
    structured, kw, hand = SOLUBILITY_CASES[case]
    ml = _model(function_tmpdir, structured=structured)
    bct = _bct(ml, solubility=True, **kw)
    got = _parse(_written(bct, function_tmpdir))
    want = _parse(hand)
    assert sorted(got[1]) == sorted(want[1])  # options, any order
    assert got[:1] + got[2:] == want[:1] + want[2:]


@SOLUBILITY_PARAMS
def test_bct_solubility_loaded(function_tmpdir, case):
    structured, kw, hand = SOLUBILITY_CASES[case]
    ml = _model(function_tmpdir, structured=structured)
    bct = MfUsgBct.load(_hand(function_tmpdir, hand), ml)
    mcomp = kw.get("mcomp", 1)
    assert bct.solubility == 1
    assert np.allclose(bct.sollim.array, np.broadcast_to(kw["sollim"], mcomp))
    assert np.allclose(bct.solslope.array, np.broadcast_to(kw["solslope"], mcomp))
    assert [_value(c) for c in bct.conc] == np.broadcast_to(
        kw.get("conc", CONC), mcomp
    ).tolist()
    got = _parse(_written(bct, function_tmpdir))
    want = _parse(hand)
    assert sorted(got[1]) == sorted(want[1])
    assert got[:1] + got[2:] == want[:1] + want[2:]


def test_bct_without_solubility_unchanged(function_tmpdir):
    """Control: no SOLUBILITY, no SOLLIM or SOLSLOPE record."""
    ml = _model(function_tmpdir)
    bct = _bct(ml)
    assert bct.solubility is None
    assert _parse(_written(bct, function_tmpdir)) == _parse(
        _head() + _consts(PRS, CONC)
    )


# ------------------------------------------- BULKD: structured grid with heat

HAND_HEAT_DIS = _head(iheat=1) + HEAT_LINE + _consts(PRS, HTCAPS, HTCONDS, CONC, TEMP)
HAND_HEAT_DISU = (
    _head(iheat=1) + HEAT_LINE + _consts(PRS, BULK, HTCAPS, HTCONDS, CONC, TEMP)
)
HAND_HEAT_DIS_ADS = (
    _head(iheat=1, iadsorb=1)
    + HEAT_LINE
    + _consts(PRS, BULK, HTCAPS, HTCONDS, ADS, CONC, TEMP)
)


@pytest.mark.parametrize("target", ["file", "handle"])
def test_bct_structured_heat_without_adsorption_write_guard(function_tmpdir, target):
    """Guard (not support): rejected before the file or handle gets anything."""
    ml = _model(function_tmpdir, structured=True)
    bct = _bct(ml, bulkd=BULK, **HEAT)
    bct.fn_path = str(function_tmpdir / "keep.bct")
    Path(bct.fn_path).write_text(SENTINEL)
    handle = io.StringIO() if target == "handle" else None
    with pytest.raises(NotImplementedError, match="structured"):
        bct.write_file(f=handle)
    assert Path(bct.fn_path).read_text() == SENTINEL
    if handle is not None:
        assert handle.getvalue() == ""


def test_bct_structured_heat_without_adsorption_load_guard(function_tmpdir):
    """Hand-written structured file as the engine reads it (no BULKD)."""
    ml = _model(function_tmpdir, structured=True)
    before = _state(ml)
    with pytest.raises(NotImplementedError, match="structured"):
        MfUsgBct.load(_hand(function_tmpdir, HAND_HEAT_DIS), ml)
    assert _state(ml) == before


@pytest.mark.parametrize(
    "structured, kw, hand",
    [
        (False, {}, HAND_HEAT_DISU),
        (True, {"iadsorb": 1, "adsorb": ADS}, HAND_HEAT_DIS_ADS),
    ],
    ids=["disu", "dis_adsorb"],
)
def test_bct_heat_bulkd_still_supported(function_tmpdir, structured, kw, hand):
    """Controls: DISU with heat, and structured with heat and adsorption."""
    ml = _model(function_tmpdir, structured=structured)
    bct = _bct(ml, bulkd=BULK, **HEAT, **kw)
    assert _parse(_written(bct, function_tmpdir)) == _parse(hand)
    ml2 = _model(function_tmpdir, name="r", structured=structured)
    bct2 = MfUsgBct.load(_hand(function_tmpdir, hand), ml2)
    assert (_value(bct2.bulkd), _value(bct2.heat)) == (BULK, TEMP)
    assert _parse(_written(bct2, function_tmpdir, "rewritten.bct")) == _parse(hand)


# ------------------------------------------------------- VARIABLE_SOLUBILITY

# per species: SOLLIMVAR, SOLSLOPEVAR, then CONC
HAND_VARIABLE_SOLUBILITY = _head(opts=" VARIABLE_SOLUBILITY") + _consts(
    PRS, 5.0, 0.5, CONC
)


def test_bct_variable_solubility_load_guard(function_tmpdir):
    """Guard (not support): raised from item 1, before any array is read."""
    ml = _model(function_tmpdir)
    before = _state(ml)
    with pytest.raises(NotImplementedError, match="VARIABLE_SOLUBILITY"):
        MfUsgBct.load(_hand(function_tmpdir, HAND_VARIABLE_SOLUBILITY), ml)
    assert _state(ml) == before


@pytest.mark.parametrize(
    "kw", [{"variable_solubility": True}, {"solubility": 2}], ids=["keyword", "value2"]
)
def test_bct_variable_solubility_construction_guard(function_tmpdir, kw):
    ml = _model(function_tmpdir)
    before = _state(ml)
    with pytest.raises(NotImplementedError, match="VARIABLE_SOLUBILITY"):
        _bct(ml, sollim=5.0, solslope=0.5, **kw)
    assert _state(ml) == before


# ------------------------------------------------------------ engine (requires_usgt)

ENG_GRID = _consts(0.2, 0.0, 0.0, 0.0)  # PRSITY, ANGLEX, DL, DT


def _eng_head(iadsorb=0, izod=0, ifod=0, opts=""):
    return (
        f" 1 0 1 1 1 {iadsorb} 0 -999.0 1.0e-6 1 0 0.0 {izod} {ifod} 0 0 0 0 0{opts}\n"
    )


def _eng_run(ws, bct_kw=None, hand=None):
    """Run a 5-node DISU model; BCT from FloPy (bct_kw) or hand-written.

    Returns (normal termination, bytes of the concentration file or None).
    """
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
    kw = {"prsity": 0.2, "idisp": 1, "dl": 0.0, "dt": 0.0, "conc": 10.0}
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
    ok = success and any("normal termination" in b.lower() for b in buff)
    return ok, (Path(ml.model_ws) / "m.con").read_bytes() if ok else None


def _conc(data):
    """Values of the last record of an unstructured single-precision CON file."""
    out, off = None, 0
    while off < len(data):
        n1, n2 = struct.unpack_from("<ii", data, off + 32)
        n = n2 - n1 + 1
        out = np.array(struct.unpack_from(f"<{n}f", data, off + 44))
        off += 44 + 4 * n
    return out


def _eng_solubility(sollim=5.0, solslope=2.0):
    kw = {"ifod": 1, "fodrw": 0.05, "solubility": True}
    return {**kw, "sollim": sollim, "solslope": solslope}


@requires_usgt
def test_bct_engine_solubility_flopy_equals_hand(function_tmpdir):
    hand = (
        _eng_head(ifod=1, opts=SOL)
        + ENG_GRID
        + _consts(5.0, 2.0, 0.05, 10.0)  # SOLLIM, SOLSLOPE, FODRW, CONC
    )
    fa = _eng_run(function_tmpdir / "flopy", _eng_solubility())
    fb = _eng_run(function_tmpdir / "hand", hand=hand)
    assert fb[0], "hand-written file did not run"
    assert fa[0] and fa[1] == fb[1]


@requires_usgt
def test_bct_engine_solubility_sensitivity(function_tmpdir):
    """SOLLIM and SOLSLOPE each change the result (the solute decays from 10)."""
    runs = {
        "ref": _eng_solubility(),
        "limit": _eng_solubility(sollim=9.9),  # crossed during the run
        "slope": _eng_solubility(solslope=0.5),
        "off": {"ifod": 1, "fodrw": 0.05},
    }
    out = {k: _eng_run(function_tmpdir / k, kw) for k, kw in runs.items()}
    assert all(ok for ok, _ in out.values())
    c = {k: _conc(v[1]) for k, v in out.items()}
    for k in ("limit", "slope", "off"):
        assert np.all(np.abs(c[k] - c["ref"]) > 1e-3), k
    # precipitated mass above the limit buffers the decay of the solute
    assert np.all(c["ref"] > c["off"])


@requires_usgt
def test_bct_engine_izod4_without_aw_flopy_equals_hand(function_tmpdir):
    """IZOD=4 without A-W_ADSORB: no decay array is read (as IZOD=0)."""
    hand = _eng_head(izod=4, ifod=1) + ENG_GRID + _consts(0.05, 10.0)  # FODRW, CONC
    fa = _eng_run(
        function_tmpdir / "flopy", {"izod": 4, "zodrw": 0.5, "ifod": 1, "fodrw": 0.05}
    )
    fb = _eng_run(function_tmpdir / "hand", hand=hand)
    fc = _eng_run(function_tmpdir / "izod0", {"ifod": 1, "fodrw": 0.05})
    assert fb[0] and fc[0], "hand-written or IZOD=0 file did not run"
    assert fb[1] == fc[1]
    assert fa[0] and fa[1] == fb[1]


@requires_usgt
def test_bct_engine_decay_sensitivity(function_tmpdir):
    """The results depend on the decay data and on the place of each array."""
    out = [
        _eng_run(function_tmpdir / f"z{i}", {"izod": 1, "zodrw": v})
        for i, v in enumerate((0.0, 0.05))
    ]
    assert out[0][0] and out[1][0]
    assert np.allclose(_conc(out[0][1]), 10.0)
    assert np.all(_conc(out[1][1]) < 10.0 - 1e-3)
    # IZOD=4 without A-W_ADSORB with one array too many (ZODRW before FODRW)
    good = _eng_head(izod=4, ifod=1) + ENG_GRID + _consts(0.05, 10.0)
    extra = _eng_head(izod=4, ifod=1) + ENG_GRID + _consts(0.5, 0.05, 10.0)
    fg = _eng_run(function_tmpdir / "good", hand=good)
    fx = _eng_run(function_tmpdir / "extra", hand=extra)
    assert fg[0]
    assert not fx[0] or fx[1] != fg[1]


@requires_usgt
def test_bct_engine_izod4_with_aw_flopy_equals_hand(function_tmpdir):
    """IZOD=4 with A-W_ADSORB and IADSORB=0: ZODRW, ZODRS and ZODRAW are read."""
    aw = {"aw_adsorb": True, "iarea_fn": 1, "ikawi_fn": 1, "awamax": 100.0}
    aw.update(alangaw=0.01, blangaw=0.1)
    decay = {"izod": 4, "zodrw": 0.011, "zodrs": 0.012, "zodraw": 0.013}
    head = _eng_head(izod=4, ifod=1, opts=" A-W_ADSORB 1 1") + ENG_GRID
    # AWAMAX; ALANGAW, BLANGAW; ZODRW, ZODRS, ZODRAW; FODRW; CONC
    hand = head + _consts(100.0, 0.01, 0.1, 0.011, 0.012, 0.013, 0.05, 10.0)
    short = head + _consts(100.0, 0.01, 0.1, 0.011, 0.013, 0.05, 10.0)  # no ZODRS
    fa = _eng_run(function_tmpdir / "flopy", {**aw, **decay, "ifod": 1, "fodrw": 0.05})
    fb = _eng_run(function_tmpdir / "hand", hand=hand)
    fs = _eng_run(function_tmpdir / "short", hand=short)
    assert fb[0], "hand-written file did not run"
    assert fa[0] and fa[1] == fb[1]
    assert not fs[0] or fs[1] != fb[1]
