"""MfUsgMdt in a fixed-format model (BAS6 without FREE).

The expected files are written by hand from GWT2MDTU1AR (USG-T 2.8.1
gwt2mdtu1.for). With IFREFM=0 the header is ``READ(IN,'(I10)') IMDTCB,
IMDTCF`` (lines 76-78, 85): by format reversion IMDTCF is the I10 field of the
next record, comment lines are not skipped, and the options (FRAHK, FRADARCY,
TSHIFTMD, SEPARATE_AI2, MULTIFILE_MD) are never read (lines 98-127 parse a
LINE that this branch does not read). Engine tests (requires_usgt) compare a
FloPy file with a hand-written file on a one-layer DISU model.
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
    MfUsgMdt,
    MfUsgOc,
    MfUsgSms,
)

# one value per array, in the GWT2MDTU1AR order (IDPF=0, one component):
# MDFLAG, VOLFRACMD, PORMD, RHOBMD, DIFFLENMD, TORTMD, KDMD, DECAYMD,
# YIELDMD, DIFFMD
SEQ = [2, 0.41, 0.21, 1.61, 0.51, 0.71, 1.1, 2.1, 3.1, 4.1]
KW = {
    "mdflag": 2,
    "volfracmd": 0.41,
    "pormd": 0.21,
    "rhobmd": 1.61,
    "difflenmd": 0.51,
    "tortmd": 0.71,
    "kdmd": [1.1],
    "decaymd": [2.1],
    "yieldmd": [3.1],
    "diffmd": [4.1],
}
SENTINEL = "SENTINEL MDT FILE\n"
OPTIONS = {
    "frahk": {"frahk": True},
    "fradarcy": {"fradarcy": True},
    "tshiftmd": {"tshiftmd": 2.0, "aiold1md": 5.1, "aiold2md": 6.1},
    "separate_ai2": {"iunitAI2": 60},
    "multifile_md": {"crootname": "md", "imdtcf": 61},
}


def _i10(value):
    return f"{value:10d}\n"


def _body(values=SEQ):
    return "".join(f"CONSTANT {v}\n" for v in values)


def _model(ws, name, free=False, idpf=0):
    m = MfUsg(modelname=name, model_ws=str(ws), structured=False)
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
    MfUsgBas(m, ifrefm=free)
    m.itrnsp = 1
    m.mcomp, m.iheat, m.imcomp = 1, 0, 0
    m.idpf = idpf
    return m


def _values(lines):
    vals = []
    for ln in lines:
        t = ln.split()
        assert t[0].upper() == "CONSTANT", ln
        vals.append(float(t[1]))
    return vals


def _loaded_values(mdt):
    arrays = [mdt.mdflag, mdt.volfracmd, mdt.pormd, mdt.rhobmd, mdt.difflenmd]
    arrays += [mdt.tortmd, mdt.kdmd[0], mdt.decaymd[0], mdt.yieldmd[0], mdt.diffmd[0]]
    return [float(a.array.ravel()[0]) for a in arrays]


# ------------------------------------------------------------------ write


def test_mdt_fixed_writes_two_i10_records(function_tmpdir):
    mdt = MfUsgMdt(_model(function_tmpdir, "w"), ipakcb=78, imdtcf=77, **KW)
    mdt.write_file()
    lines = Path(mdt.fn_path).read_text().splitlines()
    assert lines[:2] == ["        78", "        77"]
    assert not any(ln.startswith("#") for ln in lines)
    assert np.allclose(_values(lines[2:]), SEQ, rtol=1e-6)


def test_mdt_fixed_flopy_equals_hand(function_tmpdir):
    hand = _i10(78) + _i10(77) + _body()
    mdt = MfUsgMdt(_model(function_tmpdir, "f"), ipakcb=78, imdtcf=77, **KW)
    mdt.write_file()
    text = Path(mdt.fn_path).read_text()
    lines, hlines = text.splitlines(), hand.splitlines()
    assert lines[:2] == hlines[:2]
    assert len(lines) == len(hlines)
    assert np.allclose(_values(lines[2:]), _values(hlines[2:]), rtol=1e-6)
    # load the FloPy file and write it again: same records and values (the
    # array-name comments of loaded arrays differ, as in free format)
    re = MfUsgMdt.load(mdt.fn_path, _model(function_tmpdir, "f2"))
    re.fn_path = str(function_tmpdir / "f2.mdt")
    re.write_file()
    relines = Path(re.fn_path).read_text().splitlines()
    assert relines[:2] == lines[:2] and len(relines) == len(lines)
    assert _values(relines[2:]) == _values(lines[2:])


@pytest.mark.parametrize("name", list(OPTIONS))
def test_mdt_fixed_option_rejected_at_construction(function_tmpdir, name):
    with pytest.raises(ValueError, match="fixed format"):
        MfUsgMdt(_model(function_tmpdir, "c"), **{**KW, **OPTIONS[name]})


@pytest.mark.parametrize("name", list(OPTIONS))
def test_mdt_fixed_option_rejected_at_write(function_tmpdir, name):
    m = _model(function_tmpdir, "r", free=True)
    mdt = MfUsgMdt(m, **{**KW, **OPTIONS[name]})
    m.get_package("BAS6").ifrefm = False  # fixed format after construction
    p = Path(mdt.fn_path)
    p.write_text(SENTINEL)
    with pytest.raises(ValueError, match="fixed format"):
        mdt.write_file()
    assert p.read_text() == SENTINEL
    buf = io.StringIO("PREFIX\n")
    buf.seek(0, io.SEEK_END)
    pos = buf.tell()
    with pytest.raises(ValueError, match="fixed format"):
        mdt.write_file(f=buf)
    assert not buf.closed and buf.getvalue() == "PREFIX\n" and buf.tell() == pos


# ------------------------------------------------------------------- load


def test_mdt_fixed_hand_file_loads(function_tmpdir):
    p = function_tmpdir / "hand.mdt"
    p.write_text(_i10(78) + _i10(77) + _body())
    mdt = MfUsgMdt.load(str(p), _model(function_tmpdir, "l"))
    assert (mdt.ipakcb, mdt.imdtcf) == (78, 77)
    assert np.allclose(_loaded_values(mdt), SEQ, rtol=1e-6)


@pytest.mark.parametrize(
    "header",
    ["# MDT heading\n" + f"{0:9d} {0:9d}\n", f"{0:9d} {0:9d}\n"],
    ids=["heading", "one_line"],
)
def test_mdt_fixed_load_rejects_free_header(function_tmpdir, header):
    # the header FloPy wrote before (with and without the heading): USG-T
    # stops at gwt2mdtu1.for:77 on both in a fixed-format model
    p = function_tmpdir / "free.mdt"
    p.write_text(header + _body())
    with pytest.raises(ValueError, match="fixed format"):
        MfUsgMdt.load(str(p), _model(function_tmpdir, "x"))


def test_mdt_fixed_load_ignores_text_after_column_10(function_tmpdir):
    # the I10 field is columns 1-10; USG-T does not read the rest of the
    # record, so TSHIFTMD there is not set and no AIOLD arrays are read
    p = function_tmpdir / "opt.mdt"
    p.write_text("{:10d} TSHIFTMD 2.0\n".format(0) + _i10(0) + _body())
    mdt = MfUsgMdt.load(str(p), _model(function_tmpdir, "o"))
    assert (mdt.ipakcb, mdt.imdtcf, mdt.tshiftmd) == (0, 0, 0.0)
    assert np.allclose(_loaded_values(mdt), SEQ, rtol=1e-6)


@pytest.mark.parametrize(
    "header, expected",
    [
        ("\n" + _i10(77), (0, 77)),
        (" " * 10 + "\n" + _i10(77), (0, 77)),
        (_i10(78) + "\n", (78, 0)),
        (_i10(78) + " " * 10 + "\n", (78, 0)),
        ("\n\n", (0, 0)),
    ],
    ids=["rec1_empty", "rec1_spaces", "rec2_empty", "rec2_spaces", "both_empty"],
)
def test_mdt_fixed_load_blank_field_is_zero(function_tmpdir, header, expected):
    # a blank I10 field (empty record or spaces) is read as 0 by USG-T
    p = function_tmpdir / "blank.mdt"
    p.write_text(header + _body())
    mdt = MfUsgMdt.load(str(p), _model(function_tmpdir, "b"))
    assert (mdt.ipakcb, mdt.imdtcf) == expected
    assert np.allclose(_loaded_values(mdt), SEQ, rtol=1e-6)


@pytest.mark.parametrize(
    "header",
    ["  7  8    \n" + _i10(77), _i10(78) + "  7  7    \n", "78\n77\n"],
    ids=["rec1_inner_blanks", "rec2_inner_blanks", "left_justified"],
)
def test_mdt_fixed_load_blanks_in_field_ignored(function_tmpdir, header):
    # USG-T ignores the blanks inside an I10 field ("  7  8" is 78) and pads
    # a short record with blanks
    p = function_tmpdir / "inner.mdt"
    p.write_text(header + _body())
    mdt = MfUsgMdt.load(str(p), _model(function_tmpdir, "i"))
    assert (mdt.ipakcb, mdt.imdtcf) == (78, 77)
    assert np.allclose(_loaded_values(mdt), SEQ, rtol=1e-6)


@pytest.mark.parametrize(
    "header", ["abc\n" + _i10(0), _i10(78)], ids=["text", "missing_record"]
)
def test_mdt_fixed_load_rejects_bad_field(function_tmpdir, header):
    # text in the I10 field stops USG-T at gwt2mdtu1.for:77; the second
    # record is missing at the end of the file
    p = function_tmpdir / "bad.mdt"
    p.write_text(header)
    with pytest.raises(ValueError, match="fixed format"):
        MfUsgMdt.load(str(p), _model(function_tmpdir, "e"))


def test_mdt_fixed_idpf_skips_volfracmd(function_tmpdir):
    # with IDPF != 0 the header is the same two I10 records and VOLFRACMD is
    # not read (gwt2mdtu1.for 2.8.1:198-211)
    seq = SEQ[:1] + SEQ[2:]
    mdt = MfUsgMdt(_model(function_tmpdir, "d", idpf=1), ipakcb=78, imdtcf=77, **KW)
    mdt.write_file()
    lines = Path(mdt.fn_path).read_text().splitlines()
    assert lines[:2] == ["        78", "        77"]
    assert np.allclose(_values(lines[2:]), seq, rtol=1e-6)
    hand = function_tmpdir / "hand.mdt"
    hand.write_text(_i10(78) + _i10(77) + _body(seq))
    re = MfUsgMdt.load(str(hand), _model(function_tmpdir, "d2", idpf=1))
    assert (re.ipakcb, re.imdtcf) == (78, 77)
    re.fn_path = str(function_tmpdir / "d2.mdt")
    re.write_file()
    relines = Path(re.fn_path).read_text().splitlines()
    assert relines[:2] == lines[:2] and len(relines) == len(lines)
    assert _values(relines[2:]) == _values(lines[2:])


# ------------------------------------------------------- engine (requires_usgt)

ENG = {
    "mdflag": 2,
    "volfracmd": 0.5,
    "pormd": 0.3,
    "rhobmd": 1.6,
    "difflenmd": 0.5,
    "tortmd": 0.5,
    "kdmd": [0.1],
    "decaymd": [1.0e-3],
    "yieldmd": [1.0],
    "diffmd": [1.0e-2],
}
ENG_SEQ = [2, 0.5, 0.3, 1.6, 0.5, 0.5, 0.1, 1.0e-3, 1.0, 1.0e-2]


def _eng_model(ws, free):
    """5-node one-layer DISU line, one solute, no boundaries."""
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
    MfUsgBas(ml, ibound=1, strt=9.5, ifrefm=free)
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
    assert success and normal
    ws = Path(ml.model_ws)
    return _records((ws / "m.con").read_bytes()), (ws / "m.hds").read_bytes()


def _flopy_run(ws, free=False, **mdt_kw):
    ml = _eng_model(ws, free)
    MfUsgMdt(ml, **{**ENG, **mdt_kw})
    ml.write_input()
    return _run(ml)


@requires_usgt
def test_mdt_engine_fixed_flopy_equals_hand(function_tmpdir):
    con, hds = _flopy_run(function_tmpdir / "flopy")
    lines = (function_tmpdir / "flopy" / "m.mdt").read_text().splitlines()
    assert lines[:2] == ["         0", "         0"]
    # hand-written two-record file on the same fixed-format model
    ml = _eng_model(function_tmpdir / "hand", free=False)
    mdt = MfUsgMdt(ml)
    ml.write_input(
        SelPackList=[p.name[0] for p in ml.packagelist if p.name[0] != "MDT"]
    )
    Path(mdt.fn_path).write_text(_i10(0) + _i10(0) + _body(ENG_SEQ))
    hand, hand_hds = _run(ml)
    assert len(con) == len(hand) == 1
    assert np.array_equal(con[0], hand[0]) and hds == hand_hds
    # the same data in a free-format model gives the same concentrations
    free = _flopy_run(function_tmpdir / "free", free=True)[0]
    assert np.array_equal(con[0], free[0])
    assert np.all(con[0] < 10.0)  # mass diffuses into the matrix


@requires_usgt
def test_mdt_engine_fixed_sensitivity(function_tmpdir):
    base = _flopy_run(function_tmpdir / "base")[0]
    changed = _flopy_run(function_tmpdir / "changed", diffmd=[5.0e-2])[0]
    assert np.all(changed[0] < base[0])  # more diffusion into the matrix
