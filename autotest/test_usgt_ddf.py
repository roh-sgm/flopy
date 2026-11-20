"""MfUsgDdf: DDF1AR / SHARPINT1AR contract (USG-T 2.8.1 and 2.7.0 density.f).

DDF1AR reads, from the DDF unit (BAS FREE required):
  1. RHOFRESH RHOSTD CSTD (REAL) ITHICKAV IMPHDD ISHARP (INTEGER)
     [NONLINEAR NDDFTABROWS] on one record;
  2. if NONLINEAR: NDDFTABROWS list-directed records (CONC DENSITY);
  3. if ISHARP != 0 (SHARPINT1AR): one record HSEA (REAL) ISWIUN (INTEGER),
     then, when CLN is active, ZETASWI for the CLN nodes (U1DREL8).
ISWIUN > 0 saves the interface elevation (text "INTERFACE ELEV.") with the
head OC flags. DDF-B1 implements 1-3 without CLN; ISHARP with CLN stays
explicitly unsupported (engine ZETANODE CLN branch uses an undefined index).
"""

import builtins
from pathlib import Path

import numpy as np
import pytest

from autotest.usgt_exe import USGT_EXE, requires_usgt
from flopy.mfusg import (
    MfUsg,
    MfUsgBas,
    MfUsgBct,
    MfUsgCln,
    MfUsgDdf,
    MfUsgDisU,
    MfUsgLpf,
    MfUsgOc,
    MfUsgSms,
)
from flopy.modflow import ModflowChd, ModflowWel

TABLE = [(0.0, 1000.0), (17.123456, 1012.345678), (35.0, 1025.0)]


def _model(ws, bct=False, free=True, exe=None):
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
        cl12=50.0,
        fahl=100.0,
        top=10.0,
        bot=-40.0,
        area=10000.0,
        nper=1,
        perlen=1.0,
        nstp=1,
        steady=True,
    )
    MfUsgBas(ml, ibound=1, strt=[0.0, 0.3, 0.6, 0.9, 1.2], ifrefm=free)
    if bct:
        MfUsgBct(
            ml,
            itrnsp=1,
            mcomp=1,
            prsity=0.2,
            idisp=1,
            dl=0.0,
            dt=0.0,
            diffnc=0.0,
            conc=[0.0357, 0.03, 0.02, 0.01, 0.0],
        )
    return ml


def _data_lines(path):
    return [
        ln.split()
        for ln in Path(path).read_text().splitlines()
        if ln.strip() and not ln.lstrip().startswith("#")
    ]


def _cln(ml):
    return MfUsgCln(
        ml,
        ncln=1,
        iclnnds=-1,
        nndcln=1,
        nclngwc=1,
        node_prop=[[1, 1, 0, 10.0, -30.0, 1.57, 0, 0]],
        cln_gwc=[[1, 1, 50, 50, 0, 0, 10.0, 1.0, 0]],
        cln_circ=[[1, 0.5, 3.23e10]],
    )


# ---------------------------------------------------------------- write


def test_ddf_default_single_record_control(function_tmpdir):
    """Control: isharp=0 writes only the header record."""
    ml = _model(function_tmpdir)
    d = MfUsgDdf(ml, ithickav=0)
    d.write_file()
    lines = _data_lines(d.fn_path)
    assert len(lines) == 1
    assert lines[0][3:] == ["0", "0", "0"]


@pytest.mark.parametrize("isharp", [1, 2], ids=["isharp1", "isharp2"])
def test_ddf_sharp_writes_hsea_iswiun_record(function_tmpdir, isharp):
    ml = _model(function_tmpdir)
    d = MfUsgDdf(ml, ithickav=0, isharp=isharp, hsea=0.25, iswiun=60)
    d.write_file()
    lines = _data_lines(d.fn_path)
    assert len(lines) == 2
    assert int(lines[0][5]) == isharp
    assert float(lines[1][0]) == 0.25 and int(lines[1][1]) == 60


def test_ddf_sharp_record_follows_nonlinear_table(function_tmpdir):
    ml = _model(function_tmpdir)
    d = MfUsgDdf(
        ml,
        ithickav=0,
        isharp=1,
        hsea=0.0,
        iswiun=0,
        nonlinear=True,
        density_table=TABLE,
    )
    d.write_file()
    lines = _data_lines(d.fn_path)
    assert lines[0][6:] == ["NONLINEAR", "3"]
    assert [tuple(map(float, r)) for r in lines[1:4]] == TABLE
    assert [float(lines[4][0]), int(lines[4][1])] == [0.0, 0]


def test_ddf_sharp_requires_hsea(function_tmpdir):
    ml = _model(function_tmpdir)
    with pytest.raises(ValueError, match="HSEA"):
        MfUsgDdf(ml, isharp=1)


@pytest.mark.parametrize("kw", [{"hsea": 0.1}, {"iswiun": 60}], ids=["hsea", "iswiun"])
def test_ddf_sharp_values_without_isharp_rejected(function_tmpdir, kw):
    ml = _model(function_tmpdir)
    with pytest.raises(ValueError, match="ISHARP"):
        MfUsgDdf(ml, isharp=0, **kw)


@pytest.mark.parametrize("iswiun", [-1, 1.5, True], ids=["neg", "float", "bool"])
def test_ddf_iswiun_must_be_nonnegative_int(function_tmpdir, iswiun):
    ml = _model(function_tmpdir)
    with pytest.raises(ValueError, match="ISWIUN"):
        MfUsgDdf(ml, isharp=1, hsea=0.0, iswiun=iswiun)


def test_ddf_iswiun_registers_output_file(function_tmpdir):
    ml = _model(function_tmpdir)
    MfUsgDdf(ml, isharp=1, hsea=0.0, iswiun=60)
    assert 60 in ml.output_units
    assert ml.output_fnames[ml.output_units.index(60)] == "m.zta"
    ml2 = _model(function_tmpdir / "named")
    MfUsgDdf(ml2, isharp=1, hsea=0.0, iswiun=61, filenames=["a.ddf", "a.zeta"])
    assert ml2.output_fnames[ml2.output_units.index(61)] == "a.zeta"


def test_ddf_header_full_precision(function_tmpdir):
    ml = _model(function_tmpdir)
    d = MfUsgDdf(ml, rhofresh=999.973, rhostd=1024.996, cstd=0.0357, ithickav=0)
    d.write_file()
    assert [float(v) for v in _data_lines(d.fn_path)[0][:3]] == [
        999.973,
        1024.996,
        0.0357,
    ]


def test_ddf_table_full_precision(function_tmpdir):
    ml = _model(function_tmpdir)
    d = MfUsgDdf(ml, ithickav=0, nonlinear=True, density_table=TABLE)
    d.write_file()
    ml2 = _model(function_tmpdir / "re")
    assert MfUsgDdf.load(d.fn_path, ml2).density_table == TABLE


@pytest.mark.parametrize(
    "kw,match",
    [
        ({"nonlinear": True, "density_table": None}, "NONLINEAR"),
        ({"nonlinear": False, "density_table": TABLE}, "NONLINEAR"),
        (
            {"nonlinear": True, "density_table": [(0.0, 1000.0), (0.0, 1025.0)]},
            "increasing",
        ),
        (
            {"nonlinear": True, "density_table": [(35.0, 1025.0), (0.0, 1000.0)]},
            "increasing",
        ),
        ({"nonlinear": True, "density_table": [(0.0, 1000.0, 1.0)]}, "density_table"),
    ],
    ids=["empty", "without_flag", "equal_conc", "decreasing", "three_values"],
)
def test_ddf_density_table_validated(function_tmpdir, kw, match):
    ml = _model(function_tmpdir)
    with pytest.raises(ValueError, match=match):
        MfUsgDdf(ml, ithickav=0, **kw)


@pytest.mark.parametrize(
    "mutate,match",
    [
        (lambda d: setattr(d, "isharp", 1), "HSEA"),
        (lambda d: setattr(d, "density_table", TABLE), "NONLINEAR"),
        (
            lambda d: (
                setattr(d, "nonlinear", True),
                setattr(d, "density_table", [(0.0, 1000.0), ("x", 1.0)]),
            ),
            "density_table",
        ),
        (lambda d: setattr(d, "iswiun", 60), "ISHARP"),
    ],
    ids=[
        "isharp_without_hsea",
        "table_without_flag",
        "bad_row",
        "iswiun_without_sharp",
    ],
)
def test_ddf_mutation_validated_message(function_tmpdir, mutate, match):
    ml = _model(function_tmpdir)
    d = MfUsgDdf(ml, ithickav=0)
    mutate(d)
    with pytest.raises(ValueError, match=match):
        d.write_file()


MUTATIONS = {
    "isharp_without_hsea": lambda d: setattr(d, "isharp", 1),
    "table_without_flag": lambda d: setattr(d, "density_table", TABLE),
    "bad_row": lambda d: (
        setattr(d, "nonlinear", True),
        setattr(d, "density_table", [(0.0, 1000.0), ("x", 1.0)]),
    ),
    "iswiun_without_sharp": lambda d: setattr(d, "iswiun", 60),
    "cstd_zero": lambda d: setattr(d, "cstd", 0.0),
    "rhofresh_nan": lambda d: setattr(d, "rhofresh", float("nan")),
    "ithickav_fraction": lambda d: setattr(d, "ithickav", 0.5),
}


@pytest.mark.parametrize("name", list(MUTATIONS))
def test_ddf_mutation_keeps_destination(function_tmpdir, name):
    """Any rejected mutation leaves the previous file byte for byte."""
    ml = _model(function_tmpdir)
    d = MfUsgDdf(ml, ithickav=0)
    d.write_file()
    before = Path(d.fn_path).read_bytes()
    MUTATIONS[name](d)
    with pytest.raises(ValueError):
        d.write_file()
    assert Path(d.fn_path).read_bytes() == before


def test_ddf_fixed_format_not_supported(function_tmpdir):
    """DDF1AR reads 5F10.4 into integer items when BAS is not FREE (engine
    error 'Expected REAL for item 4'); FloPy refuses to write that file."""
    ml = _model(function_tmpdir, free=False)
    d = MfUsgDdf(ml, ithickav=0)
    with pytest.raises(NotImplementedError, match="FREE"):
        d.write_file()
    assert not Path(d.fn_path).exists()


def test_ddf_sharp_with_cln_not_supported(function_tmpdir):
    ml = _model(function_tmpdir)
    _cln(ml)
    d = MfUsgDdf(ml, isharp=1, hsea=0.0)
    with pytest.raises(NotImplementedError, match="CLN"):
        d.write_file()
    assert not Path(d.fn_path).exists()


def test_ddf_nonsharp_with_cln_writes(function_tmpdir):
    """Control: CLN only matters for the sharp-interface block."""
    ml = _model(function_tmpdir)
    _cln(ml)
    d = MfUsgDdf(ml, ithickav=0)
    d.write_file()
    assert len(_data_lines(d.fn_path)) == 1


# ---------------------------------------------------------------- load


@pytest.mark.parametrize(
    "text,expect",
    [
        (" 1000.0 1025.0 35.0 0 0 1\n 0.25 0\n", (1, 0.25, 0)),
        (" 1000.0 1025.0 35.0 0 0 2\n 0.0 60\n", (2, 0.0, 60)),
        (" 1000.0 1025.0 35.0 0 0 1\n# note\n -0.5 0\n", (1, -0.5, 0)),
        (" 1000.0 1025.0 35.0 0 0 1\n 0.1\n", (1, 0.1, 0)),
        (
            " 1000.0 1025.0 35.0 0 0 1 nonlinear 2\n 0 1000\n 35 1025\n 0.2 0\n",
            (1, 0.2, 0),
        ),
    ],
    ids=["iswiun0", "iswiun60", "comment", "iswiun_blank", "after_table"],
)
def test_ddf_sharp_load_hand_file(function_tmpdir, text, expect):
    src = function_tmpdir / "hand.ddf"
    src.write_text(text)
    ml = _model(function_tmpdir)
    d = MfUsgDdf.load(str(src), ml)
    assert (d.isharp, d.hsea, d.iswiun) == expect
    d.fn_path = str(function_tmpdir / "out.ddf")
    d.write_file()
    lines = _data_lines(d.fn_path)
    assert float(lines[-1][0]) == expect[1] and int(lines[-1][1]) == expect[2]


@pytest.mark.parametrize(
    "text,match",
    [
        (" 1000.0 1025.0 35.0 0 0 1\n", "HSEA"),
        (" 1000.0 1025.0 35.0 0 0 NONLINEAR 2\n 0 1000\n 35 1025\n", "ISHARP"),
        (" 1000.0 1025.0 35.0 0 0 0 NONLINEAR 2\n 0 1000\n", "NONLINEAR"),
        (" 1000.0 1025.0 35.0 0 0 0 NONLINEAR 2\n 0 1000\n 0 1025\n", "increasing"),
    ],
    ids=["sharp_missing_record", "keyword_in_isharp", "short_table", "table_order"],
)
def test_ddf_load_rejected_and_closed(function_tmpdir, text, match, monkeypatch):
    import flopy.mfusg.mfusgddf as ddfmod

    src = function_tmpdir / "bad.ddf"
    src.write_text(text)
    opened = []

    def tracking_open(*a, **k):
        fh = builtins.open(*a, **k)
        opened.append(fh)
        return fh

    monkeypatch.setattr(ddfmod, "open", tracking_open, raising=False)
    ml = _model(function_tmpdir)
    with pytest.raises(ValueError, match=match):
        MfUsgDdf.load(str(src), ml)
    assert opened and all(fh.closed for fh in opened)


def test_ddf_load_handle_left_open(function_tmpdir):
    """A caller-supplied handle is not closed by load."""
    src = function_tmpdir / "h.ddf"
    src.write_text(" 1000.0 1025.0 35.0 0 0 1\n 0.3 0\n")
    ml = _model(function_tmpdir)
    with open(src) as fh:
        d = MfUsgDdf.load(fh, ml)
        assert not fh.closed
    assert d.hsea == 0.3


def test_ddf_sharp_load_with_cln_not_supported(function_tmpdir):
    src = function_tmpdir / "s.ddf"
    src.write_text(" 1000.0 1025.0 35.0 0 0 1\n 0.0 0\n 1.0\n")
    ml = _model(function_tmpdir)
    _cln(ml)
    with pytest.raises(NotImplementedError, match="CLN"):
        MfUsgDdf.load(str(src), ml)


def test_ddf_load_henry_control(example_data_path):
    """Control: Henry.ddf has 4 tokens (ITHICKAV=01); missing integers are 0."""
    ml = MfUsg(modelname="m")
    henry = example_data_path / "mfusg_transport/Ex5_Henry/Henry.ddf"
    d = MfUsgDdf.load(str(henry), ml)
    assert (d.rhofresh, d.rhostd, d.cstd, d.ithickav, d.imphdd, d.isharp) == (
        1000.0,
        1025.0,
        35.0,
        1,
        0,
        0,
    )
    assert d.hsea is None and d.iswiun == 0


def test_ddf_linear_nonlinear_roundtrip_control(function_tmpdir):
    """Control: linear header and NONLINEAR table already round-trip."""
    ml = _model(function_tmpdir)
    d = MfUsgDdf(
        ml,
        ithickav=1,
        imphdd=1,
        nonlinear=True,
        density_table=[(0.0, 1000.0), (17.5, 1012.5), (35.0, 1025.0)],
    )
    d.write_file()
    d2 = MfUsgDdf.load(d.fn_path, _model(function_tmpdir / "re"))
    assert (d2.ithickav, d2.imphdd, d2.isharp, d2.nonlinear) == (1, 1, 0, True)
    assert d2.density_table == [(0.0, 1000.0), (17.5, 1012.5), (35.0, 1025.0)]


def test_ddf_nam_load_registers_iswiun(function_tmpdir):
    ml = _model(function_tmpdir)
    MfUsgLpf(ml, laytyp=4, hk=100.0)
    MfUsgSms(ml, linmeth=1)
    MfUsgOc(ml, stress_period_data={(0, 0): ["save head"]})
    MfUsgDdf(ml, isharp=1, hsea=0.0, iswiun=60)
    ml.write_input()
    m2 = MfUsg.load(
        "m.nam",
        model_ws=str(function_tmpdir),
        check=False,
        load_only=["disu", "bas6", "lpf", "sms", "oc", "ddf"],
    )
    assert (m2.ddf.isharp, m2.ddf.hsea, m2.ddf.iswiun) == (1, 0.0, 60)
    m2.change_model_ws(str(function_tmpdir / "re"))
    m2.write_input()
    nam = (function_tmpdir / "re" / "m.nam").read_text().splitlines()
    assert sum(" 60 " in f" {ln} " and "m.zta" in ln for ln in nam) == 1


# ---------------------------------------------------------------- check


def test_ddf_check_nonsharp_requires_bct(function_tmpdir):
    """Non-sharp DDF uses BCT CONC (engine SIGSEGV without BCT)."""
    ml = _model(function_tmpdir)
    chk = MfUsgDdf(ml, ithickav=0).check(verbose=False)
    assert "BCT" in " ".join(chk.summary_array.desc)
    ml2 = _model(function_tmpdir / "bct", bct=True)
    chk2 = MfUsgDdf(ml2, ithickav=0).check(verbose=False)
    assert "BCT" not in " ".join(chk2.summary_array.desc)


def test_ddf_check_bct_added_later_and_sharp(function_tmpdir):
    """check() looks at the model when it runs; sharp DDF does not use BCT."""
    ml = _model(function_tmpdir)
    d = MfUsgDdf(ml, ithickav=0)
    MfUsgBct(ml, itrnsp=1, mcomp=1, prsity=0.2, idisp=1, conc=0.0)
    assert "BCT" not in " ".join(d.check(verbose=False).summary_array.desc)
    ml2 = _model(function_tmpdir / "s")
    chk = MfUsgDdf(ml2, isharp=1, hsea=0.0).check(verbose=False)
    assert "BCT" not in " ".join(chk.summary_array.desc)


# ---------------------------------------------------------------- API (B1)


def test_ddf_positional_arguments_control(function_tmpdir):
    """Control: the existing positional order is unchanged."""
    ml = _model(function_tmpdir)
    d = MfUsgDdf(
        ml, 999.0, 1020.0, 30.0, 0, 1, 0, True, [(0.0, 999.0)], "dd2", 160, "x.dd2"
    )
    assert (d.rhofresh, d.rhostd, d.cstd, d.ithickav, d.imphdd, d.isharp) == (
        999.0,
        1020.0,
        30.0,
        0,
        1,
        0,
    )
    assert d.nonlinear and list(d.density_table) == [(0.0, 999.0)]
    assert d.unit_number == [160] and d.file_name == ["x.dd2"]


def test_ddf_positional_hsea_iswiun_after_filenames(function_tmpdir):
    ml = _model(function_tmpdir)
    d = MfUsgDdf(
        ml,
        1000.0,
        1025.0,
        35.0,
        0,
        0,
        1,
        False,
        None,
        "ddf",
        None,
        ["s.ddf", "s.zeta"],
        0.2,
        61,
    )
    assert (d.hsea, d.iswiun, d.file_name) == (0.2, 61, ["s.ddf"])
    assert ml.output_fnames[ml.output_units.index(61)] == "s.zeta"


NUMERIC_BAD = {
    "rhofresh_zero": ({"rhofresh": 0.0}, "RHOFRESH"),
    "rhofresh_nan": ({"rhofresh": float("nan")}, "RHOFRESH"),
    "rhofresh_text": ({"rhofresh": "1000"}, "RHOFRESH"),
    "rhostd_inf": ({"rhostd": float("inf")}, "RHOSTD"),
    "cstd_zero_linear": ({"cstd": 0.0}, "CSTD"),
    "cstd_zero_nonlinear": (
        {"cstd": 0.0, "nonlinear": True, "density_table": TABLE},
        "CSTD",
    ),
    "sharp_equal_densities": ({"isharp": 1, "hsea": 0.0, "rhostd": 1000.0}, "RHOSTD"),
    "hsea_nan": ({"isharp": 1, "hsea": float("nan")}, "HSEA"),
    "hsea_text": ({"isharp": 1, "hsea": "0.0"}, "HSEA"),
    "ithickav_float": ({"ithickav": 1.0}, "ITHICKAV"),
    "imphdd_text": ({"imphdd": "1"}, "IMPHDD"),
    "isharp_fraction": ({"isharp": 1.5, "hsea": 0.0}, "ISHARP"),
    "table_nan": (
        {"nonlinear": True, "density_table": [(0.0, float("nan"))]},
        "density_table",
    ),
    "table_bool": (
        {"nonlinear": True, "density_table": [(True, 1000.0)]},
        "density_table",
    ),
    "table_scalar_row": (
        {"nonlinear": True, "density_table": [1000.0]},
        "density_table",
    ),
}


@pytest.mark.parametrize("name", list(NUMERIC_BAD))
def test_ddf_numeric_values_validated(function_tmpdir, name):
    """Rejected before the package touches the model (no package, unit, output)."""
    kw, match = NUMERIC_BAD[name]
    ml = _model(function_tmpdir)
    state = (len(ml.packagelist), list(ml.package_units), list(ml.output_units))
    with pytest.raises(ValueError, match=match):
        MfUsgDdf(ml, **{"ithickav": 0, **kw})
    assert (len(ml.packagelist), ml.package_units, ml.output_units) == state


@pytest.mark.parametrize(
    "kw",
    [
        {"isharp": 1, "hsea": 0.0, "cstd": 0.0},
        {"rhostd": 1000.0},
        {"nonlinear": True, "density_table": [(0.0, 1000.0)]},
        {"isharp": -1, "hsea": -2.5, "iswiun": 62},
        {"ithickav": True, "imphdd": np.int64(1)},
    ],
    ids=["sharp_cstd0", "linear_equal_densities", "one_row", "isharp_neg", "int_types"],
)
def test_ddf_values_accepted_roundtrip(function_tmpdir, kw):
    """Values the engine reads without a zero denominator are not rejected."""
    ml = _model(function_tmpdir)
    d = MfUsgDdf(ml, **{"ithickav": 0, **kw})
    d.write_file()
    d2 = MfUsgDdf.load(d.fn_path, _model(function_tmpdir / "re"))
    names = ("rhofresh", "rhostd", "cstd", "ithickav", "imphdd", "isharp")
    names += ("nonlinear", "density_table", "hsea", "iswiun")
    assert [getattr(d2, n) for n in names] == [getattr(d, n) for n in names]


def test_ddf_iswiun_mutation_rejected(function_tmpdir):
    ml = _model(function_tmpdir)
    d = MfUsgDdf(ml, isharp=1, hsea=0.0, iswiun=60)
    d.write_file()
    before = Path(d.fn_path).read_bytes()
    d.iswiun = 61
    with pytest.raises(ValueError, match="ISWIUN changed"):
        d.write_file()
    assert 61 not in ml.output_units and 60 in ml.output_units
    d.iswiun = 60
    ml.remove_output(unit=60)
    with pytest.raises(ValueError, match="no longer"):
        d.write_file()
    assert Path(d.fn_path).read_bytes() == before


@pytest.mark.parametrize(
    "case", ["own_unit", "package_unit", "foreign_output", "unowned_output", "fname"]
)
def test_ddf_iswiun_collision_rejected(function_tmpdir, case):
    ml = _model(function_tmpdir)
    ml.add_output_file(70, fname="other.bin", package="XYZ")
    ml.add_output_file(71, fname="free.bin")
    kw = {"isharp": 1, "hsea": 0.0}
    kw.update(
        {
            "own_unit": {"iswiun": 152},
            "package_unit": {"iswiun": ml.bas6.unit_number[0]},
            "foreign_output": {"iswiun": 70},
            "unowned_output": {"iswiun": 71},
            "fname": {"iswiun": 60, "filenames": ["m.ddf", "other.bin"]},
        }[case]
    )
    state = (len(ml.packagelist), list(ml.output_units), list(ml.output_fnames))
    with pytest.raises(ValueError, match=r"ISWIUN|output file"):
        MfUsgDdf(ml, **kw)
    assert (len(ml.packagelist), ml.output_units, ml.output_fnames) == state


def test_ddf_output_name_kept_on_load_ws_change_and_replace(function_tmpdir):
    ml = _model(function_tmpdir)
    MfUsgLpf(ml, laytyp=4, hk=100.0)
    MfUsgSms(ml, linmeth=1)
    MfUsgOc(ml, stress_period_data={(0, 0): ["save head"]})
    MfUsgDdf(ml, isharp=1, hsea=0.0, iswiun=61, filenames=["a.ddf", "a.zeta"])
    ml.write_input()
    m2 = MfUsg.load(
        "m.nam",
        model_ws=str(function_tmpdir),
        check=False,
        load_only=["disu", "bas6", "lpf", "sms", "oc", "ddf"],
    )
    assert m2.ddf.file_name == ["a.ddf"] and m2.output_units.count(61) == 1
    assert m2.output_fnames[m2.output_units.index(61)] == "a.zeta"
    assert 61 not in m2.external_units
    with pytest.warns(UserWarning):
        MfUsgDdf(m2, isharp=1, hsea=0.1, iswiun=61, filenames=["a.ddf", "a.zeta"])
    assert m2.output_units.count(61) == 1 and m2.ddf.hsea == 0.1
    m2.change_model_ws(str(function_tmpdir / "re"))
    m2.write_input()
    nam = (function_tmpdir / "re" / "m.nam").read_text().splitlines()
    assert [ln.split()[1:] for ln in nam if ln.split()[1:2] == ["61"]] == [
        ["61", "a.zeta"]
    ]


def test_ddf_unsupported_after_write_keeps_destination(function_tmpdir):
    """CLN added after the DDF package, or BAS switched to fixed format."""
    ml = _model(function_tmpdir)
    d = MfUsgDdf(ml, isharp=1, hsea=0.0)
    d.write_file()
    before = Path(d.fn_path).read_bytes()
    _cln(ml)
    with pytest.raises(NotImplementedError, match="CLN"):
        d.write_file()
    assert Path(d.fn_path).read_bytes() == before
    ml2 = _model(function_tmpdir / "fixed")
    d2 = MfUsgDdf(ml2, ithickav=0)
    d2.write_file()
    before2 = Path(d2.fn_path).read_bytes()
    ml2.bas6.ifrefm = False
    with pytest.raises(NotImplementedError, match="FREE"):
        d2.write_file()
    assert Path(d2.fn_path).read_bytes() == before2


def _cln_model_files(ws, ddf_text):
    ml = _model(ws)
    _cln(ml)
    MfUsgDdf(ml, ithickav=0)
    ml.write_input()
    (Path(ws) / "m.ddf").write_text(ddf_text)


def test_ddf_sharp_cln_rejected_with_load_only(function_tmpdir):
    """CLN declared in the name file but not loaded (load_only)."""
    _cln_model_files(function_tmpdir, " 1000.0 1025.0 35.0 0 0 1\n 0.0 0\n 1.0\n")
    with pytest.raises(NotImplementedError, match="CLN"):
        MfUsg.load(
            "m.nam", model_ws=str(function_tmpdir), check=False, load_only=["ddf"]
        )


def test_ddf_nonsharp_cln_load_only_and_sharp_write(function_tmpdir):
    """Non-sharp DDF with CLN loads; turning it sharp is refused on write
    because the skipped CLN entry stays in the name file."""
    _cln_model_files(function_tmpdir, " 1000.0 1025.0 35.0 0 0 0\n")
    m2 = MfUsg.load(
        "m.nam", model_ws=str(function_tmpdir), check=False, load_only=["ddf"]
    )
    assert m2.ddf.isharp == 0 and m2.get_package("CLN") is None
    m2.change_model_ws(str(function_tmpdir / "re"))
    m2.ddf.isharp, m2.ddf.hsea = 1, 0.0
    with pytest.raises(NotImplementedError, match="CLN"):
        m2.ddf.write_file()
    assert not (function_tmpdir / "re" / "m.ddf").exists()


def test_ddf_load_continued_table_and_cursor(function_tmpdir):
    """List-directed rows continue across lines and blank lines; the rest of
    the last line is discarded; a caller handle is left after the record."""
    src = function_tmpdir / "c.ddf"
    src.write_text(
        "# heading\n"
        " 1000 1025 0.0357 0 0 1 nonlinear 3\n"
        " 0.0\n\n 1000.0 ignored text\n"
        " 0.01,1010.0, 99 more\n"
        " 0.0357 ,\n 1025.0\n"
        "# before HSEA\n"
        " 0.25 0 trailing\n"
        "NEXT RECORD\n"
    )
    ml = _model(function_tmpdir)
    with open(src) as fh:
        d = MfUsgDdf.load(fh, ml)
        assert fh.readline() == "NEXT RECORD\n"
    assert d.density_table == [(0.0, 1000.0), (0.01, 1010.0), (0.0357, 1025.0)]
    assert (d.isharp, d.hsea, d.iswiun) == (1, 0.25, 0)


def test_ddf_load_d_exponent_commas_tabs(function_tmpdir):
    src = function_tmpdir / "d.ddf"
    src.write_text(" 1.0D3, 1.025d3,\t35.0E0 , 01 ,0, 1\n 1.5D-1,61\n")
    d = MfUsgDdf.load(str(src), _model(function_tmpdir))
    assert (d.rhofresh, d.rhostd, d.cstd, d.ithickav, d.imphdd, d.isharp) == (
        1000.0,
        1025.0,
        35.0,
        1,
        0,
        1,
    )
    assert (d.hsea, d.iswiun) == (0.15, 61)


H0 = " 1000 1025 35 0 0 0"
LOAD_BAD = {
    "repeat_count": (H0 + " NONLINEAR 1\n 2*1000\n", "unsupported"),
    "slash": (H0 + " NONLINEAR 1\n 0.0 /\n", "unsupported"),
    "quoted": (H0 + " NONLINEAR 1\n '0.0' 1000\n", "unsupported"),
    "null_value": (H0 + " NONLINEAR 1\n 0.0,,1000\n", "null"),
    "leading_comma": (H0 + " NONLINEAR 1\n ,0.0 1000\n", "null"),
    "comment_in_table": (H0 + " NONLINEAR 1\n# x\n 0 1000\n", "density table"),
    "fraction_int": (" 1000 1025 35 0.0 0 0\n", "ITHICKAV"),
    "underscore": (" 1_000 1025 35 0 0 0\n", "RHOFRESH"),
    "nan_token": (" nan 1025 35 0 0 0\n", "RHOFRESH"),
    "nonlinear_in_imphdd": (" 1000 1025 35 0 NONLINEAR 1\n 0 1000\n", "IMPHDD"),
    "nonlinear_late": (H0 + " 7 NONLINEAR 1\n 0 1000\n", "NONLINEAR"),
    "nonlinear_no_rows": (H0 + " NONLINEAR\n", "NDDFTABROWS"),
    "nonlinear_zero_rows": (H0 + " NONLINEAR 0\n", "NDDFTABROWS"),
    "short_header": (" 1000 1025\n", "CSTD"),
    "blank_header": ("\n 1000 1025 35 0 0 0\n", "blank"),
    "blank_hsea": (" 1000 1025 35 0 0 1\n\n 0.2 0\n", "blank"),
    "indented_hash": ("  # not a comment\n" + H0 + "\n", "RHOFRESH"),
    "iswiun_fraction": (" 1000 1025 35 0 0 1\n 0.2 6.0\n", "ISWIUN"),
    "empty_file": ("", "end of file"),
    "comments_only": ("# a\n# b\n", "end of file"),
    "eof_in_row": (H0 + " NONLINEAR 2\n 0 1000\n 35\n", "end of file"),
    "cstd_zero": (" 1000 1025 0.0 0 0 0\n", "CSTD"),
}


@pytest.mark.parametrize("name", list(LOAD_BAD))
def test_ddf_load_syntax_rejected_and_closed(function_tmpdir, name, monkeypatch):
    import flopy.mfusg.mfusgddf as ddfmod

    text, match = LOAD_BAD[name]
    src = function_tmpdir / "bad.ddf"
    src.write_text(text)
    opened = []

    def tracking_open(*a, **k):
        fh = builtins.open(*a, **k)
        opened.append(fh)
        return fh

    monkeypatch.setattr(ddfmod, "open", tracking_open, raising=False)
    ml = _model(function_tmpdir)
    with pytest.raises(ValueError, match=match):
        MfUsgDdf.load(str(src), ml)
    assert opened and all(fh.closed for fh in opened)
    assert ml.get_package("DDF") is None


def test_ddf_load_errors_leave_caller_handle_open(function_tmpdir, monkeypatch):
    import flopy.mfusg.mfusgddf as ddfmod

    src = function_tmpdir / "e.ddf"
    src.write_text(" 1000.0 1025.0 35.0 0 0 1\n")
    with open(src) as fh:
        with pytest.raises(ValueError, match="HSEA"):
            MfUsgDdf.load(fh, _model(function_tmpdir))
        assert not fh.closed
    sharp = function_tmpdir / "s.ddf"
    sharp.write_text(" 1000.0 1025.0 35.0 0 0 1\n 0.0 0\n 1.0\n")
    opened = []

    def tracking_open(*a, **k):
        fh = builtins.open(*a, **k)
        opened.append(fh)
        return fh

    monkeypatch.setattr(ddfmod, "open", tracking_open, raising=False)
    ml = _model(function_tmpdir / "cln")
    _cln(ml)
    with pytest.raises(NotImplementedError, match="CLN"):
        MfUsgDdf.load(str(sharp), ml)
    assert opened and all(fh.closed for fh in opened)


# ---------------------------------------------------------------- output registry


def _registry(ml):
    """Everything a rejected DDF must leave unchanged in the model."""
    return (
        [id(p) for p in ml.packagelist],
        list(ml.package_units),
        list(ml.output_units),
        list(ml.output_fnames),
        list(ml.output_binflag),
        [list(p) for p in ml.output_packages],
        list(ml.external_units),
        list(ml.external_fnames),
        list(ml.external_binflag),
        list(ml.external_output),
        list(ml.pop_key_list),
    )


SHARP = {"isharp": 1, "hsea": 0.0}


def _name_collision_setup(ml, case):
    if case == "basename":
        ml.add_output_file(61, fname="shared.zta", binflag=True, package="OC")
        return {"iswiun": 60, "filenames": ["m.ddf", "sub/shared.zta"]}
    if case == "default_name":
        ml.add_output_file(61, fname="m.zta", binflag=True, package="OC")
        ml.add_output_file(62, fname="m.DDF.zta", binflag=True, package="XYZ")
        return {"iswiun": 60}
    if case == "package_file":
        return {"iswiun": 60, "filenames": ["m.ddf", "dir/m.bas"]}
    ml.add_external("ext.bin", 70, binflag=True)
    return {"iswiun": 60, "filenames": ["m.ddf", "ext.bin"]}


@pytest.mark.parametrize(
    "case", ["basename", "default_name", "package_file", "external_name"]
)
def test_ddf_output_effective_name_collision_rejected(function_tmpdir, case):
    """The name add_output_file would register (base name, or the default
    name with its fallback) is checked before the model changes."""
    ml = _model(function_tmpdir)
    kw = _name_collision_setup(ml, case)
    before = _registry(ml)
    with pytest.raises(ValueError, match="already"):
        MfUsgDdf(ml, **SHARP, **kw)
    assert _registry(ml) == before


def test_ddf_output_effective_names_accepted(function_tmpdir):
    """Control: a path is registered by its base name; the default name falls
    back to m.DDF.zta when m.zta belongs to another unit; the other output
    is kept."""
    ml = _model(function_tmpdir)
    d = MfUsgDdf(ml, **SHARP, iswiun=60, filenames=["m.ddf", "sub/a.zta"])
    assert ml.output_fnames[ml.output_units.index(60)] == "a.zta"
    d.write_file()
    ml2 = _model(function_tmpdir / "fallback")
    ml2.add_output_file(61, fname="m.zta", binflag=True, package="OC")
    MfUsgDdf(ml2, **SHARP, iswiun=60)
    assert list(zip(ml2.output_units, ml2.output_fnames)) == [
        (61, "m.zta"),
        (60, "m.DDF.zta"),
    ]
    assert ml2.output_packages == [["OC"], ["DDF"]]


def _change_output(ml, case):
    idx = ml.output_units.index(60)
    if case == "owner_text":
        ml.remove_output(unit=60)
        ml.add_output_file(60, fname="other.bin", binflag=False, package="OC")
    elif case == "owner_binary":
        ml.remove_output(unit=60)
        ml.add_output_file(60, fname="other.bin", binflag=True, package="OC")
    elif case == "shared_owner":
        ml.add_output_file(60, binflag=True, package="OC")
    elif case == "renamed":
        ml.output_fnames[idx] = "renamed.zta"
    elif case == "text_flag":
        ml.output_binflag[idx] = False


@pytest.mark.parametrize(
    "case", ["owner_text", "owner_binary", "shared_owner", "renamed", "text_flag"]
)
def test_ddf_output_changed_before_write_rejected(function_tmpdir, case):
    ml = _model(function_tmpdir)
    d = MfUsgDdf(ml, **SHARP, iswiun=60)
    d.write_file()
    before = Path(d.fn_path).read_bytes()
    _change_output(ml, case)
    d.hsea = 0.2
    state = _registry(ml)
    with pytest.raises(ValueError, match=r"ISWIUN 60|output unit 60"):
        d.write_file()
    assert Path(d.fn_path).read_bytes() == before
    assert _registry(ml) == state


def test_ddf_write_after_change_model_ws_control(function_tmpdir):
    """Control: changing the workspace keeps the output and writing works."""
    ml = _model(function_tmpdir)
    d = MfUsgDdf(ml, **SHARP, iswiun=60, filenames=["m.ddf", "z/a.zta"])
    ml.change_model_ws(str(function_tmpdir / "moved"))
    d.write_file()
    assert (function_tmpdir / "moved" / "m.ddf").exists()
    assert ml.output_fnames[ml.output_units.index(60)] == "a.zta"


@pytest.mark.parametrize("binflag", [False, True], ids=["text", "binary"])
def test_ddf_external_unit_rejected_when_authoring(function_tmpdir, binflag):
    ml = _model(function_tmpdir)
    ml.add_external("input.dat", 60, binflag=binflag)
    before = _registry(ml)
    with pytest.raises(ValueError, match="external file"):
        MfUsgDdf(ml, **SHARP, iswiun=60)
    assert _registry(ml) == before


def test_ddf_external_unit_rejected_on_direct_load(function_tmpdir):
    src = function_tmpdir / "s.ddf"
    src.write_text(" 1000.0 1025.0 35.0 0 0 1\n 0.0 60\n")
    ml = _model(function_tmpdir)
    ml.add_external("input.dat", 60, binflag=True)
    before = _registry(ml)
    with pytest.raises(ValueError, match="external file"):
        MfUsgDdf.load(str(src), ml)
    assert _registry(ml) == before


def _sharp_nam_model(ws):
    ml = _model(ws)
    MfUsgLpf(ml, laytyp=4, hk=100.0)
    MfUsgSms(ml, linmeth=1)
    MfUsgOc(ml, stress_period_data={(0, 0): ["save head"]})
    MfUsgDdf(ml, **SHARP, iswiun=61, filenames=["m.ddf", "a.zeta"])
    ml.write_input()
    return Path(ws) / "m.nam"


LOAD_ONLY = ["disu", "bas6", "lpf", "sms", "oc", "ddf"]


def test_ddf_nam_data_entry_before_ddf_adopted(function_tmpdir):
    """The DATA(BINARY) entry of ISWIUN listed before DDF (it is first read as
    an external file by MfUsg.load) becomes the DDF output once."""
    nam = _sharp_nam_model(function_tmpdir)
    lines = nam.read_text().splitlines()
    data = [ln for ln in lines if ln.split()[1:2] == ["61"]]
    rest = [ln for ln in lines if ln not in data]
    nam.write_text("\n".join(rest[:3] + data + rest[3:]) + "\n")
    m2 = MfUsg.load(
        "m.nam", model_ws=str(function_tmpdir), check=False, load_only=LOAD_ONLY
    )
    assert m2.output_units.count(61) == 1 and 61 not in m2.external_units
    idx = m2.output_units.index(61)
    assert (m2.output_fnames[idx], m2.output_binflag[idx]) == ("a.zeta", True)
    assert m2.output_packages[idx] == ["DDF"]


def test_ddf_nam_iswiun_entry_not_binary_rejected(function_tmpdir):
    entry = "DATA"
    nam = _sharp_nam_model(function_tmpdir)
    text = nam.read_text().replace("DATA(BINARY)      61", f"{entry:18s}61")
    nam.write_text(text)
    with pytest.raises(ValueError, match=r"DATA\(BINARY\)"):
        MfUsg.load(
            "m.nam", model_ws=str(function_tmpdir), check=False, load_only=["ddf"]
        )


@pytest.mark.parametrize("new", [{"iswiun": 61}, {"isharp": 0, "hsea": None}])
def test_ddf_replace_with_other_iswiun_rejected(function_tmpdir, new):
    ml = _model(function_tmpdir)
    old = MfUsgDdf(ml, **SHARP, iswiun=60)
    before = _registry(ml)
    with pytest.raises(ValueError, match="still registered for DDF"):
        MfUsgDdf(ml, **{**SHARP, **new})
    assert _registry(ml) == before and ml.get_package("DDF") is old


def test_ddf_replace_after_remove_output(function_tmpdir):
    """Control: once the old output is removed, a new ISWIUN is registered
    and no stale DDF output remains."""
    ml = _model(function_tmpdir)
    MfUsgDdf(ml, **SHARP, iswiun=60)
    ml.remove_output(unit=60)
    with pytest.warns(UserWarning):
        d = MfUsgDdf(ml, **SHARP, iswiun=61)
    assert [u for u, p in zip(ml.output_units, ml.output_packages) if "DDF" in p] == [
        61
    ]
    d.write_file()


def test_ddf_load_fixed_bas_not_supported(function_tmpdir, monkeypatch):
    import flopy.mfusg.mfusgddf as ddfmod

    src = function_tmpdir / "f.ddf"
    src.write_text(" 1000.0 1025.0 35.0 0 0 1\n 0.0 60\n")
    opened = []

    def tracking_open(*a, **k):
        fh = builtins.open(*a, **k)
        opened.append(fh)
        return fh

    monkeypatch.setattr(ddfmod, "open", tracking_open, raising=False)
    ml = _model(function_tmpdir, free=False)
    before = _registry(ml)
    with pytest.raises(NotImplementedError, match="FREE"):
        MfUsgDdf.load(str(src), ml)
    assert opened == [] and _registry(ml) == before
    with open(src) as fh:
        with pytest.raises(NotImplementedError, match="FREE"):
            MfUsgDdf.load(fh, ml)
        assert not fh.closed and fh.tell() == 0
    assert _registry(ml) == before


def test_ddf_load_without_bas_not_inferred_control(function_tmpdir):
    """Control: without a BAS package the format is not inferred."""
    src = function_tmpdir / "f.ddf"
    src.write_text(" 1000.0 1025.0 35.0 0 0 0\n")
    d = MfUsgDdf.load(str(src), MfUsg(modelname="m", model_ws=str(function_tmpdir)))
    assert d.cstd == 35.0


REAL32_BAD = {
    "sharp_equal_in_real": (
        {**SHARP, "rhofresh": 1000.0, "rhostd": 1000.0 + 1e-7},
        "RHOSTD",
    ),
    "rhofresh_underflow": ({"rhofresh": 1e-50}, "RHOFRESH"),
    "cstd_underflow": ({"cstd": 1e-50}, "CSTD"),
    "rhostd_overflow": ({"rhostd": 1e39}, "RHOSTD"),
    "cstd_overflow": ({"cstd": -1e39}, "CSTD"),
    "hsea_overflow": ({"isharp": 1, "hsea": 1e39}, "HSEA"),
}


@pytest.mark.parametrize("name", list(REAL32_BAD))
def test_ddf_real32_values_validated(function_tmpdir, name):
    """Denominators and ranges on the REAL values USG-T reads (protection;
    these cases were not run with USG-T)."""
    kw, match = REAL32_BAD[name]
    ml = _model(function_tmpdir)
    before = _registry(ml)
    with pytest.raises(ValueError, match=match):
        MfUsgDdf(ml, **{"ithickav": 0, **kw})
    assert _registry(ml) == before


def test_ddf_real32_values_written_unquantized(function_tmpdir):
    """Distinct REAL values are accepted and written as given."""
    ml = _model(function_tmpdir)
    d = MfUsgDdf(ml, **SHARP, rhofresh=1000.0, rhostd=1000.001, cstd=1e-30)
    d.write_file()
    assert [float(v) for v in _data_lines(d.fn_path)[0][:3]] == [
        1000.0,
        1000.001,
        1e-30,
    ]


def test_ddf_real32_collapse_rejected_on_write_and_load(function_tmpdir):
    ml = _model(function_tmpdir)
    d = MfUsgDdf(ml, **SHARP)
    d.write_file()
    before = Path(d.fn_path).read_bytes()
    d.rhostd = 1000.0 + 1e-7
    with pytest.raises(ValueError, match="RHOSTD"):
        d.write_file()
    assert Path(d.fn_path).read_bytes() == before
    src = function_tmpdir / "c.ddf"
    src.write_text(" 1000.0 1000.0000001 35.0 0 0 1\n 0.0 0\n")
    ml2 = _model(function_tmpdir / "ld")
    before2 = _registry(ml2)
    with pytest.raises(ValueError, match="RHOSTD"):
        MfUsgDdf.load(str(src), ml2)
    assert _registry(ml2) == before2


# ---------------------------------------------------------------- load context


def _ext_unit_dict(ddf="m.ddf", out="shared.zta"):
    """Minimal name-file entries (DDF on 152, ISWIUN 60 as DATA(BINARY))."""
    from types import SimpleNamespace

    return {
        152: SimpleNamespace(filetype="DDF", filename=ddf),
        60: SimpleNamespace(filetype="DATA(BINARY)", filename=out),
    }


SHARP_TEXT = " 1000.0 1025.0 35.0 0 0 1\n 0.0 60\n"


def _hook_load_construction(monkeypatch, model, on_second_call):
    """Call on_second_call() at the second _output_name check for `model`:
    the first is load()'s check before the name-file units are claimed, the
    second runs while load() builds the package."""
    orig = MfUsgDdf._output_name
    calls = []

    def hooked(m, *args, **kwargs):
        if m is model:
            calls.append(1)
            if len(calls) == 2:
                on_second_call()
        return orig(m, *args, **kwargs)

    monkeypatch.setattr(MfUsgDdf, "_output_name", staticmethod(hooked))


def test_ddf_load_context_not_shared_between_models(function_tmpdir, monkeypatch):
    """While load() of model A is building its package (paused with events),
    authoring a DDF in model B must not adopt B's external file."""
    import io
    import threading

    a = _model(function_tmpdir / "a")
    b = _model(function_tmpdir / "b")
    a.add_external("shared.zta", 60, binflag=True)
    b.add_external("shared.zta", 60, binflag=True)
    ready, release = threading.Event(), threading.Event()

    def pause():
        ready.set()
        assert release.wait(10)

    _hook_load_construction(monkeypatch, a, pause)
    result = {}

    def load_a():
        try:
            result["a"] = MfUsgDdf.load(io.StringIO(SHARP_TEXT), a, _ext_unit_dict())
        except Exception as e:  # noqa: BLE001
            result["a"] = e

    thread = threading.Thread(target=load_a)
    thread.start()
    try:
        assert ready.wait(10)
        before = _registry(b)
        with pytest.raises(ValueError, match="external file"):
            MfUsgDdf(b, **SHARP, iswiun=60, filenames=["b.ddf", "shared.zta"])
        assert _registry(b) == before
    finally:
        release.set()
        thread.join(10)
    assert not thread.is_alive()
    assert isinstance(result["a"], MfUsgDdf), result["a"]
    idx = a.output_units.index(60)
    assert (a.output_fnames[idx], a.output_binflag[idx]) == ("shared.zta", True)
    assert a.output_packages[idx] == ["DDF"] and 60 not in a.external_units
    assert "_nam_output" not in vars(MfUsgDdf)


def test_ddf_failed_load_leaves_no_context(function_tmpdir, monkeypatch):
    import io

    a = _model(function_tmpdir / "a")
    b = _model(function_tmpdir / "b")
    a.add_external("shared.zta", 60, binflag=True)
    b.add_external("shared.zta", 60, binflag=True)

    def fail():
        raise RuntimeError("construction failed")

    _hook_load_construction(monkeypatch, a, fail)
    with pytest.raises(RuntimeError, match="construction failed"):
        MfUsgDdf.load(io.StringIO(SHARP_TEXT), a, _ext_unit_dict())
    assert "_nam_output" not in vars(MfUsgDdf)
    assert 60 in a.external_units and 60 not in a.output_units
    before = _registry(b)
    with pytest.raises(ValueError, match="external file"):
        MfUsgDdf(b, **SHARP, iswiun=60, filenames=["b.ddf", "shared.zta"])
    assert _registry(b) == before


def test_ddf_load_context_bound_to_its_model(function_tmpdir):
    """A load context built for model A does not let model B adopt a file."""
    from flopy.mfusg.mfusgddf import _NamOutput

    a = _model(function_tmpdir / "a")
    b = _model(function_tmpdir / "b")
    b.add_external("shared.zta", 60, binflag=True)
    before = _registry(b)
    d = MfUsgDdf.__new__(MfUsgDdf)
    with pytest.raises(ValueError, match="external file"):
        d._setup(
            b,
            1000.0,
            1025.0,
            35.0,
            0,
            0,
            1,
            False,
            None,
            "ddf",
            None,
            ["b.ddf", "shared.zta"],
            0.0,
            60,
            nam_output=_NamOutput(a, 60, "shared.zta"),
        )
    assert _registry(b) == before


OWN_FILE = {
    "explicit": {"filenames": ["same.ddf", "same.ddf"]},
    "subfolder": {"filenames": ["m.ddf", "sub/m.ddf"]},
    "default_output": {"filenames": ["m.zta"]},
    "default_both": {"extension": "zta"},
}


@pytest.mark.parametrize("case", list(OWN_FILE))
def test_ddf_output_is_own_ddf_file_rejected(function_tmpdir, case):
    """The interface output may not be the DDF file being created; files
    already on disk and the model registry stay as they were."""
    ws = Path(function_tmpdir)
    for name in ("same.ddf", "m.ddf", "m.zta"):
        (ws / name).write_text(f"pre-existing {name}\n")
    files = {f.name: f.read_bytes() for f in ws.iterdir() if f.is_file()}
    ml = _model(ws)
    before = _registry(ml)
    with pytest.raises(ValueError, match="DDF file itself"):
        MfUsgDdf(ml, **SHARP, iswiun=60, **OWN_FILE[case])
    assert _registry(ml) == before
    assert {f.name: f.read_bytes() for f in ws.iterdir() if f.is_file()} == files


def test_ddf_load_output_is_own_ddf_file_rejected(function_tmpdir):
    import io

    ml = _model(function_tmpdir)
    before = _registry(ml)
    with pytest.raises(ValueError, match="DDF file itself"):
        MfUsgDdf.load(
            io.StringIO(SHARP_TEXT), ml, _ext_unit_dict("same.ddf", "sub/same.ddf")
        )
    assert _registry(ml) == before


# ---------------------------------------------------------------- engine


def _run(ws, ddf_kw=None, hand=None, iswiun=0):
    ml = _model(ws, exe=USGT_EXE)
    MfUsgLpf(ml, laytyp=4, hk=100.0, ss=1.0e-5, sy=0.2)
    MfUsgSms(ml, linmeth=1, mxiter=1000, iter1=200, nonlinmeth=1)
    ModflowChd(ml, stress_period_data={0: [[0, 0.0, 0.0]]})
    ModflowWel(ml, stress_period_data={0: [[4, 20.0]]})
    MfUsgOc(ml, stress_period_data={(0, 0): ["save head"]})
    if hand is None:
        MfUsgDdf(ml, **ddf_kw)
    else:
        MfUsgDdf(ml, filenames="manual.ddf")
        if iswiun:
            ml.add_output_file(iswiun, fname="m.zta", binflag=True, package="DDF")
    ml.write_input()
    if hand is not None:
        (Path(ws) / "manual.ddf").write_text(hand)
    ok, buff = ml.run_model(silent=True)
    assert ok, buff[-5:]
    head = _read_final(Path(ws) / "m.hds")
    zf = Path(ws) / "m.zta"
    zeta = _read_final(zf, text="interface elev.") if zf.exists() else None
    return head, zeta


def _read_final(path, text="headu"):
    """Single steady time (1.0) from a closed reader; values must be finite."""
    from flopy.utils import HeadUFile

    with HeadUFile(str(path), text=text) as hf:
        assert hf.get_times() == [1.0]
        data = np.asarray(hf.get_data()[0], float)
    assert np.isfinite(data).all()
    return data


@requires_usgt
def test_ddf_engine_sharp_flopy_equals_hand(function_tmpdir):
    h1, z1 = _run(
        function_tmpdir / "fp", {"ithickav": 0, "isharp": 1, "hsea": 0.0, "iswiun": 60}
    )
    h2, z2 = _run(
        function_tmpdir / "hand", hand=" 1000.0 1025.0 35.0 0 0 1\n 0.0 60\n", iswiun=60
    )
    np.testing.assert_array_equal(h1, h2)
    np.testing.assert_array_equal(np.asarray(z1), np.asarray(z2))


@requires_usgt
def test_ddf_engine_sharp_hsea_sensitivity(function_tmpdir):
    h0, z0 = _run(
        function_tmpdir / "h0", {"ithickav": 0, "isharp": 1, "hsea": 0.0, "iswiun": 60}
    )
    h1, z1 = _run(
        function_tmpdir / "h1", {"ithickav": 0, "isharp": 1, "hsea": 0.1, "iswiun": 60}
    )
    assert not np.allclose(h0, h1)
    for h, z, hsea in ((h0, z0, 0.0), (h1, z1, 0.1)):
        expect = np.clip((1.025 * hsea - h) / 0.025, -40.0, 10.0)
        np.testing.assert_allclose(np.asarray(z, float), expect, atol=1e-3)


@requires_usgt
def test_ddf_engine_sharp_differs_from_no_interface(function_tmpdir):
    hs, _ = _run(function_tmpdir / "s", {"ithickav": 0, "isharp": 1, "hsea": 0.0})
    ml = _model(function_tmpdir / "n", exe=USGT_EXE)
    MfUsgLpf(ml, laytyp=4, hk=100.0, ss=1.0e-5, sy=0.2)
    MfUsgSms(ml, linmeth=1, mxiter=1000, iter1=200, nonlinmeth=1)
    ModflowChd(ml, stress_period_data={0: [[0, 0.0, 0.0]]})
    ModflowWel(ml, stress_period_data={0: [[4, 20.0]]})
    MfUsgOc(ml, stress_period_data={(0, 0): ["save head"]})
    ml.write_input()
    ok, _ = ml.run_model(silent=True)
    assert ok
    hn = _read_final(function_tmpdir / "n" / "m.hds")
    assert hs[1] > hn[1] + 0.1


@requires_usgt
@pytest.mark.parametrize(
    "kw,hand",
    [
        (
            {"rhofresh": 1000.0, "rhostd": 1025.0, "cstd": 0.0357, "ithickav": 0},
            " 1000.0 1025.0 0.0357 0 0 0\n",
        ),
        (
            {
                "cstd": 0.0357,
                "ithickav": 0,
                "nonlinear": True,
                "density_table": [
                    (0.0, 1000.0),
                    (0.0123457, 1008.645678),
                    (0.0357, 1025.0),
                ],
            },
            " 1000.0 1025.0 0.0357 0 0 0 NONLINEAR 3\n 0.0 1000.0\n"
            " 0.0123457 1008.645678\n 0.0357 1025.0\n",
        ),
    ],
    ids=["header", "table_control"],
)
def test_ddf_engine_precision_flopy_equals_hand(function_tmpdir, kw, hand):
    """header: CSTD 0.0357 must not become 0.04. table_control: 6-digit rounding of
    the table is below the head output precision of this model (passes on develop)."""

    def run(ws, ddf_kw=None, text=None):
        ml = _model(ws, bct=True, exe=USGT_EXE)
        MfUsgLpf(ml, laytyp=4, hk=100.0, ss=1.0e-5, sy=0.2)
        MfUsgSms(ml, linmeth=1, mxiter=1000, iter1=200, nonlinmeth=1)
        ModflowChd(ml, stress_period_data={0: [[0, 0.0, 0.0]]})
        ModflowWel(ml, stress_period_data={0: [[4, 20.0]]})
        MfUsgOc(ml, stress_period_data={(0, 0): ["save head"]})
        MfUsgDdf(ml, **(ddf_kw or {"filenames": "manual.ddf"}))
        ml.write_input()
        if text:
            (Path(ws) / "manual.ddf").write_text(text)
        ok, buff = ml.run_model(silent=True)
        assert ok, buff[-5:]
        return _read_final(Path(ws) / "m.hds")

    np.testing.assert_array_equal(
        run(function_tmpdir / "fp", kw), run(function_tmpdir / "hand", text=hand)
    )


@requires_usgt
@pytest.mark.parametrize(
    "kw,hand",
    [
        (
            {"ithickav": 0, "isharp": 2, "hsea": 0.0, "iswiun": 60},
            "1000. 1025. 35. 0 0 2\n.0 60\n",
        ),
        (
            {
                "ithickav": 0,
                "isharp": 1,
                "hsea": 0.1,
                "iswiun": 60,
                "nonlinear": True,
                "density_table": [(0.0, 1000.0), (35.0, 1025.0)],
            },
            "# hand\n1.0E3,1025,3.5D1, 0,0,1 nonlinear 2\n0\n1000 rest\n35.,1025.\n"
            "# sea level\n1.0D-1 60\n",
        ),
    ],
    ids=["isharp2", "sharp_after_table"],
)
def test_ddf_engine_sharp_variants_flopy_equals_hand(function_tmpdir, kw, hand):
    """Hand files use other number formats, commas and a continued row; FloPy
    reads them as the engine does."""
    h1, z1 = _run(function_tmpdir / "fp", kw)
    h2, z2 = _run(function_tmpdir / "hand", hand=hand, iswiun=60)
    np.testing.assert_array_equal(h1, h2)
    np.testing.assert_array_equal(z1, z2)
    d = MfUsgDdf.load(
        str(function_tmpdir / "hand" / "manual.ddf"), _model(function_tmpdir / "ld")
    )
    for key, value in kw.items():
        assert getattr(d, key) == value, key
