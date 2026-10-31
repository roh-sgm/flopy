import os
import warnings
from io import StringIO
from struct import pack
from tempfile import TemporaryFile
from textwrap import dedent

import numpy as np
import pytest

from flopy.utils.util_array import Util2d, Util3d


def test_load_txt_free():
    a = np.ones((10,), dtype=np.float32) * 250.0
    fp = StringIO("10*250.0")
    fa = Util2d.load_txt(a.shape, fp, a.dtype, "(FREE)")
    np.testing.assert_equal(fa, a)
    assert fa.dtype == a.dtype

    a = np.arange(10, dtype=np.int32).reshape((2, 5))
    fp = StringIO(
        dedent(
            """\
        0 1,2,3, 4
        5 6, 7,  8 9
    """
        )
    )
    fa = Util2d.load_txt(a.shape, fp, a.dtype, "(FREE)")
    np.testing.assert_equal(fa, a)
    assert fa.dtype == a.dtype

    a = np.ones((2, 5), dtype=np.float32)
    a[1, 0] = 2.2
    fp = StringIO(
        dedent(
            """\
        5*1.0
        2.2 2*1.0, +1E-00 1.0
    """
        )
    )
    fa = Util2d.load_txt(a.shape, fp, a.dtype, "(FREE)")
    np.testing.assert_equal(fa, a)
    assert fa.dtype == a.dtype


def test_load_txt_fixed():
    a = np.arange(10, dtype=np.int32).reshape((2, 5))
    fp = StringIO(
        dedent(
            """\
        01234X
        56789
    """
        )
    )
    fa = Util2d.load_txt(a.shape, fp, a.dtype, "(5I1)")
    np.testing.assert_equal(fa, a)
    assert fa.dtype == a.dtype

    fp = StringIO(
        dedent(
            """\
        0123X
        4
        5678
        9
    """
        )
    )
    fa = Util2d.load_txt(a.shape, fp, a.dtype, "(4I1)")
    np.testing.assert_equal(fa, a)
    assert fa.dtype == a.dtype

    a = np.array([[-1, 1, -2, 2, -3], [3, -4, 4, -5, 5]], np.int32)
    fp = StringIO(
        dedent(
            """\
        -1 1-2 2-3
        3 -44 -55
    """
        )
    )
    fa = Util2d.load_txt(a.shape, fp, a.dtype, "(5I2)")
    np.testing.assert_equal(fa, a)
    assert fa.dtype == a.dtype


def test_load_block():
    a = np.ones((2, 5), dtype=np.int32) * 4
    fp = StringIO(
        dedent(
            """\
        1
        1 2 1 5 4
    """
        )
    )
    fa = Util2d.load_block(a.shape, fp, a.dtype)
    np.testing.assert_equal(fa, a)
    assert fa.dtype == a.dtype

    a = np.ones((2, 5), dtype=np.float32) * 4
    a[0:2, 1:2] = 9.0
    a[0, 2:4] = 6.0
    fp = StringIO(
        dedent(
            """\
        3
        1 2 1 5 4.0
        1 2 2 2 9.0
        1 1 3 4 6.0
    """
        )
    )
    fa = Util2d.load_block(a.shape, fp, a.dtype)
    np.testing.assert_equal(fa, a)
    assert fa.dtype == a.dtype

    a = np.zeros((2, 5), dtype=np.int32)
    a[0, 2:4] = 8
    fp = StringIO(
        dedent(
            """\
        1
        1 1 3 4 8
    """
        )
    )
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        fa = Util2d.load_block(a.shape, fp, a.dtype)
        assert len(w) == 1
        assert "blocks do not cover full array" in str(w[-1].message)
        np.testing.assert_equal(fa, a)
        assert fa.dtype == a.dtype


def test_load_bin(function_tmpdir):
    def temp_file(data):
        # writable file that is destroyed as soon as it is closed
        f = TemporaryFile(dir=function_tmpdir)
        f.write(data)
        f.seek(0)
        return f

    # INTEGER
    a = np.arange(3 * 4, dtype=np.int32).reshape((3, 4)) - 1
    fp = temp_file(a.tobytes())
    fh, fa = Util2d.load_bin((3, 4), fp, np.int32)
    assert fh is None  # no header_dtype
    np.testing.assert_equal(fa, a)
    assert fa.dtype == a.dtype

    # check warning if wrong integer type is used to read 4-byte integers
    # e.g. on platforms where int -> int64
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        fp.seek(0)
        fh, fa = Util2d.load_bin((3, 4), fp, np.int64)
        fp.close()
        assert len(w) == 1
        assert a.dtype == np.int32
        assert fh is None  # no header_dtype
        np.testing.assert_equal(fa, a)

    # REAL
    real_header_fmt = "2i2f16s3i"
    header_data = (1, 2, 3.5, 4.5, b"Hello", 6, 7, 8)
    real_header = pack(real_header_fmt, *header_data)
    assert len(real_header) == 44

    a = np.arange(10).reshape((2, 5))
    fp = temp_file(real_header + pack("10f", *list(range(10))))
    fh, fa = Util2d.load_bin((2, 5), fp, np.float32, "Head")
    fp.close()
    for h1, h2 in zip(fh[0], header_data):
        assert h1 == h2
    np.testing.assert_equal(a.astype(np.float32), fa)
    assert fa.dtype == np.float32

    # DOUBLE PRECISION
    dbl_header_fmt = "2i2d16s3i"
    dbl_header = pack(dbl_header_fmt, *header_data)
    assert len(dbl_header) == 52

    fp = temp_file(real_header + pack("10d", *list(range(10))))
    fh, fa = Util2d.load_bin((2, 5), fp, np.float64, "Head")
    fp.close()
    for h1, h2 in zip(fh[0], header_data):
        assert h1 == h2
    np.testing.assert_equal(a.astype(np.float64), fa)
    assert fa.dtype == np.float64


# double-precision arrays (U1DREL8/U2DREL8): values are REAL*8, while the
# control-record multiplier (CNSTNT) and CONSTANT values are REAL

DOUBLES = np.array([[23.1234567890123, 23.1234567890124, 22.75]])  # distinct
NOT_REAL = 23.1234567890123  # not exactly representable as float32


def _model(ws, external_path=None):
    from flopy.modflow import Modflow

    return Modflow(model_ws=str(ws), external_path=external_path)


def _open_fds():
    return len(os.listdir("/dev/fd"))


def _bin_record(values, double=True, ilay=1):
    real = "d" if double else "f"
    header = pack(
        f"<ii{real}{real}16siii",
        1,
        1,
        0.0,
        0.0,
        b"  IMMOBILE HEAD ",
        1,
        len(values),
        ilay,
    )
    return header + pack(f"<{len(values)}{real}", *values)


def test_util2d_float64_keeps_values_and_writes_17_digits(function_tmpdir):
    u2d = Util2d(_model(function_tmpdir), (1, 3), np.float64, DOUBLES, "h")
    assert u2d.array.dtype == np.float64
    np.testing.assert_array_equal(u2d.array, DOUBLES)
    assert u2d.format.fortran == "(3E24.16)"
    entry = u2d.get_file_entry()
    assert entry.startswith("INTERNAL")
    values = [float(v) for v in entry.splitlines()[1].split()]
    assert values == DOUBLES[0].tolist()  # 17 significant digits round-trip
    # float32 arrays keep their defaults
    u32 = Util2d(_model(function_tmpdir), (1, 3), np.float32, DOUBLES, "h")
    assert u32.format.fortran == "(3E15.6)"
    assert u32.array[0, 0] == u32.array[0, 1]


def test_util2d_float64_custom_format_is_honored(function_tmpdir):
    u2d = Util2d(
        _model(function_tmpdir), (1, 3), np.float64, DOUBLES, "h", fmtin="(3E15.6)"
    )
    values = [float(v) for v in u2d.get_file_entry().splitlines()[1].split()]
    assert values[0] == values[1]  # fewer digits round the values written


def test_util2d_float64_constant_only_when_exact_real(function_tmpdir):
    m = _model(function_tmpdir)
    exact = Util2d(m, (1, 3), np.float64, 22.75, "h")
    assert exact.how == "constant"
    assert exact.get_file_entry().split()[:1] == ["CONSTANT"]
    not_exact = Util2d(m, (1, 3), np.float64, NOT_REAL, "h")
    assert not_exact.how == "internal"
    values = [float(v) for v in not_exact.get_file_entry().splitlines()[1].split()]
    assert values == [NOT_REAL] * 3
    not_exact.how = "constant"
    with pytest.raises(ValueError, match="REAL"):
        not_exact.check_real_fields()
    with pytest.raises(ValueError, match="REAL"):
        not_exact.get_file_entry()
    # float32 behaviour is unchanged
    assert Util2d(m, (1, 3), np.float32, NOT_REAL, "h").how == "constant"


@pytest.mark.parametrize("cnstnt, factor", [(0.0, 1.0), (1.0, 1.0), (0.1, None)])
def test_util2d_float64_multiplier_is_real(function_tmpdir, cnstnt, factor):
    """A(J)=A(J)*CNSTNT with CNSTNT REAL; CNSTNT 0 means no multiplication."""
    factor = np.float64(np.float32(cnstnt)) if factor is None else factor
    m = _model(function_tmpdir)
    built = Util2d(m, (1, 3), np.float64, DOUBLES, "h", cnstnt=cnstnt)
    assigned = Util2d(m, (1, 3), np.float64, DOUBLES, "h")
    assigned.cnstnt = cnstnt
    u3d = Util3d(m, (2, 1, 3), np.float64, np.stack([DOUBLES] * 2), "h")
    u3d.cnstnt = cnstnt
    for u2d in [built, assigned] + u3d.util_2ds:
        assert isinstance(u2d.cnstnt, np.float32)
        np.testing.assert_array_equal(u2d.array, DOUBLES * factor)
        cr = u2d.get_file_entry().splitlines()[0].split()
        assert np.float32(float(cr[1])) == np.float32(cnstnt)
    np.testing.assert_array_equal(u3d.array, np.stack([DOUBLES] * 2) * factor)


def test_util2d_float64_fixed_record_fields(function_tmpdir):
    m = _model(function_tmpdir)
    m.array_free_format = False
    u2d = Util2d(m, (1, 3), np.float64, DOUBLES, "h", cnstnt=0.1, locat=55)
    line = u2d.get_file_entry().splitlines()[0]
    # I10, F10.0, A20, I10 (U1DREL8/U2DREL8 fixed control record)
    assert int(line[0:10]) == 55
    assert np.float32(float(line[10:20])) == np.float32(0.1)
    assert line[20:40].strip() == "(3E24.16)"
    values = [float(v) for v in u2d.get_file_entry().splitlines()[1].split()]
    assert values == DOUBLES[0].tolist()
    const = Util2d(m, (1, 3), np.float64, 22.75, "h")
    line = const.get_file_entry()
    assert int(line[0:10]) == 0 and float(line[10:20]) == 22.75
    too_long = Util2d(
        m, (1, 3), np.float64, DOUBLES, "h", cnstnt=1.2345678e-05, locat=55
    )
    with pytest.raises(ValueError, match="does not fit"):
        too_long.check_real_fields()


def test_util2d_float64_binary_output_rejected(function_tmpdir):
    m = _model(function_tmpdir, external_path="ext")
    u2d = Util2d(m, (1, 3), np.float64, DOUBLES, "h")
    u2d.format.binary = True
    with pytest.raises(NotImplementedError, match="binary"):
        u2d.get_file_entry()
    assert list((function_tmpdir / "ext").iterdir()) == []  # nothing written


def test_util3d_float64_layers_structured_and_disu(function_tmpdir):
    m = _model(function_tmpdir)
    u3d = Util3d(m, (2, 1, 3), np.float64, np.stack([DOUBLES] * 2), "h")
    u3d[1] = DOUBLES[::-1] + 1.0
    assert u3d.dtype == np.float64 and u3d[1].dtype == np.float64
    np.testing.assert_array_equal(u3d.array[1], DOUBLES[::-1] + 1.0)
    from flopy.mfusg import MfUsg

    mu = MfUsg(model_ws=str(function_tmpdir), structured=False)
    flat = np.array([23.1234567890123, 23.1234567890124, 22.75, 1.1, 1.2])
    jag = Util3d(mu, (2, None, np.array([3, 2])), np.float64, flat, "h")
    assert [u.shape for u in jag.util_2ds] == [(3,), (2,)]
    np.testing.assert_array_equal(jag.array, flat)


def test_load_bin_double_header_and_shared_cursor(function_tmpdir):
    data = _bin_record(DOUBLES[0], ilay=1) + _bin_record(DOUBLES[0] + 1, ilay=2)
    path = function_tmpdir / "shared.bin"
    path.write_bytes(data + b"TRAILING")
    with open(path, "rb") as f:
        h1, a1 = Util2d.load_bin(
            (1, 3), f, np.float64, "Head", header_precision="double"
        )
        assert f.tell() == 52 + 24
        h2, a2 = Util2d.load_bin(
            (1, 3), f, np.float64, "Head", header_precision="double"
        )
        assert f.read() == b"TRAILING"  # later data is left for the caller
        assert not f.closed  # caller's handle stays open
    assert h1["ilay"][0] == 1 and h2["ilay"][0] == 2
    np.testing.assert_array_equal(a1, DOUBLES)
    np.testing.assert_array_equal(a2, DOUBLES + 1)


@pytest.mark.parametrize("cut", [8, 40])
def test_load_bin_double_truncated(function_tmpdir, cut):
    """cut=8 truncates the data, cut=40 the header; own files are closed."""
    path = function_tmpdir / "short.bin"
    path.write_bytes(_bin_record(DOUBLES[0])[: -cut if cut == 8 else 52 - cut])
    before = _open_fds()
    with pytest.raises(ValueError, match=r"truncated|expected array size"):
        Util2d.load_bin(
            (1, 3), str(path), np.float64, "Head", header_precision="double"
        )
    assert _open_fds() == before
    with open(path, "rb") as f:
        with pytest.raises(ValueError):
            Util2d.load_bin((1, 3), f, np.float64, "Head", header_precision="double")
        assert not f.closed


@pytest.mark.parametrize(
    "line, dtype, cnstnt",
    [
        ("EXTERNAL -50 2 (FREE) 3", np.int32, 2),
        ("EXTERNAL -50 2.0 (FREE) 3", np.float32, 2.0),
        ("EXTERNAL -50 2.0 (FREE) 3", np.float64, 2.0),
        ("EXTERNAL -50 2.0", np.float64, 2.0),
        ("EXTERNAL -50", np.float32, 0.0),
        ("EXTERNAL -50 1.0 #comment", np.float64, 1.0),
    ],
)
def test_parse_control_record_negative_external_unit(line, dtype, cnstnt):
    """U1DREL*/U2DREL*: LOCAT<0 is binary and FMTIN/IPRN are not read."""
    cr = Util2d.parse_control_record(line, dtype=dtype)
    assert (cr["type"], cr["nunit"], cr["fmtin"], cr["negative_unit"]) == (
        "external",
        50,
        "(BINARY)",
        True,
    )
    assert cr["cnstnt"] == cnstnt


@pytest.mark.parametrize("dtype", [np.int32, np.float32, np.float64])
def test_parse_control_record_positive_units_unchanged(dtype):
    one = 1 if dtype == np.int32 else 1.0
    cr = Util2d.parse_control_record(f"EXTERNAL 50 {one} (BINARY) -1", dtype=dtype)
    assert (cr["type"], cr["nunit"], cr["fmtin"], cr["negative_unit"]) == (
        "external",
        50,
        "(BINARY)",
        False,
    )
    cr = Util2d.parse_control_record(f"EXTERNAL 50 {one} (FREE) 0", dtype=dtype)
    assert (cr["fmtin"], cr["iprn"], cr["negative_unit"]) == ("(FREE)", 0, False)
    cr = Util2d.parse_control_record(f"INTERNAL {one} (FREE) -1", dtype=dtype)
    assert (cr["type"], cr["cnstnt"]) == ("internal", one)
    fixed = f"{-50:10d}{one!s:>10}{'':>20}{0:10d}"
    cr = Util2d.parse_control_record(fixed, dtype=dtype)
    assert (cr["type"], cr["nunit"], cr["negative_unit"]) == ("external", 50, True)


def _load_ext(ws, line, data, dtype, **kwargs):
    from flopy.utils.mfreadnam import NamData

    (ws / "a.bin").write_bytes(data)
    fh = open(ws / "a.bin", "rb")
    try:
        ext = {50: NamData("DATA(BINARY)", str(ws / "a.bin"), fh, {})}
        rec = StringIO(line + "\n")
        rec.name = "rec.txt"
        return Util2d.load(
            rec, _model(ws), (1, 3), dtype, "h", ext_unit_dict=ext, **kwargs
        )
    finally:
        fh.close()


@pytest.mark.parametrize("double", [False, True])
def test_util2d_load_binary_multiplier_default_and_1d(function_tmpdir, double):
    dtype = np.float64 if double else np.float32
    data = _bin_record([1.0, 2.0, 3.0], double=double)
    # default (2-D readers, e.g. U2DREL8) applies the multiplier
    u2d = _load_ext(function_tmpdir, "EXTERNAL -50 2.0 (BINARY) -1", data, dtype)
    np.testing.assert_array_equal(u2d.array, [[2.0, 4.0, 6.0]])
    # the 1-D readers return before applying it
    with pytest.warns(UserWarning, match="not applied"):
        u1d = _load_ext(
            function_tmpdir,
            "EXTERNAL -50 2.0 (BINARY) -1",
            data,
            dtype,
            binary_multiplier=False,
        )
    np.testing.assert_array_equal(u1d.array, [[1.0, 2.0, 3.0]])
    assert u1d.dtype == dtype


def test_util2d_load_binary_needs_negative_unit_only_when_asked(function_tmpdir):
    data = _bin_record([1.0, 2.0, 3.0])
    u2d = _load_ext(function_tmpdir, "EXTERNAL 50 1.0 (BINARY) -1", data, np.float64)
    np.testing.assert_array_equal(u2d.array, [[1.0, 2.0, 3.0]])
    with pytest.raises(ValueError, match="negative unit"):
        _load_ext(
            function_tmpdir,
            "EXTERNAL 50 1.0 (BINARY) -1",
            data,
            np.float64,
            binary_negative_unit_only=True,
        )
    (function_tmpdir / "oc.bin").write_bytes(data)
    rec = StringIO("OPEN/CLOSE oc.bin 1.0 (BINARY) -1\n")
    rec.name = "rec.txt"
    m = _model(function_tmpdir)
    u2d = Util2d.load(rec, m, (1, 3), np.float64, "h")
    np.testing.assert_array_equal(u2d.array, [[1.0, 2.0, 3.0]])
    rec.seek(0)
    before = _open_fds()
    with pytest.raises(ValueError, match="OPEN/CLOSE"):
        Util2d.load(rec, m, (1, 3), np.float64, "h", binary_negative_unit_only=True)
    assert _open_fds() == before


def test_util2d_load_constant_float64_is_real(function_tmpdir):
    rec = StringIO(f"CONSTANT {NOT_REAL!r}\n")
    rec.name = "rec.txt"
    u2d = Util2d.load(rec, _model(function_tmpdir), (1, 3), np.float64, "h")
    assert u2d.dtype == np.float64
    np.testing.assert_array_equal(u2d.array, np.float64(np.float32(NOT_REAL)))
