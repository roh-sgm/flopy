"""
Autotests for MfUsgTvm (Time-Variant Materials package).

Tests cover:
  1. Programmatic construction → write → load round-trip
  2. Empty SP boundaries (all counts zero)
  3. Transport vs non-transport field counts in file
  4. Log base and interpolation-control round-trip
  5. Node numbering: 0-based internal ↔ 1-based file
  6. Mixed property changes across multiple boundaries
  7. Load tolerance for files with fewer boundaries than nper+1
"""

from __future__ import annotations

import io
import textwrap

import numpy as np
import pytest

import flopy
from flopy.mfusg.mfusgtvm import MfUsgTvm, _DTYPE, _PROPS


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_model(nper: int = 3, has_bct: bool = False):
    """Return a bare MfUsg stub sufficient for TVM construction."""
    import unittest.mock as mock

    # No spec so we can add arbitrary attributes without AttributeError
    m = mock.MagicMock()
    m.structured = False
    m.verbose = False
    m.nper = nper
    m.modelname = "test"
    m.namefile_ext = "nam"
    m.model_ws = "."
    m.external_path = None
    m.free_format_input = True
    m.array_free_format = True
    m.version = "mfusg"
    m.get_package.return_value = object() if has_bct else None  # truthy/falsy
    m.add_package.return_value = None
    m.pop_package.return_value = None
    # Make it pass isinstance(model, MfUsg)
    m.__class__ = flopy.mfusg.MfUsg
    return m


def _make_sp(nodes_values: dict[str, list[tuple[int, float]]]) -> dict:
    """Build a stress-period dict from {prop: [(node, value), ...]}."""
    out = {}
    for prop, pairs in nodes_values.items():
        out[prop] = np.array(pairs, dtype=_DTYPE).view(np.recarray)
    return out


def _write_and_reload(tvm: MfUsgTvm, tmp_path) -> MfUsgTvm:
    """Write TVM to a temp file and load it back using the same model."""
    path = tmp_path / "test.tvm"
    tvm.fn_path = str(path)
    tvm.write_file()
    loaded = MfUsgTvm.load(str(path), tvm.parent, nper=tvm.parent.nper)
    return loaded, path.read_text()


# ---------------------------------------------------------------------------
# 1. Programmatic round-trip — HK and SY changes
# ---------------------------------------------------------------------------

class TestTvmRoundTrip:
    def test_hk_sy_round_trip(self, tmp_path):
        """HK and SY records survive write → load with correct values."""
        nper = 3
        m = _make_model(nper=nper, has_bct=False)

        dtype = MfUsgTvm.dtype()
        sp_data = {
            0: _make_sp({"hk": [(0, 1e-4), (1, 2e-4)],
                         "sy": [(0, 0.30)]}),
            1: _make_sp({"hk": [(0, 3e-4), (1, 4e-4)],
                         "sy": [(0, 0.25)]}),
            2: _make_sp({"hk": [(0, 5e-4)]}),
            # boundary 3 (end of SP3) intentionally absent → all zeros
        }
        tvm = MfUsgTvm(
            m,
            tvmlogbasehk=10.0,
            tvmlogbasesy=0.0,
            stress_period_data=sp_data,
        )

        loaded, text = _write_and_reload(tvm, tmp_path)

        # Check global header values
        assert loaded.tvmlogbasehk == pytest.approx(10.0)
        assert loaded.tvmlogbasesy == pytest.approx(0.0)

        # Check boundary 0 HK records
        hk0 = loaded.stress_period_data[0]["hk"]
        assert len(hk0) == 2
        assert hk0[0]["node"] == 0
        assert hk0[0]["value"] == pytest.approx(1e-4)
        assert hk0[1]["node"] == 1
        assert hk0[1]["value"] == pytest.approx(2e-4)

        # Check boundary 0 SY records
        sy0 = loaded.stress_period_data[0]["sy"]
        assert len(sy0) == 1
        assert sy0[0]["node"] == 0
        assert sy0[0]["value"] == pytest.approx(0.30)

        # Boundary 3 (absent) → empty sp dict loaded
        sp3 = loaded.stress_period_data.get(3, {})
        for p in _PROPS[:-1]:  # no transport → no 'por'
            assert len(sp3.get(p, [])) == 0

    def test_all_properties_round_trip(self, tmp_path):
        """All 6 properties (transport model) survive round-trip."""
        nper = 2
        m = _make_model(nper=nper, has_bct=True)

        sp_data = {
            1: _make_sp({
                "hk":    [(0, 1.0)],
                "vka":   [(1, 2.0)],
                "ss":    [(2, 1e-5)],
                "sy":    [(3, 0.20)],
                "ddftr": [(4, 0.01)],
                "por":   [(5, 0.35)],
            }),
        }
        tvm = MfUsgTvm(
            m,
            itvmprint=1,
            tvmlogbasehk=10.0,
            tvmlogbasepor=0.0,
            stress_period_data=sp_data,
        )

        loaded, text = _write_and_reload(tvm, tmp_path)

        sp1 = loaded.stress_period_data[1]
        assert sp1["hk"][0]["value"]    == pytest.approx(1.0)
        assert sp1["vka"][0]["value"]   == pytest.approx(2.0)
        assert sp1["ss"][0]["value"]    == pytest.approx(1e-5)
        assert sp1["sy"][0]["value"]    == pytest.approx(0.20)
        assert sp1["ddftr"][0]["value"] == pytest.approx(0.01)
        assert sp1["por"][0]["value"]   == pytest.approx(0.35)


# ---------------------------------------------------------------------------
# 2. Empty SP boundaries
# ---------------------------------------------------------------------------

class TestEmptyBoundaries:
    def test_all_zero_boundaries(self, tmp_path):
        """Model with no TVM changes still writes nper+1 boundary blocks."""
        nper = 4
        m = _make_model(nper=nper, has_bct=False)
        tvm = MfUsgTvm(m)

        path = tmp_path / "empty.tvm"
        tvm.fn_path = str(path)
        tvm.write_file()
        text = path.read_text()

        # Should have exactly nper+1 = 5 SP-header lines (all-zero counts)
        sp_lines = [ln for ln in text.splitlines()
                    if "Stress period" in ln]
        assert len(sp_lines) == nper + 1

        # All counts should be zero
        for line in sp_lines:
            counts = [int(t) for t in line.split()[:5]]  # 5 fields (no transport)
            assert all(c == 0 for c in counts)

    def test_label_format(self, tmp_path):
        """Boundary labels follow the Vistas convention."""
        m = _make_model(nper=3, has_bct=False)
        tvm = MfUsgTvm(m)
        path = tmp_path / "labels.tvm"
        tvm.fn_path = str(path)
        tvm.write_file()
        lines = path.read_text().splitlines()

        sp_lines = [ln for ln in lines if "Stress period" in ln]
        assert "1 start" in sp_lines[0]
        assert "1 end"   in sp_lines[1]
        assert "2"       in sp_lines[2]
        assert "3"       in sp_lines[3]


# ---------------------------------------------------------------------------
# 3. Transport vs non-transport field counts
# ---------------------------------------------------------------------------

class TestTransportFieldCounts:
    def test_no_transport_writes_5_sp_fields(self, tmp_path):
        """Without BCT: SP header has 5 integers, global header has 6 values."""
        m = _make_model(nper=2, has_bct=False)
        tvm = MfUsgTvm(m)
        path = tmp_path / "notransport.tvm"
        tvm.fn_path = str(path)
        tvm.write_file()
        lines = path.read_text().splitlines()

        # Global header line: first non-comment, non-SP line
        global_hdr = next(ln for ln in lines if not ln.startswith("#")
                         and "Stress" not in ln)
        global_tokens = global_hdr.split()
        assert len(global_tokens) == 6  # ITVMPRINT + 5 log bases

        # SP header lines: count only tokens BEFORE the "Stress" keyword
        sp_lines = [ln for ln in lines if "Stress period" in ln]
        for ln in sp_lines:
            before_stress = ln[:ln.lower().index("stress")].split()
            assert len(before_stress) == 5

    def test_with_transport_writes_6_sp_fields(self, tmp_path):
        """With BCT: SP header has 6 integers, global header has 7 values."""
        m = _make_model(nper=2, has_bct=True)
        tvm = MfUsgTvm(m)
        path = tmp_path / "transport.tvm"
        tvm.fn_path = str(path)
        tvm.write_file()
        lines = path.read_text().splitlines()

        global_hdr = next(ln for ln in lines if not ln.startswith("#")
                         and "Stress" not in ln)
        global_tokens = global_hdr.split()
        assert len(global_tokens) == 7  # ITVMPRINT + 6 log bases

        sp_lines = [ln for ln in lines if "Stress period" in ln]
        for ln in sp_lines:
            before_stress = ln[:ln.lower().index("stress")].split()
            assert len(before_stress) == 6


# ---------------------------------------------------------------------------
# 4. Log base and interpolation-control round-trip
# ---------------------------------------------------------------------------

class TestLogBases:
    @pytest.mark.parametrize("base,expected", [
        (0.0,  0.0),    # linear
        (10.0, 10.0),   # log base 10
        (-1.0, -1.0),   # step (HK only — USG-T rejects for Ss/Sy)
        (2.0,  2.0),    # log base 2
    ])
    def test_hk_log_base_survives(self, tmp_path, base, expected):
        m = _make_model(nper=1, has_bct=False)
        tvm = MfUsgTvm(m, tvmlogbasehk=base)
        loaded, _ = _write_and_reload(tvm, tmp_path)
        assert loaded.tvmlogbasehk == pytest.approx(expected)

    def test_itvmprint_survives(self, tmp_path):
        for v in (0, 1, 2):
            m = _make_model(nper=1, has_bct=False)
            tvm = MfUsgTvm(m, itvmprint=v)
            loaded, _ = _write_and_reload(tvm, tmp_path)
            assert loaded.itvmprint == v


# ---------------------------------------------------------------------------
# 5. Node numbering
# ---------------------------------------------------------------------------

class TestNodeNumbering:
    def test_0based_stored_1based_in_file(self, tmp_path):
        """FloPy stores 0-based nodes; file contains 1-based."""
        m = _make_model(nper=2, has_bct=False)
        sp_data = {0: _make_sp({"hk": [(0, 0.01), (99, 0.02)]})}
        tvm = MfUsgTvm(m, stress_period_data=sp_data)

        path = tmp_path / "nodes.tvm"
        tvm.fn_path = str(path)
        tvm.write_file()

        # Parse only data record lines (exactly 2 tokens: node value)
        # Skip comment, SP header ("Stress"), and global header (≥6 tokens)
        data_nodes = []
        for ln in path.read_text().splitlines():
            if "Stress" in ln or ln.startswith("#"):
                continue
            toks = ln.split()
            if len(toks) != 2:   # global header has 6 tokens; data records have 2
                continue
            try:
                node = int(toks[0])
                float(toks[1])
                data_nodes.append(node)
            except ValueError:
                pass

        # Nodes 0 and 99 (0-based) → must appear as 1 and 100 (1-based)
        assert 1   in data_nodes
        assert 100 in data_nodes
        # 0-based node 0 must NOT appear as 0 in any data record
        assert 0 not in data_nodes

    def test_load_converts_to_0based(self, tmp_path):
        """Loaded data has 0-based nodes regardless of file content."""
        m = _make_model(nper=2, has_bct=False)
        sp_data = {1: _make_sp({"sy": [(5, 0.3), (10, 0.4)]})}
        tvm = MfUsgTvm(m, stress_period_data=sp_data)
        loaded, _ = _write_and_reload(tvm, tmp_path)

        sy = loaded.stress_period_data[1]["sy"]
        assert sy[0]["node"] == 5
        assert sy[1]["node"] == 10


# ---------------------------------------------------------------------------
# 6. Mixed properties across multiple boundaries
# ---------------------------------------------------------------------------

class TestMixedProperties:
    def test_different_props_per_boundary(self, tmp_path):
        """Each boundary can change a different property independently."""
        nper = 3
        m = _make_model(nper=nper, has_bct=False)
        sp_data = {
            0: _make_sp({"hk": [(0, 1.0)]}),
            1: _make_sp({"vka": [(1, 0.1)]}),
            2: _make_sp({"ss": [(2, 1e-5)]}),
            3: _make_sp({"sy": [(3, 0.25)]}),
        }
        tvm = MfUsgTvm(m, stress_period_data=sp_data)
        loaded, _ = _write_and_reload(tvm, tmp_path)

        assert loaded.stress_period_data[0]["hk"][0]["value"]  == pytest.approx(1.0)
        assert loaded.stress_period_data[1]["vka"][0]["value"] == pytest.approx(0.1)
        assert loaded.stress_period_data[2]["ss"][0]["value"]  == pytest.approx(1e-5)
        assert loaded.stress_period_data[3]["sy"][0]["value"]  == pytest.approx(0.25)


# ---------------------------------------------------------------------------
# 7. Load from string — real-world snippet (PUCOBRE format)
# ---------------------------------------------------------------------------

class TestLoadFromString:
    # Minimal TVM snippet matching the Vistas/PUCOBRE format (no transport)
    _TVM_TEXT = textwrap.dedent("""\
        # TVM package (test fixture)
            0   10.0 10.0 10.0 10.0 10.0 10.0                ITVMPRINT TVMLOGBASEHK ...
                 0         0         0         0         0         0         Stress period number 1 start
                 0         0         0         0         0         0         Stress period number 1 end
                 2         0         0         1         0         0         Stress period number 2
               101     0.00855   ITVMHK  HKNEW
               202     0.00456   ITVMHK  HKNEW
               305     0.30000   ITVMSY  SYNEW
                 0         0         0         0         0         0         Stress period number 3
    """)

    def test_load_from_string_transport(self, tmp_path):
        """Load a 7-value global header (transport format) correctly."""
        path = tmp_path / "snippet.tvm"
        path.write_text(self._TVM_TEXT)
        m = _make_model(nper=3, has_bct=True)  # transport → 7 fields expected
        tvm = MfUsgTvm.load(str(path), m, nper=3)

        assert tvm.tvmlogbasehk == pytest.approx(10.0)
        # Boundary 2: 2 HK changes
        hk2 = tvm.stress_period_data[2].get("hk", [])
        assert len(hk2) == 2
        assert hk2[0]["node"]  == 100   # 101 - 1
        assert hk2[0]["value"] == pytest.approx(0.00855)
        assert hk2[1]["node"]  == 201   # 202 - 1
        # Boundary 2: 1 SY change
        sy2 = tvm.stress_period_data[2].get("sy", [])
        assert len(sy2) == 1
        assert sy2[0]["node"]  == 304   # 305 - 1
        assert sy2[0]["value"] == pytest.approx(0.30)

    def test_load_no_transport(self, tmp_path):
        """A 6-value global header (no transport) loads without 'por'."""
        text_no_transport = textwrap.dedent("""\
                 0    0.0    0.0    0.0    0.0    0.0
                 0         0         0         0         0    Stress period number 1 start
                 0         0         0         0         0    Stress period number 1 end
        """)
        path = tmp_path / "notransport.tvm"
        path.write_text(text_no_transport)
        m = _make_model(nper=1, has_bct=False)
        tvm = MfUsgTvm.load(str(path), m, nper=1)
        assert tvm.tvmlogbasepor == pytest.approx(0.0)
        # 'por' key should be absent from all boundaries
        for ibnd, sp in tvm.stress_period_data.items():
            assert "por" not in sp

    def test_load_tolerates_fewer_boundaries(self, tmp_path):
        """File with fewer than nper+1 boundaries loads without error."""
        text = textwrap.dedent("""\
                 0    0.0    0.0    0.0    0.0    0.0
                 1         0         0         0         0    Stress period number 1 start
               999     0.001
        """)
        path = tmp_path / "short.tvm"
        path.write_text(text)
        m = _make_model(nper=5, has_bct=False)
        tvm = MfUsgTvm.load(str(path), m, nper=5)
        # Only boundary 0 should be present
        assert 0 in tvm.stress_period_data
        assert tvm.stress_period_data[0]["hk"][0]["node"] == 998  # 999-1


# ---------------------------------------------------------------------------
# 8.  _normalise_sp — accepts lists of tuples as well as recarrays
# ---------------------------------------------------------------------------

class TestNormalise:
    def test_list_of_tuples_accepted(self):
        sp_in = {"hk": [(0, 1.0), (1, 2.0)]}
        sp_out = MfUsgTvm._normalise_sp(sp_in)
        assert isinstance(sp_out["hk"], np.recarray)
        assert sp_out["hk"][0]["node"] == 0
        assert sp_out["hk"][1]["value"] == pytest.approx(2.0)

    def test_plain_array_accepted(self):
        arr = np.array([(3, 0.5)], dtype=_DTYPE)
        sp_in = {"sy": arr}
        sp_out = MfUsgTvm._normalise_sp(sp_in)
        assert isinstance(sp_out["sy"], np.recarray)

    def test_unknown_keys_ignored(self):
        sp_in = {"hk": [(0, 1.0)], "unknown_prop": [(0, 9.9)]}
        sp_out = MfUsgTvm._normalise_sp(sp_in)
        assert "unknown_prop" not in sp_out
        assert "hk" in sp_out
