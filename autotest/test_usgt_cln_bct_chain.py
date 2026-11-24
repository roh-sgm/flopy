"""CLN and BCT chain decay together (CHAINDECAY with SPATIALREACT), USG-T 2.8.1.

The whole model is written by FloPy (``write_input``), BCT and CLN included,
and compared with the same model whose BCT and CLN files are written by hand
from the USG-T 2.8.1 readers (glo2btnu1.f, cln2basu1.f, cln2props1.f). Two
species on a one-layer DISU grid with IDISP=1; species 2 is the daughter of
species 1. The generation factor of the daughter is SPTLRCT: per layer in BCT
(GW nodes) and per CLN node in CLN.
"""

from pathlib import Path

import numpy as np
import pytest

from autotest.test_usgt_cln_transport import (
    BCT_KW,
    _cln,
    _model,
    _records,
    _same,
    _transport,
)
from autotest.usgt_exe import USGT_EXE, requires_usgt
from flopy.mfusg import MfUsg, MfUsgBct, MfUsgLpf, MfUsgOc, MfUsgSms
from flopy.modflow import ModflowChd
from flopy.utils import Util2d

BCT_SPTLRCT = 0.5  # GW generation factor of species 2 from species 1
CLN_SPTLRCT = [0.8, 0.6]  # CLN generation factor, per CLN node


def _bct_kw(sptlrct=BCT_SPTLRCT):
    return {
        **BCT_KW,
        "mcomp": 2,
        "ifod": 1,
        "fodrw": [0.1, 0.0],
        "conc": [1.0, 0.0],
        "chaindecay": True,
        "spatialreact": True,
        "nparent": [0, 1],
        "jparent": [0, [1]],
        "sptlrct": [0, [sptlrct]],
    }


def _cln_kw(sptlrct=CLN_SPTLRCT):
    # species 1 has no parent; species 2: parent position 1 is species 1
    return {
        "sptlrct": [0.0, [sptlrct, 0.0]],
        "fodrw": [0.1, 0.0],
        "conc": [5.0, 0.0],
    }


def _build(ws, bct_kw, cln_kw):
    """Model with BCT and CLN built by FloPy (None: package left out)."""
    ml = _model(ws, exe=USGT_EXE)
    MfUsgLpf(ml, laytyp=0, hk=1.0, ss=1.0e-3)
    MfUsgSms(ml, linmeth=1)
    ModflowChd(ml, stress_period_data={0: [[0, 10.0, 10.0], [4, 9.0, 9.0]]})
    MfUsgBct(ml, **bct_kw)
    cln = _cln(ml, **cln_kw)
    MfUsgOc(ml, stress_period_data={(0, 4): ["save head", "save concentration"]})
    return ml, cln


def _run(ml):
    success, buff = ml.run_model(silent=True, report=True)
    assert success and any("normal termination" in str(b).lower() for b in buff)
    ws = Path(ml.model_ws)
    return (
        _records((ws / "m.clncn").read_bytes(), cln=True),
        _records((ws / "m.con").read_bytes(), cln=False),
    )


def _flopy_run(ws, bct_sptlrct=BCT_SPTLRCT, cln_sptlrct=CLN_SPTLRCT):
    ml, _ = _build(ws, _bct_kw(bct_sptlrct), _cln_kw(cln_sptlrct))
    ml.write_input()
    return _run(ml)


# BCT by hand: SGWT2BCT1G (glo2btnu1.f 2.8.1). Item 1 with the options on the
# same line (:146-170, options :64-141), PRSITY (:1562), ANGLEX (:1594), DL
# and DT (:1626, :1637), then per species (:1732): NPARENT (:1744), with
# NPARENT > 0 the JPARENT values (:1751) and, with SPATIALREACT, one SPTLRCT
# per parent and layer instead of STOTIO (:1757-1778), FODRW (IFOD=1,
# :1918-1923) and CONC.
HAND_BCT = """\
# BCT written by hand
 1 0 2 1 1 0 0 -999.0 1.0e-6 1 0 0.0 0 1 0 0 0 0 0 CHAINDECAY SPATIALREACT
CONSTANT 0.2   PRSITY
CONSTANT 0.0   ANGLEX
CONSTANT 0.0   DL
CONSTANT 0.0   DT
0   NPARENT species 1
CONSTANT 0.1   FODRW species 1
CONSTANT 1.0   CONC species 1
1   NPARENT species 2
CONSTANT 1   JPARENT species 2
CONSTANT 0.5   SPTLRCT species 2 parent 1 layer 1
CONSTANT 0.0   FODRW species 2
CONSTANT 0.0   CONC species 2
"""

# CLN by hand, flow part: item 0 NCLN ICLNNDS ICLNCB ICLNHD ICLNDD ICLNIB
# NCLNGWC NCONDUITYP and the SAVECLNCON option (cln2basu1.f :93-100, :198-200),
# NNDCLN (:237), one line per CLN node IFNO IFTYP IFDIR FLENG FELEV FANGLE
# IFLIN ICCWADI (:366-373), one line per CLN-GW connection IFNO IFNOD IFCON
# FSKIN FLENG FANISO ICGWADI (:636-642), the conduit IFNO FRAD CONDUITK
# (cln2props1.f :46-48), IBOUND (:688) and STRT (:692).
# Transport part: SGWT2BCT1CLN (glo2btnu1.f :2034-2165), DLX and ATXY (IDISP
# not 0, :2091-2096), then per species MCOMP SPTLRCT arrays (:2107-2118),
# FODRW (IFOD=1, :2134) and CONC (:2144-2158).
HAND_CLN = """\
# CLN written by hand
1 -1 0 0 0 0 2 1 SAVECLNCON 37
CONSTANT 2   NNDCLN
1 1 0 10.0 0.5 1.57 0 0
2 1 0 10.0 0.5 1.57 0 0
1 1 0 0.0 10.0 1.0 0
2 5 0 0.0 10.0 1.0 0
1 0.5 3.23e10
CONSTANT 1   IBOUND
CONSTANT 9.5   STRT
CONSTANT 0.0   DLX
CONSTANT 0.0   ATXY
CONSTANT 0.0   SPTLRCT species 1 parent position 1
CONSTANT 0.0   SPTLRCT species 1 parent position 2
CONSTANT 0.1   FODRW species 1
CONSTANT 5.0   CONC species 1
INTERNAL 1.0 (FREE) -1   SPTLRCT species 2 parent position 1
0.8 0.6
CONSTANT 0.0   SPTLRCT species 2 parent position 2
CONSTANT 0.0   FODRW species 2
CONSTANT 0.0   CONC species 2
"""


def _hand_run(ws):
    """Same model; the BCT and CLN files are the hand-written ones."""
    ml, _ = _build(ws, _bct_kw(), _cln_kw())
    ml.write_input(
        SelPackList=[
            p.name[0] for p in ml.packagelist if p.name[0] not in ("BCT", "CLN")
        ]
    )
    (Path(ws) / "m.bct").write_text(HAND_BCT)
    (Path(ws) / "m.cln").write_text(HAND_CLN)
    return _run(ml)


def test_cln_bct_chain_written_by_flopy(function_tmpdir):
    """Without engine: the BCT and CLN chain data FloPy writes."""
    ml, cln = _build(function_tmpdir, _bct_kw(), _cln_kw())
    ml.write_input()
    bct_text = (Path(function_tmpdir) / "m.bct").read_text().splitlines()
    assert "CHAINDECAY" in bct_text[1] and "SPATIALREACT" in bct_text[1]
    body = [ln.split("#")[0].split() for ln in bct_text[6:]]
    assert body == [
        ["0"],
        ["CONSTANT", "1.000000E-01"],
        ["CONSTANT", "1.000000E+00"],
        ["1"],
        ["INTERNAL", "1", "(1I10)", "-1"],
        ["1"],
        ["CONSTANT", "5.000000E-01"],
        ["CONSTANT", "0.000000E+00"],
        ["CONSTANT", "0.000000E+00"],
    ]
    got = _transport((Path(function_tmpdir) / "m.cln").read_text())
    assert got == [
        ("dll", [0.0, 0.0]),
        ("dlm", [0.0, 0.0]),
        ("sptlrct", [0.0, 0.0]),
        ("sptlrct", [0.0, 0.0]),
        ("fodrw", [0.1, 0.1]),
        ("conc", [5.0, 5.0]),
        ("sptlrct", [pytest.approx(v) for v in CLN_SPTLRCT]),
        ("sptlrct", [0.0, 0.0]),
        ("fodrw", [0.0, 0.0]),
        ("conc", [0.0, 0.0]),
    ]


@requires_usgt
def test_cln_bct_chain_flopy_equals_hand(function_tmpdir):
    flopy_out = _flopy_run(function_tmpdir / "flopy")
    hand_out = _hand_run(function_tmpdir / "hand")
    assert [x[0] for x in flopy_out[0]] == ["CONC01", "CONC02"]
    assert [x[0] for x in flopy_out[1]] == ["CONC01", "CONC02"]
    assert _same(flopy_out[0], hand_out[0]) and _same(flopy_out[1], hand_out[1])
    # the daughter is generated in the GW and in the CLN nodes
    assert np.all(flopy_out[1][1][1] > 0.0) and np.all(flopy_out[0][1][1] > 0.0)


@requires_usgt
def test_cln_bct_chain_sensitivity(function_tmpdir):
    base = _flopy_run(function_tmpdir / "base")
    # CLN SPTLRCT: the daughter in the CLN nodes changes, the parent does not
    cln0 = _flopy_run(function_tmpdir / "cln0", cln_sptlrct=0.0)
    assert np.array_equal(base[0][0][1], cln0[0][0][1])
    assert np.all(base[0][1][1] != cln0[0][1][1])
    # BCT SPTLRCT: the daughter in the matrix changes, the parent does not
    bct0 = _flopy_run(function_tmpdir / "bct0", bct_sptlrct=0.0)
    assert np.array_equal(base[1][0][1], bct0[1][0][1])
    assert np.all(base[1][1][1] != bct0[1][1][1])


def _chain_data(ml):
    bct, cln = ml.get_package("BCT"), ml.get_package("CLN")
    return (
        list(bct.nparent),
        bct.jparent[1].array.tolist(),
        # BCT: per parent, one Util2d per layer
        [[[u.array.ravel().tolist() for u in p] for p in s] for s in bct.sptlrct[1:]],
        [[u.array.ravel().tolist() for u in s] for s in cln.sptlrct],
        [u.array.ravel().tolist() for u in cln.conc],
    )


def _load(ws):
    return MfUsg.load("m.nam", model_ws=str(ws), exe_name=USGT_EXE, check=False)


def test_cln_bct_chain_load_and_rewrite(function_tmpdir):
    """Without engine: load the whole model and write it again."""
    ws = Path(function_tmpdir)
    ml, _ = _build(ws / "flopy", _bct_kw(), _cln_kw())
    ml.write_input()
    loaded = _load(ws / "flopy")
    assert _chain_data(loaded) == (
        [0, 1],
        [1],
        [[[[BCT_SPTLRCT] * 5]]],
        [
            [[0.0, 0.0], [0.0, 0.0]],
            [[pytest.approx(v) for v in CLN_SPTLRCT], [0.0, 0.0]],
        ],
        [[5.0, 5.0], [0.0, 0.0]],
    )
    assert all(
        isinstance(u, Util2d) for s in loaded.get_package("CLN").sptlrct for u in s
    )
    loaded.change_model_ws(str(ws / "rewritten"))
    loaded.write_input()
    for name in ("m.cln",):
        assert _transport((ws / "rewritten" / name).read_text()) == _transport(
            (ws / "flopy" / name).read_text()
        )
    assert _chain_data(_load(ws / "rewritten")) == _chain_data(loaded)


@requires_usgt
def test_cln_bct_chain_load_rewrite_runs_the_same(function_tmpdir):
    ws = Path(function_tmpdir)
    flopy_out = _flopy_run(ws / "flopy")
    loaded = _load(ws / "flopy")
    loaded.change_model_ws(str(ws / "rewritten"))
    loaded.write_input()
    out = _run(loaded)
    assert _same(flopy_out[0], out[0]) and _same(flopy_out[1], out[1])
