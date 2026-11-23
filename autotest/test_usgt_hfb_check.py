"""MfUsgHfb under the package and model checks.

``MfUsg.load`` runs ``model.check()``, and ``Package.check`` reads
``stress_period_data`` as an MfList whenever the attribute is set. MfUsgHfb
keeps a plain dict there (empty for static barriers, one recarray per stress
period with TRANSIENT_HFB), so the check must leave it alone: as with
ModflowHfb there is no HFB-specific check and ``check`` returns None.

The model files are the HFB example of ``mfusg_test/rch_evt_tests`` (copied,
because the checks write next to the model files) and small hand-built
packages. No engine is needed.
"""

import shutil
from pathlib import Path

import numpy as np
import pytest

from flopy.mfusg import MfUsg, MfUsgHfb
from flopy.modflow import ModflowDis

_NAM = "usg_rch_evt_hfb.nam"
_CASES = ["static", "transient", "parameters"]
# text written after the heading line, taken from the code before the change
_TEXT = {
    "static": (
        "         0         0         2\n"
        "         1         2       0.0001\n"
        "         3         4       0.0002\n"
    ),
    "transient": (
        "         0         0         2  TRANSIENT_HFB\n"
        "1\n"
        "         1         2       0.0001\n"
        "         3         4       0.0002\n"
        "1\n"
        "         5         6       0.0003\n"
        "         7         8       0.0004\n"
        "-1\n"
    ),
    "parameters": (
        "         1         2         1\n"
        "hp HFB 2.0 2\n"
        "         1         2          0.5\n"
        "         3         4         0.25\n"
        "         5         6        0.001\n"
        "         1\n"
        "hp\n"
    ),
}


def _example(example_data_path, ws):
    shutil.copytree(example_data_path / "mfusg_test" / "rch_evt_tests", ws)
    return ws


def _model(ws, name, nper=3):
    ml = MfUsg(structured=False, model_ws=str(ws), modelname=name)
    ModflowDis(ml, nlay=1, nrow=1, ncol=1, nper=nper)
    return ml


def _rows(*rows):
    dtype = MfUsgHfb.get_default_dtype(structured=False)
    return np.array(list(rows), dtype=dtype).view(np.recarray)


def _hfb(ml, case, unitnumber=None):
    """Add to ml an MfUsgHfb with static, transient or parameter barriers."""
    if case == "static":
        return MfUsgHfb(ml, hfb_data=_rows((0, 1, 1.0e-4), (2, 3, 2.0e-4)))
    if case == "transient":
        return MfUsgHfb(
            ml,
            hfb_data=_rows((0, 1, 1.0e-4), (2, 3, 2.0e-4)),
            transient=True,
            stress_period_data={1: _rows((4, 5, 3.0e-4), (6, 7, 4.0e-4))},
            unitnumber=unitnumber,
        )
    return MfUsgHfb(
        ml,
        hfb_data=_rows((4, 5, 1.0e-3)),
        parameters={"hp": {"parval": "2.0", "data": [(0, 1, 0.5), (2, 3, 0.25)]}},
        acthfb_names=["hp"],
        unitnumber=unitnumber,
    )


def test_hfb_example_load_runs_model_check(function_tmpdir, example_data_path):
    """Loading the HFB example with the default check=True does not fail."""
    ws = _example(example_data_path, function_tmpdir / "model")
    m = MfUsg.load(_NAM, model_ws=ws)
    hfb = m.get_package("HFB6")
    assert isinstance(hfb, MfUsgHfb) and not hfb.transient
    assert hfb.stress_period_data == {}
    assert list(hfb.hfb_data["node1"]) == [21, 22]
    assert (ws / "usg_rch_evt_hfb.chk").is_file()
    assert hfb.check(f=None, verbose=False) is None
    assert m.check(f=None, verbose=False) is not None


@pytest.mark.parametrize("case", ["transient", "parameters"])
def test_hfb_example_variant_runs_model_check(function_tmpdir, example_data_path, case):
    """The example with transient or parameter barriers: the model check
    works on the built model and on the model loaded back with check=True."""
    ws = _example(example_data_path, function_tmpdir / "model")
    m = MfUsg.load(_NAM, model_ws=ws, check=False)
    m.remove_package("HFB6")
    _hfb(m, case, unitnumber=26)  # the unit of the example; 29 is its DISU
    assert m.check(f=None, verbose=False) is not None
    m.write_input()

    hfb = MfUsg.load(_NAM, model_ws=ws).get_package("HFB6")
    if case == "transient":
        assert hfb.transient and set(hfb.stress_period_data) == {0, 1}
        assert list(hfb.stress_period_data[1]["node1"]) == [4, 6]
    else:
        assert hfb.acthfb_names == ["hp"]
        assert list(hfb.parameters["hp"]["data"]["node1"]) == [0, 2]
    assert hfb.check(f=None, verbose=False) is None


@pytest.mark.parametrize("case", _CASES)
def test_hfb_package_check(function_tmpdir, case):
    """check() returns None for a built and for a loaded package and, given a
    file name, writes the base-class notice (the path ModflowHfb takes)."""
    hfb = _hfb(_model(function_tmpdir, case), case)
    assert hfb.check(f=None, verbose=False) is None
    hfb.write_file()
    loaded = MfUsgHfb.load(hfb.fn_path, _model(function_tmpdir, f"{case}2"), nper=3)
    assert loaded.check(f=None, verbose=False) is None
    assert loaded.check(f="hfb.chk", verbose=False) is None
    assert (function_tmpdir / "hfb.chk").read_text() == (
        "check method not implemented for HFB6 Package."
    )


@pytest.mark.parametrize("case", _CASES)
def test_hfb_stress_period_data_is_still_a_dict(function_tmpdir, case):
    spd = _hfb(_model(function_tmpdir, case), case).stress_period_data
    assert type(spd) is dict
    if case == "transient":
        assert set(spd) == {1} and isinstance(spd[1], np.recarray)
    else:
        assert spd == {}


@pytest.mark.parametrize("case", _CASES)
def test_hfb_written_text_unchanged(function_tmpdir, case):
    hfb = _hfb(_model(function_tmpdir, case), case)
    hfb.write_file()
    assert Path(hfb.fn_path).read_text() == f"{hfb.heading}\n" + _TEXT[case]
