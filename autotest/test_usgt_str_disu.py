"""STR on an unstructured (DISU) MfUsg model.

MfUsg registers the base ``ModflowStr`` for STR, which reads and writes the
node-based reach list of an unstructured grid. With free format this is the
record USG-T reads (gwf2str7u1.f:916-925 in 2.8.1 and 2.7.0: NODE, SEG, REACH
and five reals). SUB and SWT keep their guards.

The model files are the STR example of ``mfusg_test/rch_evt_tests`` (copied,
because the load writes the check summary next to the model files). No engine
is needed.
"""

import shutil

import pytest

from flopy.mfusg import MfUsg, MfUsgSub, MfUsgSwt
from flopy.modflow import ModflowDis, ModflowStr

_NAM = "usg_rch_evt_str.nam"
# text after the heading line that upstream FloPy (develop 996687f3) writes
# when it loads and rewrites this example
_TEXT = (
    "3 1 0 0 1 86400.0 53 50 \n"
    "         3         0         0  # stress period 1\n"
    "8 1 1 1.0 0.0 0.10000000149011612 2.0 2.5 \n"
    "5 1 2 0.0 0.0 0.10000000149011612 1.899999976158142 2.4000000953674316 \n"
    "2 1 3 0.0 0.0 0.10000000149011612 1.7999999523162842 2.299999952316284 \n"
    "1.0 0.0010000000474974513 0.10000000149011612 \n"
    "1.0 0.0010000000474974513 0.10000000149011612 \n"
    "1.0 0.0010000000474974513 0.10000000149011612 \n"
    "        -1         0         0  # stress period 2\n"
)


def _example(example_data_path, ws):
    shutil.copytree(example_data_path / "mfusg_test" / "rch_evt_tests", ws)
    return ws


def test_str_disu_example_loads(function_tmpdir, example_data_path):
    """The STR example loads on its DISU grid with the base class."""
    ws = _example(example_data_path, function_tmpdir / "model")
    m = MfUsg.load(_NAM, model_ws=ws)
    assert not m.structured
    pkg = m.get_package("STR")
    assert type(pkg) is ModflowStr
    reaches = pkg.stress_period_data[0]
    assert reaches.dtype.names[:3] == ("node", "segment", "reach")
    assert list(reaches["node"]) == [7, 4, 1]  # 0-based; 8, 5, 2 in the file
    assert list(reaches["segment"]) == [1, 1, 1]
    assert list(reaches["reach"]) == [1, 2, 3]


def test_str_disu_rewrite_matches_upstream(function_tmpdir, example_data_path):
    """Load and rewrite: the STR file is the text the base class writes, and
    it loads back to the same reaches."""
    ws = _example(example_data_path, function_tmpdir / "model")
    m = MfUsg.load(_NAM, model_ws=ws)
    pkg = m.get_package("STR")
    m.write_input()
    assert (ws / "usg_rch_evt.str").read_text() == f"{pkg.heading}\n" + _TEXT

    again = MfUsg.load(_NAM, model_ws=ws).get_package("STR")
    for name in pkg.stress_period_data[0].dtype.names:
        assert list(again.stress_period_data[0][name]) == list(
            pkg.stress_period_data[0][name]
        ), name


@pytest.mark.parametrize("structured", [False, True])
def test_str_registry_is_base_class(structured):
    assert MfUsg(structured=structured).mfnam_packages["str"] is ModflowStr


@pytest.mark.parametrize("cls, label", [(MfUsgSub, "SUB"), (MfUsgSwt, "SWT")])
def test_sub_swt_still_guarded(function_tmpdir, cls, label):
    ml = MfUsg(structured=False, model_ws=str(function_tmpdir))
    ModflowDis(ml, nlay=1, nrow=1, ncol=1)
    assert ml.mfnam_packages[label.lower()] is cls
    with pytest.raises(NotImplementedError, match=f"{label}.*compatibility"):
        cls(ml)
