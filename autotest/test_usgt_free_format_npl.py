"""free_format_npl of an unstructured (DISU) model loaded by MfUsg.

``MfUsg.load`` leaves ``free_format_npl`` unset, so every ``(FREE)`` array
record is written on one line: a layer of HK, or the whole JA, CL12 and FAHL
arrays of DISU. Setting it after the load wraps those arrays at that many
values per line. Arrays written with a Fortran format, such as the DISU top
and bottom ``(<nodes of the layer>E15.6)``, keep one layer per line either way.

USG-T reads ``(FREE)`` arrays with list-directed input, which takes the
values from as many lines as needed: ``READ(LOCAT,*) (A(J),J=1,JJ)`` in
utl7u1.f, U1DREL :754-755, U1DREL8 :891-892, U1DINT :1009-1010, U2DINT
:1145-1146, U2DREL :1375-1376 and U2DREL8 :1535-1536 (the file is the same in
2.8.1 and 2.7.0). The engine runs write each example both ways and need
USGT_EXE.
"""

from pathlib import Path

import pytest

from autotest.test_usg_transport_exe import _nam_listing, _usgt_run_report
from autotest.usgt_exe import USGT_EXE, requires_usgt
from flopy.mbase import run_model
from flopy.mfusg import MfUsg

_EX8 = ("mfusg_transport/Ex8_Lake", "lak_usg_01.nam", 289)
# name: (example folder, name file, nodes per layer, output files compared)
_MODELS = {
    "ex8-lake-289": (*_EX8, ("lak_usg_01.hds", "lak_usg_01.con")),
    "ex7-multispecies-2000": (
        "mfusg_transport/Ex7_Matrix_Diffusion/Multispecies",
        "USG_Multispecies.nam",
        2000,
        ("USG_Multispecies.hds", "USG_Multispecies.con"),
    ),
    "ex3-nest-11716": (
        "mfusg_transport/Ex3_CLN_Conduit/Nest",
        "Conduit_Nest.nam",
        11716,
        (
            "Conduit_Nest.hds",
            "Conduit_Nest.con",
            "Conduit_Nest_CLN.HDS",
            "Conduit_Nest_CLN.CON",
        ),
    ),
}


def _write(example_data_path, ws, npl, model=_EX8, set_npl=True):
    """Load an example without touching it and write it into ws."""
    subdir, namefile, _ = model
    ml = MfUsg.load(namefile, model_ws=str(example_data_path / subdir), check=False)
    if set_npl:
        ml.free_format_npl = npl
    ml.change_model_ws(str(ws))
    ml.write_input()
    return ml


def _longest_line(ws):
    """(characters, values, file name) of the longest line of the files in ws."""
    best = (0, 0, "")
    for path in sorted(Path(ws).iterdir()):
        if path.is_file():
            for line in path.read_text(errors="replace").splitlines():
                best = max(best, (len(line), len(line.split()), path.name))
    return best


def _free_array_values(path):
    """Most values on one line among the (FREE) arrays of a file."""
    most, free = 0, False
    for line in Path(path).read_text().splitlines():
        if line.split()[:1] in (["INTERNAL"], ["CONSTANT"]):  # control record
            free = "(FREE)" in line
        elif free:
            most = max(most, len(line.split()))
    return most


def test_load_disu_leaves_free_format_npl_unset(example_data_path):
    subdir, namefile, _ = _EX8
    ml = MfUsg.load(namefile, model_ws=str(example_data_path / subdir), check=False)
    assert not ml.structured
    assert ml.free_format_npl is None


def test_rewrite_disu_default_writes_each_free_array_on_one_line(
    function_tmpdir, example_data_path
):
    ml = _write(example_data_path, function_tmpdir, None, set_npl=False)
    assert _free_array_values(ml.disu.fn_path) == ml.disu.njag


def test_free_format_npl_set_after_load_is_honored(function_tmpdir, example_data_path):
    ml = _write(example_data_path, function_tmpdir, 10)
    assert _free_array_values(ml.disu.fn_path) == 10


@requires_usgt
@pytest.mark.parametrize("name", list(_MODELS))
def test_usgt_exe_long_free_format_lines(function_tmpdir, example_data_path, name):
    """The engine reads an example written with each (FREE) array on one
    line, and gives the same heads and concentrations as with 10 values per
    line (the longest lines are then the layers written with a format)."""
    *model, outputs = _MODELS[name]
    namefile, nodelay = model[1], model[2]
    results = {}
    for npl in (None, 10):
        ws = function_tmpdir / f"npl_{npl}"
        _write(example_data_path, ws, npl, model)
        chars, values, fname = _longest_line(ws)
        print(
            f"\n{name} npl={npl}: longest line {chars} chars, {values} values, {fname}"
        )
        if npl is None:
            assert values > nodelay
        else:
            assert values <= nodelay
        success, buff = run_model(
            USGT_EXE, namefile, model_ws=str(ws), silent=True, report=True
        )
        assert success and any("normal termination" in b.lower() for b in buff), (
            _usgt_run_report(
                f"{name}, free_format_npl={npl}: no normal termination",
                ws,
                success,
                buff,
                _nam_listing(ws, namefile),
            )
        )
        results[npl] = {out: (ws / out).read_bytes() for out in outputs}
    for out in outputs:
        assert len(results[None][out]) > 0, out
        assert results[None][out] == results[10][out], out
