# USG-Transport improvements on this fork

This fork of [modflowpy/flopy](https://github.com/modflowpy/flopy) carries
USG-Transport additions on `develop`, aimed at improving how
`MfUsg.load(...).write_input()` round-trips existing USG-T models and
extending post-processing utilities for transport output. Nothing else is
changed.

Tested against **USG-T 1.8** and **USG-T 2.7** binaries on both structured and
unstructured grids.

These are incremental improvements on top of the USG-T work already in
upstream flopy. They live here while testing continues.

## What's added

### Round-trip fixes (original five)

| Branch | What it does |
|---|---|
| `usgt/rch-transport-fix` | In `MfUsgRch.write_file`, `"# Stress period {kper + 1}"` becomes an f-string so the comment substitutes per SP, and `INRECH` is followed by `INIRCH` only when `NRCHOP == 2` (matching the read path). Both details were harmless at runtime but left the written RCH file not round-trippable. |
| `usgt/cln-load-none-unit` | In `MfUsgCln.__init__`, treat `None` entries in `unitnumber` as 0 instead of calling `int(None)`. This comes up when a CLN-declared output unit is not declared in the NAM's `ext_unit_dict`. |
| `usgt/tib-package` | Adds a new `MfUsgTib` class for the Transient Ibound package. Minimal text-preserving round-trip — enough to load + write an existing TIB. |
| `usgt/bas-preserve-unstructured` | `MfUsgBas.write_file` re-emits the `UNSTRUCTURED` keyword when `parent.structured is False`. The load side already reads the token; adding it to write closes the round-trip. |
| `usgt/nam-rebase-output-paths` | `BaseModel._reset_external` stores the basename of output files on `change_model_ws`. Previously the absolute paths from the original workspace were preserved, so a relocated model wrote its binary output back into the source folder. On one real test this inflated runtime roughly 15× because of remote-folder I/O. |

### New packages and utilities

| What | File | What it does |
|---|---|---|
| `MfUsgChd` | `flopy/mfusg/mfusgchd.py` | CHD package for unstructured USG-T grids. Node-based (replaces k/i/j), supports AUX concentration variables. Full `load` and `write_file` for the USG-T format (`NACT    Stress Period N` headers, `-1` reuse). |
| `MfUsgRiv` | `flopy/mfusg/mfusgriv.py` | RIV package for unstructured USG-T grids. Supports AUX concentration and a trailing reach-ID column (`irch`) that is written positionally without being declared as AUX. `irch` is auto-detected from the first data row when loading. |
| `MfUsgEts` | `flopy/mfusg/mfusgets.py` | Segmented Evapotranspiration (ETS) package for USG-T 2.7. Distinct from the simpler EVT: supports `NETSEG > 1` with per-SP `PXDP`/`PETM` segment arrays, the `IESFACTOR` transport flag, and named parameters. Full `load` and `write_file`. |
| `MfusgTransportListBudget` | `flopy/utils/mflistfile.py` | Reads transport species budget from a USG-T listing file for a single species. Handles both **old** USG-T format (transport blocks use `VOLUMETRIC BUDGET`, same keyword as flow) and **new** format (transport blocks use `MASS BUDGET`). Instantiate once per species: `MfusgTransportListBudget("model.lst", species=2)`. Returns the same recarrays / DataFrames as `MfusgListBudget`. |

### Minor fixes

| File | Fix |
|---|---|
| `flopy/mfusg/mfusgbct.py` | Removed a stray `print()` left from development. Fixed file-handle management in `write_file`: the file is now closed only when opened internally (`close_on_exit` flag), so callers that pass an open handle are not surprised. |

## Install

```bash
conda create -n usgt -y python=3.12 numpy pandas matplotlib jupyter shapely
conda activate usgt
git clone --branch develop https://github.com/roh-sgm/flopy.git
cd flopy
pip install -e .
```

## Quick check

```python
import flopy
print(flopy.__version__)                              # 3.11.0.dev0 (or newer)
print('MfUsgTib' in dir(flopy.mfusg))                # True
print('MfUsgChd' in dir(flopy.mfusg))                # True
print('MfUsgRiv' in dir(flopy.mfusg))                # True
print('MfUsgEts' in dir(flopy.mfusg))                # True
print('MfusgTransportListBudget' in dir(flopy.utils)) # True

# Multi-species transport budget reader
from flopy.utils import MfusgTransportListBudget
lst_s1 = MfusgTransportListBudget("model.lst", species=1)
lst_s2 = MfusgTransportListBudget("model.lst", species=2)
inc, cum = lst_s1.get_budget()
df_inc, df_cum = lst_s2.get_dataframes(start_datetime="2000-01-01")
```

All other flopy APIs behave exactly as upstream; the patches only change the
text that `MfUsg.load(...).write_input()` produces for USG-T models.

## Validation

Exercised end-to-end against an unstructured USG-T model with
BCT + CLN + TIB transport (≈ 112 k nodes, 493 stress periods). Running the
flopy-rewritten model through the same USG-T binary as the reference run
produces a listing file that matches bit-for-bit on every budget column —
flow and transport — (`max|diff| = 0` on all stress periods, across 20 flow
and 16 mass components). Wall-clock runtime is equivalent to the reference.

Tested binaries:
- USG-T 2.7.0 (ARM) — round-trip sanity on a small reference model
- USG-T 1.8 (ARM and x86) — full validation on the transport model above

## Upstream

Held on the fork while testing continues. No PRs to `modflowpy/flopy` have
been opened yet.
