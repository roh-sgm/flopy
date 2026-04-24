# USG-Transport improvements on this fork

This fork of [modflowpy/flopy](https://github.com/modflowpy/flopy) carries
five small USG-Transport additions on `develop`, aimed at improving how
`MfUsg.load(...).write_input()` round-trips existing USG-T models. Nothing
else is changed.

Tested against **USG-T 1.8** and **USG-T 2.7** binaries on both structured and
unstructured grids.

These are incremental improvements on top of the USG-T work already in
upstream flopy. They live here while testing continues.

## What's added

| Branch | What it does |
|---|---|
| `usgt/rch-transport-fix` | In `MfUsgRch.write_file`, `"# Stress period {kper + 1}"` becomes an f-string so the comment substitutes per SP, and `INRECH` is followed by `INIRCH` only when `NRCHOP == 2` (matching the read path). Both details were harmless at runtime but left the written RCH file not round-trippable. |
| `usgt/cln-load-none-unit` | In `MfUsgCln.__init__`, treat `None` entries in `unitnumber` as 0 instead of calling `int(None)`. This comes up when a CLN-declared output unit is not declared in the NAM's `ext_unit_dict`. |
| `usgt/tib-package` | Adds a new `MfUsgTib` class for the Transient Ibound package. Minimal text-preserving round-trip — enough to load + write an existing TIB. |
| `usgt/bas-preserve-unstructured` | `MfUsgBas.write_file` re-emits the `UNSTRUCTURED` keyword when `parent.structured is False`. The load side already reads the token; adding it to write closes the round-trip. |
| `usgt/nam-rebase-output-paths` | `BaseModel._reset_external` stores the basename of output files on `change_model_ws`. Previously the absolute paths from the original workspace were preserved, so a relocated model wrote its binary output back into the source folder. On one real test this inflated runtime roughly 15× because of remote-folder I/O. |

All five are merged onto `develop`. The feature branches remain for
cherry-picking or possible upstream PRs later.

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
print(flopy.__version__)                   # 3.11.0.dev0 (or newer)
print('MfUsgTib' in dir(flopy.mfusg))      # True
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
