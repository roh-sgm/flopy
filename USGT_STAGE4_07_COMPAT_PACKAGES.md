# Stage 4.7 - Compatibility-Only Packages

Goal: decide whether SFR, STR, GAGE, FHB, SUB, and SWT should remain
compatibility-only or receive USG-T-specific semantic implementations.

## Current State

These packages use base FloPy classes or compatibility behavior:

- `SFR`
- `STR`
- `GAGE`
- `FHB`
- `SUB`
- `SWT`

Current project preference is CLN-based coupling, so these are lower priority.

## Outcome — Stage 4.6A (STR/SUB/SWT guard; SFR kept) — executed

The final gap audit (`USGT_STAGE4_09_FINAL_GAP_AUDIT.md`) flagged these six
base-class packages as the only place a USG-T model could be written silently
wrong. Stage 4.6A closes that with a **targeted guard**:

- **STR, SUB, SWT — guarded.** `flopy/mfusg/mfusgcompat.py` defines thin
  subclasses (`MfUsgStr`/`MfUsgSub`/`MfUsgSwt`) registered in
  `MfUsg.mfnam_packages` in place of the base classes. On an **unstructured**
  `MfUsg` model (`version == "mfusg"` and `not structured`) both **load** and
  **authoring/`__init__`** raise `NotImplementedError` before any file is read or
  written (no partial file). On a *structured* `MfUsg` model they delegate to the
  base class, and plain `flopy.modflow.Modflow` use is unaffected (it keeps the
  base classes). Message: compatibility-only; USG-T unstructured/DISU not
  implemented; use CLN / implement an MfUsg-specific package / keep as raw
  external.
- **SFR — NOT guarded (kept as compatibility/base support).** The original
  "guard all four" plan was too aggressive: `examples/data/freyberg_usg` is a
  DISU model that uses SFR and is loaded → written → **run** via base
  `ModflowSfr2` in `autotest/test_usg.py` (`test_freyburg_usg`/`_external`). So
  unstructured SFR is empirically supported by the base class; guarding it would
  break passing tests and contradict the repo. SFR is documented as
  compatibility/base support, **DISU-validated via `freyberg_usg`** — *not* a
  Full USG-T transport implementation.
- **FHB, GAGE — unchanged.** They round-trip via their base classes in the `Ex8`
  real model; left as compatibility.

Tests: `test_mfusg_compat_str_sub_swt_load_guard`,
`test_mfusg_compat_str_sub_swt_write_guard`,
`test_mfusg_compat_guard_noop_and_sfr_untouched` (`-k compat`, 3 passed); the
`freyberg_usg` SFR tests in `test_usg.py` still pass; transport focused **241**,
combined **245** (ARM).

## Work Order

Audit one package at a time:

1. SFR
2. STR
3. GAGE
4. FHB
5. SUB
6. SWT

For each:

- identify whether USG-T 2.7 has format differences from base MODFLOW,
- check whether any target/example model uses it,
- decide: implement, raw/text support, or compatibility-only,
- add tests if implementation changes.

## Promotion Criteria

Do not promote a package beyond compatibility-only unless it has:

- Fortran-derived spec,
- semantic load/write,
- from-scratch authoring test,
- round-trip test,
- explicit unsupported-mode behavior.

## Agent Prompt

```text
Goal: audit one compatibility-only package for USG-T-specific support, starting with SFR.

Read first:
- USGT_STAGE4_MASTER_PLAN.md
- USGT_STAGE4_07_COMPAT_PACKAGES.md
- USGT_roadmap.md
- the relevant base FloPy package class
- the corresponding USG-T 2.7 Fortran source

Task:
Audit SFR first. Determine whether USG-T 2.7 input differs from the base FloPy class enough to require `MfUsgSfr` work. If yes, propose/implement a minimal semantic package with tests. If no, strengthen documentation and tests justifying compatibility-only.

Validation:
python -m pytest autotest/test_usg_transport.py -q
python -m pytest autotest/test_usg_transport_exe.py -q
USGT_EXE=/Users/roh.sgm/Documents/GitHub_Projects/GW-Software-Compiled-via-Claude/gfortran/usgt_2.7/usgt_270_arm python -m pytest autotest/test_usg_transport.py autotest/test_usg_transport_exe.py -q
git diff --check
git status --short

Constraints:
- Push only to origin, never upstream.
- Do one package per card unless the result is documentation-only.
```
