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
