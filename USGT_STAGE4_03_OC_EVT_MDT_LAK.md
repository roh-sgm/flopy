# Stage 4.3 - OC / EVT / MDT / LAK Fullness Pass

Goal: promote the remaining project-relevant `OK but not Full` packages only
where Fortran audit, semantic authoring, round-trip tests, and unsupported-mode
behavior justify it.

## Package Order

1. `OC`
2. `EVT`
3. `MDT`
4. `LAK`

This order keeps the smaller output/control packages ahead of large LAK.

## Current State

### OC

Status: implemented and tested, not Full.

Known gap:

- full output-block round-trip and authoring coverage is not exhaustive,
- `BOOTSTRAPPING` is parsed on load but not fully protected by from-scratch
  authoring tests.

### EVT

Status: implemented and tested, not Full.

Known gap:

- transport `IETFACTOR`/`ETFACTOR` authoring exists,
- full per-`NEVTOP` array authoring combinations are not exhaustively tested.

### MDT

Status: implemented and round-trip validated via Ex7, not Full.

Known gap:

- field order is not independently Fortran-audited,
- no from-scratch authoring test.

### LAK

Status: implemented enough for Ex8 load/write, not Full.

Known gap:

- from-scratch authoring of `TABLEINPUT`,
- `TRANSPORTBOUNDARY`,
- complete lake connectivity,
- bathymetry/stage tables,
- transport boundary concentrations.

## Work Plan

Treat each package as its own mini-card. Do not try to promote all four in one
large commit unless the changes are documentation-only.

### Card A - OC Fullness

Required:

- audit `mfusgoc.py` against USG-T OC behavior,
- add from-scratch `BOOTSTRAPPING` test if supported,
- add output-block coverage for common save/print concentration/head/budget
  combinations,
- decide whether compact/non-compact or rare modes are in scope.

### Card B - EVT Fullness

Required:

- audit `gwf2evt8u1.f`,
- add `NEVTOP` combination tests,
- add authoring tests with and without transport,
- verify 0-based/1-based behavior where node/layer indices apply,
- protect `IETFACTOR`/`ETFACTOR` reload.

### Card C - MDT Fullness

Required:

- audit `gwt2mdtu1.for`,
- document exact dataset order and conditional reads,
- add from-scratch MDT authoring test,
- reload written MDT and compare all semantic fields,
- include multi-species if Fortran supports/uses it.

### Card D - LAK Fullness

Required:

- audit `gwf2lak7u1.f`,
- define a minimal valid synthetic lake model,
- add from-scratch authoring for `TABLEINPUT` and `TRANSPORTBOUNDARY`,
- test load/write/reload,
- only add executable validation if the synthetic LAK model is stable and cheap.

## Documentation Updates

After each card, update:

- `USGT_roadmap.md`,
- `USGT_PACKAGE_BACKLOG.md`,
- `USGT_improvements.md`.

Only promote a package to `Full` after its own card is complete.

## Validation

```bash
python -m pytest autotest/test_usg_transport.py -q
python -m pytest autotest/test_usg_transport_exe.py -q
USGT_EXE=/Users/roh.sgm/Documents/GitHub_Projects/GW-Software-Compiled-via-Claude/gfortran/usgt_2.7/usgt_270_arm python -m pytest autotest/test_usg_transport.py autotest/test_usg_transport_exe.py -q
git diff --check
git status --short
```

## First Agent Prompt For This Theme

```text
Goal: complete the OC package fullness pass for USG-T 2.7 without touching EVT/MDT/LAK yet.

Read first:
- USGT_STAGE4_MASTER_PLAN.md
- USGT_STAGE4_03_OC_EVT_MDT_LAK.md
- USGT_roadmap.md
- flopy/mfusg/mfusgoc.py
- autotest/test_usg_transport.py
- USG-T 2.7 Fortran source files relevant to OC/output control.

Task:
Audit OC behavior, then add focused semantic authoring/round-trip tests for the missing OC modes, especially BOOTSTRAPPING and output-block combinations. Promote OC only if the tests justify it; otherwise keep the status honest and document exactly what remains.

Validation:
python -m pytest autotest/test_usg_transport.py -q
python -m pytest autotest/test_usg_transport_exe.py -q
USGT_EXE=/Users/roh.sgm/Documents/GitHub_Projects/GW-Software-Compiled-via-Claude/gfortran/usgt_2.7/usgt_270_arm python -m pytest autotest/test_usg_transport.py autotest/test_usg_transport_exe.py -q
git diff --check
git status --short

Constraints:
- Push only to origin, never upstream.
- Do not touch MF6-TID.
- Do not bundle OC, EVT, MDT, and LAK changes in one commit unless docs only.
```
