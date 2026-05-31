# Stage 4.5 - QRT/DRT Rare Controls

Goal: close the remaining rare-control gaps in QRT and DRT: QRT `TRANSIENTQ`
and recipient-node `U1DINT` `EXTERNAL`/`OPEN/CLOSE` controls.

## Current State

QRT and DRT are strong for common authoring:

- node-based records,
- return-flow recipients,
- AUX concentrations,
- `CHANGEC`,
- `ITMP/-1` reuse,
- main-list `SFAC`/`OPEN-CLOSE`/`EXTERNAL` load with expanded valid write.

Remaining gaps:

- `MfUsgQrt` rejects `TRANSIENTQ`.
- `_usgt_returnflow` rejects recipient `EXTERNAL` / `OPEN/CLOSE`.
- Parameter preservation is handled separately in
  `USGT_STAGE4_04_PARAMETERS.md`.

## Fortran Sources

- `gwf2QRT8u.f`
- `gwf2drt8u.f`
- `utl7u1.f` for `U1DINT`

## Required Work

1. Audit `TRANSIENTQ` in QRT:
   - records,
   - time-series semantics,
   - stress-period interaction,
   - write format.
2. Audit recipient-node `U1DINT` external controls:
   - `INTERNAL`,
   - `CONSTANT`,
   - `EXTERNAL`,
   - `OPEN/CLOSE`.
3. Implement support where reasonable.
4. If preserving external recipient files is too much, at least load and expand
   them with an honest `Expanded valid write` label.

## Required Tests

- QRT `TRANSIENTQ` explicit supported case or explicit documented deferral.
- DRT recipient list via `OPEN/CLOSE`.
- QRT recipient list via `EXTERNAL` resolved through `ext_unit_dict`.
- Expanded write/reload if controls are not preserved.
- No regression in existing inline/constant recipient tests.

## Documentation Updates

Update:

- `USGT_roadmap.md`,
- `USGT_PACKAGE_BACKLOG.md`,
- `USGT_improvements.md`.

## Agent Prompt

```text
Goal: implement the remaining QRT/DRT rare list controls for USG-T 2.7, starting with recipient U1DINT EXTERNAL/OPEN/CLOSE.

Read first:
- USGT_STAGE4_MASTER_PLAN.md
- USGT_STAGE4_05_QRT_DRT_RARE_CONTROLS.md
- flopy/mfusg/mfusgqrt.py
- flopy/mfusg/mfusgdrt.py
- flopy/mfusg/_usgt_returnflow.py
- /Users/roh.sgm/Documents/SGM Co/GDM/01_Software/GSI/USG-Transport 2.7.0 (2026.03.21)/USGT_V_2-7-0_Source_Code/gwf2QRT8u.f
- /Users/roh.sgm/Documents/SGM Co/GDM/01_Software/GSI/USG-Transport 2.7.0 (2026.03.21)/USGT_V_2-7-0_Source_Code/gwf2drt8u.f

Task:
Audit and implement recipient-node U1DINT EXTERNAL/OPEN/CLOSE support for QRT/DRT. Then audit QRT TRANSIENTQ and either implement it or leave a precise explicit unsupported decision with tests. Prefer load-and-expanded-write if full preservation is too large.

Validation:
python -m pytest autotest/test_usg_transport.py -q
python -m pytest autotest/test_usg_transport_exe.py -q
USGT_EXE=/Users/roh.sgm/Documents/GitHub_Projects/GW-Software-Compiled-via-Claude/gfortran/usgt_2.7/usgt_270_arm python -m pytest autotest/test_usg_transport.py autotest/test_usg_transport_exe.py -q
git diff --check
git status --short

Constraints:
- Push only to origin, never upstream.
- Do not mix parameter preservation into this task.
```
