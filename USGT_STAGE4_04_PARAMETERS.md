# Stage 4.4 - MODFLOW Parameter Preservation

Goal: implement or formally close full MODFLOW-style parameter preservation for
USG-T packages that currently expand parameters or fail explicitly.

## Current State

- `ETS`: parameterized input loads and writes expanded non-parametric arrays with
  `NPETS=0`.
- `HFB`: `NPHFB>0` fails explicitly.
- `SGB`: `NPSGB>0` fails explicitly.
- `QRT`: `NPQRT>0` fails explicitly.
- `DRT`: `NPDRT>0` fails explicitly.

This is honest and safe, but not literal full support.

## Fortran Sources

Primary shared source:

- `parutl7.f`

Package sources:

- `gwf2ets8u1.f`
- `gwf2hfb7u1.f`
- `glo2sgbu1.f`
- `gwf2QRT8u.f`
- `gwf2drt8u.f`

## Design Requirements

- Build one shared design for parameters before package-specific edits.
- Avoid ad hoc string preservation if structured parameter objects are feasible.
- Preserve load/write behavior for:
  - parameter definitions,
  - active parameter instances per stress period,
  - cluster records,
  - parameter values and multipliers,
  - non-parameterized records mixed with parameterized ones.
- If a package writes expanded arrays by design, keep the label `Expanded valid
  write`; do not call it Full.

## Recommended Implementation Order

1. ETS array parameters, because partial expansion already exists.
2. HFB named parameters.
3. SGB list parameters.
4. DRT list parameters.
5. QRT list parameters.

## Required Tests

For each package:

- parameterized load,
- write preserving parameter syntax,
- reload preserved file,
- active parameter record with at least one non-parameterized record where
  Fortran allows it,
- explicit failure for unsupported parameter variants if any remain.

## Documentation Updates

Update:

- `USGT_roadmap.md`,
- `USGT_PACKAGE_BACKLOG.md`,
- `USGT_improvements.md`.

Promote only the package modes that preserve parameter syntax. If ETS continues
to expand, leave it as `Expanded valid write`.

## Validation

```bash
python -m pytest autotest/test_usg_transport.py -q
python -m pytest autotest/test_usg_transport_exe.py -q
USGT_EXE=/Users/roh.sgm/Documents/GitHub_Projects/GW-Software-Compiled-via-Claude/gfortran/usgt_2.7/usgt_270_arm python -m pytest autotest/test_usg_transport.py autotest/test_usg_transport_exe.py -q
git diff --check
git status --short
```

## Agent Prompt

```text
Goal: design and begin implementing MODFLOW-style parameter preservation for USG-T packages, starting with ETS.

Read first:
- USGT_STAGE4_MASTER_PLAN.md
- USGT_STAGE4_04_PARAMETERS.md
- USGT_roadmap.md
- flopy/mfusg/mfusgets.py
- flopy/mfusg/mfusghfb.py
- flopy/mfusg/mfusgsgb.py
- flopy/mfusg/mfusgdrt.py
- flopy/mfusg/mfusgqrt.py
- /Users/roh.sgm/Documents/SGM Co/GDM/01_Software/GSI/USG-Transport 2.7.0 (2026.03.21)/USGT_V_2-7-0_Source_Code/parutl7.f

Task:
Audit the shared parameter machinery in the Fortran and propose/implement a shared FloPy design. Start with ETS parameter preservation unless the audit shows a better first package. Add tests that prove parameter syntax is preserved on write and reload. Do not degrade the existing expanded-valid-write behavior unless replacing it with true preservation.

Validation:
python -m pytest autotest/test_usg_transport.py -q
python -m pytest autotest/test_usg_transport_exe.py -q
USGT_EXE=/Users/roh.sgm/Documents/GitHub_Projects/GW-Software-Compiled-via-Claude/gfortran/usgt_2.7/usgt_270_arm python -m pytest autotest/test_usg_transport.py autotest/test_usg_transport_exe.py -q
git diff --check
git status --short

Constraints:
- Push only to origin, never upstream.
- Do not call any parameterized package Full unless preservation is implemented and tested.
- Prefer a shared parameter abstraction over five unrelated parsers.
```
