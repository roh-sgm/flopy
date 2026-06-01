# Stage 4.4 - MODFLOW Parameter Preservation

Goal: implement or formally close full MODFLOW-style parameter preservation for
USG-T packages that currently expand parameters or fail explicitly.

## Current State

- `ETS`: **parameter-preserving (Stage 4.4A, executed)** — array parameters for
  the ETSR rate now load → write → reload with their parameter syntax intact
  (`NPETS>0` in item 2a, definition blocks, per-period activation records,
  including `INSTANCES`). An opt-in expanded fallback (`expand_parameters=True`)
  keeps the old `NPETS=0` behavior; from-scratch parameter *authoring* still
  raises `NotImplementedError`. See `USGT_STAGE4_04_PARAMETERS_ETS.md`.
- `HFB`: `NPHFB>0` fails explicitly.
- `SGB`: `NPSGB>0` fails explicitly.
- `QRT`: `NPQRT>0` fails explicitly.
- `DRT`: `NPDRT>0` fails explicitly.

This is honest and safe. ETS array-parameter preservation is now implemented;
the list-parameter packages remain explicit-fail pending the list-parameter
write path.

## Shared abstraction (Stage 4.4A)

`flopy/mfusg/_usgt_parameters.py` reuses the already-shared
`flopy.modflow.ModflowParBc` *parser* and adds the missing **write** side for
the **array-parameter** form (`UPARARRRP`/`UPARARRSUB2`):
`write_array_parameter_defs`, `read_active_array_parameters`,
`write_active_array_parameters`. The parsed `ModflowParBc.bc_parms` dict is the
structured representation preserved on the package — no ad-hoc string
preservation, no second parser.

The analogous **list-parameter** write path (`UPARLSTRP`/`UPARLSTSUB`, for
SGB/DRT/QRT and HFB-style barriers) is the next package step. It would: parse
the `PARAMETER NP MXL` header, store the per-parameter `NLST` list rows plus
multiplier (`PARVAL`) and instances, and on write emit the definition blocks
plus per-period activation records (scaling rows by the parameter value on
expansion, per `UPARLSTSUB`). Until implemented, those packages keep their
explicit `NotImplementedError`.

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

1. ETS array parameters, because partial expansion already exists. **DONE
   (Stage 4.4A).**
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

Promote only the package modes that preserve parameter syntax. ETS now
**preserves** ETSR array parameters (Stage 4.4A) and is documented as
parameter-preserving for that path; it is still not `Full` because from-scratch
parameter authoring is unsupported. The list-parameter packages stay
explicit-fail until their write path lands.

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
