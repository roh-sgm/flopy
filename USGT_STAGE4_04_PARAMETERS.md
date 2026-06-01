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
- `HFB`: **parameter-preserving (Stage 4.4B, executed)** — list parameters
  (`UPARLSTRP`/`UPARLSTSUB`) now load → write → reload with their syntax intact
  (`NPHFB>0`, definition blocks with `NLST` barrier rows, `NACTHFB` + active
  names). Barrier lists (parameter `NLST` rows and non-parametric `NHFBNP` rows)
  also support leading `SFAC`/`OPEN/CLOSE`/`EXTERNAL` list controls via the shared
  `_usgt_list.begin_list_block` (list-control follow-up). From-scratch parameter
  authoring and `TRANSIENT_HFB`+`NPHFB>0` raise `NotImplementedError`;
  `INSTANCES` are unsupported (matches the Fortran). See
  `USGT_STAGE4_04_PARAMETERS_HFB.md`.
- `SGB`: `NPSGB>0` fails explicitly.
- `QRT`: `NPQRT>0` fails explicitly.
- `DRT`: `NPDRT>0` fails explicitly.

This is honest and safe. ETS (array) and HFB (list) parameter preservation are
implemented; the remaining list-parameter packages (SGB/DRT/QRT) stay
explicit-fail pending their own preservation.

## Shared abstraction (Stage 4.4A/4.4B)

`flopy/mfusg/_usgt_parameters.py` holds two parameter-preservation paths:

- **Array parameters** (`UPARARRRP`/`UPARARRSUB2`, ETS): reuses the already-shared
  `flopy.modflow.ModflowParBc` *parser* and adds the missing **write** side —
  `write_array_parameter_defs`, `read_active_array_parameters`,
  `write_active_array_parameters`. The parsed `ModflowParBc.bc_parms` dict is the
  structured representation preserved on the package.
- **List parameters** (`UPARLSTRP`/`UPARLSTSUB`, HFB) — Stage 4.4B:
  `read_list_parameter_header` / `write_list_parameter_header` (the
  `PARNAM PARTYP PARVAL NLST [INSTANCES n]` header) and
  `read_active_list_parameters` / `write_active_list_parameters` (the activation
  names). The per-parameter `NLST` rows are package-specific, so the row
  reader/writer stays in the package (HFB keeps `_read_hfb_rows`/`_write_hfb_rows`
  with 0-based↔1-based conversion).

No ad-hoc string preservation, no second parser. The list-parameter helpers are
reusable by SGB/DRT/QRT when their preservation lands; until then those packages
keep their explicit `NotImplementedError`.

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
2. HFB named parameters. **DONE (Stage 4.4B).**
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
**preserves** ETSR array parameters (Stage 4.4A) and HFB now **preserves** list
parameters (Stage 4.4B); both are documented as parameter-preserving for those
paths but stay not `Full` because from-scratch parameter authoring is
unsupported (and, for HFB, `TRANSIENT_HFB`+`NPHFB>0` is unsupported). The
remaining list-parameter packages (SGB/DRT/QRT) stay explicit-fail until their
write path lands.

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
