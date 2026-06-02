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
- `SGB`: **parameter-definition-preserving only (Stage 4.4C + review
  follow-up)** — the `PARAMETER NPSGB MXS` record and the per-parameter
  definitions (`UPARLSTRP` header + `NLST` `NODE GRADIENT [aux]` rows) load →
  write → reload with their syntax intact; SFAC inert on the gradient (Fortran
  ISCLOC=2), AUX preserved, `ITMP<0` reuse kept. **Active SGB parameters are
  unsupported**: USG-T 2.7 defines them as `PARTYP='SGB'` (`UPARLSTRP`,
  glo2sgbu1.f:97) but activates them as `PTYP='G'` (`UPARLSTSUB`,
  glo2sgbu1.f:185), a type conflict (`parutl7.f:684/800`) that aborts the run.
  So `NP>0` and `INSTANCES` raise `NotImplementedError`; the writer validates
  `MXS` and inconsistent state raises `ValueError` (no partial file). See
  `USGT_STAGE4_04_PARAMETERS_SGB.md`.
- `DRT`: **parameter-preserving (Stage 4.4D, executed)** — list parameters
  (`UPARLSTRP`/`SGWF2DRT8LS`) load → write → reload with their syntax intact
  (`NPDRT`/`MXL` in item 1, definitions with `NLST` drain rows **and their
  RETURNFLOW recipients / spreading blocks**, per-SP `ITMP NP` + active names).
  DRT is internally consistent — definition and activation both use
  `PARTYP='DRT'` (`gwf2drt8u.f:135` / `:1108`), so active parameters are valid
  (unlike SGB). Parameter value scales `COND` (Fortran `IPVL=5`); AUX/CHANGEC
  preserved; `ITMP<0` reuse kept. `INSTANCES` and from-scratch authoring raise
  `NotImplementedError`; `MXL`/consistency validated (`ValueError`, no partial
  file). **Review follow-up:** `MXADRT` is sized to the active total per period
  (non-parametric + active-parameter rows; the Fortran aborts if `NDRTCL >
  MXADRT`); and an activated SPREAD (`NR<0`) parameter is structural-round-trip
  only, not execution-guaranteed (USG-T copies `DRTF` but not `NodDRT` on
  activation). See `USGT_STAGE4_04_PARAMETERS_DRT.md`.
- `QRT`: `NPQRT>0` fails explicitly.

This is honest and safe. ETS (array) and HFB + SGB + DRT (list) parameter
preservation are implemented; the remaining list-parameter package (QRT) stays
explicit-fail pending its own preservation (and a `PARTYP` audit like SGB/DRT).

## Shared abstraction (Stage 4.4A/4.4B/4.4C)

`flopy/mfusg/_usgt_parameters.py` holds two parameter-preservation paths:

- **Array parameters** (`UPARARRRP`/`UPARARRSUB2`, ETS): reuses the already-shared
  `flopy.modflow.ModflowParBc` *parser* and adds the missing **write** side —
  `write_array_parameter_defs`, `read_active_array_parameters`,
  `write_active_array_parameters`. The parsed `ModflowParBc.bc_parms` dict is the
  structured representation preserved on the package.
- **List parameters** (`UPARLSTAL`/`UPARLSTRP`/`UPARLSTSUB`, HFB 4.4B + SGB 4.4C):
  `read_list_parameter_count` / `write_list_parameter_count` (the leading
  `PARAMETER NP MXL` record, for packages that carry it — SGB/DRT/QRT);
  `read_list_parameter_header` / `write_list_parameter_header` (the
  `PARNAM PARTYP PARVAL NLST [INSTANCES n]` definition header); and
  `read_active_list_parameters` / `write_active_list_parameters` (the activation
  names). The per-parameter `NLST` rows are package-specific, so the row
  reader/writer stays in the package (HFB `_read_hfb_rows`/`_write_hfb_rows`; SGB
  `_read_sgb_rows`/`_write_sgb_rows`; DRT `_read_drain_rows`/`_write_drain_line`,
  the last carrying per-row RETURNFLOW recipients), each with 0-based↔1-based
  conversion.

No ad-hoc string preservation, no second parser. The list-parameter helpers are
reusable by QRT when its preservation lands; until then it keeps its explicit
`NotImplementedError`.

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
3. SGB list parameters. **DONE (Stage 4.4C; definition-preserving only —
   active SGB params abort the Fortran).**
4. DRT list parameters. **DONE (Stage 4.4D; full active-parameter preservation —
   DRT is consistent).**
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
**preserves** ETSR array parameters (Stage 4.4A); HFB **preserves** list
parameters (Stage 4.4B); SGB is **definition-preserving only** (Stage 4.4C +
review follow-up — USG-T 2.7 rejects an SGB activation with a `PARTYP='SGB'` vs
`'G'` type conflict); DRT **preserves** list parameters including activations
and recipients (Stage 4.4D — DRT is internally consistent, both `'DRT'`). All
stay not `Full` because from-scratch parameter authoring is unsupported (and,
per package: HFB `TRANSIENT_HFB`+`NPHFB>0`; SGB active parameters; DRT/HFB/SGB
`INSTANCES`). The remaining list-parameter package (QRT) stays explicit-fail
until its write path lands (and a `PARTYP` audit like SGB/DRT).

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
