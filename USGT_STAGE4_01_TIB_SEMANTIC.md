# Stage 4.1 - TIB Semantic Support

Goal: replace `MfUsgTib` raw/text-only support with semantic TIB authoring,
load, write, and tests, while preserving a safe raw round-trip fallback if
needed. The most important outcome is that users can create and write a TIB
package from Python data without first loading an existing `.tib` file.

## Current State

`MfUsgTib` currently:

- stores `raw_body` from loaded files,
- writes that text back,
- can accept hand-built raw `blocks`,
- does not parse TIB stress-period data into semantic arrays/lists,
- does not offer a clean from-scratch constructor for project workflows.

Roadmap status: `Raw/text round-trip`, not authoring-ready.

## Fortran Source

Primary file:

`glo2basu1.f`

Primary subroutine:

`GWF2TIB1RP`

The agent must audit the Fortran before implementation. Do not rely only on the
existing docstring.

## Expected TIB Concepts To Model

From the current notes, each stress period contains:

- a header: `NIB0 NIB1 NIBM1 [NICB0 NICB1 NICBM1]`,
- `NIB0` node numbers read by `U1DINT`,
- `NIB1` records for inactive-to-active or active boundary behavior,
- `NIBM1` records for active-to-inactive or related mode,
- optional transport concentration blocks when BCT/transport is active:
  `NICB0`, `NICB1`, `NICBM1`.

The agent must verify exact meaning, record order, reuse behavior, and whether
nodes are GWF nodes, CLN nodes, or both.

## Design Requirements

- Internal nodes must be 0-based.
- File I/O must be 1-based.
- Preserve USG-T valid `U1DINT` behavior for node lists.
- Support from-scratch Python authoring without requiring an existing TIB file;
  this is the primary acceptance criterion.
- Provide clear constructor arguments and/or helper data classes/functions that
  write valid TIB stress-period records from structured Python inputs.
- Preserve existing round-trip behavior for files that fit the parser.
- If an uncommon syntax is not supported, fail explicitly or preserve via a
  documented raw fallback; do not silently write partial data.
- Do not break the existing raw round-trip tests.

## Suggested API Shape

The agent may revise this after reading the Fortran, but should aim for a clear
semantic model such as:

```python
MfUsgTib(
    model,
    stress_period_data={
        0: {
            "ib0_nodes": [...],
            "ib1_records": [...],
            "ibm1_records": [...],
            "conc_ib0": ...,
            "conc_ib1": ...,
            "conc_ibm1": ...,
        }
    },
)
```

If a different structure better matches Fortran, document it in the class
docstring and tests.

## Required Tests

Add focused tests in `autotest/test_usg_transport.py`:

- from-scratch authoring of a minimal non-transport TIB package, written without
  loading any `.tib` file first,
- from-scratch authoring of a transport-enabled TIB package with concentration
  blocks, written without loading any `.tib` file first,
- semantic load of a minimal non-transport TIB file,
- semantic write and reload from scratch,
- `U1DINT` node list with multiple nodes per line,
- 0-based internal nodes and 1-based file output,
- raw/text fallback still round-trips the existing fixture,
- unsupported syntax fails explicitly or remains raw by documented design.

If possible, include one tiny executable smoke case using USG-T 2.7 only if it is
stable and cheap. Otherwise keep executable validation out of this card.

## Documentation Updates

Update:

- `USGT_roadmap.md`
- `USGT_PACKAGE_BACKLOG.md`
- `USGT_improvements.md`

Promote TIB only if semantic constructor + parser + writer + tests are complete.
Otherwise change status from raw/text to partial with exact unsupported modes.

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
Goal: implement from-scratch semantic TIB authoring/writing for USG-T 2.7 in FloPy, while preserving safe round-trip behavior.

Read first:
- USGT_STAGE4_MASTER_PLAN.md
- USGT_STAGE4_01_TIB_SEMANTIC.md
- USGT_roadmap.md
- USGT_PACKAGE_BACKLOG.md
- flopy/mfusg/mfusgtib.py
- autotest/test_usg_transport.py
- /Users/roh.sgm/Documents/SGM Co/GDM/01_Software/GSI/USG-Transport 2.7.0 (2026.03.21)/USGT_V_2-7-0_Source_Code/glo2basu1.f

Task:
Audit `GWF2TIB1RP` in the Fortran and implement semantic `MfUsgTib` from-scratch authoring and valid writing for TIB stress-period records. Loading existing files and round-tripping are required, but the central requirement is that a user can build a TIB package directly from Python/numpy inputs and write a valid USG-T 2.7 `.tib` file. Preserve 0-based internal node indexing and 1-based file I/O. Keep the current raw/text round-trip path only as a documented fallback if needed; do not silently write incomplete semantic data.

Required tests:
- from-scratch non-transport authoring without loading a `.tib` first,
- from-scratch transport-enabled authoring with concentration blocks,
- semantic non-transport load/write/reload,
- U1DINT node lists with multiple nodes per line,
- transport concentration blocks,
- 0-based internal / 1-based file indexing,
- existing raw round-trip behavior remains protected.

Update docs:
- USGT_roadmap.md
- USGT_PACKAGE_BACKLOG.md
- USGT_improvements.md

Validation:
python -m pytest autotest/test_usg_transport.py -q
python -m pytest autotest/test_usg_transport_exe.py -q
USGT_EXE=/Users/roh.sgm/Documents/GitHub_Projects/GW-Software-Compiled-via-Claude/gfortran/usgt_2.7/usgt_270_arm python -m pytest autotest/test_usg_transport.py autotest/test_usg_transport_exe.py -q
git diff --check
git status --short

Constraints:
- Push only to origin, never upstream.
- Do not touch MF6-TID or unrelated packages.
- Do not promote TIB to Full unless semantic parser, writer, authoring tests, and round-trip tests are all present.
```

## Stage 4.1 review follow-up — resolved

Two P2 findings from the Stage 4.1 review were fixed (TIB status unchanged,
still `Full (authoring)`):

1. **Ambiguous input modes.** `MfUsgTib.__init__` now treats
   `stress_period_data`, `blocks`, and `raw_body` as mutually exclusive: zero
   non-empty modes is allowed (no-op all-zero TIB), but two or more raise
   `ValueError`, so `write_file` never silently prefers `raw_body`.
2. **Premature EOF under `parse=True`.** `_parse_semantic` now raises if the
   file ends before all `nper` stress-period headers are read (instead of
   returning partial data). `load(parse=True)` catches it and falls back to the
   raw round-trip, so a short file is never rewritten with synthetic no-op
   stress periods.

Tests added in `autotest/test_usg_transport.py`:
`test_mfusgtib_rejects_mixed_input_modes`,
`test_mfusgtib_parse_rejects_truncated_file`.

Validation: focused **109 passed**, exe **3 passed**, combined **112 passed**
under the USG-T 2.7 ARM binary; `git diff --check` clean.

### Polish — explicit-empty modes closed

The first fix counted modes with a truthiness test (`if val`) while `write_file`
branches on `is not None`. That let an explicitly-empty mode slip past
validation and still steer the writer into a silent branch — e.g.
`stress_period_data={..}, raw_body=""` (writer would take the empty raw branch)
or `stress_period_data={}, blocks={..}` (writer would take the empty semantic
branch). Validation now counts a mode as supplied by **explicit presence**
(`is not None`), matching `write_file`, so a single explicitly-empty mode is
valid but any mix raises `ValueError`. Regression cases added to
`test_mfusgtib_rejects_mixed_input_modes`. TIB stays `Full (authoring)`;
focused **109 passed**, exe **3 passed**, combined **112 passed** under the ARM
binary.
