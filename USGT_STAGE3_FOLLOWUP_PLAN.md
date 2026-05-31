# USG-T Stage 3 Follow-Up Plan

Date: 2026-05-31

Basis: `USGT_STAGE3_REREVIEW.md`

> **Status: COMPLETED — executed by commit `d3725bd8`** (follow-up range
> `6618181c..d3725bd8` on `develop`). All re-review findings are Fixed (see the
> Resolution section of `USGT_STAGE3_REREVIEW.md`). The work order below is
> retained verbatim as the historical plan; it has already been carried out.

Historical context (at the time this plan was written): after fetching
`origin`, `develop` and `origin/develop` both pointed to `6618181c` and no
post-review agent fixes were visible yet. This plan converted the then-open
re-review findings into an explicit work order, which is now done.

Local macOS USG-T 2.7 executable available for validation:

```text
/Users/roh.sgm/Documents/GitHub_Projects/GW-Software-Compiled-via-Claude/gfortran/usgt_2.7/usgt_270_arm
```

Use it through `USGT_EXE`; do not hardcode this absolute path into tests.

## Goal

Close the Stage 3 re-review findings so the branch is internally consistent and
ready to use as the next upstream-readiness baseline.

This pass should not add new USG-T package features. It is a cleanup of
test-contracts, documentation status, and verification wording.

## Scope

Files expected to change:

- `autotest/test_usg_transport.py`
- `autotest/test_usg_transport_exe.py`
- `USGT_roadmap.md`
- `USGT_STAGE3_COMPLETION_PLAN.md`
- `USGT_PACKAGE_BACKLOG.md` if executable-tier wording is still inaccurate
- `USGT_UPSTREAM_INFRA.md` if executable-tier wording is still inaccurate
- `USGT_improvements.md` only if Stage 3 summary wording needs alignment
- `USGT_STAGE3_REREVIEW.md` if the agent records the resolution table there

Do not modify package implementation files unless a test-contract cleanup
reveals a real bug.

## Work Items

### 1. Make executable selection consistent

Problem: `autotest/test_usg_transport_exe.py` supports `USGT_EXE`, but the
real-model `Ex1..Ex9` tests in `autotest/test_usg_transport.py` still hardcode
`mfusg_gsi`.

Tasks:

- Add near the imports in `autotest/test_usg_transport.py`:

```python
USGT_EXE = os.environ.get("USGT_EXE", "mfusg_gsi")
```

- Replace every executable-gated real-model decorator:

```python
@requires_exe("mfusg_gsi")
```

with:

```python
@requires_exe(USGT_EXE)
```

- Replace every real-model load call:

```python
exe_name="mfusg_gsi"
```

with:

```python
exe_name=USGT_EXE
```

Acceptance:

- Default behavior is unchanged when `USGT_EXE` is unset.
- A non-default executable name or absolute path works through `USGT_EXE`.
- Docs no longer overstate behavior.

### 2. Reconcile executable-tier documentation

Problem: docs describe the whole executable tier as `USGT_EXE`-gated, but that
is only true after Work Item 1.

Tasks:

- Review and update:
  - `USGT_PACKAGE_BACKLOG.md`
  - `USGT_UPSTREAM_INFRA.md`
  - `USGT_improvements.md`
  - `autotest/test_usg_transport_exe.py` module docstring
- Make the wording precise:
  - default focused tests are still `autotest/test_usg_transport.py`;
  - real-model `Ex1..Ex9` executable tests and from-scratch executable smoke
    tests both use `USGT_EXE` after this fix;
  - the executable tier remains optional and skips when the executable cannot be
    resolved.

Acceptance:

- There is no doc sentence claiming `USGT_EXE` behavior that the tests do not
  actually implement.

### 3. Convert roadmap audit from planning state to final state

Problem: `USGT_roadmap.md` has a Stage 3 audit table that still reads like a
pre-completion plan.

Tasks:

- Rewrite the Stage 3 audit table as final-state documentation.
- Use final package states from the main coverage table:
  - `DIS`, `DISU`, `CLN`, `DPF`, `PCB` are `✅ Full`.
  - `OC`, `EVT`, `MDT`, `LAK`, `GNC` remain plain `✅` with their reason.
  - `TIB`, `GSF` are `Raw/text round-trip`.
  - `ETS`, `HFB`, `DPT` remain `Partial` / `Expanded valid write` as documented.
  - `SFR`, `STR`, `GAGE`, `FHB`, `SUB`, `SWT` remain `Compatibility-only`.
- Remove or rewrite phrases like "until their Card lands" because Stage 3 cards
  have already landed.

Acceptance:

- A reviewer cannot find contradictory status claims for `OC`, `EVT`, `MDT`,
  or the packages promoted in Stage 3 Card 8.

### 4. Fix Stage 3 baseline and whitespace

Problem: `USGT_STAGE3_COMPLETION_PLAN.md` says Stage 3 is 101 passed at the top
but still says 95 passed in `Current Baseline`. It also has a blank line at EOF.

Tasks:

- Update `Current Baseline` to say the current focused suite is **101 passed**.
- If `95 passed` remains anywhere in Stage 3 context, label it as the historical
  Phase 2 + polish baseline.
- Remove the final blank line at EOF.

Acceptance:

```bash
git diff --check -- USGT_STAGE3_COMPLETION_PLAN.md
```

passes.

### 5. Align transport executable test wording with assertions

Problem: `test_usgt_exe_minimal_transport_from_scratch` says it checks
"physically bounded values", but it only checks run success, `.con` existence,
and mass-budget closure.

Tasks, choose one:

- Preferred conservative fix: change the docstring to say it produces a
  concentration file and closes the species-isolated transport budget.
- Optional stronger fix: if reading `tran.con` is reliable, add an assertion that
  concentration values are finite and physically bounded for the smoke model.

Acceptance:

- The test description matches what the test actually verifies.

### 6. Record resolution

Tasks:

- Add a short "Resolution" section at the top of `USGT_STAGE3_REREVIEW.md` with
  each finding marked `Fixed`.
- Include the final verification results.
- Do not delete the original findings; preserve them as historical review
  context.

Acceptance:

- The doc reads top-to-bottom as a closed re-review, similar to the Phase 2
  review docs.

## Verification Commands

Run after edits:

```bash
python -m pytest autotest/test_usg_transport.py -q
python -m pytest autotest/test_usg_transport_exe.py -q
git diff --check 17f96dd0..HEAD -- .gitignore USGT_PACKAGE_BACKLOG.md USGT_STAGE3_COMPLETION_PLAN.md USGT_UPSTREAM_INFRA.md USGT_improvements.md USGT_roadmap.md USGT_STAGE3_REREVIEW.md autotest/test_usg_transport.py autotest/test_usg_transport_exe.py
```

Optional targeted executable-name check:

```bash
USGT_EXE=mfusg_gsi python -m pytest autotest/test_usg_transport.py autotest/test_usg_transport_exe.py -q
```

Recommended local macOS USG-T 2.7 check:

```bash
USGT_EXE=/Users/roh.sgm/Documents/GitHub_Projects/GW-Software-Compiled-via-Claude/gfortran/usgt_2.7/usgt_270_arm python -m pytest autotest/test_usg_transport.py autotest/test_usg_transport_exe.py -q
```

Expected:

- Focused suite: **101 passed** or higher if tests are added.
- Executable smoke suite: **2 passed** when `USGT_EXE` resolves, clean skip
  otherwise.
- `git diff --check`: clean.

## Done Criteria

- All findings in `USGT_STAGE3_REREVIEW.md` are resolved.
- No package status overclaims remain in the docs.
- `USGT_EXE` is the single executable-selection contract for all USG-T
  executable tests.
- The branch remains focused on Stage 3 cleanup; no unrelated MF6-TID or
  upstream branch changes are included.
