# USG-T Stage 3 Re-Review

Date: 2026-05-31

- Original review target (pre-fix): `17f96dd0..6618181c`
- Follow-up range (the fixes): `6618181c..d3725bd8`
- Resolution commit: `d3725bd8`

This document is **closed**: the Resolution section below records the
follow-up that fixed every finding (commit `d3725bd8`); the original review
text is retained underneath as historical context. It is not a pre-fix-only
document.

Scope: Stage 3 completion pass after `USGT_STAGE3_COMPLETION_PLAN.md`. This
review focuses on consistency between the new tests, the roadmap/backlog docs,
and the stated upstream/CI objective. Stage 3 did not modify package source
files, so this is primarily a documentation and test-contract review rather than
a fresh Fortran semantic audit.

## Resolution (2026-05-31) — all findings Fixed

| Finding | Status | What changed |
|---|---|---|
| P1 — `USGT_EXE` does not control all executable tests | **Fixed** | `test_usg_transport.py` now defines `USGT_EXE = os.environ.get("USGT_EXE", "mfusg_gsi")` once; the 16 real-model `Ex*` run tests use `@requires_exe(USGT_EXE)` and `exe_name=USGT_EXE`. Same single contract as `test_usg_transport_exe.py`. Verified: `USGT_EXE=<abs path to usgt_270_arm>` runs the real-model tests (no `mfusg_gsi` name required); default still resolves `mfusg_gsi`. |
| P2 — roadmap audit still contains pre-completion targets | **Fixed** | The Stage 3 audit table in `USGT_roadmap.md` is rewritten as final state: `DIS/DISU/CLN/DPF/PCB = ✅ Full`; `OC/EVT/MDT/LAK = ✅ (intentionally not Full)` with reasons; `TIB/GSF = Raw/text`; `ETS/HFB/DPT` Partial/Expanded; base classes Compatibility-only. "until their Card lands" removed. |
| P2 — completion plan stale baseline | **Fixed** | `USGT_STAGE3_COMPLETION_PLAN.md` `Current Baseline` now says **101 passed** (95 labeled the historical Phase 2 + polish baseline); trailing EOF blank line removed (`git diff --check` clean). |
| P3 — transport executable test overclaims concentration validation | **Fixed** | `test_usgt_exe_minimal_transport_from_scratch` docstring now states it produces a `.con` and closes the species-isolated transport mass budget (no concentration-value claim); the `.con` is USG node-based binary. |

Verification after the fixes (commit `d3725bd8`), by executable availability —
counts depend on whether `USGT_EXE` resolves, because the real-model `Ex*` and
from-scratch executable tests are `@requires_exe`-gated:

- **USG-T 2.7 resolves via `USGT_EXE`** (e.g. the local ARM binary
  `USGT_EXE=.../usgt_2.7/usgt_270_arm`), both suites together:
  **103 passed**.
- **`mfusg_gsi` available on PATH by default** (no `USGT_EXE` set):
  `autotest/test_usg_transport.py` → **101 passed**;
  `autotest/test_usg_transport_exe.py` → **2 passed**.
- **No resolvable executable** (e.g. `USGT_EXE` unset and `mfusg_gsi` absent,
  or `USGT_EXE` pointing at a missing path), both suites together:
  **85 passed, 18 skipped** (the 16 `Ex*` + 2 from-scratch executable tests
  skip cleanly).
- `git diff --check 17f96dd0..d3725bd8`: **clean**.

CI note: an environment without a USG-T executable should expect
**85 passed, 18 skipped**, not a literal `101 passed`; the executable tier is
optional by design.

Everything below this line is the **original re-review**, preserved for history;
its findings are resolved per the table above.

## Verification Run (original)

Commands run locally:

```bash
python -m pytest autotest/test_usg_transport.py -q
python -m pytest autotest/test_usg_transport_exe.py -q
git diff --check 17f96dd0..6618181c -- .gitignore USGT_PACKAGE_BACKLOG.md USGT_STAGE3_COMPLETION_PLAN.md USGT_UPSTREAM_INFRA.md USGT_improvements.md USGT_roadmap.md autotest/test_usg_transport.py autotest/test_usg_transport_exe.py
```

Results:

- `autotest/test_usg_transport.py`: **101 passed**.
- `autotest/test_usg_transport_exe.py`: **2 passed**.
- `git diff --check`: **failed** on one whitespace issue:
  `USGT_STAGE3_COMPLETION_PLAN.md:424: new blank line at EOF`.

## High-Level Summary

The Stage 3 work is directionally good: the base suite is green, an opt-in
executable smoke tier now exists, and the docs make a much more honest
distinction between Full semantic support, raw/text round-trip, expanded valid
write, and compatibility-only packages.

The remaining issues are not deep package bugs. They are mostly contract
problems between what the docs claim and what the tests currently do. These
should be fixed before handing the branch to another agent as an upstream-ready
baseline, because they affect reproducibility and make the roadmap look more
complete than it actually is.

## Findings

### P1 - `USGT_EXE` does not control all executable tests

Files:

- `autotest/test_usg_transport.py`
- `autotest/test_usg_transport_exe.py`
- `USGT_PACKAGE_BACKLOG.md`
- `USGT_UPSTREAM_INFRA.md`

The new from-scratch executable tests use:

```python
USGT_EXE = os.environ.get("USGT_EXE", "mfusg_gsi")
@requires_exe(USGT_EXE)
```

But the existing real-model `Ex1..Ex9` tests in
`autotest/test_usg_transport.py` still hardcode:

```python
@requires_exe("mfusg_gsi")
MfUsg.load(..., exe_name="mfusg_gsi", ...)
```

This contradicts the docs, which describe the whole executable tier as
`USGT_EXE`-gated. It also weakens the upstream/CI story: a tagged release,
absolute executable path, or CI-installed binary with another name cannot run
the real-model tests without renaming the executable to `mfusg_gsi`.

Fix:

- Define `USGT_EXE = os.environ.get("USGT_EXE", "mfusg_gsi")` once in
  `autotest/test_usg_transport.py`.
- Replace all `@requires_exe("mfusg_gsi")` decorators on `Ex1..Ex9` tests with
  `@requires_exe(USGT_EXE)`.
- Replace all `exe_name="mfusg_gsi"` in those tests with `exe_name=USGT_EXE`.
- Keep the default as `mfusg_gsi` so current local behavior is preserved.
- Update docs only after the tests actually share the same executable-selection
  contract.

Acceptance:

- `USGT_EXE=/absolute/path/to/usgt python -m pytest autotest/test_usg_transport.py -q`
  can run the real-model executable tests without requiring the binary name
  `mfusg_gsi`.
- Default behavior still skips/runs based on `mfusg_gsi` when `USGT_EXE` is
  unset.

### P2 - Stage 3 roadmap audit still contains pre-completion targets

File: `USGT_roadmap.md`

The package coverage table now says:

- `OC`: kept `✅`, not Full.
- `EVT`: kept `✅`, not Full.
- `MDT`: kept `✅`, not Full.
- `CLN`, `DPF`, `DIS`, `DISU`, `PCB`: promoted to Full.

But the Stage 3 status audit table still says things like:

- `OC | ✅ | finish now | Full semantic ...`
- `EVT | ✅ | finish now | Full semantic ...`
- `MDT | ✅ (unverified) | verify or demote | Full semantic or honest demotion`

Those were useful planning targets, but after Stage 3 they read like unresolved
or overclaimed work. This undermines the "honest status" criterion.

Fix:

- Convert the Stage 3 audit table from planning language to final disposition.
- Mark promoted packages as `✅ Full`.
- Mark `OC`, `EVT`, and `MDT` as intentionally retained at plain `✅`, with the
  exact reason from the coverage table.
- Remove or rewrite the sentence that says plain `✅` rows remain until their
  card lands, because the cards have landed.

Acceptance:

- A reviewer can read `USGT_roadmap.md` top-to-bottom without seeing conflicting
  statuses for the same package.

### P2 - Stage 3 completion plan has stale baseline text

File: `USGT_STAGE3_COMPLETION_PLAN.md`

The top of the file says the Stage 3 suite is now **101 passed**, but the
`Current Baseline` section still says:

```text
Current expected result after Phase 2 + polish: 95 passed.
```

That number is now historical, not current. The same file also has a trailing
blank line at EOF that fails `git diff --check`.

Fix:

- Change `Current Baseline` to **101 passed**.
- If `95 passed` is kept, explicitly label it as the Phase 2 historical
  baseline.
- Remove the final blank line so `git diff --check` passes.
- Check whether `USGT_PACKAGE_BACKLOG.md` and `USGT_improvements.md` use `95
  passed` only for the Phase 2 section, not Stage 3.

Acceptance:

- `git diff --check 17f96dd0..HEAD -- USGT_STAGE3_COMPLETION_PLAN.md` passes.
- Stage 3 docs consistently identify the focused suite as 101 tests.

### P3 - Transport executable test overclaims concentration validation

File: `autotest/test_usg_transport_exe.py`

The transport smoke-test docstring says it "produces a concentration file with
physically bounded values", but the assertions only check:

- run success,
- `.con` file exists,
- transport list budget exists,
- `PRESCRIBED_CONCS_IN` is present,
- final mass-budget discrepancy is small.

That is still a useful smoke test, but it does not check concentration values.

Fix, choose one:

- Add a lightweight concentration-read assertion if the `.con` file can be read
  reliably for this structured smoke model, for example finite values and a
  simple physical bound.
- Or change the docstring to say the test verifies that concentration output is
  produced and the transport budget closes.

Acceptance:

- The test description matches the assertions.

## Agent Work Order

1. Fix executable selection consistency in `autotest/test_usg_transport.py` by
   making `Ex1..Ex9` use the same `USGT_EXE` contract as
   `autotest/test_usg_transport_exe.py`.
2. Clean up `USGT_roadmap.md` so the Stage 3 status audit is final-state
   documentation, not stale planning language.
3. Clean up `USGT_STAGE3_COMPLETION_PLAN.md`: current suite is 101 passed,
   Phase 2 was 95 passed, and `git diff --check` must be clean.
4. Align the transport executable test docstring with what it really asserts, or
   add the missing concentration-value assertion.
5. Re-run:

```bash
python -m pytest autotest/test_usg_transport.py -q
python -m pytest autotest/test_usg_transport_exe.py -q
git diff --check 17f96dd0..HEAD -- .gitignore USGT_PACKAGE_BACKLOG.md USGT_STAGE3_COMPLETION_PLAN.md USGT_UPSTREAM_INFRA.md USGT_improvements.md USGT_roadmap.md autotest/test_usg_transport.py autotest/test_usg_transport_exe.py
```

Expected result:

- Focused suite: **101 passed** or higher if tests are added.
- Executable smoke suite: **2 passed** when the executable is available, cleanly
  skipped otherwise.
- `git diff --check`: clean.
