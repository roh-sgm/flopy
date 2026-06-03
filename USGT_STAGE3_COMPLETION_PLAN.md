# USG-T 2.7 FloPy Completion Plan

Date: 2026-05-31

Purpose: guide the next agents through the remaining work needed to move this
fork from strong project support to an upstream-ready USG-Transport 2.7 support
story.

Primary FloPy references:

- `USGT_roadmap.md`
- `USGT_PACKAGE_BACKLOG.md`
- `USGT_PHASE2_REVIEW.md`
- `USGT_PHASE2_REREVIEW.md`
- `USGT_improvements.md`

Primary USG-T source reference:

`/Users/roh.sgm/Documents/SGM Co/GDM/01_Software/GSI/USG-Transport 2.7.0 (2026.03.21)/USGT_V_2-7-0_Source_Code`

## Stage 3 status: COMPLETE (2026-05-31)

All ten cards are done (commits on `develop`, pushed to `origin`). Default
focused suite `python -m pytest autotest/test_usg_transport.py -q` →
**101 passed** with a USG-T executable, or **85 passed, 16 skipped** without one
(see Current Baseline below); the opt-in `autotest/test_usg_transport_exe.py`
runs end-to-end under `USGT_EXE` (default `mfusg_gsi`). Outcomes:

- **C1 audit:** per-package decision table in `USGT_roadmap.md`; TIB/GSF
  relabeled Raw/text.
- **C2 parameters:** Expanded valid write (ETS) / explicit fail (NP* list
  params); ETS expand test added.
- **C3 DPT A-W_ADSORBIM:** spec'd; kept explicit-unsupported (rare).
- **C4 TIB:** grammar recorded; kept Raw/text v1.
- **C5 LAK:** kept `✅` not Full (Ex8-validated; authoring deferred).
- **C6 base-class:** SFR/STR/GAGE/FHB/SUB/SWT = Compatibility-only (FHB/GAGE
  exercised via Ex8).
- **C7 recipient U1DINT:** INTERNAL/CONSTANT supported; external = explicit fail.
- **C8 plain-`✅`:** promoted CLN/DPF/DIS/DISU/PCB to Full; OC/EVT/MDT honest `✅`.
- **C9 executable:** `USGT_EXE`-gated from-scratch flow+transport e2e tests.
- **C10 upstream:** `USGT_UPSTREAM_INFRA.md` proposal.

## Current Baseline

The main USG-T authoring path is now strong. The default focused suite is:

```bash
python -m pytest autotest/test_usg_transport.py -q
```

Current expected result depends on whether a USG-T executable resolves (the
`Ex1..Ex9` run tests are `@requires_exe`-gated):

- with a USG-T executable available (`mfusg_gsi` on PATH, or `USGT_EXE` set):
  **101 passed**;
- with no resolvable executable: **85 passed, 16 skipped** (the 16 `Ex*` run
  tests skip cleanly).

The opt-in `autotest/test_usg_transport_exe.py` adds 2 from-scratch executable
tests; both suites together are **103 passed** under `USGT_EXE` (USG-T 2.7) or
**85 passed, 18 skipped** with no executable. The earlier **95 passed** figure
was the Phase 2 + polish baseline and is now historical.

Known local repo state: generated USG-T transport `.CBB` outputs are ignored
for the `Ex3_CLN_Conduit/Dispersion` validation folder because they are too
large to commit.

## Definition of Done

Do not call a package or feature `Full semantic` unless all are true:

- It has a Fortran-derived I/O spec in the agent notes or docs.
- It can be authored from scratch with Python/numpy inputs.
- It can load existing files.
- It writes valid USG-T 2.7 input.
- It has separate authoring and round-trip tests.
- Unsupported modes fail explicitly with actionable messages.
- `USGT_roadmap.md`, `USGT_PACKAGE_BACKLOG.md`, and `USGT_improvements.md` are
  updated.

Valid non-full labels:

- `Expanded valid write`: richer input syntax is read but output is normalized
  to expanded valid USG-T input.
- `Raw/text round-trip`: existing text is preserved, but no semantic API is
  claimed.
- `Compatibility-only`: base FloPy class is used and USG-T-specific records are
  not independently validated.

## Stage 3 Work Order

### 0. Keep The Baseline Clean

Goal: keep the branch stable while agents add deeper support.

Tasks:

- Confirm `.gitignore` keeps large generated `.CBB` files out of `git status`.
- Run `python -m pytest autotest/test_usg_transport.py -q` before and after each
  package card.
- Avoid touching MF6-TID work.
- Keep all pushes to `origin`, not upstream.
- Preserve existing user changes; do not rewrite unrelated files.

Acceptance:

- `git status --short` shows only intended source/doc changes.
- The USG-T focused suite remains green after every card (101 passed after
  Stage 3; 95 was the Phase 2 + polish baseline).

### 1. Package Status Audit

Goal: produce a fresh, honest map of what is complete, partial, raw, or
compatibility-only.

Tasks:

- Re-read `USGT_roadmap.md` and mark every non-`Full` row with one of:
  `finish now`, `explicitly defer`, or `compatibility-only`.
- Check packages currently marked plain `✅` but not `Full`: `DIS`, `DISU`, `OC`,
  `EVT`, `LAK`, `CLN`, `PCB`, `MDT`, `DPF`, `TIB`, `GSF`.
- For each plain `✅`, decide whether it should become:
  - `Full semantic`,
  - `Expanded valid write`,
  - `Raw/text round-trip`,
  - or remain `Compatibility-only`.
- Update the roadmap before implementation if status labels overclaim.

Acceptance:

- There is no ambiguous `✅` status for a package that matters to USG-T 2.7
  authoring.
- The roadmap reflects authoring reality, not only load/write success.

### 2. Parameter Strategy Card

Packages: `ETS`, `HFB`, `SGB`, `QRT`, `DRT`

Goal: decide and implement one consistent v1 policy for MODFLOW-style
parameters.

Current state (at the start of Stage 3 — **now historical**; the Stage 4.4/4.6
parameter work has since preserved + authored these. Authoritative state:
`USGT_roadmap.md` / `USGT_STAGE4_09_FINAL_GAP_AUDIT.md`):

- `ETS`: parameterized input loads and writes expanded `NPETS=0`.
- `HFB`: `NPHFB>0` fails explicitly.
- `SGB/QRT/DRT`: named parameters fail explicitly.

Decision required:

- Option A: keep `Expanded valid write` / explicit unsupported behavior for all
  parameterized inputs.
- Option B: implement parameter preservation and active-parameter records.

Recommended v1 path:

- Keep expanded/explicit behavior for Stage 3 unless a real model requires
  parameter preservation.
- Strengthen tests and docs so upstream reviewers see that the behavior is
  intentional.

Tasks:

- Audit `parutl7.f` plus package-specific Fortran parameter calls.
- Add or verify negative tests for `NPSGB>0`, `NPQRT>0`, `NPDRT>0`, `NPHFB>0`.
- Add one expanded-write test for ETS parameterized input if not already strong
  enough.
- Update docs to say parameter preservation is deferred by design.

Acceptance:

- No parameterized package writes incomplete parameter syntax.
- Every unsupported parameter path raises `NotImplementedError`.
- Status labels say `Expanded valid write` or `Partial`, not `Full`, where
  preservation is missing.

### 3. DPT Immobile Air-Water Adsorption

Package: `DPT`

Fortran: `gwt2dptu1.f`, `dpt2aw_adsorb.f`

Goal: either implement `A-W_ADSORBIM` semantically or keep it explicitly
unsupported with full documentation.

Tasks:

- Write a Fortran-derived spec for `A-W_ADSORBIM`: option tokens, extra function
  indices, array order, conditional reads, and interaction with `DPF`,
  `IDISPIM`, and BCT `A-W_ADSORB`.
- Decide product scope:
  - implement full semantic authoring/load/write; or
  - keep explicit unsupported failure because this mode is rare.
- If implementing:
  - add constructor fields,
  - add load/write support,
  - add authoring and round-trip tests.
- If deferring:
  - add a focused test proving the option fails before any shifted reads,
  - document why it remains Partial.

Acceptance:

- DPT array layout cannot silently shift under any supported or unsupported
  adsorption option.

### 4. TIB Semantic Authoring Decision

Package: `TIB`

Current state: raw/text round-trip.

Goal: decide whether TIB needs semantic authoring for upstream readiness.

Tasks:

- Audit the TIB Fortran read path and write a small grammar/spec.
- Identify the minimal useful semantic API: node lists, `U1DINT` controls,
  stress-period behavior, and output normalization.
- Decide:
  - keep raw/text round-trip for v1; or
  - implement semantic constructor plus writer.
- If implementing semantic mode, preserve raw-body fallback for unknown files.

Required tests if semantic:

- Minimal from-scratch TIB.
- Multi-node-per-line `U1DINT`.
- Round-trip existing raw file.
- Mixed raw/semantic failure behavior if unsupported combinations exist.

Acceptance:

- TIB status is honest. It should not be presented as authoring-ready unless a
  semantic constructor and tests exist.

### 5. LAK Authoring Completion

Package: `LAK`

Current state: Ex8 real-model load/write validates `TABLEINPUT` and
`TRANSPORTBOUNDARY`; from-scratch authoring deferred.

Goal: close LAK if it matters for upstream completeness.

Tasks:

- Audit `gwf2lak7u1.f` and USG-T lake transport extensions.
- Define constructor arguments for `TABLEINPUT`, `TRANSPORTBOUNDARY`, and lake
  transport concentrations.
- Add a minimal synthetic LAK model that writes those options from scratch.
- Add a round-trip test preserving Ex8 behavior.

Acceptance:

- LAK can move from plain `✅` to `Full semantic`, or docs explicitly keep it
  `✅ not Full` with a rationale.

### 6. Base-Class Compatibility Packages

Packages: `SFR`, `STR`, `GAGE`, `FHB`, `SUB`, `SWT`

Current state: compatibility-only / Partial.

Goal: decide whether these are truly out of scope or need USG-T semantics.

Recommended path:

- Keep them compatibility-only unless a real USG-T 2.7 model in the target use
  cases requires them.
- Do not spend time implementing surface-water packages if CLN is the intended
  coupling route for this project.

Tasks:

- For each package, inspect Fortran for USG-T-specific node/unstructured record
  differences.
- Search real validation models for package usage.
- Classify each as:
  - `out of scope`,
  - `round-trip compatibility needed`,
  - or `semantic authoring needed`.
- If in scope, create a separate card per package. Do not bundle all six into
  one PR.

Acceptance:

- Roadmap says clearly which packages are compatibility-only and why.
- No base-class package is implied to support USG-T-specific records unless it
  has tests.

### 7. Recipient `U1DINT` External Controls

Packages: `QRT`, `DRT`

Current state: main lists handle `SFAC`, `OPEN/CLOSE`, `EXTERNAL`; recipient
`U1DINT` blocks support `INTERNAL` and `CONSTANT`, while external recipient
lists fail explicitly.

Goal: decide whether recipient-node `U1DINT` external controls need support.

Tasks:

- Verify in Fortran whether recipient-node `U1DINT` blocks accept
  `EXTERNAL`/`OPEN/CLOSE` in the same way as standard arrays.
- Search real USG-T models for recipient lists using external controls.
- If needed, extend `_usgt_returnflow.read_u1dint_list` to resolve external
  controls using the same path policy as `_usgt_list`.
- Add tests for `INTERNAL`, `CONSTANT`, and any external modes implemented.

Acceptance:

- Recipient list behavior is either supported or explicitly documented as an
  unsupported rare syntax.

### 8. Plain `✅` Package Hardening

Packages: `OC`, `EVT`, `PCB`, `MDT`, `DPF`, `CLN`, `DISU`

Goal: convert important plain `✅` packages to `Full semantic` where justified.

Tasks:

- `OC`: add/verify authoring tests for ATS adaptive time stepping and
  BOOTSTRAPPING.
- `EVT`: add/verify transport concentration array authoring.
- `PCB`: verify against Fortran source, not only field order.
- `MDT`: independently verify against Fortran source or demote status.
- `DPF`: decide whether current focused tests justify `Full semantic`; if yes,
  update roadmap, otherwise list missing fields.
- `CLN`: decide whether current connectivity/geometries tests justify
  `Full semantic`; if yes, update roadmap.
- `DISU`: keep large-grid formatting tests protected.

Acceptance:

- Important packages are either `Full semantic` with tests or clearly marked
  less than Full.

### 9. Executable End-To-End Validation

Goal: prove that FloPy-authored USG-T inputs run with USG-T 2.7 and produce
stable outputs.

Tasks:

- Create a slow/optional test marker for USG-T executable tests.
- Define how the executable path is supplied locally and in CI:
  - env var such as `USGT_EXE`,
  - local fixture path,
  - or downloaded tagged release once available.
- Select a small matrix:
  - one minimal model per major package family,
  - one real Ex model load/write/run,
  - one multi-species transport budget check.
- Compare outputs at the right tolerance:
  - `.lst` budget summaries,
  - `.hds`,
  - `.con`,
  - `.cbb/.cbc`,
  - species-isolated budget records.
- Keep large generated binary outputs out of git.

Acceptance:

- Local optional suite can run USG-T 2.7 end-to-end.
- Default CI remains light, but there is a documented slow validation tier.

### 10. Upstream Infrastructure Track

Goal: satisfy the concerns in the November 2024 email before asking USGS to
accept deeper USG-T support.

This track is partly outside FloPy and should be handled in parallel only after
the package-level story is stable.

Tasks:

- Verify current upstream reality before implementation:
  - FloPy GitHub workflows,
  - `get-modflow`,
  - `modflowpy/install-modflow-action`,
  - `pymake`,
  - any current MODFLOW-USG / USG-T release source.
- Prepare or improve a USG-Transport source repository:
  - tagged releases,
  - source archive assets,
  - build instructions,
  - CI compile on Linux/macOS/Windows if feasible,
  - smoke tests that run the official executable.
- Add USG-T executable autotests to that source repo before asking FloPy to
  depend on it.
- Only after source-repo CI exists:
  - propose `pymake` changes if needed,
  - propose `get-modflow` / install action changes if needed,
  - wire optional FloPy USG-T executable tests to that tagged release.

Acceptance:

- There is a reproducible tagged USG-T release source.
- There is CI proving the executable builds and passes smoke tests.
- FloPy does not depend on an ad-hoc zip with unknown test status.

## Suggested Agent Cards

Use one agent per card unless the card is purely documentation.

1. `stage3/status-audit`: update roadmap statuses and produce a package-by-
   package decision table.
2. `stage3/parameters`: ETS/HFB/SGB/QRT/DRT parameter policy and tests.
3. `stage3/dpt-aw-adsorbim`: implement or explicitly close DPT
   `A-W_ADSORBIM`.
4. `stage3/tib-authoring-decision`: semantic TIB spec and implementation
   decision.
5. `stage3/lak-authoring`: TABLEINPUT/TRANSPORTBOUNDARY from-scratch authoring.
6. `stage3/base-compat`: SFR/STR/GAGE/FHB/SUB/SWT scope decisions.
7. `stage3/recipient-u1dint`: QRT/DRT recipient-list external controls.
8. `stage3/plain-checkmarks`: harden OC/EVT/PCB/MDT/DPF/CLN/DISU status.
9. `stage3/executable-validation`: optional USG-T 2.7 run suite.
10. `stage3/upstream-infra`: tagged release / pymake / get-modflow path.

## Priority Recommendation

Recommended order:

1. Status audit.
2. Parameter strategy.
3. DPT explicit/semantic decision.
4. TIB and LAK authoring decisions.
5. Plain `✅` hardening.
6. Executable validation.
7. Upstream infrastructure.
8. Compatibility-only packages only if a real model requires them.

Rationale: upstream reviewers will care first that status claims are honest,
tests are clear, and unsupported modes fail predictably. Full support for rare
packages is less persuasive than a clean, reproducible, well-tested core.
