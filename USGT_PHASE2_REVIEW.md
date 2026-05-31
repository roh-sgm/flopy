# USG-T Phase 2 Review

Date: 2026-05-30

Review target: `b0de152f..31a4075f` on `develop`.

Ignore for this review: `8fddfe9a` (`MF6-TID` work-in-progress).

Primary reference:

`/Users/roh.sgm/Documents/SGM Co/GDM/01_Software/GSI/USG-Transport 2.7.0 (2026.03.21)/USGT_V_2-7-0_Source_Code`

## Resolution status (2026-05-30) — all findings addressed

| Finding | Status | What changed |
|---|---|---|
| P0 — Structured DRT authoring broken | **Fixed** | `MfUsgDrt.__init__` now raises `NotImplementedError` on structured models (use `ModflowDrt`); dead structured `write_file` branch removed; `load` still delegates structured files to `ModflowDrt.load`. Tests: structured construction raises, classic `ModflowDrt` structured authoring writes a valid file, structured load delegates to base. |
| P0 — SGB/QRT/DRT main-list controls | **Fixed (Expanded valid write)** | New `flopy/mfusg/_usgt_list.py::begin_list_block` consumes `SFAC` (scales Q for QRT, COND for DRT; inert on SGB gradient per Fortran ISCLOC), `OPEN/CLOSE`, and `EXTERNAL` (via `ext_unit_dict`) before row parsing; unresolved `EXTERNAL` raises `NotImplementedError`. Writes emit expanded inline rows. Tests: SFAC (SGB/QRT/DRT incl. recipients + spreading), OPEN/CLOSE, EXTERNAL negative. |
| P1 — BAS bare `IHM` | **Fixed** | `IHM` parsed with an optional integer: bare `IHM` or `IHM <option>` → `iuihm=0`; `IHM <int>` → that unit. No `IndexError`/`ValueError`. Tests: `FREE IHM`, `FREE IHM 123`, `FREE IHM SY-ALL`, programmatic `iuihm=0/99`. |
| P1 — QRT/DRT recipient metadata divergence | **Fixed** | `write_file` validates `recipient_nodes` against the stress list when `RETURNFLOW` is active (exactly one list per record; omitted ⇒ all-zero; mismatch raises `ValueError`). Tests: too-short/too-long fail for QRT and DRT; omitted ⇒ zero recipients works. |
| P2 — TABRICH node-count contract | **Fixed** | `_tabrich.node_count` now: structured → `nlay*nrow*ncol`; unstructured+DISU → `disu.nodes`; unstructured+DIS → grid-product fallback; otherwise a clear `ValueError`. Test covers all four. |

Verification: `python -m pytest autotest/test_usg_transport.py -q` → **92 passed**
(80 prior + 12 Phase 2 regression tests). No raw `ValueError`/`IndexError` for
valid USG-T list-control / option syntax; unsupported syntax raises
`NotImplementedError`.

## High-Level Summary

The implementation is directionally strong for the main USG-T authoring goal:
new packages can be created from Python inputs, basic 0-based internal indexing
is handled, QRT/DRT recipient ordering matches the Fortran in the common inline
cases, and TABRICH is written in the correct BCF/LPF location.

However, several "Full semantic" claims are too strong until Phase 2 closes the
gaps below. The current tests cover the happy-path inline formats well, but they
do not yet exercise standard MODFLOW list-control syntax (`SFAC`, `EXTERNAL`,
`OPEN/CLOSE`), structured DRT authoring through the new `MfUsgDrt` registration,
or optional BAS `IHM` syntax. These are not cosmetic issues: they affect
round-trip coverage for existing files and correctness when users author
packages from scratch.

Phase 2 should focus on correctness and contract clarity, not new package
surface area. The immediate goal is to make the implementation match the status
claims, or demote the affected package status honestly until support is added.

## Verification Run

Command run:

```bash
python -m pytest autotest/test_usg_transport.py -q
```

Result:

```text
80 passed, 4 warnings in 41.58s
```

Additional probes found failures not covered by the current suite:

- Programmatic structured `MfUsgDrt(...).write_file()` raises
  `AttributeError: 'dict' object has no attribute 'data'`.
- BAS load with a bare `IHM` option raises `IndexError`.
- SGB/QRT/DRT load of a stress-period list beginning with `SFAC` raises
  `ValueError` while the Fortran accepts that list-control syntax.

## Phase 2 Findings

### P0 - Structured DRT authoring is broken

Files:

- `flopy/mfusg/mfusg.py`
- `flopy/mfusg/mfusgdrt.py`

Evidence:

- `MfUsg` now registers `drt` to `MfUsgDrt`.
- `MfUsgDrt.__init__` builds USG-style `dict` stress-period storage directly.
- `MfUsgDrt.write_file()` delegates to `ModflowDrt.write_file(self)` when
  `self.parent.structured` is true.
- Base `ModflowDrt.write_file()` expects base-class storage (`MfList`), not the
  USG-style `dict`, so structured programmatic authoring crashes.

Risk:

This violates the "programmatic authoring from scratch" requirement. It also
makes the handoff statement "structured delegates to base" misleading: load may
delegate, but constructor/write do not form a valid base-class object.

Required fix:

- Choose one explicit structured contract:
  - keep structured `DRT` registration on `ModflowDrt`, and use `MfUsgDrt` only
    for unstructured USG-T DRT8; or
  - make `MfUsgDrt` call/compose the base constructor for structured models and
    keep base `MfList` storage in that path.
- Add a structured authoring regression test that writes a valid classic DRT
  file from scratch.
- Add an unstructured DRT authoring test after the fix to ensure the USG path
  was not regressed.

Acceptance:

- `MfUsgDrt(MfUsg(structured=True), ...).write_file()` either succeeds with
  valid base DRT syntax or is no longer the path exposed by the package
  registry.

### P0 - SGB/QRT/DRT do not support or explicitly reject main-list controls

Files:

- `flopy/mfusg/mfusgsgb.py`
- `flopy/mfusg/mfusgqrt.py`
- `flopy/mfusg/mfusgdrt.py`

Evidence:

- Fortran list readers support standard controls such as `SFAC`, `EXTERNAL`,
  and `OPEN/CLOSE` for these stress-period lists.
- The Python loaders read the next line as a data row immediately:
  - SGB parses row tokens directly in `MfUsgSgb.load`.
  - QRT parses row tokens directly in `_parse_sink_tokens`.
  - DRT parses row tokens directly in `_parse_drain_tokens`.
- When a valid Fortran file uses `SFAC`, the current code raises a raw
  `ValueError` while trying to parse `SFAC` as a node id.

Risk:

This is a round-trip blocker for existing USG-T input files that use normal
MODFLOW list-control syntax. It also breaks the design rule that unsupported
syntax must fail explicitly.

Required fix:

- Implement a shared USG-T list-control reader for these main package lists, or
  use an existing FloPy utility if it can preserve the USG-T row semantics.
- Minimum supported syntax for Phase 2 should be:
  - inline/free rows,
  - `SFAC`,
  - `EXTERNAL`,
  - `OPEN/CLOSE`.
- If any control remains unsupported, detect it before row parsing and raise
  `NotImplementedError` with package name and unsupported control.
- Preserve internal 0-based `node` ids after applying the control syntax.
- Decide whether writes always emit expanded inline rows or preserve controls.
  Expanded inline writes are acceptable if documented as `Expanded valid write`.

Required tests:

- SGB load/write with `SFAC`.
- QRT load/write with `SFAC`, including `RETURNFLOW` recipient blocks.
- DRT load/write with `SFAC`, including single recipient and spreading
  recipient blocks.
- At least one `OPEN/CLOSE` or `EXTERNAL` test per package family, using
  temporary files and `ext_unit_dict` where appropriate.
- Negative tests for any intentionally unsupported control.

Acceptance:

- Valid Fortran list-control input no longer crashes with a raw parser error.
- Package status is updated honestly:
  - `Full semantic` if supported semantically with tests.
  - `Expanded valid write` if read controls are expanded on output.
  - `Partial` if controls are still explicitly unsupported.

### P1 - BAS bare `IHM` option fails to load

File:

- `flopy/mfusg/mfusgbas.py`

Evidence:

- `glo2basu1.f` handles `IHM` by reading an optional `IUIHM` unit and then
  checking whether `IUIHM > 0`.
- `MfUsgBas.load` currently does `int(opts[opts.index("IHM") + 1])`, which
  fails when `IHM` is the last token and can also misinterpret the next option
  as a unit.

Risk:

This blocks round-trip loading of valid BAS option lines and makes
`RICHARDS_HP/IHM` support less robust than the Fortran.

Required fix:

- Parse `IHM` as an option with optional integer argument.
- If no integer follows, set `iuihm=0`.
- If the next token is another option, leave it as an option and set `iuihm=0`.
- Decide whether writing should preserve bare `IHM` or normalize to `IHM 0`.
  Normalizing is acceptable if documented.

Required tests:

- Load BAS with `FREE IHM`.
- Load BAS with `FREE IHM 123`.
- Load BAS with `FREE IHM SY-ALL`.
- Programmatic write with `ihm=True, iuihm=0`.
- Programmatic write with `ihm=True, iuihm>0`.

Acceptance:

- Bare `IHM` loads successfully and produces `ihm=True, iuihm=0`.

### P1 - QRT/DRT recipient metadata can silently diverge from stress data

Files:

- `flopy/mfusg/mfusgqrt.py`
- `flopy/mfusg/mfusgdrt.py`

Evidence:

- Writers default missing recipient lists to empty lists.
- If `recipient_nodes[kper]` is shorter than `stress_period_data[kper]`, later
  records are written with zero recipients instead of failing.
- If recipient lists are longer than the stress-period records, extra metadata
  is ignored.

Risk:

This is a programmatic authoring bug: a user can build an in-memory package that
looks like it has return-flow recipients, but the written USG-T file silently
loses some of them.

Required fix:

- Validate `recipient_nodes` whenever `RETURNFLOW` is active.
- For every explicit stress period, require exactly one recipient list per
  stress record.
- Permit omitted `recipient_nodes` only if the package clearly means all
  records have zero recipients, and document that behavior.
- Reject extra recipient lists.

Required tests:

- QRT mismatch too short fails.
- QRT mismatch too long fails.
- DRT mismatch too short fails.
- DRT mismatch too long fails.
- Valid zero-recipient authoring still works if intentionally supported.

Acceptance:

- Return-flow recipient metadata cannot be silently dropped or shifted during
  write.

### P2 - TABRICH node counting needs a clearer unstructured contract

File:

- `flopy/mfusg/_tabrich.py`

Evidence:

- `_tabrich.node_count()` assumes `model.structured=False` always means a DISU
  package exists.
- Several test patterns in the suite use `MfUsg(structured=False)` with a
  classic `DIS` object for lightweight package tests.

Risk:

TABRICH authoring/loading can fail with `AttributeError` in lightweight models
or mixed compatibility setups. That may be acceptable if TABRICH is declared
DISU-only for unstructured models, but the failure should be explicit.

Required fix:

- Decide the supported contract:
  - require DISU when `structured=False`, with a clear `ValueError`; or
  - fall back to `DIS` dimensions when a classic DIS package is present.
- Add one regression test for the chosen behavior.

Acceptance:

- The failure mode is intentional and documented, or the fallback works.

## Recommended Phase 2 Work Order

1. Update package status labels before coding:
   - mark DRT structured support as reopened,
   - mark SGB/QRT/DRT list-control support as incomplete,
   - mark BAS `IHM` support as incomplete until optional parsing is fixed.
2. Fix structured DRT authoring or remove the structured path from `MfUsgDrt`.
3. Add a shared main-list reader for SGB/QRT/DRT controls.
4. Fix BAS optional `IHM` parsing.
5. Add recipient-list validation for QRT and DRT.
6. Clarify TABRICH unstructured node-count behavior.
7. Expand `autotest/test_usg_transport.py` with the required tests above.
8. Update `USGT_roadmap.md`, `USGT_improvements.md`, and
   `USGT_PACKAGE_BACKLOG.md` after the tests pass.

## Agent Acceptance Checklist

Before Phase 2 can be called complete:

- The exact pytest command above still passes.
- New regression tests fail before the fixes and pass after them.
- No raw `ValueError`/`IndexError` is raised for valid USG-T input syntax.
- Unsupported syntax raises `NotImplementedError` with an actionable message.
- Package status labels do not overclaim beyond implemented and tested support.
- Programmatic authoring tests exist for every fixed behavior.
- Round-trip or expanded-write tests exist for every fixed reader behavior.

