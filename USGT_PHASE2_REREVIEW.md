# USG-T Phase 2 Re-Review

Date: 2026-05-31

Review target (re-review): `31a4075f..248fbdc8` on `develop`. This re-review is
resolved by the polish pass `248fbdc8..476a7297` (see Resolution status below).

Scope note: ignore `8fddfe9a` and the MF6-TID files that appear in the raw diff
for the re-review range, because that commit sits between those two endpoints.

Primary reference:

`/Users/roh.sgm/Documents/SGM Co/GDM/01_Software/GSI/USG-Transport 2.7.0 (2026.03.21)/USGT_V_2-7-0_Source_Code`

## Resolution status (2026-05-31) — all re-review findings addressed

| Finding | Status | What changed |
|---|---|---|
| P1 — Quoted `OPEN/CLOSE` filenames | **Fixed** | `_usgt_list._open_close_filename` parses the filename quote-aware from the raw line: bare token, single/double-quoted, and single-quoted names containing spaces. Tests: SGB `'rows.dat'` and `"rows.dat"`, DRT `'rows with spaces.dat'`; expanded write still emits inline rows (no `OPEN/CLOSE`). |
| P2 — Positive `EXTERNAL` pytest coverage | **Fixed** | Added a positive `EXTERNAL` test for SGB, QRT, and DRT via a minimal `ext_unit_dict` (`NamData`); the SGB and QRT external files begin with `SFAC`, covering source-switch + scale together. |
| P2 — DRT omitted-recipient zero path | **Fixed** | Added a DRT `RETURNFLOW` test with multiple records and omitted `recipient_nodes`; write succeeds and reload asserts `recipient_nodes[0] == [[], []]` (mirrors the QRT test). |
| P3 — `USGT_PHASE2_REVIEW.md` mixed open/closed text | **Fixed** | Added a "Current status" note; relabeled the prose as "(original, pre-fix)" / "Original Findings (pre-fix)"; labeled the `80 passed` run as the pre-fix baseline; removed the trailing blank line at EOF. |

Verification: `python -m pytest autotest/test_usg_transport.py -q` → **95 passed**
(92 + 3 polish-pass tests). No raw filesystem/parse error for quoted
`OPEN/CLOSE`; `EXTERNAL` via `ext_unit_dict` is covered in default CI.

## Current status

The polish pass is **complete** and this re-review is **closed** (suite: 95
passed; resolved in `248fbdc8..476a7297`). Everything below this line is the
**original re-review (pre-polish)**, preserved for history; its findings are
resolved per the resolution table above and must not be read as still open.

---

## High-Level Summary (original re-review, pre-polish)

The five Phase 2 findings are substantially closed. The main authoring/round-trip
contracts are better now: structured `MfUsgDrt` no longer crashes silently, the
SGB/QRT/DRT list-control path handles the common `SFAC`, `OPEN/CLOSE`, and
`EXTERNAL` cases, BAS `IHM` no longer raises on bare syntax, and QRT/DRT
recipient metadata is validated before write.

The remaining issues are narrower hardening tasks. I do not see a new P0 design
blocker, but one valid `OPEN/CLOSE` syntax still fails, and a few test/doc gaps
make the "all closed" story slightly too broad. Treat this as a polish pass
before calling Phase 2 fully upstream-ready.

## Verification Run (original re-review, pre-polish baseline)

> Pre-polish baseline only. After the polish pass the suite passes **95** — see
> "Current status" above.

Command run:

```bash
python -m pytest autotest/test_usg_transport.py -q
```

Result:

```text
92 passed, 4 warnings in 42.39s
```

Additional manual probes (pre-polish):

- `EXTERNAL` with a resolvable `ext_unit_dict` was manually checked for SGB,
  QRT, and DRT; all three loaded and applied expected node/scale semantics.
- `OPEN/CLOSE rows.dat` works.
- `OPEN/CLOSE 'rows.dat'` failed with `FileNotFoundError` because the helper did
  not strip quotes before joining the path. *(Fixed in the polish pass: quoted
  names — including single-quoted names with spaces — are now parsed.)*

## Original Findings (pre-polish)

> All findings below are **resolved** — see the resolution table at the top of
> this file. They are retained verbatim for historical context.

### P1 - Quoted `OPEN/CLOSE` filenames are not parsed like Fortran/FloPy

File:

- `flopy/mfusg/_usgt_list.py`

Evidence:

- `begin_list_block` takes `tok[1]` directly for `OPEN/CLOSE` and joins it with
  `model.model_ws`.
- The Fortran `URWORD` parser supports single-quoted words, including file
  names with spaces or commas.
- Existing FloPy utility parsing strips single and double quote characters from
  `OPEN/CLOSE` filenames.
- Manual probe: `OPEN/CLOSE rows.dat` succeeds, but `OPEN/CLOSE 'rows.dat'`
  tries to open a literal path containing the quote characters and raises
  `FileNotFoundError`.

Risk:

Valid USG-T input can still fail in the newly added list-control path. This is
especially relevant for model folders or external list files with spaces in
their names.

Required fix:

- Parse `OPEN/CLOSE` filenames with the same quote handling used by existing
  FloPy utilities, or reuse a shared parser if possible.
- At minimum, strip surrounding single and double quotes.
- Prefer supporting single-quoted names with spaces because `URWORD` explicitly
  allows them.

Required tests:

- SGB `OPEN/CLOSE 'rows.dat'`.
- QRT or DRT `OPEN/CLOSE 'rows with spaces.dat'`.
- Confirm expanded write still emits inline rows and does not preserve
  `OPEN/CLOSE`.

Acceptance:

- Quoted `OPEN/CLOSE` filenames that USG-T Fortran accepts load without a raw
  filesystem error.

### P2 - Positive `EXTERNAL` coverage is missing from the pytest suite

Files:

- `autotest/test_usg_transport.py`
- `flopy/mfusg/_usgt_list.py`

Evidence:

- The suite tests unresolved `EXTERNAL` as an explicit `NotImplementedError`.
- It does not test successful `EXTERNAL` resolution via `ext_unit_dict`.
- Manual probes with minimal `ext_unit_dict` entries succeeded for SGB, QRT, and
  DRT, including `SFAC` inside the external file.

Risk:

The code path works in a manual probe, but future edits could break it without
CI catching the regression.

Required fix:

- Add one positive `EXTERNAL` test for each package family, or one parameterized
  test covering SGB/QRT/DRT.
- Include at least one external file that begins with `SFAC` so source-switching
  and scale application are covered together.

Acceptance:

- The claim "`EXTERNAL` via `ext_unit_dict` is supported" is backed by default
  pytest coverage.

### P2 - DRT omitted-recipient zero path lacks a direct regression test

Files:

- `autotest/test_usg_transport.py`
- `flopy/mfusg/mfusgdrt.py`

Evidence:

- QRT has a regression test for `RETURNFLOW` with omitted `recipient_nodes`,
  which should write all-zero recipients.
- DRT uses the same validation contract but does not have the equivalent direct
  test.

Risk:

Low implementation risk because the code is symmetric, but the documented
contract should be covered for both packages.

Required fix:

- Add a DRT authoring test with `RETURNFLOW`, multiple records, and no
  `recipient_nodes`.
- Reload and assert `recipient_nodes[0] == [[], ...]`.

Acceptance:

- Both QRT and DRT have explicit zero-recipient omission tests.

### P3 - `USGT_PHASE2_REVIEW.md` still mixes resolved and original-open text

File:

- `USGT_PHASE2_REVIEW.md`

Evidence:

- The top table says all findings are addressed.
- Later sections still say the gaps remain and show the old `80 passed` result.
- `git diff --check 31a4075f..248fbdc8` also reports a blank-line-at-EOF issue
  in this file. The same command reports trailing whitespace in MF6-TID docs,
  but those files are explicitly out of scope for this USG-T review.

Risk:

This is not a runtime issue, but it makes the review artifact ambiguous for the
next agent. The document should read as a resolved review with an archived
"original findings" section, not as both open and closed at the same time.

Required fix:

- Rename the old open sections to "Original Findings".
- Add a short "Current Status" section after the resolution table.
- Update the old verification result or clearly label it as the pre-fix result.
- Remove the extra blank line at EOF.

Acceptance:

- A new agent can read the file top-to-bottom without wondering whether Phase 2
  is open or closed.

## Recommended Work Order

1. Fix quoted `OPEN/CLOSE` parsing in `_usgt_list.begin_list_block`.
2. Add positive `EXTERNAL` tests for SGB/QRT/DRT.
3. Add the DRT omitted-recipient regression test.
4. Clean up `USGT_PHASE2_REVIEW.md` so it is internally consistent.
5. Run `python -m pytest autotest/test_usg_transport.py -q`.
