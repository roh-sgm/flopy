# Stage 4.4D — MODFLOW Parameter Preservation: DRT (executed)

Scope: implement load → write → reload **preservation** of DRT list parameters
(`NPDRT > 0`), **only after** auditing the Fortran for the kind of `PARTYP`
mismatch that made SGB active parameters non-executable. Only `MfUsgDrt` and the
shared parameter helper are touched; SGB/HFB/ETS/QRT are not (beyond shared
docs).

## Fortran audit (mandatory, post-SGB)

### `gwf2drt8u.f`

- **Item 1 (AR, `GWF2DRT8U1AR`)**: `MXADRT IDRTCB NPDRT MXL [options]` read
  directly (`gwf2drt8u.f:46/50-53`) — `NPDRT`/`MXL` are in item 1, there is **no**
  separate `PARAMETER` line (unlike SGB). Options: `AUX`, `RETURNFLOW`
  (`IDRTFL=4`), `NOPRINT`, `SPREAD <MXSPREADNDS>`, `CHANGEC` (`IDCHANGEC`).
  `NDRTVL = 5 + NAUX + 2 + IDRTFL`.
- **Definitions (items 2-3)**: for `K=1..NPDRT`,
  `UPARLSTRP(...,'DRT','DRT',1,NUMINST)` (`gwf2drt8u.f:135`) — `UPARLSTRP`
  enforces `PARTYP == PTYPX`, here `'DRT'`. Each definition's `NLST` rows are read
  by `SGWF2DRT8LR` (the **same** list reader used for non-parametric rows), so a
  parameter row carries its `NODE EL COND NR [PROP] [IDCHNGTYP] [aux]` and its
  RETURNFLOW recipients / spreading `U1DINT` block.
- **Per stress period (RP, `GWF2DRT8U1RP`)**: `ITMP NP` when `NPDRT>0`, else just
  `ITMP` (`gwf2drt8u.f:186-199`). `ITMP<0` reuses non-parametric drains; `ITMP>0`
  reads them. Then `PRESET('DRT')` and, if `NP>0`, `NP` activations via
  `SGWF2DRT8LS` (`gwf2drt8u.f:240`).

### Critical check: definition vs activation `PARTYP` (the SGB question)

- definition: `UPARLSTRP(...,'DRT','DRT',...)` → `PARTYP='DRT'` (`gwf2drt8u.f:135`);
- activation: `SGWF2DRT8LS` sets `PTYP='DRT'` (`gwf2drt8u.f:1108`) and enforces
  `PARTYP(IP) == PTYP` (`gwf2drt8u.f:1133`).

**Both are `'DRT'` — consistent.** Unlike SGB (`'SGB'` definition vs `'G'`
activation, which aborts), an active DRT parameter is valid and executable in
USG-T 2.7. So DRT active-parameter preservation is implemented (not just
definitions).

### Other interactions

- **Scaling**: `SGWF2DRT8LS` uses `IPVL1=IPVL2=5`, so the parameter value scales
  `COND` (column 5). `SFAC` (in `SGWF2DRT8LR`, `ISCLOC=5`) also scales `COND`.
- **RETURNFLOW / recipient_nodes / SPREAD**: parameter rows are read by the same
  `SGWF2DRT8LR` as non-parametric rows, so they carry recipients exactly the same
  way — single inline recipient (`NR>0`) or a `-NR` spreading `U1DINT` block.
  Preserving a parameter therefore requires preserving its per-row recipients.
  **Execution caveat (review):** activation copies only `DRTF`
  (`SGWF2DRT8LS`), **not** `NodDRT`, and the spreading consumer
  (`GWF2DRT8U1AD` / RP) walks `NodDRT` sequentially (`gwf2drt8u.f:262-276`). An
  inline recipient (`NR>0`) lives in `DRTF` and is copied, so an *activated*
  inline-recipient parameter is fine; an *activated* SPREAD (`NR<0`) parameter is
  **structural round-trip only — not execution-guaranteed**. FloPy preserves the
  recipients faithfully on load/write but does not claim every activated-SPREAD
  case runs in USG-T 2.7.
- **MXADRT (item 1)**: `IDRTPB=MXADRT+1`, `MXDRT=MXADRT+MXL`; in RP `MXADRT` is
  the cap on `NDRTCL`, and `SGWF2DRT8LS` does `NDRTCL = NDRTCL + NLST` and aborts
  if `NDRTCL > MXADRT` (`gwf2drt8u.f:1173-1174`). So `MXADRT` must cover the
  *active* total per period (non-parametric drains + the `NLST` rows of every
  active parameter), not just the non-parametric count.
- **CHANGEC / IDCHNGTYP** and **AUX**: present in the parameter rows (same dtype),
  preserved.
- **INSTANCES**: `UPARLSTRP`/`UINSRP` support `NUMINST>0` (instance blocks).
- **ITMP<0 reuse**: applies to the non-parametric drains; a reuse period can still
  carry `NP>0` activations.

## Design — shared helper + per-row recipients

The shared `flopy/mfusg/_usgt_parameters.py` helpers are reused unchanged:
`read_list_parameter_header` / `write_list_parameter_header` (the
`PARNAM PARTYP PARVAL NLST` header) and `read_active_list_parameters` /
`write_active_list_parameters` (activation names). DRT does **not** use the
`PARAMETER NP MXL` helpers (its `NPDRT`/`MXL` are in item 1).

The `NLST` drain rows + recipients stay in DRT, reusing the existing row
machinery: a new `_read_drain_rows` (factored out of `load`, shared by the
definitions and the non-parametric rows) returns `(recarray, recipient_lists)`
via `_parse_drain_tokens` + `begin_list_block`; `_write_drain_line` writes each
row + its spreading block. So a preserved parameter keeps **per-row
recipient_nodes**.

## DRT implementation (`MfUsgDrt`)

- `load`: parse `NPDRT`/`MXL` from item 1; when `NPDRT>0`, read each definition
  (`read_list_parameter_header` + `NLST` rows via `_read_drain_rows`, 0-based)
  into `self.parameters = {name: {"partyp","parval","nlst","data",
  "recipient_nodes"}}` and keep `MXL` as `self.mxl`. Per stress period, read
  `ITMP NP` (`NP` only when `NPDRT>0`), the non-parametric drains (`ITMP<0`
  reuse), then `NP` activation names into `self.active_params = {kper: [name,
  ...]}`.
- `write_file`: when preserving, write `NPDRT MXL` in item 1, the definition
  blocks (header + rows + recipients), then per SP `ITMP NP` + the
  non-parametric drains + the activation names. Validate first (before opening
  the file). Non-parametric DRT (authoring, `ITMP<0` reuse, RETURNFLOW, SPREAD,
  CHANGEC, AUX, SFAC, structured rejection/delegation) is unchanged.
- **MXADRT (item-1 field 1)** is `_max_active_drains()` — the max over stress
  periods of `len(non-parametric rows) + sum(nlst of active parameters)` (review
  follow-up). Previously it counted only the non-parametric rows, which would
  under-size `MXADRT` whenever a period activates parameters and the Fortran
  would abort (`NDRTCL > MXADRT`). `ITMP<0` reuse carries the previous period's
  non-parametric count into the sum.
- Indexing: internal 0-based, file 1-based, for both definition rows (and their
  recipients) and non-parametric rows.

### Validation / explicit failures

- **From-scratch parameter authoring is now supported (Stage 4.6B)** — pass
  `parameters={name: {...}}` + `active_params={kper: [...]}`; see
  `USGT_STAGE4_10_DRT_PARAMETER_AUTHORING.md`. (Previously this raised
  `NotImplementedError`.)
- **`active_params` with no `parameters` defined** → `ValueError`.
- **`INSTANCES` (`NUMINST>0`)** → `NotImplementedError` on load (instances
  combined with per-row recipient lists are not implemented yet); the authoring
  API has no instance field.
- **`MXL`** → auto-computed as the total definition rows when omitted;
  `ValueError` when given and `<` that total.
- **Inconsistent definitions** → `ValueError` (`partyp` not `DRT`, missing
  `parval`, empty `data`, `nlst != len(data)`, `recipient_nodes` length, or
  recipients without `RETURNFLOW`); **`active_params`** → `ValueError` for an
  undefined name, a duplicate activation in a period, or a stress period outside
  `0..nper-1`.
- **Contract hardening (Stage 4.6B review follow-up)** → `ValueError`, all
  before the file is opened: a definition/active **name** that is not a single
  whitespace-free token or exceeds 10 chars (Fortran `CHARACTER*10` PARNAM,
  upper-cased); a **duplicate definition name** case-insensitively (`dp`/`DP`);
  a **`parval`** that is a multi-token string (e.g. `"1 2"`); a **negative**
  parametric `data["node"]` or `recipient_nodes` entry (nodes are 0-based
  non-negative integers).
- All validation runs before the file is opened — no partial file.

## Decision: honest status

DRT is **`✅ Full (authoring)`** for the non-parametric package, is
**parameter-preserving for `NPDRT` (definitions + activations + recipients)**,
and (Stage 4.6B) supports **from-scratch `NPDRT>0` parameter authoring** — DRT is
internally consistent in the Fortran, so active parameters are real. Remaining
caveats on the parametric axis:

- parameter `INSTANCES` are unsupported (`NotImplementedError`), although the
  Fortran supports them; the from-scratch authoring API has no instance field;
- the `ITMP<0` reuse of non-parametric drains is written **expanded** (rows
  re-emitted, not re-written as `-1`), matching the existing non-parametric DRT
  behavior; data round-trips, the `-1` syntax does not;
- **activated SPREAD (`NR<0`) parameter recipients are structural round-trip
  only, not execution-guaranteed** (review): USG-T 2.7 copies `DRTF` but not
  `NodDRT` on activation, so FloPy preserves the recipients on load/write but
  does not claim every activated-SPREAD case runs. Inline (`NR>0`) recipients
  travel in `DRTF` and are execution-safe.

Roadmap label: **`✅ Full (authoring) / Parameter-preserving (NPDRT, incl.
activations + recipients) / Expanded valid write (list controls)`**.

## Review follow-up (executed)

Three fixes from review:

1. **P1 — MXADRT under-counted active parameters.** `write_file` set `MXADRT` to
   the max non-parametric row count only; with active parameters the Fortran
   would abort (`NDRTCL = NDRTNP + Σ nlst > MXADRT`). Now `MXADRT =
   _max_active_drains()` includes the active-parameter rows per period (with
   `ITMP<0` reuse carrying the previous non-parametric count). The roundtrip test
   now asserts `MXADRT=2` (1 non-parametric + 1 active `NLST=1`), and a new test
   covers two active params (`NLST` 2 + 1) ⇒ `MXADRT=4`.
2. **P2 — SPREAD-in-parameter is structural-only.** Documented above and in the
   test docstring: activated `NR<0` recipients round-trip but are not
   execution-guaranteed (no `NodDRT` copy on activation).
3. **P3 — stale roadmap.** `USGT_roadmap.md` Gap §7 (and the DRT row) no longer
   say `NPDRT>0` fails explicitly.

Second review follow-up (Stage 4.4E review): active-parameter names are now
resolved **case-insensitively** in `_max_active_drains` (via the shared
`resolve_list_parameter` helper), matching the validation and the Fortran
(`SGWF2DRT8LS` UPCASEs both names). Previously a definition `dp` activated as
`DP` passed validation but its `NLST` was not added to `MXADRT`, under-sizing the
header. Also, activating the same parameter more than once in a stress period
(case-insensitively) now raises `ValueError` ("already activated" in the
Fortran), before the file is opened.

## Tests (`-k mfusgdrt`, 24 passed)

New (Stage 4.4D + review follow-up):

- `test_mfusgdrt_parameterized_roundtrip` — item-1 `NPDRT/MXL` + definition with
  an inline RETURNFLOW recipient + per-SP activation; 0-based / 1-based; asserts
  `MXADRT=2` (non-parametric + active rows).
- `test_mfusgdrt_parameter_mxadrt_counts_active_rows` — two active params
  (`NLST` 2 + 1) + 1 non-parametric drain ⇒ written `MXADRT=4` (review P1).
- `test_mfusgdrt_parameter_spread_recipients` — a definition row keeps SPREAD
  (multi-node `U1DINT`) recipients + CHANGEC + AUX; docstring states the
  structural-round-trip-only contract for activated SPREAD (review P2).
- `test_mfusgdrt_parameter_mixed_and_active` — `ITMP NP` with non-parametric
  drains and an active parameter in one period.
- `test_mfusgdrt_parameter_reuse_with_active` — `ITMP<0` reuse + active parameter.
- `test_mfusgdrt_parameter_sfac_scales_cond` — SFAC in a parameter block scales
  COND.
- `test_mfusgdrt_parameter_instances_unsupported` — `INSTANCES` →
  `NotImplementedError`.
- `test_mfusgdrt_parameter_active_without_defs_fails` — active params without
  definitions → `ValueError`, no partial file. (Repurposed from the old
  `_from_scratch_fails`: from-scratch authoring *with* definitions is now
  supported — Stage 4.6B; see `USGT_STAGE4_10_DRT_PARAMETER_AUTHORING.md`.)
- `test_mfusgdrt_parameter_mxl_too_small_fails` — `MXL <` total →
  `ValueError`, no partial file.
- `test_mfusgdrt_parameter_inconsistent_fails` — `len(data) != nlst` →
  `ValueError`, no partial file.
- `test_mfusgdrt_parameter_active_undefined_fails` — active name not defined →
  `ValueError`, no partial file.

Second review follow-up:

- `test_mfusgdrt_parameter_active_case_insensitive` — definition `dp`, activation
  `DP` ⇒ `MXADRT=2`.
- `test_mfusgdrt_parameter_duplicate_active_fails` — `dp` + `DP` in one period →
  `ValueError`, no partial file.

Unchanged (non-parametric regression): the existing DRT tests
(`inline_single_recipient`, `changec_concentration`, `spread_multi_node`,
`stress_period_reuse`, `nam_registry`, `structured_construction_raises`,
`load_structured_delegates_to_base`, `sfac_single_and_spread`,
`recipient_count_mismatch_fails`, `zero_recipients_when_omitted`,
`recipient_constant_u1dint`). The replaced `test_mfusgdrt_parameters_fail_explicitly`
is superseded (DRT no longer fails on parameterized load; the explicit-failure
contract is in the negative tests above).

No from-scratch parametric **executable** smoke was added (a convergent
parameterized DRT model is not cheap); the writer is audited against the Fortran
and round-trips in FloPy. Executable validation remains a documented manual tier.

## Validation

```bash
python -m pytest autotest/test_usg_transport.py -k mfusgdrt -q
python -m pytest autotest/test_usg_transport.py -k "mfusgdrt or usgt_list" -q
python -m pytest autotest/test_usg_transport.py -q
python -m pytest autotest/test_usg_transport_exe.py -q
USGT_EXE=/Users/roh.sgm/Documents/GitHub_Projects/GW-Software-Compiled-via-Claude/gfortran/usgt_2.7/usgt_270_arm python -m pytest autotest/test_usg_transport.py autotest/test_usg_transport_exe.py -q
git diff --check
git status --short
```

Results (after the second review follow-up): `-k mfusgdrt` **24 passed**;
`-k "mfusgqrt or mfusgdrt or usgt_list"` **47 passed**; focused **217 passed**;
exe **4 passed**; combined **221 passed** under the USG-T 2.7 ARM binary.
