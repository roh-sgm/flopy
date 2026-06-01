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
- Indexing: internal 0-based, file 1-based, for both definition rows (and their
  recipients) and non-parametric rows.

### Validation / explicit failures

- **From-scratch authoring** (active parameters with no loaded definitions) →
  `NotImplementedError`.
- **`INSTANCES` (`NUMINST>0`)** → `NotImplementedError` on load (instances
  combined with per-row recipient lists are not implemented yet).
- **`MXL` hardening** → `ValueError` when definitions are present and `MXL<=0` or
  `MXL <` the total definition rows.
- **Inconsistent definitions** → `ValueError` (missing keys, `len(data)!=nlst`,
  `len(recipient_nodes)!=nlst`, or an active name not defined).
- All validation runs before the file is opened — no partial file.

## Decision: honest status

DRT stays **`✅ Full (authoring)`** for the non-parametric package and is now
**parameter-preserving for `NPDRT` (definitions + activations + recipients)** —
DRT is internally consistent in the Fortran, so active parameters are real. It is
**not** "Full" on the parametric axis because:

- from-scratch parameter authoring is unsupported (`NotImplementedError`);
- parameter `INSTANCES` are unsupported (`NotImplementedError`), although the
  Fortran supports them;
- the `ITMP<0` reuse of non-parametric drains is written **expanded** (rows
  re-emitted, not re-written as `-1`), matching the existing non-parametric DRT
  behavior; data round-trips, the `-1` syntax does not.

Roadmap label: **`✅ Full (authoring) / Parameter-preserving (NPDRT, incl.
activations + recipients) / Expanded valid write (list controls)`**.

## Tests (`-k mfusgdrt`, 21 passed)

New (Stage 4.4D):

- `test_mfusgdrt_parameterized_roundtrip` — item-1 `NPDRT/MXL` + definition with
  an inline RETURNFLOW recipient + per-SP activation; 0-based / 1-based.
- `test_mfusgdrt_parameter_spread_recipients` — a definition row keeps SPREAD
  (multi-node `U1DINT`) recipients + CHANGEC + AUX.
- `test_mfusgdrt_parameter_mixed_and_active` — `ITMP NP` with non-parametric
  drains and an active parameter in one period.
- `test_mfusgdrt_parameter_reuse_with_active` — `ITMP<0` reuse + active parameter.
- `test_mfusgdrt_parameter_sfac_scales_cond` — SFAC in a parameter block scales
  COND.
- `test_mfusgdrt_parameter_instances_unsupported` — `INSTANCES` →
  `NotImplementedError`.
- `test_mfusgdrt_parameter_from_scratch_fails` — active params without
  definitions → `NotImplementedError`, no partial file.
- `test_mfusgdrt_parameter_mxl_too_small_fails` — `MXL <` total →
  `ValueError`, no partial file.
- `test_mfusgdrt_parameter_inconsistent_fails` — `len(data) != nlst` →
  `ValueError`, no partial file.
- `test_mfusgdrt_parameter_active_undefined_fails` — active name not defined →
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

Results: `-k mfusgdrt` **21 passed**; `-k "mfusgdrt or usgt_list"` **24 passed**;
focused **202 passed**; exe **4 passed**; combined **206 passed** under the USG-T
2.7 ARM binary.
