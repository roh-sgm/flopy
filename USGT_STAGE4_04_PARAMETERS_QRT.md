# Stage 4.4E — MODFLOW Parameter Preservation: QRT (executed)

Scope: implement load → write → reload **preservation** of QRT list parameters
(`NPQRT > 0`), **only after** auditing the Fortran (after SGB's `PARTYP` mismatch
and DRT's `MXADRT`/SPREAD lessons). Only `MfUsgQrt` and shared docs are touched;
DRT/SGB/HFB/ETS are not.

## Fortran audit (mandatory)

### `gwf2QRT8u.f`

- **Item 1 (AR, `GWF2QRT8U1AR`)**: `MXAQRT MXRTCELLS IQRTCB NPQRT MXL [options]`
  read directly (`:56`/`:60-61`) — `NPQRT`/`MXL` are in item 1, no separate
  `PARAMETER` line. Options: `AUX`, `RETURNFLOW` (`IQRTFL=2`), `AUTOFLOWREDUCE`,
  `CHANGEC` (`IQCHANGEC`), `NOPRINT`, `IUNIT_AFR_QRT`, `TRANSIENTQ`.
  `NQRTVL = 5 + NAUX + IQRTFL`.
- **Definitions (items 2-3)**: for `K=1..NPQRT`,
  `UPARLSTRP(...,'QRT','QRT',1,NUMINST)` (`:183`). Each definition's `NLST` rows
  are read by `SGWF2QRT8LR`, the same reader as non-parametric rows; recipient
  `U1DINT` blocks follow **after all the rows** of the block (`:1050-1070`), in
  sink order — not interleaved per-row like DRT.
- **Per stress period (RP, `GWF2QRT8U1RP`)**: `ITMP NP` when `NPQRT>0`, else just
  `ITMP`. `ITMP<0` reuses non-parametric sinks; `ITMP>0` reads them. Then, if
  `NP>0`, `NP` activations via `SGWF2QRT8LS` (`:330`).

### Critical checks (the SGB/DRT questions)

1. **`PARTYP` definition vs activation:** definition `UPARLSTRP(...,'QRT','QRT')`
   (`:183`) and activation `SGWF2QRT8LS` sets `PTYP='QRT'` (`:1118`) and enforces
   `PARTYP(IP)==PTYP` (`:1143`). **Both `'QRT'` — consistent.** So an active QRT
   parameter is *type-valid* (like DRT, unlike SGB).

2. **`PARVAL` scaling — Fortran bug.** `QRTF` layout (with `RETURNFLOW`):
   `QRTF(1)=node`, `QRTF(4)=Q`, `QRTF(5)=NumRT` (recipient count, set at `:854`),
   `QRTF(6)=Rfprop`. `SGWF2QRT8LS` scales with `IPVL1=IPVL2=5`
   (`:1119-1120`), i.e. it multiplies `QRTF(5)=NumRT` by the parameter value —
   **not `Q`**. `SFAC`, by contrast, uses `ISCLOC=4` (`:746-747`) and correctly
   scales `Q`. So the QRT parameter value scales the *recipient count*, which is
   a bug: with `PARVAL != 1.0` it corrupts `NumRT`. **FloPy does not apply the
   parameter value to `Q`.**

3. **`MXAQRT`** = `IQRTPB-1`; `SGWF2QRT8LS` does `NQRTCL = NQRTCL + NLST` and
   aborts if `NQRTCL > MXAQRT` (`:1183-1184`). So `MXAQRT` must cover the active
   total per period (non-parametric + active-parameter rows), like DRT's
   `MXADRT`.

4. **`MXRTCELLS` / `NodQRT`.** Recipients are read into `NodQRT` after the rows.
   `SGWF2QRT8LS` (`:1196-1240`) copies `QRTF` on activation but **not** `NodQRT`,
   so an *activated* parameter's RETURNFLOW recipients are not resolved at run
   time (same as DRT).

5. **`INSTANCES`**: `UPARLSTRP`/`UINSRP` support `NUMINST>0`.

6. **`TRANSIENTQ`**: a separate transient-flow time-series option, out of scope.

## Decision: honest status — structural-preserving, not execution-guaranteed

QRT stays **`✅ Full (authoring)`** for the non-parametric package. On the
parametric axis it is **structurally preserved (load → write → reload), not
execution-guaranteed**: `NPQRT>0` definitions, their recipient blocks, and the
per-SP activations all round-trip faithfully, but FloPy does **not** claim an
activated QRT parameter runs in USG-T 2.7, because of two Fortran issues
(check 2 and check 4 above): the parameter value scales `NumRT` not `Q`, and
`NodQRT` is not copied on activation. It is **not** "Full" on the parametric
axis. Other gaps: from-scratch authoring and `INSTANCES` raise
`NotImplementedError`; the `ITMP<0` reuse of non-parametric sinks is written
expanded; `TRANSIENTQ` still fails explicitly.

Roadmap label: **`✅ Full (authoring) / Parameter structural-preserving (NPQRT,
not execution-guaranteed) / Expanded valid write (list controls)`**.

## QRT implementation (`MfUsgQrt`)

- `load`: parse `NPQRT`/`MXL` from item 1; when `NPQRT>0`, read each definition
  (`read_list_parameter_header` + `NLST` rows + recipient blocks via the new
  `_read_sink_rows`, 0-based) into `self.parameters = {name: {"partyp","parval",
  "nlst","data","recipient_nodes"}}` and keep `MXL` as `self.mxl`. Per SP read
  `ITMP NP` (`NP` only when `NPQRT>0`), the non-parametric sinks (`ITMP<0`
  reuse), then `NP` activation names into `self.active_params`.
- `write_file`: when preserving, write `NPQRT MXL` in item 1, the definition
  blocks (header + rows + recipient blocks), then per SP `ITMP NP` + the
  non-parametric sinks + the activation names. `MXAQRT = _max_active_sinks()`
  (non-parametric + active-parameter rows per period); `MXRTCELLS =
  _max_rt_cells()` (max non-parametric recipients per period and the
  parameter-definition recipients). Validate first (no partial file).
  Non-parametric QRT (authoring, `ITMP<0` reuse, RETURNFLOW, CHANGEC, AUX, SFAC,
  AUTOFLOWREDUCE) is unchanged.
- Indexing: internal 0-based, file 1-based, for definition rows (and recipients)
  and non-parametric rows.

### Validation / explicit failures

- **From-scratch authoring** (active parameters with no loaded definitions) →
  `NotImplementedError`.
- **`INSTANCES` (`NUMINST>0`)** → `NotImplementedError` on load.
- **`MXL<=0` or `MXL <` total definition rows** → `ValueError`.
- **Inconsistent definitions** → `ValueError` (missing keys, `len(data)!=nlst`,
  `len(recipient_nodes)!=nlst`, or an active name not defined).
- **`TRANSIENTQ`** → `NotImplementedError`.
- All validation runs before the file is opened — no partial file.

## Tests (`-k mfusgqrt`, 18 passed)

New (Stage 4.4E):

- `test_mfusgqrt_parameterized_roundtrip` — item-1 `NPQRT/MXL` + definition (no
  RETURNFLOW) + per-SP activation; `MXAQRT=2`; 0-based / 1-based.
- `test_mfusgqrt_parameter_returnflow_recipients` — definition row keeps its
  RETURNFLOW recipients (`U1DINT`) + CHANGEC + AUX; `MXRTCELLS` reflects them.
- `test_mfusgqrt_parameter_mxaqrt_counts_active_rows` — two active params
  (`NLST` 2 + 1) + 1 non-parametric ⇒ `MXAQRT=4`.
- `test_mfusgqrt_parameter_reuse_with_active` — `ITMP<0` reuse + active parameter.
- `test_mfusgqrt_parameter_sfac_scales_q` — SFAC scales `Q` (ISCLOC=4), unlike
  the parameter value.
- `test_mfusgqrt_parameter_instances_unsupported` — `INSTANCES` →
  `NotImplementedError`.
- `test_mfusgqrt_parameter_from_scratch_fails` — active params without defs →
  `NotImplementedError`, no partial file.
- `test_mfusgqrt_parameter_mxl_too_small_fails` — `MXL <` total → `ValueError`.
- `test_mfusgqrt_parameter_inconsistent_fails` — `len(data)!=nlst` and
  `len(recipient_nodes)!=nlst` → `ValueError`, no partial file.
- `test_mfusgqrt_parameter_active_undefined_fails` — active name not defined →
  `ValueError`.

Unchanged (non-parametric regression): `minimal_authoring`,
`returnflow_concentration`, `multi_recipient_and_pure_sink`, `nam_registry`,
`sfac_with_recipients_expanded_write`, `recipient_count_mismatch_fails`,
`zero_recipients_when_omitted`, and `unsupported_modes_fail_explicitly` (now
asserts only `TRANSIENTQ`; the `NPQRT>0`-fails case is superseded).

No from-scratch parametric **executable** smoke was added — and, given the
`PARVAL→NumRT` Fortran bug and the un-copied `NodQRT`, an executable test of an
*activated* QRT parameter would not be meaningful. The writer is audited against
the Fortran and round-trips in FloPy.

## Validation

```bash
python -m pytest autotest/test_usg_transport.py -k mfusgqrt -q
python -m pytest autotest/test_usg_transport.py -k "mfusgqrt or usgt_list" -q
python -m pytest autotest/test_usg_transport.py -q
python -m pytest autotest/test_usg_transport_exe.py -q
USGT_EXE=/Users/roh.sgm/Documents/GitHub_Projects/GW-Software-Compiled-via-Claude/gfortran/usgt_2.7/usgt_270_arm python -m pytest autotest/test_usg_transport.py autotest/test_usg_transport_exe.py -q
git diff --check
git status --short
```

Results: `-k mfusgqrt` **18 passed**; `-k "mfusgqrt or usgt_list"` **21 passed**;
focused **213 passed**; exe **4 passed**; combined **217 passed** under the USG-T
2.7 ARM binary.
