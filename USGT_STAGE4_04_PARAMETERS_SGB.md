# Stage 4.4C — MODFLOW Parameter Preservation: SGB (executed)

Scope: implement load → write → reload **preservation** of SGB list parameters
(`NPSGB > 0`), audited against `glo2sgbu1.f` and `parutl7.f`, without breaking the
existing non-parametric authoring. Only SGB and the shared parameter helper are
touched; ETS/HFB and DRT/QRT are not.

## Fortran audit

### `glo2sgbu1.f`

- **Item 1a (AR, `GLO2SGBU1AR`)**: `URDCOM` reads comments, then
  `UPARLSTAL(IN,IOUT,LINE,NPSGB,MXS)` consumes an optional leading
  `PARAMETER NPSGB MXS` record (and reads the next line as the header). `NPSGB`
  is the number of named parameters; `MXS` the max parameter list entries.
- **Item 1 (header)**: `MXACTS ISGBCB [AUX <name> ...] [NOPRINT]`.
  `NSGBVL = 5 + NAUX` (node + gradient + a CBC slot + 2 structured dummy slots,
  plus AUX).
- **Items 2-3 (parameter definitions)**: for `K=1..NPSGB`,
  `UPARLSTRP(...,'SGB','SGB',1,NUMINST)` reads
  `PARNAM PARTYP PARVAL NLST [INSTANCES n]` (`ITERP=1`). SGB **does support
  `INSTANCES`** (`NUMINST>0` reads `NUMINST` instance blocks, each
  `UINSRP`-named with `NLST/NUMINST` rows). The `NLST` rows are read with
  `ULSTRDU` (`NODE GRADIENT [aux ...]`), which also accepts leading `SFAC` /
  `OPEN/CLOSE` / `EXTERNAL` list controls.
- **Per stress period (RP, `GLO2SGBU1RP`)**: `ITMP NP`. `ITMP<0` reuses the
  previous period's non-parametric rows; `ITMP>0` reads `NNPSGB=ITMP`
  non-parametric rows (`ULSTRDU`). Then `PRESET('Q')` and, if `NP>0`, `NP`
  `UPARLSTSUB` activations (each a parameter name, plus an instance name when the
  parameter is time-varying).
- **Scaling (`ISCLOC1=ISCLOC2=2`)**: both `SFAC` (in `ULSTRDU`) and the parameter
  value (in `UPARLSTSUB`) scale **internal column 2**, which for SGB is a *dummy*
  column — **not** the gradient (`SGB(4)`). So SFAC and the parameter value are
  *inert on the gradient* that FloPy stores; this matches the existing
  non-parametric SFAC behavior.

## Design — shared helper

`flopy/mfusg/_usgt_parameters.py` gains the leading-`PARAMETER` record helpers
(`UPARLSTAL`): `read_list_parameter_count` / `write_list_parameter_count`. The
existing list-parameter helpers are reused unchanged:
`read_list_parameter_header` / `write_list_parameter_header` (the
`PARNAM PARTYP PARVAL NLST` definition header) and `read_active_list_parameters`
/ `write_active_list_parameters` (the activation names). The `NLST`
`NODE GRADIENT [aux]` rows stay in SGB (`_read_sgb_rows` / `_write_sgb_rows`, with
0-based↔1-based conversion and SFAC consumed via the shared
`_usgt_list.begin_list_block`).

## SGB implementation (`MfUsgSgb`)

- `load`: consume the optional `PARAMETER NPSGB MXS`; when `NPSGB>0`, read each
  definition (`read_list_parameter_header` + `NLST` rows via `_read_sgb_rows`,
  0-based) into `self.parameters = {name: {"partyp","parval","nlst","data"}}` and
  keep `MXS` as `self.mxs`. Per stress period, keep the non-parametric rows
  (`stress_period_data`, `ITMP<0` reuse) and read `NP` activation names into
  `self.active_params = {kper: [name, ...]}`. AUX preserved.
- `write_file`: when preserving, emit `PARAMETER NPSGB MXS`, the header, the
  definitions + rows, then per SP `ITMP NP` + the non-parametric rows + the
  activation names. Validate the parameter state up front (before opening the
  file). Non-parametric SGB (authoring, `ITMP<0` reuse, AUX, list controls) is
  unchanged.
- Indexing: internal 0-based, file 1-based, for both definition rows and
  non-parametric rows.

### Validation / explicit failures

- `INSTANCES` (`NUMINST>0`): the Fortran supports them, but FloPy does not yet —
  `load` raises `NotImplementedError`.
- From-scratch parameter authoring (active parameters with no loaded
  definitions): `write_file` raises `NotImplementedError`.
- Inconsistent preserved state (`load`: `NP>0` without definitions; `write_file`:
  a definition missing `partyp`/`parval`/`nlst`/`data`, `len(data) != nlst`, or
  an active name not defined case-insensitively): `ValueError`, with no partial
  file written (validation runs before the file is opened).

## Decision: honest status

SGB stays **`✅ Full (authoring)` for the non-parametric package** and is now
**parameter-preserving for `NPSGB`** (load → write → reload). It is **not**
promoted to "Full" on the parametric axis, because:

- from-scratch parameter authoring is unsupported (`NotImplementedError`);
- parameter `INSTANCES` are unsupported (`NotImplementedError`), although the
  Fortran supports them;
- the `ITMP<0` reuse of non-parametric rows is written **expanded** (the rows are
  re-emitted, not re-written as `-1`), matching the existing non-parametric SGB
  behavior; data round-trips, the `-1` syntax does not.

Roadmap label: **`✅ Full (authoring) / Parameter-preserving (NPSGB) / Expanded
valid write (list controls)`**.

## Tests (`-k mfusgsgb`, 14 passed)

New (Stage 4.4C):

- `test_mfusgsgb_parameterized_roundtrip` — PARAMETER record + definition + per-SP
  activation; 0-based internal / 1-based file.
- `test_mfusgsgb_parameter_mixed_with_nonparam_and_active` — `ITMP NP` with both
  non-parametric rows and an active parameter in one stress period.
- `test_mfusgsgb_parameter_reuse_with_active` — `ITMP<0` reuse + active parameter.
- `test_mfusgsgb_parameter_aux` — AUX on parameter definition rows.
- `test_mfusgsgb_parameter_sfac_inert_on_gradient` — SFAC consumed, gradient
  unscaled.
- `test_mfusgsgb_parameter_open_close` — definition rows via `OPEN/CLOSE`.
- `test_mfusgsgb_parameter_instances_unsupported` — `INSTANCES` →
  `NotImplementedError`.
- `test_mfusgsgb_parameter_from_scratch_fails` — active params without
  definitions → `NotImplementedError`, no partial file.
- `test_mfusgsgb_parameter_inconsistent_fails` — `len(data) != nlst` →
  `ValueError`, no partial file.

Unchanged (non-parametric regression): `test_mfusgsgb_roundtrip` (incl. `ITMP<0`
reuse), `test_mfusgsgb_aux_roundtrip`, `test_mfusgsgb_programmatic_authoring`,
`test_mfusgsgb_nam_registry`, `test_mfusgsgb_sfac_and_open_close`. The replaced
`test_mfusgsgb_parameters_fail_explicitly` is superseded (SGB no longer fails on
parameterized load; the explicit-failure contract is in the negative tests
above).

No from-scratch parametric **executable** smoke was added (a convergent
parameterized SGB model is not cheap); the writer is audited line-by-line against
the Fortran and round-trips in FloPy. Executable validation remains a documented
manual tier.

## Validation

```bash
python -m pytest autotest/test_usg_transport.py -k mfusgsgb -q
python -m pytest autotest/test_usg_transport.py -k "mfusgsgb or usgt_list" -q
python -m pytest autotest/test_usg_transport.py -q
python -m pytest autotest/test_usg_transport_exe.py -q
USGT_EXE=/Users/roh.sgm/Documents/GitHub_Projects/GW-Software-Compiled-via-Claude/gfortran/usgt_2.7/usgt_270_arm python -m pytest autotest/test_usg_transport.py autotest/test_usg_transport_exe.py -q
git diff --check
git status --short
```

Results: `-k mfusgsgb` **14 passed**; `-k "mfusgsgb or usgt_list"` **17 passed**;
focused **192 passed**; exe **4 passed**; combined **196 passed** under the USG-T
2.7 ARM binary.
