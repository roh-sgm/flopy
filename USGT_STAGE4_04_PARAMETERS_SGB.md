# Stage 4.4C — MODFLOW Parameter Preservation: SGB (executed; review follow-up)

Scope: implement load → write → reload **preservation** of SGB list parameters
(`NPSGB > 0`), audited against `glo2sgbu1.f` and `parutl7.f`, without breaking the
existing non-parametric authoring. Only SGB and the shared parameter helper are
touched; ETS/HFB and DRT/QRT are not.

**Review follow-up outcome:** the review found that USG-T 2.7 cannot actually
*activate* an SGB parameter (a `PARTYP='SGB'` definition vs a `PTYP='G'`
activation is a type conflict that aborts the run). SGB is therefore
**definition-preserving only** — parameter definitions round-trip, but active
parameters (`NP>0`) and `INSTANCES` raise `NotImplementedError`, and the writer
is hardened (`MXS` validated; no partial files). See the critical finding below.

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

### Critical finding (review): active SGB parameters abort USG-T 2.7

The parameter **type** required by the definition and by the activation are
**contradictory**:

- the definition is read with `UPARLSTRP(...,'SGB','SGB',1,...)`
  (`glo2sgbu1.f:97`), and `UPARLSTRP` enforces `PARTYP == PTYPX`
  (`parutl7.f:684`, "Parameter type must be: SGB") → a valid definition needs
  `PARTYP='SGB'`;
- the activation is read with `UPARLSTSUB(IN,'SGB',IOUTU,'G',...)`
  (`glo2sgbu1.f:185`), and `UPARLSTSUB` enforces `PARTYP == PTYP`
  (`parutl7.f:800`, "Parameter type conflict") → activating that parameter needs
  `PARTYP='G'`.

A parameter has a single `PARTYP`, so it cannot be both `'SGB'` and `'G'`. Any
**active** SGB parameter (`NP>0` in a stress period) therefore trips the
type-conflict guard in `UPARLSTSUB` and **aborts the run** (`USTOP`). Defining
SGB parameters (`NPSGB>0`) is fine as long as they are never activated; an
activation is not executable in USG-T 2.7. This is a bug in the official USG-T
2.7 source, not a FloPy limitation.

## Design — shared helper

`flopy/mfusg/_usgt_parameters.py` gains the leading-`PARAMETER` record helpers
(`UPARLSTAL`): `read_list_parameter_count` / `write_list_parameter_count`. The
existing `read_list_parameter_header` / `write_list_parameter_header` (the
`PARNAM PARTYP PARVAL NLST` definition header) are reused unchanged. The `NLST`
`NODE GRADIENT [aux]` rows stay in SGB (`_read_sgb_rows` / `_write_sgb_rows`, with
0-based↔1-based conversion and SFAC consumed via the shared
`_usgt_list.begin_list_block`). The activation helpers
(`read_active_list_parameters` / `write_active_list_parameters`) are **not** used
by SGB, because SGB activations are unsupported (see the critical finding above).

## SGB implementation (`MfUsgSgb`) — definition-preserving only

- `load`: consume the optional `PARAMETER NPSGB MXS`; when `NPSGB>0`, read each
  definition (`read_list_parameter_header` + `NLST` rows via `_read_sgb_rows`,
  0-based) into `self.parameters = {name: {"partyp","parval","nlst","data"}}` and
  keep `MXS` as `self.mxs`. Per stress period, keep the non-parametric rows
  (`stress_period_data`, `ITMP<0` reuse). A per-SP `NP>0` (active parameter)
  raises `NotImplementedError` (the type-conflict above). AUX preserved.
- `write_file`: when definitions are present, emit `PARAMETER NPSGB MXS`, the
  header, the definitions + rows, then per SP `ITMP 0` + the non-parametric rows
  (NP is always 0; active parameters are never written). Validate the parameter
  state up front (before opening the file). Non-parametric SGB (authoring,
  `ITMP<0` reuse, AUX, list controls) is unchanged.
- Indexing: internal 0-based, file 1-based, for both definition rows and
  non-parametric rows.

### Validation / explicit failures

- **Active parameters (`NP>0`)** → `NotImplementedError` on load and on write
  (the `PARTYP='SGB'` vs `'G'` mismatch aborts USG-T 2.7). The message cites
  `glo2sgbu1.f:97/185` and `parutl7.f:684/800`.
- **`INSTANCES` (`NUMINST>0`)** → `NotImplementedError` on load (a consequence:
  instanced definitions could only ever be activated, which is unsupported).
- **`MXS` hardening** → `ValueError` when definitions are present and `MXS<=0`
  (e.g. a from-scratch `PARAMETER 1 0`) or `MXS <` the total definition rows.
- **Inconsistent definitions** → `ValueError` (a definition missing
  `partyp`/`parval`/`nlst`/`data`, or `len(data) != nlst`).
- All validation runs **before** the file is opened, so no partial file is
  written.

## Decision: honest status — definition-preserving only

SGB stays **`✅ Full (authoring)` for the non-parametric package**. On the
parametric axis it is **definition-preserving only**: `NPSGB>0` parameter
definitions (with **no** activations) round-trip (load → write → reload), but
**active SGB parameters are unsupported** because USG-T 2.7 itself rejects them
(`PARTYP='SGB'` definition vs `PTYP='G'` activation → "Parameter type conflict",
`USTOP`). It is **not** "Full" on the parametric axis:

- active parameters (`NP>0`) and `INSTANCES` raise `NotImplementedError`;
- the `ITMP<0` reuse of non-parametric rows is written **expanded** (rows
  re-emitted, not re-written as `-1`), matching the existing non-parametric SGB
  behavior; data round-trips, the `-1` syntax does not.

Roadmap label: **`✅ Full (authoring) / Parameter-definition-preserving (NPSGB,
no activations) / Expanded valid write (list controls)`**.

## Tests (`-k mfusgsgb`, 15 passed)

Definition-preserving (new):

- `test_mfusgsgb_definitions_preserved_roundtrip` — `PARAMETER` record +
  definition rows; 0-based internal / 1-based file; no activations.
- `test_mfusgsgb_definitions_aux` — AUX on parameter definition rows.
- `test_mfusgsgb_definitions_sfac_inert_on_gradient` — SFAC consumed, gradient
  unscaled.
- `test_mfusgsgb_definitions_open_close` — definition rows via `OPEN/CLOSE`.

Explicit failures (new):

- `test_mfusgsgb_active_parameters_unsupported` — a file with `NP>0` →
  `NotImplementedError` (the `SGB`/`G` mismatch).
- `test_mfusgsgb_parameter_instances_unsupported` — `INSTANCES` →
  `NotImplementedError`.
- `test_mfusgsgb_parameter_mxs_zero_fails` — definitions + `MXS=0`
  (`PARAMETER 1 0`) → `ValueError`, no partial file.
- `test_mfusgsgb_parameter_mxs_too_small_fails` — `MXS <` total rows →
  `ValueError`.
- `test_mfusgsgb_parameter_inconsistent_fails` — `len(data) != nlst` →
  `ValueError`, no partial file.
- `test_mfusgsgb_active_params_constructed_fails` — manually constructed
  activations → `NotImplementedError`, no partial file.

Unchanged (non-parametric regression): `test_mfusgsgb_roundtrip` (incl. `ITMP<0`
reuse), `test_mfusgsgb_aux_roundtrip`, `test_mfusgsgb_programmatic_authoring`,
`test_mfusgsgb_nam_registry`, `test_mfusgsgb_sfac_and_open_close`. The earlier
Stage 4.4C tests that asserted active-parameter preservation are superseded.

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

Results (after the review follow-up): `-k mfusgsgb` **15 passed**;
`-k "mfusgsgb or usgt_list"` **18 passed**; focused **193 passed**; exe
**4 passed**; combined **197 passed** under the USG-T 2.7 ARM binary.
