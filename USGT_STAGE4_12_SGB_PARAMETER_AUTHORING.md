# Stage 4.6C-B — SGB from-scratch parameter-*definition* authoring (`NPSGB>0`)

Date: 2026-06-02 · Base: `develop` @ `690f79a7`

Goal: let a user build SGB parameter **definitions** from Python (`NPSGB>0`) and
write a valid file, with no prior load — reusing the DRT 4.6B / HFB 4.6C-A
pattern. **Active SGB parameters remain unsupported** (USG-T 2.7 reads SGB
definitions with `PARTYP='SGB'` but activates them with `PTYP='G'`, a type
conflict that aborts the run), so this is *definition authoring only*: every
stress period is written with `NP=0` and any activation raises
`NotImplementedError`. QRT/DRT/HFB/ETS are not touched.

## Audit

The on-file format for SGB parameter definitions is identical whether they came
from `load` or the user — the writers (`write_list_parameter_count` →
`PARAMETER NPSGB MXS`; `write_list_parameter_header` → `UPARLSTRP`; the row
writer) were already Fortran-audited and round-tripped (Stage 4.4C). So this is
a FloPy-side ergonomics + validation change in `mfusgsgb.py`.

Fortran (`glo2sgbu1.f`): the optional leading `PARAMETER NPSGB MXS` (`UPARLSTAL`)
dimensions the parameter list; definitions are read with `UPARLSTRP`
(`PARTYP='SGB'`, glo2sgbu1.f:97). Activation (`UPARLSTSUB`, glo2sgbu1.f:185) uses
`PTYP='G'`; the mismatch (`parutl7.f:684/800`) aborts the run — hence active SGB
params are out of scope by design, not by FloPy limitation. `MXS` must be ≥ the
total of all definition `NLST`, so `mxs = Σ nlst` is the minimal valid value.

## What changed (code)

Only `flopy/mfusg/mfusgsgb.py` (+ the shared `_usgt_parameters` helpers and a
docstring cleanup there):

- new **`_normalize_param(name, pdef)`** — validates + canonicalizes one
  definition (`partyp` defaults to / must be `SGB`, original casing preserved;
  `data` recarray or array-like built with the active SGB dtype, incl. AUX
  columns when present; nodes non-negative 0-based; `nlst` computed/validated)
  via the shared `check_parameter_name` / `check_parval`;
- **`_validate_parameter_write`** now normalizes every definition, rejects
  activations (`NotImplementedError`), checks case-insensitive duplicate names,
  auto-computes `MXS = Σ nlst` (validated `≥` total if given), and returns
  `(params, mxs)`;
- `write_file` uses the canonical `(params, mxs)` and writes `NP=0` every period.
  An empty period now writes `ITMP=0` (zero rows) until some non-param rows have
  appeared (then `ITMP=-1` reuse) — a definition-only file therefore reloads to a
  non-empty stress-period dict instead of crashing `MfList`.

Internal nodes stay 0-based; the file is 1-based.

## API

```python
sgb = MfUsgSgb(
    mfusg_model,
    parameters={
        "gp": {
            "parval": "1.5",                 # scales the SGB gradient slot
            "data": [(4, 0.01), (9, 0.02)],  # 0-based (node, gradient[, aux...])
            # "partyp": "SGB" (default), "nlst": 2 (auto from data)
        },
    },
    # "mxs": auto = sum(nlst) when omitted
)
sgb.write_file()    # PARAMETER NPSGB MXS + UPARLSTRP defs; NP=0 every period
```

AUX columns: pass a `dtype` that carries them (e.g.
`Package.add_to_dtype(MfUsgSgb.get_default_dtype(), ["C01"], np.float64)`).

## Validations (all before the file is opened — no partial file)

| Condition | Result |
|---|---|
| `active_params` / per-SP `NP>0` / loaded activation | `NotImplementedError` (Fortran `SGB`/`G` conflict) |
| name not a single whitespace-free token, or > 10 chars | `ValueError` |
| duplicate definition name (case-insensitive) | `ValueError` |
| `partyp` not `SGB` | `ValueError` |
| `parval` missing/blank or a multi-token string (`"1 2"`) | `ValueError` |
| `data` missing/empty | `ValueError` |
| `nlst` given ≠ `len(data)` | `ValueError` |
| explicit `mxs` `<` total definition rows | `ValueError` |
| negative node | `ValueError` |
| parameter `INSTANCES` | `NotImplementedError` (load) |

## Tests (`autotest/test_usg_transport.py`)

- `test_mfusgsgb_parameter_authoring_from_scratch` — ergonomic input, auto
  `NPSGB`/`MXS`, per-SP `NP=0`, 0↔1-based, reload.
- `test_mfusgsgb_parameter_authoring_with_aux` — AUX concentration column.
- `test_mfusgsgb_parameter_authoring_mixed_with_nonparam` — defs + non-param rows.
- `test_mfusgsgb_parameter_authoring_auto_counts` — `NPSGB`/`MXS` auto.
- Negatives: `_duplicate_names_fails`, `_bad_name_fails`, `_bad_parval_fails`,
  `_empty_data_fails`, `_negative_node_fails` — each `ValueError`, no partial file.
- `test_mfusgsgb_parameter_mxs_zero_auto_computes` — repurposed from the old
  `_mxs_zero_fails` (MXS=0 now auto-computes). The active-param negatives
  (`_active_parameters_unsupported`, `_active_params_constructed_fails`,
  `_parameter_instances_unsupported`), `_mxs_too_small_fails`, and
  `_inconsistent_fails` stay green.

Results: `-k mfusgsgb` **24**, `-k "mfusgsgb or mfusghfb or mfusgdrt"` **93**,
focused **272**, exe **4**, combined **276** (USG-T 2.7 ARM).

## Status after Stage 4.6C-B

SGB: **Full (non-parametric authoring) + parameter-definition preservation **and**
from-scratch definition authoring** + Expanded valid write (list controls).
**Active SGB parameters remain unsupported** (Fortran `PARTYP='SGB'` vs `'G'`
conflict) — `NotImplementedError`; `INSTANCES` likewise. Verification: FloPy
load → write → reload + the Fortran format/sizing audit; the from-scratch file
uses the same writers as the exe-exercised preservation path.

## Validation

```bash
python -m pytest autotest/test_usg_transport.py -k mfusgsgb -q   # 24 passed
python -m pytest autotest/test_usg_transport.py -q               # 272 passed
python -m pytest autotest/test_usg_transport_exe.py -q           # 4 passed
USGT_EXE=.../usgt_270_arm python -m pytest autotest/test_usg_transport.py autotest/test_usg_transport_exe.py -q   # 276 passed
git diff --check
```
