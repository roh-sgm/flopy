# Stage 4.6C-A — HFB from-scratch MODFLOW list-parameter authoring (`NPHFB>0`)

Date: 2026-06-02 · Base: `develop` @ `ad80ea7f`

Goal: promote `MfUsgHfb` from *parameter-preserving only* to **from-scratch
parameter authoring** for **non-transient** HFB (`NPHFB>0`), reusing the DRT
4.6B pattern, while keeping the existing load → write → reload preservation. A
user can now build a parameterized `MfUsgHfb` entirely in Python (no prior
load), write a valid `NPHFB>0` file, reload it with the same semantics, and
activate parameters. SGB/QRT are **not** touched (Stage 4.6C-B).

## Audit

The on-file format for a parameterized HFB is identical whether the definitions
came from `load` or from the user — "from-scratch" is purely a FloPy-side
distinction about where `self.parameters` / `self.acthfb_names` came from. The
write side already used the shared, Fortran-audited writers
(`write_list_parameter_header` → `UPARLSTRP`; `write_active_list_parameters`;
`_write_hfb_rows`), and the preservation path round-trips (Stage 4.4B). So this
stage does **not** change the file format; it makes building `self.parameters`
from Python ergonomic and validated.

Fortran sizing (`gwf2hfb7u1.f`): item 1 is `NPHFB MXFBP NHFBNP`; the parameter
rows are stored from `IHFBPB` into `MXHFB = (NHFBNP + MXFBP) + MXFBP`
(lines 34–36, 65–67), and `UPARLSTRP` reads them with `PARTYP='HFB '`. So
**`MXFBP` must be ≥ the total of all definition `NLST`** — `mxfb = Σ nlst` is the
minimal valid value. Names/`parval` follow the shared `parutl7.f` `UPARLSTRP`
rules (`PARNAM` is one `URWORD` word in a `CHARACTER*10` buffer, upper-cased;
`PARVAL` is one numeric token).

## What changed (code)

- **`flopy/mfusg/_usgt_parameters.py`**: factored the name/`parval` validation
  shared with DRT into `check_parameter_name(name, what, prefix)` and
  `check_parval(parval, name, prefix)` (the latter accepts any `numbers.Real`
  — `int`/`float`/`np.integer`/`np.floating`, excluding `bool` — or a single
  string token).
- **`flopy/mfusg/mfusgdrt.py`**: its `_check_param_name` / `_check_parval` now
  delegate to the shared helpers (behavior unchanged; DRT's 4.6B tests are the
  regression).
- **`flopy/mfusg/mfusghfb.py`**:
  - new **`_normalize_param(name, pdef, structured)`** — validates + canonicalizes
    one definition (`partyp` defaults to / must be `HFB`; `data` recarray *or*
    array-like built with the active barrier dtype; `nlst` computed/validated;
    barrier indices non-negative; preserves the original `partyp` text so a
    loaded file's casing round-trips);
  - **`_validate_parameter_write`** now normalizes every definition, auto-computes
    the counts (`NPHFB = len(params)`, `MXFBP = Σ nlst`, `NACTHFB =
    len(acthfb_names)` when omitted/0; validated if given), checks
    case-insensitive duplicate definition names and active-name validity, and
    returns `(params, nphfb, mxfb, nacthfb)`;
  - `write_file` triggers the parameter path on **any** parameter intent
    (`parameters`, `acthfb_names`, or `nphfb>0`) and writes the auto-computed
    counts; the old `NotImplementedError` "cannot author from scratch" is gone —
    `NPHFB>0`/active-names with no definitions now raises a clear `ValueError`.

Internal barrier indices stay 0-based (`node1`/`node2` or `k`/`irow*`/`icol*`);
the file is 1-based. Works for structured and unstructured grids.

## API

```python
hfb = MfUsgHfb(
    mfusg_model,
    parameters={
        "hp": {
            "parval": "2.0",                       # scales HYDCHR
            "data": [(0, 1, 0.5), (2, 3, 0.25)],   # 0-based; (node1, node2, hydchr)
            # "partyp": "HFB" (default), "nlst": 2 (auto from data)
        },
    },
    acthfb_names=["hp"],                            # nphfb/mxfb/nacthfb auto
)
hfb.write_file()                                    # valid NPHFB>0 file
```

Structured barrier rows are `(k, irow1, icol1, irow2, icol2, hydchr)`.

## Validations (all before the file is opened — no partial file)

| Condition | Result |
|---|---|
| `NPHFB>0` / `acthfb_names` but no `parameters` | `ValueError` |
| name not a single whitespace-free token, or > 10 chars | `ValueError` |
| duplicate definition name (case-insensitive) | `ValueError` |
| `partyp` not `HFB` | `ValueError` |
| `parval` missing/blank or a multi-token string (`"1 2"`) | `ValueError` |
| `data` missing/empty | `ValueError` |
| `nlst` given ≠ `len(data)` | `ValueError` |
| explicit `nphfb` ≠ `len(parameters)` | `ValueError` |
| explicit `mxfb` `<` total parameter rows | `ValueError` |
| negative barrier index | `ValueError` |
| active name undefined / duplicate (case-insensitive) | `ValueError` |
| `TRANSIENT_HFB` + parameters | `NotImplementedError` |
| parameter `INSTANCES` | `NotImplementedError` (load) |

## Tests (`autotest/test_usg_transport.py`)

- `test_mfusghfb_parameter_authoring_from_scratch` — unstructured, ergonomic
  input, write → reload, NPHFB / auto-MXFBP / defs / active names, 0↔1-based.
- `test_mfusghfb_parameter_authoring_structured` — structured `k/irow/icol`.
- `test_mfusghfb_parameter_authoring_mixed_nonparam` — `NHFBNP` + parameter.
- `test_mfusghfb_parameter_authoring_auto_counts` — `nphfb`/`mxfb`/`nacthfb` auto.
- Negatives: `_duplicate_definition_names_fails`, `_bad_name_fails`,
  `_bad_parval_fails`, `_empty_data_fails`, `_mxfb_too_small_fails`,
  `_negative_index_fails` — each `ValueError`, no partial file.
- `test_mfusghfb_parameter_active_without_defs_fails` — repurposed from the old
  "from-scratch fails" test (now `ValueError`). The existing preservation /
  count-/nlst-/undefined-/duplicate-active / TRANSIENT+params tests stay green.

Results: `-k mfusghfb` **30**, `-k "mfusghfb or mfusgdrt"` **69**, focused
**263**, exe **4**, combined **267** (USG-T 2.7 ARM).

## Status after Stage 4.6C-A

HFB: **Parameter-preserving + from-scratch parameter authoring (non-transient)**.
Still `⚠️ Partial` / not `Full`: `TRANSIENT_HFB` + `NPHFB>0` and parameter
`INSTANCES` remain `NotImplementedError`. Verification: FloPy load → write →
reload (structured + unstructured) + the line-by-line Fortran format/sizing
audit; the from-scratch file uses the same writers as the preservation path.

## Validation

```bash
python -m pytest autotest/test_usg_transport.py -k mfusghfb -q   # 30 passed
python -m pytest autotest/test_usg_transport.py -q               # 263 passed
python -m pytest autotest/test_usg_transport_exe.py -q           # 4 passed
USGT_EXE=.../usgt_270_arm python -m pytest autotest/test_usg_transport.py autotest/test_usg_transport_exe.py -q   # 267 passed
git diff --check
```
