# Stage 4.6C-C — QRT from-scratch list-parameter authoring (`NPQRT>0`, structural)

Date: 2026-06-02 · Base: `develop` @ `1d0d5401`

Goal: let a user build a parameterized `MfUsgQrt` from Python (`NPQRT>0`), write
a valid file, reload it with the same semantics, and activate parameters per
stress period — reusing the DRT 4.6B / HFB 4.6C-A / SGB 4.6C-B pattern. This is
**structural authoring**: like the Stage 4.4E preservation path, an *activated*
QRT parameter round-trips but is **not execution-guaranteed** (see the caveat
below). DRT/HFB/SGB/ETS are not touched. This closes the list-parameter
from-scratch family (DRT/HFB/SGB/QRT).

## Audit

The on-file format for a parameterized QRT is identical whether definitions came
from `load` or the user; the writers (`write_list_parameter_header` →
`UPARLSTRP`, the sink-block + recipient `U1DINT` writer, per-SP `ITMP NP` +
`write_active_list_parameters`) were Fortran-audited and round-tripped in Stage
4.4E. `preserve = bool(self.parameters) or any(active_params)` already triggered
the param path, so a fully-canonical dict already wrote; the gaps were ergonomic
(normalize `data`→recarray, auto `nlst`/`MXL`/`recipient_nodes`, `partyp`
default, name/`parval`/node validation, `kper` range) and the misleading
from-scratch `NotImplementedError`.

Fortran (`gwf2QRT8u.f`): item 1 is `MXAQRT MXRTCELLS IQRTCB NPQRT MXL`;
definitions are read with `UPARLSTRP` (`PARTYP='QRT'`, :183) and dimensioned into
`MXQRT = MXAQRT + MXL`, so **`MXL` ≥ the total of all definition `NLST`**
(`mxl = Σ nlst` minimal). `MXAQRT` must cover the active total per period
(`NQRTNP + Σ active NLST`; `SGWF2QRT8LS` aborts if `NQRTCL > MXAQRT`).

**Execution caveat (unchanged, Stage 4.4E):** QRT is type-consistent
(`PARTYP='QRT'` both sides, :183/:1118), so an active QRT parameter is
type-*valid*, but two Fortran issues mean activation is **not
execution-guaranteed**: `SGWF2QRT8LS` uses `IPVL=5`, i.e. the parameter value
scales `QRTF(5)=NumRT` (recipient count) **not** `Q` (a bug — `SFAC` scales `Q`
at `ISCLOC=4`); and `NodQRT` is not copied on activation. FloPy therefore
authors/preserves the *structure* faithfully but does **not** apply the
parameter value to `Q`.

## What changed (code)

Only `flopy/mfusg/mfusgqrt.py`:

- new **`_normalize_param(name, pdef)`** — validates + canonicalizes one
  definition (`partyp` defaults to / must be `QRT`, casing preserved; `data`
  recarray or array-like built with the active QRT dtype; `nlst` and
  `recipient_nodes` computed/validated; nodes and recipient nodes non-negative;
  recipients require `RETURNFLOW`) via the shared `check_parameter_name` /
  `check_parval`;
- **`_validate_parameter_write`** now normalizes every definition, checks
  case-insensitive duplicate names, auto-computes `MXL = Σ nlst` (validated `≥`
  total if given), validates `active_params` (stress period `0..nper-1`, names
  defined case-insensitively, no duplicate activation in a period), and returns
  `(params, mxl)`; the old from-scratch `NotImplementedError` is gone —
  `active_params` with no definitions now raises `ValueError`;
- `write_file`, `_max_active_sinks(params)`, `_max_rt_cells(params)` consume the
  canonical `(params, mxl)`; `NPQRT = len(params)`.

Internal nodes/recipients stay 0-based; the file is 1-based.

## API

```python
qrt = MfUsgQrt(
    mfusg_model,
    options=["RETURNFLOW"],
    parameters={
        "qp": {
            "parval": "1.0",                        # scales NumRT in the Fortran (caveat)
            "data": [(0, -100.0, 0.75), (4, -30.0, 0.0)],  # 0-based (node, q, rfprop)
            "recipient_nodes": [[8, 9], []],        # RETURNFLOW only; 0-based
            # "partyp": "QRT" (default), "nlst": 2 (auto)
        },
    },
    # "mxl": auto = sum(nlst)
    active_params={0: ["qp"]},                      # per zero-based stress period
)
qrt.write_file()                                    # valid NPQRT>0 file
```

## Validations (all before the file is opened — no partial file)

| Condition | Result |
|---|---|
| `active_params` but no `parameters` defined | `ValueError` |
| name not a single whitespace-free token, or > 10 chars | `ValueError` |
| duplicate definition name (case-insensitive) | `ValueError` |
| `partyp` not `QRT` | `ValueError` |
| `parval` missing/blank or a multi-token string (`"1 2"`) | `ValueError` |
| `data` missing/empty | `ValueError` |
| `nlst` given ≠ `len(data)` | `ValueError` |
| explicit `mxl` `<` total parameter rows | `ValueError` |
| `recipient_nodes` length ≠ `nlst`, or recipients without `RETURNFLOW` | `ValueError` |
| negative node or recipient | `ValueError` |
| active name undefined / duplicate (case-insensitive) | `ValueError` |
| `active_params` stress period outside `0..nper-1` | `ValueError` |
| `INSTANCES`; `TRANSIENTQ` + `NPQRT>0`; external-unit recipients | `NotImplementedError` |

## Tests (`autotest/test_usg_transport.py`)

`test_mfusgqrt_parameter_authoring_from_scratch`,
`_returnflow_recipients`, `_mixed_nonparam_active` (`MXAQRT` counts both),
`_auto_mxl`, `_active_case_insensitive`; negatives `_partyp_invalid_fails`,
`_duplicate_definition_names_fails`, `_bad_name_or_parval_fails`,
`_kper_out_of_range_fails`, `_negative_node_or_recipient_fails`;
`test_mfusgqrt_parameter_active_without_defs_fails` (repurposed from the old
`_from_scratch_fails`, now `ValueError`). The Stage 4.4E preservation /
4.5A TRANSIENTQ / 4.5B recipient-control tests stay green.

Results: `-k mfusgqrt` **47**, `-k "mfusgqrt or mfusgdrt or mfusgsgb"` **110**,
focused **282**, exe **4**, combined **286** (USG-T 2.7 ARM).

## Status after Stage 4.6C-C

QRT: **Full (non-parametric authoring) + structural parameter
authoring/preservation** + `TRANSIENTQ` (4.5A) + recipient `U1DINT`
EXTERNAL/OPEN-CLOSE (4.5B). Active QRT parameters are **structural, not
execution-guaranteed** (Fortran scales `NumRT` not `Q`; `NodQRT` not copied).
`INSTANCES`, `TRANSIENTQ`+`NPQRT>0`, and external-unit recipients stay
`NotImplementedError`.

## Validation

```bash
python -m pytest autotest/test_usg_transport.py -k mfusgqrt -q   # 47 passed
python -m pytest autotest/test_usg_transport.py -q               # 282 passed
python -m pytest autotest/test_usg_transport_exe.py -q           # 4 passed
USGT_EXE=.../usgt_270_arm python -m pytest autotest/test_usg_transport.py autotest/test_usg_transport_exe.py -q   # 286 passed
git diff --check
```
