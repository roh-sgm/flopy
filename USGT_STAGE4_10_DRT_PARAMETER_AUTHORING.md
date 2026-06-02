# Stage 4.6B — DRT from-scratch MODFLOW list-parameter authoring (`NPDRT>0`)

Date: 2026-06-02 · Base: `develop` @ `4c8e5870`

Goal: promote DRT from *parameter-preserving only* to **from-scratch parameter
authoring** for `NPDRT>0`, while keeping the existing load → write → reload
preservation. A user can now build a parameterized `MfUsgDrt` entirely in Python
(no prior load), write a valid `NPDRT>0` file, reload it with the same semantics,
and activate parameters per stress period. This is the pilot for the shared
parameter-authoring path that Stage 4.6C will reuse for HFB / SGB(defs) / QRT.

## Audit

The on-file format for a parameterized DRT is identical whether the definitions
came from `load` or from the user — "from-scratch" is purely a FloPy-side
distinction about where `self.parameters` / `self.active_params` originated. The
write side already used the shared, Fortran-audited writers
(`write_list_parameter_header` → `UPARLSTRP` header `PARNAM PARTYP PARVAL NLST`;
the per-SP `ITMP NP` header + `write_active_list_parameters`; `_write_drain_line`
for the `NLST` rows with their RETURNFLOW recipients / spreading blocks), and the
preservation path round-trips (and is exercised by the Stage 4.4D tests). So
Stage 4.6B does **not** change the file format; it makes building
`self.parameters` from Python ergonomic and validated.

Fortran sizing (`gwf2drt8u.f`): item 1 is `MXADRT IDRTCB NPDRT MXL`; the
parameter rows are stored from `IDRTPB = MXADRT+1` into `MXDRT = MXADRT + MXL`
(line 114–115), so **`MXL` must be ≥ the total of all definition `NLST`**.
`MXADRT` must cover the active total per period (`NDRTNP + Σ active NLST`;
`SGWF2DRT8LS` aborts if `NDRTCL > MXADRT`). DRT is type-consistent — definition
and activation both use `PARTYP='DRT'` (`gwf2drt8u.f:135`/`:1108`) — so active
DRT parameters are valid (unlike SGB).

## What changed (code)

Only `flopy/mfusg/mfusgdrt.py`:

- **`_normalize_param(name, pdef)`** (new): validates + canonicalizes one
  definition. Accepts an ergonomic from-scratch dict *or* an already-canonical
  loaded one and returns `{partyp, parval, nlst, data (recarray),
  recipient_nodes}`:
  - `partyp` defaults to `"DRT"` and must be `"DRT"`;
  - `parval` required;
  - `data` required & non-empty; a recarray is kept, otherwise built from any
    array-like with the active `dtype` (e.g. a list of tuples), nodes 0-based;
  - `nlst` computed from `len(data)` when omitted, validated when given;
  - `recipient_nodes` defaults to empty lists when omitted; must be one list per
    row; non-empty recipients require `RETURNFLOW`.
- **`_validate_parameter_write()`** (rewritten): normalizes every definition,
  resolves `MXL` (`= Σ nlst` when omitted/0, else validated `≥` that total), and
  validates `active_params` (stress period in `0..nper-1`; names defined
  case-insensitively; no duplicate activation in a period). Returns
  `(params, mxl)`. The old `NotImplementedError` "cannot author from scratch" is
  gone; `active_params` with **no** definitions now raises a clear `ValueError`.
- `write_file` consumes the returned canonical `params`/`mxl`; `_max_active_drains`
  / `_max_spread_nodes` take the canonical `params`.

Internal nodes/recipients stay 0-based; the file is 1-based (`+1` on write).

## API

```python
drt = MfUsgDrt(
    mfusg_model,                      # unstructured MfUsg
    options=["RETURNFLOW", "CHANGEC", "AUX C01"],
    parameters={
        "dp": {
            "parval": "1.5",          # scales COND (Fortran IPVL=5)
            "data": [                 # recarray or array-like (dtype-built); 0-based nodes
                (0, 5.0, 10.0, 0.7, 2, 0.15),
                (3, 4.0, 20.0, 0.5, 1, 0.25),
            ],
            "recipient_nodes": [[8, 9], [7]],   # RETURNFLOW only; 0-based
            # "partyp": "DRT"  (default), "nlst": 2  (auto from data)
        },
    },
    # "mxl": auto = sum(nlst) when omitted
    active_params={0: ["dp"]},        # per zero-based stress period
)
drt.write_file()                      # valid NPDRT>0 file
```

## Validations (all before the file is opened — no partial file)

| Condition | Result |
|---|---|
| `partyp` not `DRT` | `ValueError` |
| `parval` missing/blank, or a multi-token string (e.g. `"1 2"`) | `ValueError` |
| `data` missing/empty | `ValueError` |
| `nlst` given ≠ `len(data)` | `ValueError` |
| `recipient_nodes` length ≠ `nlst` | `ValueError` |
| recipients present without `RETURNFLOW` | `ValueError` |
| `mxl` given `<` total definition rows | `ValueError` |
| definition/active **name** not a single whitespace-free token | `ValueError` |
| definition/active **name** > 10 chars (Fortran `CHARACTER*10` PARNAM) | `ValueError` |
| duplicate **definition** name (case-insensitive) | `ValueError` |
| `active_params` name not defined (case-insensitive) | `ValueError` |
| duplicate activation in one period (case-insensitive) | `ValueError` |
| `active_params` stress period outside `0..nper-1` | `ValueError` |
| `active_params` but no `parameters` defined | `ValueError` |
| negative parametric `data["node"]` or `recipient_nodes` entry | `ValueError` |
| parameter `INSTANCES` | unsupported (load raises; no authoring field) |

The name and `parval` rules come from `parutl7.f` `UPARLSTRP`: `PARNAM` (`PN`)
is read as one `URWORD` word into a `CHARACTER*10` buffer and upper-cased, and
`PARVAL` (`PV`) is read as one numeric token — so longer names truncate/collide
and only a single whitespace-free token is meaningful (Stage 4.6B review
follow-up; this is contract hardening, not a status change).

## Tests (`autotest/test_usg_transport.py`)

- `test_mfusgdrt_parameter_authoring_from_scratch` — ergonomic input, write →
  reload, `NPDRT`/auto-`MXL`/defs/activation, 0-based ↔ 1-based.
- `test_mfusgdrt_parameter_authoring_returnflow_changec_aux` — RETURNFLOW +
  CHANGEC + AUX + per-row `recipient_nodes`.
- `test_mfusgdrt_parameter_authoring_mixed_nonparam_active` — non-param row +
  active parameter in one SP; `MXADRT` counts both.
- `test_mfusgdrt_parameter_authoring_mxl_auto` — `MXL` auto = Σ nlst.
- `test_mfusgdrt_parameter_authoring_partyp_invalid_fails`,
  `_kper_out_of_range_fails`, `_recipients_without_returnflow_fails` — new
  validations, no partial file.
- `test_mfusgdrt_parameter_active_without_defs_fails` — repurposed from the old
  "from-scratch fails" test (now `ValueError`, since authoring *with* definitions
  is supported).

Review follow-up (contract hardening) — 5 more:
`test_mfusgdrt_parameter_authoring_duplicate_definition_names_fails` (`dp`/`DP`),
`_bad_name_fails` (blank / whitespace / multi-token / >10 chars),
`_bad_parval_fails` (`"1 2"` / blank), `_negative_node_fails`,
`_negative_recipient_fails` — each `ValueError`, no partial file.

The Stage 4.4D preservation tests (`parameterized_roundtrip`, `mxadrt`, `spread`,
`mixed`, `reuse`, `sfac`, `instances_unsupported`, `mxl_too_small`,
`inconsistent`, `active_undefined`, `case_insensitive`, `duplicate`) stay green.
After the follow-up: `-k mfusgdrt` **39**, focused **253**, combined **257** (ARM).

## Status after Stage 4.6B

DRT: **Full (authoring) — incl. from-scratch `NPDRT>0` parameters** /
parameter-preserving / Expanded valid write (list controls). Caveats unchanged:
`INSTANCES` unsupported; an activated `SPREAD` (`NR<0`) parameter round-trips
structurally but is not execution-guaranteed (USG-T copies `DRTF` but not
`NodDRT` on activation).

## Validation

```bash
python -m pytest autotest/test_usg_transport.py -k mfusgdrt -q   # 39 passed
python -m pytest autotest/test_usg_transport.py -q               # 253 passed
python -m pytest autotest/test_usg_transport_exe.py -q           # 4 passed
USGT_EXE=.../usgt_270_arm python -m pytest autotest/test_usg_transport.py autotest/test_usg_transport_exe.py -q   # 257 passed
git diff --check
```

Verification level: Fortran format audit + FloPy load → write → reload equality
(the from-scratch file uses the same writers as the exe-exercised preservation
path). An end-to-end `USGT_EXE` run of a from-scratch parametric DRT was not
added (it needs a full from-scratch DISU model); the focused/combined suites are
green.
