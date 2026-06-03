# Stage 4.7A — A5 list-parameter edge cases (HFB transient params + INSTANCES)

Date: 2026-06-03 · Base: `develop` @ `72f716b6`

Goal: close/formalize the last gap-audit item **A5** — `TRANSIENT_HFB`+`NPHFB>0`
and parameter `INSTANCES` in the list-parameter packages (HFB, DRT, QRT, SGB).
Implement only branches USG-T 2.7 supports consistently; where the Fortran
rejects them or they are execution-unsafe, keep/strengthen the explicit failure
with tests.

**Decision after the Fortran audit: A5 is closed as "deferred, explicit-fail"** —
none of the branches is implementable execution-safe. The change here is to
**strengthen the `NotImplementedError` messages with the exact Fortran
lines/aborts** and add the one missing negative test (HFB INSTANCES on load).
No behavior change.

## Fortran audit

### HFB `INSTANCES` — USG-T itself aborts

`gwf2hfb7u1.f:131-138`:

```
DO 20 K = 1,NPHFB
  CALL UPARLSTRP(LSTSUM,MXHFB,INHFB,IOUT,IP,'HFB ','HFB ',1,NUMINST)
  IF(NUMINST.GT.0) THEN
    WRITE(IOUT,*) ' INSTANCES ARE NOT SUPPORTED FOR HFB'
    CALL USTOP(' ')
  END IF
```

HFB explicitly **`USTOP`s** when a parameter has `INSTANCES`. Not implementable —
the executable rejects it. (`MfUsgHfb.load` already raised; message now quotes
this.)

### HFB `TRANSIENT_HFB` + `NPHFB>0` — duplicate-parameter abort

`TRANSIENT_HFB` (`gwf2hfb7u1.f:54`, `ITRHFB=1`) re-reads the barrier dataset every
stress period. With `NPHFB>0` the second period would re-process the parameter
definitions through `UPARLSTRP` with `ITERP=1`, which aborts on a repeated name:
`parutl7.f:604-609` — `IF(PARTYP(NP).NE.' ' .AND. ITERP.EQ.1) THEN ... 'Duplicate
parameter name' ... CALL USTOP`. So the combination is execution-unsafe; FloPy
keeps it `NotImplementedError` (message now cites the abort).

### `INSTANCES` machinery for list params (DRT/QRT) — structural-only at best

`UPARLSTRP`/`UINSRP`/`UPARLSTSUB` (`parutl7.f`) do support `INSTANCES` generically:
`UPARLSTRP` reserves `NLST*NUMINST` rows (`:641-646`), and on activation
`UPARLSTSUB` copies the selected instance's `NLST` rows into the active list
(`:812-872`, `NLST=NLST/NUMINST`, `III=I-1+IPLOC(1,IP)+(NI-1)*NLST`). **But** for
DRT/QRT the recipient-node lists are read with `U1DINT` **outside** the parameter
`RLIST`:

- **DRT**: spreading recipients (`NodDRT`, `NR<0`) are not copied on activation —
  an activated DRT parameter is already *structural-only* without instances.
- **QRT**: doubly unsafe — the parameter value scales `QRTF(5)=NumRT`, not `Q`
  (`IPVL=5`, a Fortran bug), and `NodQRT` is not copied on activation.

So per-instance DRT/QRT parameters cannot be honored execution-safe, and FloPy
does not model the instance structure in its list-parameter reader. Kept
`NotImplementedError` (messages now state the structural-only reason).

### SGB active params / `INSTANCES` — type-conflict abort (unchanged)

SGB defines parameters as `PARTYP='SGB'` (`UPARLSTRP`, `glo2sgbu1.f:97`) but
activates them as `PTYP='G'` (`UPARLSTSUB`, `glo2sgbu1.f:185`) — a type conflict
that aborts (`parutl7.f`). Active SGB params are therefore impossible, and
`INSTANCES` (which only matter when a parameter is activated) are equally
blocked. Decision unchanged from Stage 4.4C/4.6C-B.

## What changed (code — message-only, no behavior change)

`flopy/mfusg/`:

- **mfusghfb.py**: the `INSTANCES` load error now quotes `gwf2hfb7u1.f:135-137`
  ("INSTANCES ARE NOT SUPPORTED FOR HFB" + `USTOP`); the `TRANSIENT_HFB`+params
  `write_file` error now cites the `ITERP=1` "Duplicate parameter name" abort
  (`parutl7.f:604-609`).
- **mfusgdrt.py / mfusgqrt.py**: the `INSTANCES` load errors now state the
  structural-only reason (recipient `NodDRT`/`NodQRT` lists are read outside the
  parameter RLIST and not copied on activation; QRT also scales `NumRT` not `Q`).
- **mfusgsgb.py**: unchanged — its message already cites the `SGB`/`G` conflict.

## Tests

New: `test_mfusghfb_parameter_instances_unsupported` — an HFB file with a
parameter carrying `INSTANCES 2` raises `NotImplementedError` on load (the
previously-untested HFB branch). The existing negative tests cover the rest:
`test_mfusghfb_transient_with_parameters_fails`,
`test_mfusgsgb_active_parameters_unsupported`,
`test_mfusgsgb_parameter_instances_unsupported`,
`test_mfusgqrt_parameter_instances_unsupported`,
`test_mfusgdrt_parameter_instances_unsupported`. `-k "mfusghfb or mfusgdrt or
mfusgqrt or mfusgsgb"` **141**, focused **304**, exe **4**, combined **308**
(USG-T 2.7 ARM).

(`ruff`: the four package files format clean; a pre-existing E501 at
`mfusghfb.py:254`, unrelated to A5, is left untouched per the surgical-changes
policy.)

## Status — A5 closed (deferred, explicit-fail), package labels unchanged

Every A5 branch is now an explicit, Fortran-justified failure with a test:

| Branch | Fortran verdict | FloPy |
|---|---|---|
| HFB `INSTANCES` | USG-T `USTOP`s ("…NOT SUPPORTED FOR HFB") | `NotImplementedError` (load) |
| HFB `TRANSIENT_HFB`+`NPHFB>0` | `ITERP=1` duplicate-name abort | `NotImplementedError` (write) |
| DRT `INSTANCES` | structural-only (NodDRT not copied) | `NotImplementedError` (load) |
| QRT `INSTANCES` | structural-only + scales NumRT not Q | `NotImplementedError` (load) |
| SGB active params / `INSTANCES` | `PARTYP='SGB'` vs `'G'` abort | `NotImplementedError` (load) |

No package status label changes: HFB/QRT/DRT stay as in 4.6 (parameter-preserving
+ from-scratch authoring; QRT structural-only), SGB definition-preserving +
authoring. This is the **honest closeout of A5**: the remaining list-parameter
edge cases are not implementable execution-safe, so they fail explicitly with the
Fortran reason cited.

## Validation

```bash
python -m pytest autotest/test_usg_transport.py -k "mfusghfb or mfusgdrt or mfusgqrt or mfusgsgb" -q   # 141
python -m pytest autotest/test_usg_transport.py -q          # 304
python -m pytest autotest/test_usg_transport_exe.py -q      # 4
USGT_EXE=.../usgt_270_arm python -m pytest autotest/test_usg_transport.py autotest/test_usg_transport_exe.py -q   # 308
git diff --check
```
