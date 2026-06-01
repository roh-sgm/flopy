# Stage 4.4A — MODFLOW Parameter Preservation: ETS (executed)

Scope: design a shared abstraction to **preserve** MODFLOW-style named
parameters in USG-T packages, and use **ETS** as the first real implementation.
HFB/SGB/DRT/QRT are out of scope for this commit beyond design/docs (they keep
their explicit `NotImplementedError` on `NP*>0`).

## Fortran audit

### Shared parameter machinery (`parutl7.f`)

Two families of named parameters, both already *parsed* by FloPy's
`flopy.modflow.ModflowParBc`:

- **Array parameters** — `UPARARRAL` / `UPARARRRP` / `UPARARRSUB2`. Used by the
  surface array packages (ETS, EVT, RCH). A definition is
  `PARNAM PARTYP PARVAL NCLU [INSTANCES n]`, followed by (per instance) an
  instance-name line and `NCLU` cluster lines. With `ILFLG=0` (surface packages)
  a cluster line is `MLTARR ZONARR [zone ...]` (no layer field). Activation
  (`UPARARRSUB2`) reads `count` lines of `PNAME [INSTANCE] [IPF]`.
- **List parameters** — `UPARLSTAL` / `UPARLSTRP` / `UPARLSTSUB`. Used by the
  list/boundary packages (SGB, DRT, QRT, HFB-style). Each parameter owns `NLST`
  list rows scaled by the parameter value.

### ETS (`gwf2ets8u1.f`)

- **Item 2a** carries `NETSOP IETSCB NPETS NETSEG [IESFACTOR]`. ETS calls
  `UPARARRAL(-1, ...)`; because `IN<0`, `parutl7.f:14` (`IF(IN.GT.0)`) **skips**
  the `PARAMETER N` decode. So USG-T ETS has **no separate `PARAMETER` line** —
  `NPETS` comes from item 2a (unlike MODFLOW-2005 ETS).
- After item 2a (and the optional transport `ESFACTOR` record), if `NPETS>0`,
  `NPETS` parameter definitions are read with `UPARARRRP(IN,IOUT,N,0,PTYP,1,1,0)`
  (`ILFLG=0`, `ITVP=1` → instances allowed; `PTYP` must be `ETS`).
- **Per stress period**, only **ETSR** (max ET rate) is parameterized. When
  `NPETS>0` and `INETSR>=0`, `INETSR` is the **number of active parameters** that
  period (must be `>=1`), and `UPARARRSUB2` reads that many activation lines in
  place of the ETSR array. `INETSR<0` reuses the previous period.
- `ETSS` (surface), `ETSX` (extinction depth), `IETS` (layer indicator), and
  `PXDP`/`PETM` (segments) are **always plain arrays** — never parameterized. So
  a parameterized ETS file inherently **mixes** a parameterized ETSR with
  non-parameterized arrays (this is the "mix" the plan asks to support/test).

## Design — shared abstraction

`flopy/mfusg/_usgt_parameters.py` adds the missing **write** side for the
array-parameter form, reusing `ModflowParBc.loadarray` as the (already shared)
parser. No second parser is introduced. Structured objects, not ad-hoc strings:
the parsed `ModflowParBc.bc_parms` dict is the structured representation and is
preserved verbatim on the package.

Functions:

- `write_array_parameter_defs(f, pak_parms)` — serialize the definition blocks
  (`UPARARRRP` grammar, `ILFLG=0`), including `INSTANCES`.
- `read_active_array_parameters(f, count, pak_parms)` — read per-period
  activation records into `[(name, instance_or_None), ...]`, consuming an
  instance token only for time-varying parameters; a trailing print flag (`IPF`)
  is not preserved (it is listing output, not data).
- `write_active_array_parameters(f, records)` — write the activation records.

The list-parameter write path (SGB/DRT/QRT/HFB) is **designed but not
implemented** here; see `USGT_STAGE4_04_PARAMETERS.md`. Those packages keep
their explicit `NotImplementedError`.

## ETS implementation

- `MfUsgEts.load(..., expand_parameters=False)`:
  - **default (`False`)**: *preserve*. Keep `npets`, store the parsed
    `ModflowParBc` on `self.parameters`, and store per-period activation records
    on `self.evtr_parm = {kper: [(name, instance), ...]}`. The in-memory ETSR
    array is still filled by expansion (so `.evtr` is meaningful and
    `NONE`/`ALL` clusters expand without MULT/ZONE packages).
  - **`True`**: legacy *expand* — parameters become concrete ETSR arrays and the
    package writes `NPETS=0` ("Expanded valid write").
- `MfUsgEts.write_file()`:
  - If `npets>0` and `self.parameters` is set: write `NPETS` in item 2a (no
    `PARAMETER` line), the definition blocks, and per period `INETSR =
    len(records)` (or `-1` to reuse) with the activation records in place of the
    ETSR array — while `ETSS`/`ETSX`/`IETS`/`PXDP`/`PETM` stay plain arrays.
  - If `npets>0` and `self.parameters` is `None` (from-scratch authoring):
    `NotImplementedError` (authoring parameter definitions from Python is out of
    scope).
  - If `npets==0`: unchanged expanded write.

## Tests (`-k mfusgets`, 10 passed)

New:

- `test_mfusgets_parameterized_load_preserves` — default load keeps `npets=1`,
  `parameters`, and the per-SP activation record.
- `test_mfusgets_parameterized_write_preserves_syntax` — write keeps `NPETS=1`
  in item 2a, the def block + cluster, the activation record, **no** `PARAMETER`
  line, and the parameterized ETSR mixed with two plain `CONSTANT` arrays
  (ETSS+ETSX).
- `test_mfusgets_parameterized_roundtrip` — load → write → reload keeps the
  syntax and the `5.0E-4` value.
- `test_mfusgets_parameter_instances_roundtrip` — a time-varying parameter
  (`INSTANCES 2`) preserves its active instance (`etsrate spring`).

Regression (the legacy fallback, made explicit):

- `test_mfusgets_parameterized_load_expands_to_npets0` now passes
  `expand_parameters=True` and still asserts `NPETS=0` / no `PARAMETER` /
  expanded array.
- `test_mfusgets_parameterized_write_fails_explicitly` — from-scratch `npets=1`
  (no loaded defs) still raises `NotImplementedError`.

The four existing non-parameter ETS tests (construction, NETSEG>1 write,
NETSOP=2 authoring, IESFACTOR) are unchanged.

No from-scratch parametric **executable** smoke was added: building a convergent
parameterized model is not cheap, and the writer is audited line-by-line against
`gwf2ets8u1.f` + `parutl7.f` and round-trips in FloPy. Executable validation of
the preserved file remains a documented manual tier.

## Status decision

ETS stays **⚠️ Partial / not `Full`**, but the parameter dimension is upgraded
from *Expanded valid write* to **parameter-preserving for the ETSR array
parameter** (load → write → reload, instances included), with:

- an opt-in expanded fallback (`expand_parameters=True`), kept as the honest
  fallback; and
- an explicit `NotImplementedError` for from-scratch parameter *authoring*.

Not promoted to `Full` because from-scratch parameter authoring is unsupported.
(Only ETSR can be parameterized in the Fortran, so that is not a gap; ETS zonal
time-series is a separate EVT-side gap tracked elsewhere.)

## Validation

```bash
python -m pytest autotest/test_usg_transport.py -k mfusgets -q
python -m pytest autotest/test_usg_transport.py -q
python -m pytest autotest/test_usg_transport_exe.py -q
USGT_EXE=/Users/roh.sgm/Documents/GitHub_Projects/GW-Software-Compiled-via-Claude/gfortran/usgt_2.7/usgt_270_arm python -m pytest autotest/test_usg_transport.py autotest/test_usg_transport_exe.py -q
git diff --check
git status --short
```

Results: `-k mfusgets` **10 passed**; focused **169 passed**; exe **4 passed**;
combined **173 passed** under the USG-T 2.7 ARM binary.
