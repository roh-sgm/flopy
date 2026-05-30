# Spec: Transient Idomain (TID) for MODFLOW 6 in FloPy

**Status:** Draft — for implementation in a dedicated Claude Code session.  
**Repo:** `https://github.com/roh-sgm/flopy` (branch `develop`)  
**Related:** `flopy/mfusg/mfusgtib.py` (USG-T TIB reference implementation)

---

## 1. Background: What USG-T TIB Does

The USG-T `TIB` (Transient Ibound) package lets you change the `ibound` status
of individual cells each stress period without rebuilding the model. Three
counters per SP control the transitions:

| Counter | Meaning |
|---------|---------|
| `NIB0`  | Cells switching **active → inactive** (ibound → 0) |
| `NIB1`  | Cells switching **inactive → active** (ibound → 1, head assigned) |
| `NIBM1` | Cells switching **constant-head → inactive** (ibound: -1 → 0) |

When transport (BCT) is active, three matching counters follow (`NICB0`,
`NICB1`, `NICBM1`) to change the concentration ibound in parallel.

**Typical use case:** Simulating progressive cell activation — a tailings
impoundment filling over time, a lake expanding, excavation stages — where the
active domain changes each stress period.

---

## 2. Why MF6 Is Different (and Why We Need Both `Mf6Tid` and the BMI)

USG-T inherits MODFLOW-2005's `BAS` package concept where `ibound` is a
stress-period array the simulator re-reads every SP from the input file.

MF6 replaced this with `idomain` in the discretization packages (`DIS`,
`DISV`, `DISU`):

- **`idomain` in MF6 is defined once** in the DIS/DISV/DISU input file. There
  is no stress-period variant.
- There is **no native MF6 package** equivalent to TIB.
- MF6 exposes `idomain` through the **BMI (Basic Model Interface)** as a
  writable variable, enabling runtime modification from Python.

### Why two components are needed

The BMI is a low-level API — it can set `idomain`, but it has no concept of
*when* to do so or *which cells* should change. That knowledge must live
somewhere: that is what `Mf6Tid` stores.

```
Mf6Tid (schedule)  →  Mf6TidRunner  →  BMI API  →  MF6 binary
  "what / when"          "glue"         "how"        "runs"
```

- **`Mf6Tid`** = the schedule object (Python data, no binary required). Analogous
  to the TIB input file in USG-T. Can be built, inspected, serialised, and
  loaded without running MF6 at all.
- **`Mf6TidRunner`** = the runtime glue that reads the schedule and calls BMI
  before each SP. Requires MF6 shared library + `modflowapi`.

The naming uses **TID** (Transient Idomain) rather than TIB to be consistent
with MF6 terminology (`idomain`, not `ibound`) and to make the distinction from
the USG-T package explicit.

---

## 3. Recommended Architecture

### `Mf6Tid` — schedule container

Stores per-SP activation/deactivation lists. Does **not** write a separate MF6
input file (there is no such format). Instead serialises to a sidecar JSON so
the schedule survives `write_simulation()` / `load()` round-trips without
needing a live MF6 binary.

```python
from flopy.mf6.tid import Mf6Tid

tid = Mf6Tid(gwf)

# Add transitions (0-based stress period index)
# Node identifiers: flat int for DISU, (lay, row, col) tuple for DIS,
# (lay, cell2d) tuple for DISV.
tid.add(kper=5,  activate=[1002, 1003, 1004])   # idomain → 1 for these nodes
tid.add(kper=10, deactivate=[1002])              # idomain → 0 for this node

# Inspect the cumulative idomain at any SP
tid.idomain_at(kper=10)   # → np.ndarray, shape (nnodes,)

# Persist
tid.write()               # writes <model_ws>/<model_name>.tid.json
tid2 = Mf6Tid.load(gwf)   # reload from JSON
```

### `Mf6TidRunner` — BMI runner

```python
from flopy.mf6.tid import Mf6TidRunner

runner = Mf6TidRunner(sim, tid)
runner.run()   # steps through all SPs, applying idomain changes before each SP
```

Internally:
1. `mf6.initialize(namefile)` — load MF6 shared library and prepare.
2. For each SP:
   a. Resolve cumulative `idomain` from schedule (§5 below).
   b. `mf6.set_value("GWF-1/DIS/IDOMAIN", idomain_state)`.
   c. Advance all time steps in the SP with `mf6.update()`.
3. `mf6.finalize()`.

---

## 4. File Locations (proposed)

```
flopy/
  mf6/
    tid.py     ← Mf6Tid + Mf6TidRunner
```

Export from `flopy/mf6/__init__.py`:
```python
from .tid import Mf6Tid, Mf6TidRunner
```

Sidecar file: `<model_ws>/<model_name>.tid.json`  
Format:
```json
{
  "5":  {"activate": [1002, 1003, 1004], "deactivate": []},
  "10": {"activate": [],                 "deactivate": [1002]}
}
```
Keys are zero-based SP indices as strings. Nodes are stored as flat integers
(after resolving DIS/DISV tuples at `add()` time).

---

## 5. `idomain` State Machine

The runner tracks **cumulative** state — not just per-SP deltas — so that cells
activated in SP 5 remain active in SP 6 unless explicitly deactivated:

```python
idomain_state = initial_idomain.copy()   # read from DIS/DISV/DISU at init

for kper in range(nper):
    for node in tid.schedule.get(kper, {}).get("activate", []):
        idomain_state[node] = 1
    for node in tid.schedule.get(kper, {}).get("deactivate", []):
        idomain_state[node] = 0
    mf6.set_value(idomain_var, idomain_state)
    while current_sp == kper:
        mf6.update()
```

This matches USG-T TIB semantics exactly.

---

## 6. Node Addressing

| Package | `add()` accepts | Resolved internally to |
|---------|-----------------|------------------------|
| `DIS`   | `(lay, row, col)` tuples (0-based) | flat node index |
| `DISV`  | `(lay, cell2d)` tuples (0-based)   | flat node index |
| `DISU`  | flat int (0-based)                  | flat node index |

Resolution via `gwf.modelgrid.get_node()`. All nodes stored as flat ints in
the JSON sidecar.

---

## 7. Transport / Concentration Ibound

USG-T TIB has matching concentration counters (`NICB0`, `NICB1`, `NICBM1`).
In MF6, transport cells follow the flow `idomain` automatically — there is no
separate concentration ibound exposed by BMI.

**v1 decision:** Concentration ibound support is out of scope. Flow `idomain`
via BMI is sufficient.

---

## 8. Verification Criteria

1. `tid.idomain_at(kper=N)` returns the expected array after a known sequence
   of activate/deactivate calls.
2. A small DIS model (e.g. 3×3×1, 10 SPs) run with `Mf6TidRunner` produces a
   listing file where activated cells show budget contributions starting at the
   correct SP and inactive cells show zero.
3. Re-activating a previously deactivated cell restores normal budget behaviour.
4. `Mf6Tid.load(gwf)` round-trips a written `.tid.json` without data loss.
5. Node conversion (DIS tuple → flat index) agrees with
   `modelgrid.get_node()`.

---

## 9. Dependencies

- `modflowapi` (pip: `modflowapi`) — required for `Mf6TidRunner` only.
  Import-guarded:
  ```python
  try:
      import modflowapi
  except ImportError:
      raise ImportError("modflowapi is required for Mf6TidRunner. pip install modflowapi")
  ```
- `flopy.mf6` — already present.
- `Mf6Tid` itself: only `json` + `numpy`.

---

## 10. Out of Scope (v1)

- Concentration idomain (no BMI equivalent in MF6).
- Head assignment for activating cells — MF6 uses the last known head for
  newly-activated cells. If a specific start head is needed, set it via the
  BMI `X` array before activating or use the IC package.
- Constant-head → inactive transition (`NIBM1` equivalent) — CHD package in
  MF6 is separate from `idomain`; deactivating a CHD cell requires removing it
  from the CHD stress period data in the same SP (out of scope here).
- Command-line runner (non-BMI). If BMI is not available, `Mf6Tid` can still
  be used to generate a per-SP `idomain` array that the user applies manually.
- GUI / Vistas integration.

---

## 11. Reference

- USG-T Fortran source: subroutine `GWF2TIB1RP` in `glo2basu1.f`
- USG-T FloPy implementation: `flopy/mfusg/mfusgtib.py`
- MF6 writable BMI variables: `mf6.get_var_names()` — look for `*/DIS/IDOMAIN`
- `modflowapi` package: https://github.com/MODFLOW-USGS/modflowapi
