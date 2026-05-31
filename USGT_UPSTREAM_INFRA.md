# USG-Transport Upstream Infrastructure Track (Stage 3 Card 10)

Date: 2026-05-31

Goal: define the executable + release infrastructure needed before asking USGS /
`modflowpy/flopy` to accept deeper USG-Transport support. This answers the
November 2024 review criterion: **FloPy package support alone is not enough — a
reproducible executable release with its own CI and stable tests is required.**

This track is mostly *outside* FloPy and should proceed in parallel only after
the package-level story (Cards 1–9) is stable. FloPy must not depend on an
ad-hoc zip of unknown provenance.

## 1. Current upstream reality (verified 2026-05-31)

- **`flopy get-modflow`** (`flopy/utils/get_modflow.py`) downloads released
  binary bundles from GitHub: default `MODFLOW-ORG/executables`, plus the
  `modflow6` / `modflow6-nightly-build` repos. It hardcodes no program names —
  it extracts whatever assets a tagged release ships.
- The **USGS MODFLOW-USG** (`mfusg`) binary is distributed inside the
  `MODFLOW-ORG/executables` bundle (built with `pymake`).
- **USG-Transport** (the GSI Environmental transport-enabled build used here,
  invoked as `mfusg_gsi`, version 2.7.0) is a **distinct program from a
  separate source tree** (`USGT_V_2-7-0_Source_Code`). It is **not** in the
  `MODFLOW-ORG/executables` bundle and is **not** installable via `get-modflow`.
- FloPy/`modflow_devtools` `@requires_exe` resolves an executable by name on
  `PATH` **or by absolute path**. Both the from-scratch
  `test_usg_transport_exe.py` tests and the real-model `Ex1..Ex9` run tests in
  `test_usg_transport.py` share one contract — `@requires_exe(USGT_EXE)` with
  `USGT_EXE = os.environ.get("USGT_EXE", "mfusg_gsi")` — and skip cleanly unless
  the executable resolves. There is currently **no reproducible,
  CI-installable USG-T executable**; it must be supplied locally via `USGT_EXE`
  (a name on `PATH` or an absolute path; default `mfusg_gsi`), and must be
  **USG-Transport 2.7**.

## 2. The gap

| Criterion (Nov-2024) | Today | Needed |
|---|---|---|
| Reproducible release | proprietary zip / local build | tagged source release + prebuilt assets |
| Executable CI | none public | compile + smoke-test on Linux/macOS/Windows |
| Stable executable tests | only this fork's optional suite | autotests living *with* the executable source |
| FloPy dependency | none (tests skip) | optional, pinned to a tagged release |

## 3. Proposed path (in order)

1. **Stand up a USG-Transport source repository** (e.g.
   `<org>/usg-transport`) with:
   - the 2.7.0 source tree under version control,
   - **tagged releases** (`v2.7.0`, …) with source archive assets,
   - documented build instructions (gfortran + an Intel option), and a
     `Makefile`/`meson`/`pymake` build path.
2. **Executable CI in that source repo** (GitHub Actions):
   - compile on Linux, macOS, Windows,
   - run **smoke tests** (a handful of minimal + reference inputs, assert
     normal termination and compare `LST`/`HDS`/`CON`/`CBB` to stored
     references),
   - publish built binaries as release assets per platform.
3. **Executable autotests with the source** (not only in FloPy): the minimal
   from-scratch flow/transport models in `test_usg_transport_exe.py` are a good
   seed and can be mirrored/expanded in the source repo so the executable is
   validated independently of FloPy.
4. **Only after that CI exists**, wire FloPy:
   - propose a `pymake` target and/or a `get-modflow`-compatible release source
     (a dedicated `--repo`/owner entry, or inclusion in a transport bundle),
   - point the optional FloPy USG-T exe tests at the tagged release via
     `USGT_EXE` (already supported) or a `get-modflow`-installed path,
   - keep the dependency **optional**: default CI stays green without the exe.

## 4. What this fork already provides toward the criteria

- A clean, honest, well-tested FloPy package layer (a 101-test focused suite:
  synthetic authoring/round-trip plus `Ex1..Ex9` real-model loads; status audit;
  explicit failure for unsupported modes) — see `USGT_roadmap.md`.
- An **opt-in executable validation tier** (`USGT_EXE`-gated) that runs
  from-scratch and real models end-to-end and checks budgets — the seed for the
  source-repo smoke tests (Card 9).
- Documented bit-for-bit real-model validation (Model A / Model B) in
  `USGT_improvements.md`.

## 5. Acceptance

- A reproducible **tagged USG-T release source** exists.
- **CI proves** the executable builds and passes smoke tests on the major OSes.
- FloPy depends on that tagged release only **optionally**, never on an ad-hoc
  zip with unknown test status.

## Out of scope here

Building the executable CI and the external source repository is outside the
`modflowpy/flopy` fork. This document is the plan; execution happens in the
USG-Transport source project.
