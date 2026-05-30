"""Transient Idomain (TID) for MODFLOW 6.

Two classes:
  Mf6Tid       -- schedule container (pure Python, no MF6 binary required)
  Mf6TidRunner -- BMI runner that applies the schedule via modflowapi
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from .mfmodel import MFModel


def _resolve_nodes(gwf: "MFModel", nodes) -> list[int]:
    """Convert node identifiers to flat 0-based indices.

    Accepts:
      - int / list[int]         — DISU flat index, passed through
      - tuple / list[tuple]     — DIS (lay,row,col) or DISV (lay,cell2d)
    """
    if not nodes:
        return []
    mg = gwf.modelgrid
    result = []
    for n in nodes:
        if isinstance(n, (int, np.integer)):
            result.append(int(n))
        elif isinstance(n, tuple):
            result.append(int(mg.get_node(n)[0]))
        else:
            raise TypeError(
                f"Node identifier must be int or tuple, got {type(n)}"
            )
    return result


class Mf6Tid:
    """Transient Idomain schedule for MODFLOW 6.

    Parameters
    ----------
    gwf : MFModel
        The GWF model.  Used to resolve node indices and locate the workspace.

    Examples
    --------
    >>> tid = Mf6Tid(gwf)
    >>> tid.add(kper=2, activate=[(0, 1, 1)])
    >>> tid.add(kper=5, deactivate=[(0, 1, 1)])
    >>> idomain = tid.idomain_at(kper=3)
    >>> tid.write()
    >>> tid2 = Mf6Tid.load(gwf)
    """

    def __init__(self, gwf: "MFModel") -> None:
        self._gwf = gwf
        # {kper: {"activate": [node,...], "deactivate": [node,...]}}
        self.schedule: dict[int, dict[str, list[int]]] = {}

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def add(
        self,
        kper: int,
        activate=None,
        deactivate=None,
    ) -> None:
        """Register cell transitions for a stress period.

        Parameters
        ----------
        kper : int
            Zero-based stress period index.
        activate : list, optional
            Nodes to switch to active (idomain → 1).
            Each item is an int (DISU) or (lay,row,col)/(lay,cell2d) tuple.
        deactivate : list, optional
            Nodes to switch to inactive (idomain → 0).
        """
        act = _resolve_nodes(self._gwf, activate or [])
        deact = _resolve_nodes(self._gwf, deactivate or [])

        if kper not in self.schedule:
            self.schedule[kper] = {"activate": [], "deactivate": []}
        self.schedule[kper]["activate"].extend(act)
        self.schedule[kper]["deactivate"].extend(deact)

    def idomain_at(self, kper: int) -> np.ndarray:
        """Cumulative idomain state at the *start* of stress period *kper*.

        Reads the initial idomain from the DIS/DISV/DISU package and applies
        all scheduled transitions up to and including *kper*.

        Returns
        -------
        np.ndarray, shape (nnodes,), dtype int
        """
        state = self._initial_idomain().copy()
        for k in range(kper + 1):
            entry = self.schedule.get(k, {})
            for node in entry.get("activate", []):
                state[node] = 1
            for node in entry.get("deactivate", []):
                state[node] = 0
        return state

    def write(self) -> None:
        """Serialise schedule to ``<model_ws>/<model_name>.tid.json``."""
        payload = {
            str(k): {
                "activate": v["activate"],
                "deactivate": v["deactivate"],
            }
            for k, v in sorted(self.schedule.items())
        }
        with open(self._json_path(), "w") as f:
            json.dump(payload, f, indent=2)

    @classmethod
    def load(cls, gwf: "MFModel") -> "Mf6Tid":
        """Load a previously written ``.tid.json`` sidecar.

        Parameters
        ----------
        gwf : MFModel
            Must be the same model that produced the sidecar.
        """
        tid = cls(gwf)
        with open(tid._json_path()) as f:
            raw = json.load(f)
        for k_str, entry in raw.items():
            k = int(k_str)
            tid.schedule[k] = {
                "activate": [int(n) for n in entry.get("activate", [])],
                "deactivate": [int(n) for n in entry.get("deactivate", [])],
            }
        return tid

    # ------------------------------------------------------------------
    # Internals
    # ------------------------------------------------------------------

    def _json_path(self) -> Path:
        ws = Path(self._gwf.simulation.simulation_data.mfpath.get_sim_path())
        return ws / f"{self._gwf.name}.tid.json"

    def _initial_idomain(self) -> np.ndarray:
        """Flat initial idomain from the discretization package."""
        dis = self._gwf.dis
        try:
            arr = dis.idomain.array
        except AttributeError:
            arr = None

        if arr is None:
            return np.ones(self._gwf.modelgrid.nnodes, dtype=int)
        return arr.flatten().astype(int)


# ---------------------------------------------------------------------------
# BMI runner
# ---------------------------------------------------------------------------


class Mf6TidRunner:
    """Apply an Mf6Tid schedule to a MODFLOW 6 instance via modflowapi.

    Parameters
    ----------
    sim : MFSimulation
        A FloPy MF6 simulation object (already written to disk).
    tid : Mf6Tid
        The activation/deactivation schedule.
    lib_path : str or Path, optional
        Path to the MF6 shared library (``.so`` / ``.dylib`` / ``.dll``).
        Required; no auto-detection is attempted.

    Examples
    --------
    >>> runner = Mf6TidRunner(sim, tid, lib_path="/usr/local/lib/libmf6.so")
    >>> runner.run()
    """

    def __init__(self, sim, tid: Mf6Tid, lib_path=None) -> None:
        try:
            import modflowapi as _modflowapi
        except ImportError as exc:
            raise ImportError(
                "modflowapi is required for Mf6TidRunner. "
                "Install it with: pip install modflowapi"
            ) from exc

        if lib_path is None:
            raise ValueError(
                "lib_path (path to the MF6 shared library) is required. "
                "Example: Mf6TidRunner(sim, tid, lib_path='/path/to/libmf6.so')"
            )

        self._modflowapi = _modflowapi
        self._sim = sim
        self._tid = tid
        self._lib_path = Path(lib_path)

    def run(self) -> None:
        """Step through all stress periods, applying idomain before each.

        The DIS package must have idomain=1 for every cell that will ever be
        activated by the schedule, because MF6 excludes idomain=0 cells from
        its reduced node system at initialisation time and they cannot be added
        back later.  The recommended pattern is to write the model with
        idomain=1 everywhere and use Mf6Tid to control which cells start
        inactive (add them to the kper=0 deactivate list).
        """
        sim_path = Path(self._sim.simulation_data.mfpath.get_sim_path())

        # Period structure from FloPy TDIS (avoids BMI variable name fragility)
        tdis = self._sim.tdis
        pd = tdis.perioddata.get_data()  # recarray with fields perlen/nstp/tsmult

        # working_directory sets the run dir; initialize() takes no arguments
        mf6 = self._modflowapi.ModflowApi(
            str(self._lib_path), working_directory=str(sim_path)
        )
        mf6.initialize()

        gwf = self._tid._gwf
        model_name = gwf.name.upper()

        # GWF/IBOUND is the solver's live active-cell flag (rank 1, size=nnodes).
        # DIS/IDOMAIN is rank-3 and its set_value is unsupported in MF6 6.x BMI;
        # more importantly, idomain=0 cells are excluded from the reduced node
        # system at init time and cannot be reactivated by changing DIS/IDOMAIN.
        ibound_tag = mf6.get_var_address("IBOUND", model_name)
        ibound_ptr = mf6.get_value_ptr(ibound_tag)  # live view of MF6 memory

        # idomain_state mirrors ibound_ptr values as the schedule is applied.
        # Start from all-ones (DIS should have idomain=1 for all live cells).
        idomain_state = np.ones(len(ibound_ptr), dtype=np.int32)

        for kper, row in enumerate(pd):
            entry = self._tid.schedule.get(kper, {})
            for node in entry.get("activate", []):
                idomain_state[node] = 1
            for node in entry.get("deactivate", []):
                idomain_state[node] = 0

            ibound_ptr[:] = idomain_state

            nstp = int(row["nstp"])
            for _ in range(nstp):
                mf6.update()

        mf6.finalize()
