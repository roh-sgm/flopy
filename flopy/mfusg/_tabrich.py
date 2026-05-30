"""
Shared TABRICH helpers for the USG-Transport Richards-equation tabular input
used by BCF and LPF (items 1c and 1d in ``gwf2bcf-lpf-u1.f``).

When the ``TABRICH`` option is active, USG-T reads -- immediately after the
package option line and *before* LAYCON:

* Item 1c, ``IUZONTAB``: one integer retention-zone index per node, read with
  ``U1DINT`` (a standard array control record).
* Item 1d, ``RETCRVS(NUZONES, NUTABROWS, 3)``: for each zone, ``NUTABROWS``
  rows of three free-format values -- capillary head, saturation, relative
  permeability.

The Python representation of ``RETCRVS`` is an ndarray of shape
``(nuzones, nutabrows, 3)``, matching the Fortran read order (zone outer, row
middle, the three columns inner).
"""

import numpy as np

from ..utils import Util2d


def node_count(model):
    """Total number of nodes for the IUZONTAB zone-map array."""
    if model.structured:
        nrow, ncol, nlay, _ = model.nrow_ncol_nlay_nper
        return nlay * nrow * ncol
    dis = model.get_package("DISU")
    return int(dis.nodes)


def make_iuzontab(model, iuzontab):
    """Wrap a user/loaded zone map as a Util2d integer node array."""
    if isinstance(iuzontab, Util2d):
        return iuzontab
    nodes = node_count(model)
    return Util2d(model, (nodes,), np.int32, iuzontab, "iuzontab")


def validate_retcrvs(retcrvs, nuzones, nutabrows):
    """Coerce/validate RETCRVS to shape (nuzones, nutabrows, 3)."""
    arr = np.asarray(retcrvs, dtype=np.float64)
    expected = (int(nuzones), int(nutabrows), 3)
    if arr.shape != expected:
        raise ValueError(
            f"retcrvs has shape {arr.shape}, expected {expected} "
            "(nuzones, nutabrows, 3): capillary head, saturation, "
            "relative permeability."
        )
    return arr


def write_tabrich(f_obj, iuzontab, retcrvs):
    """Write items 1c (IUZONTAB) and 1d (RETCRVS)."""
    f_obj.write(iuzontab.get_file_entry())
    nuzones, nutabrows, _ = retcrvs.shape
    for izon in range(nuzones):
        for irow in range(nutabrows):
            c0, c1, c2 = retcrvs[izon, irow]
            f_obj.write(f" {c0:.6e} {c1:.6e} {c2:.6e}\n")


def read_tabrich(f_obj, model, nuzones, nutabrows, ext_unit_dict=None):
    """Read items 1c (IUZONTAB) and 1d (RETCRVS).

    Returns
    -------
    (Util2d, np.ndarray)
        The zone-map array and the ``(nuzones, nutabrows, 3)`` retention
        curves.
    """
    nodes = node_count(model)
    iuzontab = Util2d.load(
        f_obj, model, (nodes,), np.int32, "iuzontab", ext_unit_dict
    )
    nuzones = int(nuzones)
    nutabrows = int(nutabrows)
    retcrvs = np.empty((nuzones, nutabrows, 3), dtype=np.float64)
    for izon in range(nuzones):
        for irow in range(nutabrows):
            toks = f_obj.readline().split()
            retcrvs[izon, irow, :] = (
                float(toks[0]),
                float(toks[1]),
                float(toks[2]),
            )
    return iuzontab, retcrvs
