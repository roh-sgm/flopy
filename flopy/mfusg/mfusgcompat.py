"""
Compatibility guards for base MODFLOW-2005 packages under MODFLOW-USG.

STR, SUB, and SWT have no USG-T-specific FloPy implementation. Their base
MODFLOW-2005 classes assume a structured grid, and the USG-T unstructured
(node-indexed) item layouts have not been validated, so using them on an
unstructured (DISU) ``MfUsg`` model could silently read or write a non-USG-T
file. These thin subclasses are registered in ``MfUsg.mfnam_packages`` in place
of the base classes so such use **fails explicitly** (Stage 4.6A).

Scope and non-goals:

* The guard fires **only** for an unstructured ``MfUsg`` model
  (``model.version == "mfusg"`` and ``not model.structured``). On a *structured*
  ``MfUsg`` model the wrappers delegate to the base class unchanged, and plain
  ``flopy.modflow.Modflow`` use is entirely unaffected (it keeps registering the
  base classes directly).
* **SFR is intentionally not guarded.** ``examples/data/freyberg_usg`` is a
  DISU model that uses SFR and is loaded, written, and run via the base
  ``ModflowSfr2`` in ``autotest/test_usg.py`` — so unstructured SFR is
  empirically supported by the base class. SFR stays compatibility/base support
  (it is *not* promoted to a Full USG-T transport implementation).
* **FHB and GAGE are intentionally not guarded** — they round-trip via their
  base classes in the ``Ex8`` real model.

The guard raises in ``__init__`` (before the package is constructed or added to
the model) and in ``load`` (before any file is read), so no partial file is
written and no partial state is left on the model.
"""

from ..modflow.mfstr import ModflowStr
from ..modflow.mfsub import ModflowSub
from ..modflow.mfswt import ModflowSwt


def _guard_compat(model, label):
    """Raise ``NotImplementedError`` for a compat package on unstructured MfUsg.

    No-op for structured ``MfUsg`` models and for non-``MfUsg`` models (e.g.
    plain ``flopy.modflow.Modflow``), so base-class behavior is untouched there.
    """
    if getattr(model, "version", None) == "mfusg" and not getattr(
        model, "structured", True
    ):
        raise NotImplementedError(
            f"{label} is a compatibility-only package for MODFLOW-USG: USG-T "
            "unstructured/DISU support is not implemented, so loading or writing "
            f"{label} on an unstructured MfUsg model could produce a non-USG-T "
            "file. Use CLN for the coupling, implement an MfUsg-specific "
            f"{label} package, or keep the file as a raw external input if you "
            "need it."
        )


class MfUsgStr(ModflowStr):
    """STR guarded for unstructured MfUsg (see module docstring)."""

    def __init__(self, model, *args, **kwargs):
        _guard_compat(model, "STR")
        super().__init__(model, *args, **kwargs)

    @classmethod
    def load(cls, f, model, *args, **kwargs):
        _guard_compat(model, "STR")
        return ModflowStr.load(f, model, *args, **kwargs)


class MfUsgSub(ModflowSub):
    """SUB guarded for unstructured MfUsg (see module docstring)."""

    def __init__(self, model, *args, **kwargs):
        _guard_compat(model, "SUB")
        super().__init__(model, *args, **kwargs)

    @classmethod
    def load(cls, f, model, *args, **kwargs):
        _guard_compat(model, "SUB")
        return ModflowSub.load(f, model, *args, **kwargs)


class MfUsgSwt(ModflowSwt):
    """SWT guarded for unstructured MfUsg (see module docstring)."""

    def __init__(self, model, *args, **kwargs):
        _guard_compat(model, "SWT")
        super().__init__(model, *args, **kwargs)

    @classmethod
    def load(cls, f, model, *args, **kwargs):
        _guard_compat(model, "SWT")
        return ModflowSwt.load(f, model, *args, **kwargs)
