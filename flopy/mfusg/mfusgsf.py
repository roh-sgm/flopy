"""mfusgsf module. Contains the MfUsgGsf class.

Grid Specification File (GSF) — USG-Transport / MODFLOW-USG.

Format (unstructured):
    UNSTRUCTURED [GWF]
    NNODES  NLAY  [...]
    NVERTS
    X(1)  Y(1)  Z(1)
    ...
    X(NVERTS)  Y(NVERTS)  Z(NVERTS)
    NODENO  XC  YC  ZC  LAY  NVERTS  V1  V2  ...  VN  (one line per node)

This implementation is a text round-tripper: it stores the file verbatim and
writes it back unchanged. A ``to_grid()`` method delegates to
``UnstructuredGrid.from_gridspec()`` for full geometric parsing.
"""

from __future__ import annotations

from pathlib import Path

from ..pakbase import Package
from .mfusg import MfUsg


class MfUsgGsf(Package):
    """MFUSG Grid Specification File (GSF) Package.

    Parameters
    ----------
    model : flopy.mfusg.MfUsg
        Parent model.
    lines : list of str, optional
        Raw lines of the GSF file (including the UNSTRUCTURED header).
        Stored verbatim and written unchanged.
    extension : str, optional
        File extension (default "gsf").
    unitnumber : int, optional
        File unit number (default 110).
    filenames : str or list of str, optional
        Package filename override.

    Notes
    -----
    Use ``MfUsgGsf.load(path, model)`` to populate from an existing file.
    Use ``gsf.to_grid()`` to get a fully-parsed ``UnstructuredGrid`` object
    (delegates to ``flopy.discretization.UnstructuredGrid.from_gridspec``).
    """

    def __init__(
        self,
        model,
        lines=None,
        extension="gsf",
        unitnumber=None,
        filenames=None,
    ):
        msg = (
            "Model object must be of type flopy.mfusg.MfUsg\n"
            f"but received type: {type(model)}."
        )
        assert isinstance(model, MfUsg), msg

        if unitnumber is None:
            unitnumber = MfUsgGsf._defaultunit()

        filenames = self._prepare_filenames(filenames, num=1)

        super().__init__(
            model,
            extension=extension,
            name=self._ftype(),
            unit_number=unitnumber,
            filenames=filenames,
        )
        self._generate_heading()

        self.lines = list(lines) if lines else []
        self.parent.add_package(self)

    # ------------------------------------------------------------------
    # write_file
    # ------------------------------------------------------------------

    def write_file(self, check=False):
        """Write the GSF file verbatim."""
        with open(self.fn_path, "w") as f:
            for line in self.lines:
                f.write(line if line.endswith("\n") else line + "\n")

    # ------------------------------------------------------------------
    # to_grid
    # ------------------------------------------------------------------

    def to_grid(self):
        """Return an ``UnstructuredGrid`` parsed from this GSF file.

        Requires the GSF file to exist on disk (call ``write_file()`` first
        if it was created programmatically).

        Returns
        -------
        flopy.discretization.UnstructuredGrid
        """
        from ..discretization.unstructuredgrid import UnstructuredGrid

        path = self.fn_path
        if not Path(path).exists():
            self.write_file()
        return UnstructuredGrid.from_gridspec(path, split_vertices=True)

    # ------------------------------------------------------------------
    # load
    # ------------------------------------------------------------------

    @classmethod
    def load(cls, f, model, nper=None, ext_unit_dict=None, check=False):
        """Load an existing GSF file.

        Parameters
        ----------
        f : str or file-like
            Path to the .gsf file, or an open file handle.
        model : MfUsg
            Model to attach the package to.

        Returns
        -------
        MfUsgGsf
        """
        msg = (
            "Model object must be of type flopy.mfusg.MfUsg\n"
            f"but received type: {type(model)}."
        )
        assert isinstance(model, MfUsg), msg

        if model.verbose:
            print("loading gsf package file...")

        if hasattr(f, "read"):
            fh = f
        else:
            fh = open(f, "r")

        try:
            lines = fh.readlines()
        finally:
            if not hasattr(f, "read"):
                fh.close()

        unitnumber = None
        filenames = [None]
        if ext_unit_dict is not None:
            unitnumber_ext, fname_ext = model.get_ext_dict_attr(
                ext_unit_dict, filetype=cls._ftype()
            )
            if unitnumber_ext is not None:
                unitnumber = unitnumber_ext
            if fname_ext is not None:
                filenames = [fname_ext]

        return cls(
            model,
            lines=lines,
            unitnumber=unitnumber,
            filenames=filenames,
        )

    @staticmethod
    def _ftype():
        return "GSF"

    @staticmethod
    def _defaultunit():
        return 110
