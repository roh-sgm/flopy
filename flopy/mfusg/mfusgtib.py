"""mfusgtib module.  Contains the MfUsgTib class. Note that the user can access
the MfUsgTib class as `flopy.mfusg.MfUsgTib`.

Transient Ibound (TIB) Package — USG-Transport.

Reference: Fortran subroutine GWF2TIB1RP in glo2basu1.f of USG-Transport.
Per-stress-period format (free format):

    NIB0 NIB1 NIBM1 [NICB0 NICB1 NICBM1]   <- 3 ints without transport, 6 with
    [NIB0 node numbers read via U1DINT — i.e. an INTERNAL/EXTERNAL array block]
    [NIB1 lines, each: ICELL  [HEAD value | AVHEAD]]
    [NIBM1 lines, each: ICELL [HEAD value | AVHEAD]]
    [with transport: concentration blocks for NICB0, NICB1, NICBM1]

The minimal class below preserves per-SP content verbatim (text block
round-trip). Semantic parsing can be added later if needed. This is enough
for models that load, rewrite, and run through MfUsg without manual patching.
"""

from __future__ import annotations

from ..pakbase import Package
from .mfusg import MfUsg


class MfUsgTib(Package):
    """MFUSG Transient Ibound (TIB) Package Class.

    Parameters
    ----------
    model : flopy.mfusg.MfUsg
        Parent model.
    blocks : dict[int, str]
        Per-stress-period block text, keyed by zero-based kper. Each value is
        the full content for that SP starting with the header line
        (NIB0 NIB1 NIBM1 [NICB0 NICB1 NICBM1]) and ending at the last line of
        that period's data. Trailing newline on each block is fine.
    extension : str, optional
        File extension (default "tib").
    unitnumber : int, optional
        File unit number (default 102, matching common GMS/GV8 convention).
    filenames : str or list of str, optional
        Package filename override.

    Notes
    -----
    This is a faithful text round-tripper for TIB. It supports models where
    TIB content comes from an existing file (the common USG-T workflow). For
    models that need TIB authored programmatically, use `blocks` with
    hand-constructed strings per SP.
    """

    def __init__(
        self,
        model,
        blocks=None,
        raw_body=None,
        extension="tib",
        unitnumber=None,
        filenames=None,
    ):
        msg = (
            "Model object must be of type flopy.mfusg.MfUsg\n"
            f"but received type: {type(model)}."
        )
        assert isinstance(model, MfUsg), msg

        if unitnumber is None:
            unitnumber = MfUsgTib._defaultunit()

        filenames = self._prepare_filenames(filenames, num=1)

        super().__init__(
            model,
            extension=extension,
            name=self._ftype(),
            unit_number=unitnumber,
            filenames=filenames,
        )
        self._generate_heading()

        self.blocks = dict(blocks) if blocks else {}
        self.raw_body = raw_body
        self.parent.add_package(self)

    def write_file(self, check=False):
        """Write the TIB package file."""
        nper = self.parent.nper
        with open(self.fn_path, "w") as f:
            f.write(f"{self.heading}\n")
            if self.raw_body is not None:
                f.write(self.raw_body)
                if self.raw_body and not self.raw_body.endswith("\n"):
                    f.write("\n")
                return
            for kper in range(nper):
                block = self.blocks.get(kper, "")
                if not block:
                    # No TIB entry for this SP: emit a header of zeros so USG-T
                    # advances. We don't know whether the parent BCT is active,
                    # so emit 6 zeros — the 3-zero form is also valid but 6 is
                    # always safe when a BCT is present.
                    inbct = getattr(self.parent, "bct", None)
                    n = 6 if inbct is not None else 3
                    f.write((" 0" * n) + "\n")
                else:
                    if not block.endswith("\n"):
                        block = block + "\n"
                    f.write(block)

    @classmethod
    def load(cls, f, model, nper=None, ext_unit_dict=None, check=False):
        """Load an existing TIB file into `blocks`."""
        msg = (
            "Model object must be of type flopy.mfusg.MfUsg\n"
            f"but received type: {type(model)}."
        )
        assert isinstance(model, MfUsg), msg

        if model.verbose:
            print("loading tib package file...")

        if hasattr(f, "read"):
            fh = f
            filename = getattr(f, "name", None)
        else:
            filename = f
            fh = open(f, "r")

        try:
            lines = fh.readlines()
        finally:
            if not hasattr(f, "read"):
                fh.close()

        # Strip comment header lines (any leading lines starting with "#")
        i = 0
        while i < len(lines) and lines[i].lstrip().startswith("#"):
            i += 1

        raw_body = "".join(lines[i:])

        # Construct package
        # Preserve original unit + filename from NAM
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
            raw_body=raw_body,
            unitnumber=unitnumber,
            filenames=filenames,
        )

    @staticmethod
    def _ftype():
        return "TIB"

    @staticmethod
    def _defaultunit():
        return 102
