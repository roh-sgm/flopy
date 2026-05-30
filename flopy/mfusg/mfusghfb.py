"""
mfusghfb module.  Contains the MfUsgHfb class.
"""

import numpy as np

from ..modflow.mfhfb import ModflowHfb
from ..pakbase import Package
from ..utils.recarray_utils import create_empty_recarray


class MfUsgHfb(ModflowHfb):
    """MODFLOW-USG Hydraulic Flow Barrier (HFB6) Package.

    Extends :class:`ModflowHfb` with TRANSIENT_HFB support and proper
    registration for USG-T models.  The HFB applies only to the
    Groundwater Flow Process; it has no transport AUX variables.

    Parameters
    ----------
    model : MfUsg
        Model object.
    nhfbnp : int
        Number of non-parametric HFB barriers (dataset 4).
    hfb_data : np.recarray, optional
        Static barrier definitions.  Fields depend on grid type:
        - Unstructured: ``node1``, ``node2``, ``hydchr``
        - Structured:   ``k``, ``irow1``, ``icol1``, ``irow2``, ``icol2``, ``hydchr``
        Indices are **zero-based** (converted to 1-based on write).
    nacthfb : int
        Number of active HFB parameters per stress period (default 0).
    options : list of str, optional
        Extra options written to the header line (e.g. ``['NOPRINT']``).
        Do **not** include ``TRANSIENT_HFB`` here; use *transient=True*.
    transient : bool
        If True, write ``TRANSIENT_HFB`` in the header and include
        per-stress-period IHFBRD + optional updated barrier data.
    stress_period_data : dict, optional
        Only used when *transient=True*.  Keys are zero-based stress
        period indices.  Values are either:
        - ``None`` / missing key → IHFBRD = -1 (reuse previous SP data)
        - np.recarray with same dtype as *hfb_data* → IHFBRD = 1 and NHFBNP
          records are written.  Indices are zero-based.
    extension : str
        File extension (default ``'hfb'``).
    unitnumber : int, optional
    filenames : str or list of str, optional

    Notes
    -----
    Parameters (NPHFB > 0) are not yet supported; set NPHFB = 0 (default).

    Examples
    --------
    >>> import flopy, numpy as np
    >>> m = flopy.mfusg.MfUsg()
    >>> hfb = flopy.mfusg.MfUsgHfb.load("model.hfb", m, nper=12)
    """

    def __init__(
        self,
        model,
        nhfbnp=0,
        hfb_data=None,
        nacthfb=0,
        options=None,
        transient=False,
        stress_period_data=None,
        extension="hfb",
        unitnumber=None,
        filenames=None,
        add_package=True,
    ):
        if unitnumber is None:
            unitnumber = ModflowHfb._defaultunit()
        if options is None:
            options = []

        # Build hfb_data recarray
        dtype = MfUsgHfb.get_default_dtype(structured=model.structured)
        if hfb_data is None and stress_period_data:
            first = next(
                (v for _, v in sorted(stress_period_data.items()) if v is not None),
                None,
            )
            if first is not None:
                hfb_data = first
                nhfbnp = len(first)
        if hfb_data is None:
            hfb_data = create_empty_recarray(0, dtype)
            nhfbnp = 0
        elif not isinstance(hfb_data, np.recarray):
            hfb_data = np.array(hfb_data, dtype=dtype).view(np.recarray)
        if nhfbnp == 0:
            nhfbnp = len(hfb_data)

        # Parent __init__ handles Package registration, self.hfb_data, etc.
        # add_package is intercepted here so we can control it independently.
        super().__init__(
            model,
            nhfbnp=nhfbnp,
            hfb_data=hfb_data,
            nacthfb=nacthfb,
            options=options,
            extension=extension,
            unitnumber=unitnumber,
            filenames=self._prepare_filenames(filenames),
        )

        self.transient = transient
        self.stress_period_data = stress_period_data or {}

        # The parent always adds itself; remove then re-add only if requested.
        if not add_package:
            model.pop_package(self._ftype())

    # ------------------------------------------------------------------
    # write_file
    # ------------------------------------------------------------------

    def write_file(self):
        """Write HFB6 package file honouring TRANSIENT_HFB when set."""
        if self.nphfb > 0:
            raise NotImplementedError(
                "MfUsgHfb.write_file does not yet support parameterized HFB "
                "definitions (NPHFB > 0)."
            )
        structured = self.parent.structured
        nper = self.parent.nper

        with open(self.fn_path, "w") as f:
            f.write(f"{self.heading}\n")

            # Dataset 1: NPHFB MXFB NHFBNP [OPTIONS]
            header = f"{self.nphfb:10d}{self.mxfb:10d}{self.nhfbnp:10d}"
            for opt in self.options:
                header += f"  {opt}"
            if self.transient:
                header += "  TRANSIENT_HFB"
            f.write(header + "\n")

            if not self.transient:
                self._write_hfb_rows(f, self.hfb_data, structured)
                return

            # Fortran reads IHFBRD at the start of each stress period.  When
            # IHFBRD > 0 it then reads NHFBNP HFB rows; IHFBRD is a flag, not a
            # row count.
            for kper in range(nper):
                sp_data = self.stress_period_data.get(kper, None)
                if sp_data is None and kper == 0 and self.nhfbnp > 0:
                    sp_data = self.hfb_data

                if sp_data is None:
                    f.write("-1\n")
                    continue

                if not isinstance(sp_data, np.recarray):
                    sp_data = np.array(sp_data, dtype=self.hfb_data.dtype).view(np.recarray)
                if len(sp_data) != self.nhfbnp:
                    raise ValueError(
                        "Transient HFB stress_period_data entries must contain "
                        f"NHFBNP={self.nhfbnp} rows; got {len(sp_data)} for kper {kper}."
                    )
                f.write("1\n")
                self._write_hfb_rows(f, sp_data, structured)

    @staticmethod
    def _write_hfb_rows(f, data, structured):
        for row in data:
            if structured:
                vals = [
                    f"{int(row['k']) + 1:10d}",
                    f"{int(row['irow1']) + 1:10d}",
                    f"{int(row['icol1']) + 1:10d}",
                    f"{int(row['irow2']) + 1:10d}",
                    f"{int(row['icol2']) + 1:10d}",
                    f"{float(row['hydchr']):13.6g}",
                ]
            else:
                vals = [
                    f"{int(row['node1']) + 1:10d}",
                    f"{int(row['node2']) + 1:10d}",
                    f"{float(row['hydchr']):13.6g}",
                ]
            f.write("".join(vals) + "\n")

    @staticmethod
    def _read_hfb_rows(f, nrows, dtype, structured):
        data = create_empty_recarray(nrows, dtype)
        for i in range(nrows):
            t = f.readline().split()
            if structured:
                data[i] = (
                    int(t[0]) - 1,
                    int(t[1]) - 1,
                    int(t[2]) - 1,
                    int(t[3]) - 1,
                    int(t[4]) - 1,
                    float(t[5]),
                )
            else:
                data[i] = (int(t[0]) - 1, int(t[1]) - 1, float(t[2]))
        return data

    # ------------------------------------------------------------------
    # load
    # ------------------------------------------------------------------

    @classmethod
    def load(cls, f, model, nper=None, ext_unit_dict=None):
        """Load an HFB6 package from a USG-T file.

        Parameters
        ----------
        f : str or file-like
        model : MfUsg
        nper : int, optional
        ext_unit_dict : dict, optional

        Returns
        -------
        MfUsgHfb
        """
        if model.verbose:
            print("loading mfusg hfb6 package file...")

        if nper is None:
            nper = model.nper
        structured = model.structured

        openfile = not hasattr(f, "read")
        if openfile:
            f = open(f, "r")

        # Skip comments
        line = f.readline()
        while line.startswith("#"):
            line = f.readline()

        # Dataset 1: NPHFB MXFB NHFBNP [OPTIONS]
        tokens = line.split()
        nphfb = int(tokens[0])
        mxfb = int(tokens[1])
        nhfbnp = int(tokens[2])
        options = []
        transient = False
        for tok in tokens[3:]:
            if tok.upper() == "TRANSIENT_HFB":
                transient = True
            elif tok.upper() not in ("OPTIONS",):
                options.append(tok)

        if nphfb > 0:
            raise NotImplementedError(
                "MfUsgHfb.load: parametric HFBs (NPHFB > 0) not yet supported."
            )

        dtype = cls.get_default_dtype(structured=structured)
        stress_period_data = {}
        if transient:
            hfb_data = create_empty_recarray(0, dtype)
            for kper in range(nper):
                line = f.readline()
                if not line:
                    break
                ihfbrd = int(line.split()[0])
                if ihfbrd > 0:
                    sp_arr = cls._read_hfb_rows(f, nhfbnp, dtype, structured)
                    stress_period_data[kper] = sp_arr
                    if len(hfb_data) == 0:
                        hfb_data = sp_arr.copy()
        else:
            hfb_data = cls._read_hfb_rows(f, nhfbnp, dtype, structured)

        if openfile:
            f.close()

        unitnumber = None
        filenames = [None]
        if ext_unit_dict is not None:
            unitnumber, filenames[0] = model.get_ext_dict_attr(
                ext_unit_dict, filetype=cls._ftype()
            )

        return cls(
            model,
            nhfbnp=nhfbnp,
            hfb_data=hfb_data,
            nacthfb=0,
            options=options,
            transient=transient,
            stress_period_data=stress_period_data,
            extension="hfb",
            unitnumber=unitnumber,
            filenames=filenames,
        )
