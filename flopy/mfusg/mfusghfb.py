"""
mfusghfb module.  Contains the MfUsgHfb class.
"""

import numpy as np

from ..modflow.mfhfb import ModflowHfb
from ..pakbase import Package
from ..utils.recarray_utils import create_empty_recarray
from ._usgt_list import begin_list_block
from ._usgt_parameters import (
    read_active_list_parameters,
    read_list_parameter_header,
    write_active_list_parameters,
    write_list_parameter_header,
)


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
    nphfb : int
        Number of HFB barriers defined by named parameters (dataset 1). Set by
        :meth:`load` when the file is parameterized; from-scratch parameter
        authoring is not supported (see Notes).
    mxfb : int
        Maximum number of barriers defined by parameters (MXFBP in the Fortran).
        Preserved on load/write.
    nacthfb : int
        Number of active HFB parameters (dataset 5). Set by :meth:`load`.
    parameters : dict or None
        Preserved HFB list-parameter definitions, keyed by parameter name:
        ``{name: {"partyp": str, "parval": str, "nlst": int, "data": recarray}}``
        where ``data`` holds the NLST barrier rows (zero-based, same dtype as
        *hfb_data*). Set by :meth:`load` when ``NPHFB > 0``.
    acthfb_names : list of str or None
        Names of the active parameters (dataset 6), in file order. Set by
        :meth:`load`.
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
    HFB uses MODFLOW **list** parameters (``UPARLSTRP`` / ``UPARLSTSUB`` in
    ``parutl7.f``): each parameter owns NLST barrier rows and its value scales
    the barrier ``HYDCHR`` factor. As of Stage 4.4B, parameterized HFB files are
    **preserved** (load -> write -> reload keeps ``NPHFB``, the definition blocks,
    and the active-parameter records). HFB does not support parameter
    ``INSTANCES`` (the Fortran aborts), so neither does this class.

    Two parameter modes are intentionally not supported and fail explicitly with
    ``NotImplementedError``: from-scratch parameter *authoring* (``nphfb > 0``
    with no loaded definitions) and ``TRANSIENT_HFB`` combined with parameters
    (``NPHFB > 0``), which would re-read/redefine parameters each stress period.
    Non-parametric HFB (static and ``TRANSIENT_HFB``) is unchanged.

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
        nphfb=0,
        mxfb=0,
        parameters=None,
        acthfb_names=None,
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
            nphfb=nphfb,
            mxfb=mxfb,
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

        # Preserved HFB list-parameter definitions and active-parameter records
        # (set by load when NPHFB > 0).
        self.parameters = parameters
        self.acthfb_names = acthfb_names if acthfb_names is not None else []

        # The parent always adds itself; remove then re-add only if requested.
        if not add_package:
            model.pop_package(self._ftype())

    # ------------------------------------------------------------------
    # write_file
    # ------------------------------------------------------------------

    def write_file(self):
        """Write HFB6 package file honouring TRANSIENT_HFB when set."""
        preserve = self.nphfb > 0 and self.parameters is not None
        if self.nphfb > 0 and self.transient:
            raise NotImplementedError(
                "MfUsgHfb.write_file does not support TRANSIENT_HFB combined "
                "with HFB parameters (NPHFB > 0): the Fortran re-reads the "
                "parameter definitions every stress period (UPARLSTRP ITERP=1), "
                "which would redefine them. Use non-transient parameterized HFB, "
                "or NPHFB=0 for transient barriers."
            )
        if self.nphfb > 0 and self.parameters is None:
            raise NotImplementedError(
                "MfUsgHfb.write_file cannot author HFB parameter definitions "
                "from scratch (NPHFB > 0 without loaded parameter data). "
                "Parameter preservation is supported for files read by "
                "MfUsgHfb.load; for from-scratch input use NPHFB=0."
            )
        if self.nphfb > 0 and self.nacthfb != len(self.acthfb_names):
            raise ValueError(
                "MfUsgHfb.write_file: nacthfb "
                f"({self.nacthfb}) must equal the number of active parameter "
                f"names ({len(self.acthfb_names)})."
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

            # Parameterized HFB (non-transient): dataset 2-3 parameter
            # definitions, then dataset 4 non-parametric barriers, then
            # dataset 5-6 active-parameter records.
            if preserve:
                for name, pdef in self.parameters.items():
                    write_list_parameter_header(
                        f, name, pdef["partyp"], pdef["parval"], pdef["nlst"]
                    )
                    self._write_hfb_rows(f, pdef["data"], structured)
                self._write_hfb_rows(f, self.hfb_data, structured)
                f.write(f"{self.nacthfb:10d}\n")
                write_active_list_parameters(f, self.acthfb_names)
                return

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
                    sp_data = np.array(sp_data, dtype=self.hfb_data.dtype).view(
                        np.recarray
                    )
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
    def _read_hfb_rows(f, nrows, dtype, structured, model=None, ext_unit_dict=None):
        data = create_empty_recarray(nrows, dtype)
        if nrows == 0:
            return data
        # A barrier list may begin with MODFLOW list controls (SFAC, OPEN/CLOSE,
        # EXTERNAL) per SGWF2HFB7RL/RLU. begin_list_block consumes them and
        # returns the source to read the rows from, the SFAC scale, and the first
        # already-read row. SFAC scales the barrier FACTOR/HYDCHR (HFB(6,II) =
        # FACTOR*SFAC in the Fortran).
        source, sfac, line, to_close = begin_list_block(
            f, model=model, ext_unit_dict=ext_unit_dict, package="HFB"
        )
        for i in range(nrows):
            if i > 0:
                line = source.readline()
            t = line.split()
            if structured:
                data[i] = (
                    int(t[0]) - 1,
                    int(t[1]) - 1,
                    int(t[2]) - 1,
                    int(t[3]) - 1,
                    int(t[4]) - 1,
                    float(t[5]) * sfac,
                )
            else:
                data[i] = (int(t[0]) - 1, int(t[1]) - 1, float(t[2]) * sfac)
        if to_close is not None:
            to_close.close()
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

        dtype = cls.get_default_dtype(structured=structured)
        stress_period_data = {}
        parameters = None
        acthfb_names = []
        nacthfb = 0

        if nphfb > 0:
            if transient:
                raise NotImplementedError(
                    "MfUsgHfb.load does not support TRANSIENT_HFB combined with "
                    "HFB parameters (NPHFB > 0): the Fortran re-reads/redefines "
                    "the parameter definitions every stress period (UPARLSTRP "
                    "ITERP=1). This combination is out of scope."
                )
            # Dataset 2-3: NPHFB list-parameter definitions (header + NLST rows).
            parameters = {}
            for _ in range(nphfb):
                name, partyp, parval, nlst, numinst = read_list_parameter_header(
                    f.readline()
                )
                if numinst > 0:
                    raise NotImplementedError(
                        "MfUsgHfb.load: HFB parameter INSTANCES are not "
                        "supported (gwf2hfb7u1.f aborts when NUMINST>0)."
                    )
                data = cls._read_hfb_rows(
                    f, nlst, dtype, structured, model, ext_unit_dict
                )
                parameters[name] = {
                    "partyp": partyp,
                    "parval": parval,
                    "nlst": nlst,
                    "data": data,
                }
            # Dataset 4: barriers not defined by parameters.
            hfb_data = cls._read_hfb_rows(
                f, nhfbnp, dtype, structured, model, ext_unit_dict
            )
            # Dataset 5-6: number of active parameters and their names.
            nacthfb = int(f.readline().split()[0])
            acthfb_names = read_active_list_parameters(f, nacthfb)
        elif transient:
            hfb_data = create_empty_recarray(0, dtype)
            for kper in range(nper):
                line = f.readline()
                if not line:
                    break
                ihfbrd = int(line.split()[0])
                if ihfbrd > 0:
                    sp_arr = cls._read_hfb_rows(
                        f, nhfbnp, dtype, structured, model, ext_unit_dict
                    )
                    stress_period_data[kper] = sp_arr
                    if len(hfb_data) == 0:
                        hfb_data = sp_arr.copy()
        else:
            hfb_data = cls._read_hfb_rows(
                f, nhfbnp, dtype, structured, model, ext_unit_dict
            )

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
            nphfb=nphfb,
            mxfb=mxfb,
            nhfbnp=nhfbnp,
            hfb_data=hfb_data,
            nacthfb=nacthfb,
            parameters=parameters,
            acthfb_names=acthfb_names,
            options=options,
            transient=transient,
            stress_period_data=stress_period_data,
            extension="hfb",
            unitnumber=unitnumber,
            filenames=filenames,
        )
