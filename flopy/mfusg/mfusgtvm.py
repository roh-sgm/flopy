"""mfusgtvm module. Contains the MfUsgTvm class.

Time-Variant Materials (TVM2) Package — USG-Transport.

Reference: Fortran subroutines TVMU2AR / TVMU2RP / TVMU2READ in tvmu2.f
of USG-Transport 2.7.

File structure
--------------
0.  [#Text]  — optional comment lines, must be at the top

1.  ITVMPRINT  TVMLOGBASEHK  TVMLOGBASEVKA  TVMLOGBASESS  TVMLOGBASESY
    TVMDDFTR  [TVMLOGBASEPOR]
    — 7 values when BCT transport is active, 6 otherwise.

For each SP boundary (nper + 1 total — one for the *start* of SP 1, then
one for the *end* of each stress period):

2.  NTVMHK  NTVMVKA  NTVMSS  NTVMSY  NTVMDDFTR  [NTVMPOR]
    — 6 values with transport, 5 without.

3.  ITVMHK  HKNEW        (× NTVMHK)
4.  ITVMVKA  VKANEW       (× NTVMVKA)
5.  ITVMSS  SSNEW         (× NTVMSS)
6.  ITVMSY  SYNEW         (× NTVMSY)
7.  ITVMDDFTR  DDFTRNEW   (× NTVMDDFTR)
8.  ITVMPOR  PORNEW        (× NTVMPOR, transport only)

All node numbers in the file are 1-based; FloPy stores them 0-based.
"""

from __future__ import annotations

import numpy as np

from ..pakbase import Package
from ..utils.recarray_utils import create_empty_recarray
from .mfusg import MfUsg

# Property names ordered exactly as USG-T reads/writes them.
_PROPS = ("hk", "vka", "ss", "sy", "ddftr", "por")
# Dtype for per-node records (0-based node stored internally).
_DTYPE = np.dtype([("node", np.int32), ("value", np.float64)])


def _next_data_line(fobj):
    """Return next non-blank line from *fobj*, or ``None`` at EOF."""
    while True:
        line = fobj.readline()
        if not line:
            return None
        s = line.strip()
        if s:
            return s


class MfUsgTvm(Package):
    """MFUSG Time-Variant Materials (TVM2) Package.

    Parameters
    ----------
    model : MfUsg
        Parent model object.
    itvmprint : int, optional
        Listing-file verbosity (0 = quiet, 1 = summary, 2 = detailed).
        Default 0.
    tvmlogbasehk : float, optional
        HK interpolation control: 0 = linear, > 1 = logarithmic (value is
        the base), < 0 = step function. Default 0.
    tvmlogbasevka : float, optional
        VKA interpolation control. Same semantics as *tvmlogbasehk*.
    tvmlogbasess : float, optional
        Ss interpolation control. Step function (< 0) not supported by USG-T.
    tvmlogbasesy : float, optional
        Sy interpolation control. Step function not supported.
    tvmddftr : float, optional
        DDFTR (dual-domain transfer coefficient) interpolation control.
    tvmlogbasepor : float, optional
        Porosity interpolation control. Step function not supported.
        Written only when the model has an active BCT transport package.
    stress_period_data : dict, optional
        Mapping ``{ibnd: sp_dict}`` where *ibnd* is the **boundary index**
        (0 = start of SP 1, k = end of SP k, up to *nper*).  Each
        *sp_dict* is a ``dict[str, array-like]`` with any subset of the
        keys ``'hk'``, ``'vka'``, ``'ss'``, ``'sy'``, ``'ddftr'``,
        ``'por'``.  Values are recarrays (or anything coercible to one)
        with dtype ``[('node', int32), ('value', float64)]``.  Nodes are
        **0-based**; the file writer adds 1 automatically.  Missing
        boundaries default to all-zero counts (no changes that boundary).
    extension : str, optional
        File extension. Default ``'tvm'``.
    unitnumber : int, optional
        FORTRAN unit number. Default 109.
    filenames : str or list of str, optional
        Package filename override.

    Notes
    -----
    TVM2 requires the LPF package with IKCFLAG = 0 (nodal hydraulic
    conductivities).  Step interpolation (log-base < 0) is not supported
    by USG-T for Ss, Sy, or porosity.  Log bases must be > 1 or exactly
    0 (linear).

    The file contains **nper + 1** SP-boundary blocks: one for the *start*
    of stress period 1 (initial conditions) followed by one for the *end*
    of each stress period.  USG-T interpolates within each SP between the
    boundary values bracketing that period.

    When the model has no BCT package the ``'por'`` key and
    ``tvmlogbasepor`` are silently omitted from the written file.

    Examples
    --------
    >>> import numpy as np
    >>> import flopy
    >>> m = flopy.mfusg.MfUsg(modelname="test")
    >>> dtype = flopy.mfusg.MfUsgTvm.dtype()
    >>> # Three SP model; change HK at boundary 1 (end of SP 1)
    >>> sp = {1: {"hk": np.array([(0, 1e-4), (4, 2e-4)], dtype=dtype)}}
    >>> tvm = flopy.mfusg.MfUsgTvm(m, tvmlogbasehk=10.0,
    ...                             stress_period_data=sp)
    """

    _PROPS = _PROPS
    _DTYPE = _DTYPE

    def __init__(
        self,
        model,
        itvmprint=0,
        tvmlogbasehk=0.0,
        tvmlogbasevka=0.0,
        tvmlogbasess=0.0,
        tvmlogbasesy=0.0,
        tvmddftr=0.0,
        tvmlogbasepor=0.0,
        stress_period_data=None,
        extension="tvm",
        unitnumber=None,
        filenames=None,
    ):
        assert isinstance(model, MfUsg), (
            "Model object must be of type flopy.mfusg.MfUsg\n"
            f"but received type: {type(model)}."
        )

        if unitnumber is None:
            unitnumber = MfUsgTvm._defaultunit()

        filenames = self._prepare_filenames(filenames, num=1)

        super().__init__(
            model,
            extension=extension,
            name=self._ftype(),
            unit_number=unitnumber,
            filenames=filenames,
        )
        self._generate_heading()

        self.itvmprint = int(itvmprint)
        self.tvmlogbasehk = float(tvmlogbasehk)
        self.tvmlogbasevka = float(tvmlogbasevka)
        self.tvmlogbasess = float(tvmlogbasess)
        self.tvmlogbasesy = float(tvmlogbasesy)
        self.tvmddftr = float(tvmddftr)
        self.tvmlogbasepor = float(tvmlogbasepor)

        self.stress_period_data: dict[int, dict[str, np.recarray]] = {}
        for ibnd, sp in (stress_period_data or {}).items():
            self.stress_period_data[int(ibnd)] = self._normalise_sp(sp)

        self.parent.add_package(self)

    # ------------------------------------------------------------------
    # Class helpers
    # ------------------------------------------------------------------

    @classmethod
    def dtype(cls) -> np.dtype:
        """Recarray dtype for TVM node-value records."""
        return cls._DTYPE

    @classmethod
    def _normalise_sp(cls, sp_dict: dict) -> dict[str, np.recarray]:
        """Convert any list/array values in *sp_dict* to recarrays."""
        out: dict[str, np.recarray] = {}
        for prop in cls._PROPS:
            data = sp_dict.get(prop)
            if data is None:
                continue
            if isinstance(data, np.recarray):
                out[prop] = data
            else:
                out[prop] = np.array(list(data), dtype=cls._DTYPE).view(
                    np.recarray
                )
        return out

    def _has_transport(self) -> bool:
        """True when the model has an active BCT transport package."""
        return self.parent.get_package("BCT") is not None

    @staticmethod
    def _sp_label(ibnd: int) -> str:
        """Return the SP-boundary comment label used by Vistas."""
        if ibnd == 0:
            return "Stress period number 1 start"
        if ibnd == 1:
            return "Stress period number 1 end"
        return f"Stress period number {ibnd}"

    # ------------------------------------------------------------------
    # write_file
    # ------------------------------------------------------------------

    def write_file(self, check: bool = False) -> None:
        """Write the TVM2 package file."""
        nper = self.parent.nper
        transport = self._has_transport()
        active_props = self._PROPS if transport else self._PROPS[:-1]

        with open(self.fn_path, "w") as f:
            f.write(f"{self.heading}\n")

            # Item 1: global interpolation controls
            bases = [
                self.tvmlogbasehk,
                self.tvmlogbasevka,
                self.tvmlogbasess,
                self.tvmlogbasesy,
                self.tvmddftr,
            ]
            if transport:
                bases.append(self.tvmlogbasepor)
            f.write(f"{self.itvmprint:10d}")
            for b in bases:
                f.write(f"{b:10g}")
            f.write("\n")

            # nper + 1 SP boundaries
            for ibnd in range(nper + 1):
                sp = self.stress_period_data.get(ibnd, {})
                counts = [len(sp.get(p, [])) for p in active_props]
                f.write("".join(f"{n:10d}" for n in counts))
                f.write(f"    {self._sp_label(ibnd)}\n")
                for prop in active_props:
                    for row in sp.get(prop, []):
                        node_1based = int(row["node"]) + 1
                        value = float(row["value"])
                        f.write(f"{node_1based:10d}  {value:13.6g}\n")

    # ------------------------------------------------------------------
    # load
    # ------------------------------------------------------------------

    @classmethod
    def load(
        cls,
        f,
        model,
        nper: int | None = None,
        ext_unit_dict: dict | None = None,
        check: bool = False,
    ) -> "MfUsgTvm":
        """Load an existing TVM2 file.

        Parameters
        ----------
        f : str or file-like
            Path or open file object.
        model : MfUsg
            Parent model.
        nper : int, optional
            Number of stress periods. Inferred from *model* when not given.
        ext_unit_dict : dict, optional
            External unit dictionary from name-file parsing.

        Returns
        -------
        MfUsgTvm
        """
        assert isinstance(model, MfUsg), (
            "Model object must be of type flopy.mfusg.MfUsg\n"
            f"but received type: {type(model)}."
        )
        if model.verbose:
            print("loading tvm package file...")

        if nper is None:
            nper = model.nper

        openfile = not hasattr(f, "read")
        if openfile:
            f = open(f, "r")

        try:
            # Skip leading comment (#) lines
            line = f.readline()
            while line and line.lstrip().startswith("#"):
                line = f.readline()

            if not line:
                raise EOFError("TVM file is empty or contains only comments.")

            # Item 1: global header
            # Filter to contiguous numeric tokens; stop at first non-numeric.
            num_tokens: list[str] = []
            for tok in line.split():
                try:
                    float(tok)
                    num_tokens.append(tok)
                except ValueError:
                    break

            if len(num_tokens) < 6:
                raise ValueError(
                    f"TVM global header must have ≥ 6 numeric values; "
                    f"got {len(num_tokens)}:\n  {line.strip()}"
                )
            transport_in_file = len(num_tokens) >= 7

            itvmprint = int(float(num_tokens[0]))
            tvmlogbasehk = float(num_tokens[1])
            tvmlogbasevka = float(num_tokens[2])
            tvmlogbasess = float(num_tokens[3])
            tvmlogbasesy = float(num_tokens[4])
            tvmddftr = float(num_tokens[5])
            tvmlogbasepor = float(num_tokens[6]) if transport_in_file else 0.0

            active_props = cls._PROPS if transport_in_file else cls._PROPS[:-1]
            n_sp_fields = len(active_props)

            # Read nper + 1 SP boundaries
            stress_period_data: dict[int, dict] = {}
            for ibnd in range(nper + 1):
                hdr = _next_data_line(f)
                if hdr is None:
                    break  # file ends before all boundaries — tolerate

                # Parse up to n_sp_fields integers; ignore trailing label text
                counts: list[int] = []
                for tok in hdr.split()[:n_sp_fields]:
                    try:
                        counts.append(int(tok))
                    except ValueError:
                        break
                # Pad with zeros if fewer fields than expected
                counts.extend([0] * (n_sp_fields - len(counts)))

                sp: dict[str, np.recarray] = {}
                for prop, n in zip(active_props, counts):
                    if n <= 0:
                        continue
                    arr = create_empty_recarray(n, cls._DTYPE)
                    for i in range(n):
                        rec = _next_data_line(f)
                        if rec is None:
                            raise EOFError(
                                f"TVM boundary {ibnd}, property '{prop}': "
                                f"expected {n} records but file ended early."
                            )
                        toks = rec.split()
                        arr[i]["node"] = int(toks[0]) - 1  # 1-based → 0-based
                        arr[i]["value"] = float(toks[1])
                    sp[prop] = arr

                stress_period_data[ibnd] = sp

        finally:
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
            itvmprint=itvmprint,
            tvmlogbasehk=tvmlogbasehk,
            tvmlogbasevka=tvmlogbasevka,
            tvmlogbasess=tvmlogbasess,
            tvmlogbasesy=tvmlogbasesy,
            tvmddftr=tvmddftr,
            tvmlogbasepor=tvmlogbasepor,
            stress_period_data=stress_period_data,
            unitnumber=unitnumber,
            filenames=filenames,
        )

    # ------------------------------------------------------------------
    # Package metadata
    # ------------------------------------------------------------------

    @staticmethod
    def _ftype() -> str:
        return "TVM"

    @staticmethod
    def _defaultunit() -> int:
        return 109
