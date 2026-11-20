"""
mfusgddf module.

Contains the MfUsgDdf class. Note that the user can access
the MfUsgDdf class as `flopy.mfusg.MfUsgDdf`.
"""

import math
import numbers
import os
import re

import numpy as np

from ..pakbase import Package
from ..utils.utils_def import get_unitnumber_from_ext_unit_dict
from .mfusg import MfUsg

# Fortran F30.0 / list-directed real and I30 integer fields (URWORD, READ *)
_REAL_RE = re.compile(r"[+-]?(\d+\.?\d*|\.\d+)([eEdD][+-]?\d+)?")
_INT_RE = re.compile(r"[+-]?\d+")


def _words(line):
    """Split a record as URWORD does (blank, comma and tab separate words)."""
    return [w for w in re.split(r"[ ,\t]+", line.strip()) if w]


def _parse_real(token, name):
    if not _REAL_RE.fullmatch(token):
        raise ValueError(f"DDF: {name} must be a real number, got {token!r}")
    return float(token.replace("d", "e").replace("D", "e"))


def _parse_int(token, name):
    if not _INT_RE.fullmatch(token):
        raise ValueError(f"DDF: {name} must be an integer, got {token!r}")
    return int(token)


def _as_real(value, name):
    if isinstance(value, bool) or not isinstance(value, numbers.Real):
        raise ValueError(f"DDF: {name} must be a real number, got {value!r}")
    value = float(value)
    if not math.isfinite(value):
        raise ValueError(f"DDF: {name} must be finite, got {value!r}")
    return value


def _as_real32(value, name):
    """The value as USG-T stores it in a REAL (single precision) variable."""
    with np.errstate(over="ignore", under="ignore"):
        value32 = float(np.float32(value))
    if not math.isfinite(value32):
        raise ValueError(
            f"DDF: {name} {value!r} is outside the single-precision (REAL) range"
        )
    return value32


def _as_int(value, name, allow_bool=True):
    if (isinstance(value, bool) and not allow_bool) or not isinstance(
        value, numbers.Integral
    ):
        raise ValueError(f"DDF: {name} must be an integer, got {value!r}")
    return int(value)


def _next_record(f, what):
    """Return the next record that does not start with '#' (URDCOM)."""
    while True:
        line = f.readline()
        if line == "":
            raise ValueError(f"DDF: end of file before the {what}")
        if not line.startswith("#"):
            break
    if not line.strip():
        raise ValueError(f"DDF: blank {what} is not supported (USG-T would read zeros)")
    return line


def _read_table(f, nrows):
    """Read NDDFTABROWS list-directed records of two values (CONC DENSITY).

    As READ(INDDF,*): a row may continue on the next lines and blank lines
    are skipped; what follows the second value on the last line is ignored.
    Null values, '/', repeat counts and quotes are not supported.
    """
    rows = []
    for irow in range(1, nrows + 1):
        what = f"NONLINEAR density table row {irow} of {nrows}"
        vals = []
        prev_comma = True
        while len(vals) < 2:
            line = f.readline()
            if line == "":
                raise ValueError(f"DDF: end of file in the {what}")
            for tok in re.findall(r"[^\s,]+|,", line):
                if len(vals) == 2:
                    break
                if tok == ",":
                    if prev_comma:
                        raise ValueError(f"DDF: null value in the {what}")
                    prev_comma = True
                    continue
                prev_comma = False
                if any(c in tok for c in "/*'\""):
                    raise ValueError(
                        f"DDF: unsupported list-directed syntax {tok!r} in the {what}"
                    )
                vals.append(_parse_real(tok, what))
        rows.append(tuple(vals))
    return rows


class _NamOutput:
    """ISWIUN entry of a name file, validated by one MfUsgDdf.load() call for
    one model; only that construction may adopt it as the DDF output."""

    __slots__ = ("fname", "iswiun", "model")

    def __init__(self, model, iswiun, fname):
        self.model, self.iswiun, self.fname = model, iswiun, fname

    def matches(self, model, iswiun, fname):
        return self.model is model and (self.iswiun, self.fname) == (iswiun, fname)


class MfUsgDdf(Package):
    """MODFLOW-USG Transport Density Driven Flow (DDF) package class.

    Parameters
    ----------
    model : model object
        the model object (of type :class:`flopy.mfusg.MfUsg`) to which
        this package will be added.
    rhofresh : float
        density of freshwater (not zero).
    rhostd : float
        density of the standard solution (with isharp != 0 it must differ
        from rhofresh).
    cstd : float
        concentration of the standard solution, in the units of BCT CONC
        (not zero when isharp == 0, also with nonlinear).
    ithickav : int
        a flag indicating if thickness weighted averaging for density term
        0 = arithmetic averaging; 1 = thickness weighting averaging
    imphdd : int
        a flag of hydraulic head term in the density formulation
        0 = explicit (on the right-hand side vector) symmetry of the matrix
        1 = implicit (on the left-hand side matrix) asymmetric matrix
    isharp : int
        sharp-interface flag; any value other than 0 selects the sharp
        interface model (default 0).
    nonlinear : bool
        if True, density is a tabular function of concentration
        (density_table) and NONLINEAR NDDFTABROWS follows ISHARP.
    density_table : list of (float, float)
        rows (CONC, DENSITY) with strictly increasing concentrations;
        required (at least one row) when nonlinear is True.
    extension : str, optional
        file extension (default is 'ddf').
    unitnumber : int, optional
        fortran unit number for this package (default is none).
    filenames : str or list of str
        filenames to use for the package and the interface output file:
        [DDF file, interface elevation file]. if filenames=none the names
        are created from the model name and the extensions 'ddf' and
        'zta'. if a single string is passed it is the DDF file name.
        default is none.
    hsea : float, optional
        sea level (HSEA), an elevation in model length units. required when
        isharp != 0 and only allowed then.
    iswiun : int, optional
        unit number on which the interface elevation is saved (ISWIUN) when
        isharp != 0; 0 (default) means no output. a positive value registers
        a binary output file (extension 'zta') saved with the OC head flags,
        readable with HeadUFile(text="interface elev."). it cannot be
        changed after the package is created. the unit and the registered
        file name (the base name, as add_output_file registers it) must not
        belong to another package, output or external file; when loading
        from a name file its DATA(BINARY) entry is used. a model with an
        output still registered for DDF on another unit (a replaced or
        removed DDF package) is refused until that output is removed with
        model.remove_output.

    Notes
    -----
    USG-T reads the first record (RHOFRESH RHOSTD CSTD ITHICKAV IMPHDD
    ISHARP [NONLINEAR NDDFTABROWS]) in free format; missing integers are 0.
    The three densities and concentration are read into a REAL (single
    precision) temporary, so the header is effectively single precision
    in the engine even though it is written with full precision; the
    NONLINEAR table is read in double precision. With isharp != 0 the
    record HSEA ISWIUN follows the table.

    The zero-denominator checks (RHOFRESH, CSTD with isharp == 0, RHOSTD -
    RHOFRESH with isharp != 0) and the range of RHOFRESH, RHOSTD, CSTD and
    HSEA are checked on the single-precision values; the values are written
    as given.

    Not supported (NotImplementedError): BAS without FREE (USG-T reads the
    first record with 5F10.4 and cannot read the integer items), when
    writing and, if the BAS package is already in the model, when loading
    (without a BAS package the format is not inferred); and isharp != 0
    with a CLN package (USG-T then also reads ZETASWI for the CLN nodes).
    The reader accepts a subset of list-directed input for the table.
    Non-sharp DDF uses BCT concentrations; check() reports a
    missing BCT package. With a structured grid, non-sharp DDF and saved
    heads stopped USG-T 2.8.1 and 2.7.0 with a segmentation fault; the
    sharp interface was not run on a structured grid.

    Examples
    --------

    >>> import flopy
    >>> m = flopy.mfusg.MfUsg()
    >>> ddf = flopy.mfusg.MfUsgDdf(m)
    """

    def __init__(
        self,
        model,
        rhofresh=1000.0,
        rhostd=1025.0,
        cstd=35.0,
        ithickav=1,
        imphdd=0,
        isharp=0,
        nonlinear=False,
        density_table=None,
        extension="ddf",
        unitnumber=None,
        filenames=None,
        hsea=None,
        iswiun=0,
    ):
        self._setup(
            model,
            rhofresh,
            rhostd,
            cstd,
            ithickav,
            imphdd,
            isharp,
            nonlinear,
            density_table,
            extension,
            unitnumber,
            filenames,
            hsea,
            iswiun,
            nam_output=None,
        )

    def _setup(
        self,
        model,
        rhofresh,
        rhostd,
        cstd,
        ithickav,
        imphdd,
        isharp,
        nonlinear,
        density_table,
        extension,
        unitnumber,
        filenames,
        hsea,
        iswiun,
        nam_output,
    ):
        """Body of __init__. nam_output is the _NamOutput built by load() for
        this model, or None (authoring)."""
        msg = (
            "Model object must be of type flopy.mfusg.MfUsg\n"
            f"but received type: {type(model)}."
        )
        assert isinstance(model, MfUsg), msg

        # set default unit number of one is not specified
        if unitnumber is None:
            unitnumber = self._defaultunit()

        values = self._checked_values(
            rhofresh,
            rhostd,
            cstd,
            ithickav,
            imphdd,
            isharp,
            nonlinear,
            density_table,
            hsea,
            iswiun,
        )
        if isinstance(filenames, list):
            filenames = list(filenames)
        filenames = self._prepare_filenames(filenames, 2)
        fname = filenames[1] and os.path.basename(filenames[1])
        out_name = self._output_name(
            model,
            values["iswiun"],
            unitnumber,
            filenames[1],
            nam_output=isinstance(nam_output, _NamOutput)
            and nam_output.matches(model, values["iswiun"], fname),
            own_file=self._own_file_name(model, filenames[0], extension),
        )

        # call base package constructor
        super().__init__(
            model,
            extension=extension,
            name=self._ftype(),
            unit_number=unitnumber,
            filenames=filenames[0],
        )

        self._generate_heading()
        for key, value in values.items():
            setattr(self, key, value)
        self._iswiun_registered = self.iswiun
        self._iswiun_fname = out_name
        if self.iswiun > 0:
            model.add_output_file(
                self.iswiun,
                fname=out_name,
                extension="zta",
                binflag=True,
                package=self._ftype(),
            )

        self.parent.add_package(self)
        return

    @staticmethod
    def _checked_values(
        rhofresh,
        rhostd,
        cstd,
        ithickav,
        imphdd,
        isharp,
        nonlinear,
        density_table,
        hsea,
        iswiun,
    ):
        """Validate the DDF input; return the normalized values."""
        v = {
            "rhofresh": _as_real(rhofresh, "RHOFRESH"),
            "rhostd": _as_real(rhostd, "RHOSTD"),
            "cstd": _as_real(cstd, "CSTD"),
            "ithickav": _as_int(ithickav, "ITHICKAV"),
            "imphdd": _as_int(imphdd, "IMPHDD"),
            "isharp": _as_int(isharp, "ISHARP"),
            "nonlinear": bool(nonlinear),
        }
        # USG-T reads RHOFRESH RHOSTD CSTD (and HSEA) into a REAL variable:
        # the denominators are checked on the single-precision values
        r32 = {k: _as_real32(v[k], k.upper()) for k in ("rhofresh", "rhostd", "cstd")}
        if r32["rhofresh"] == 0.0:
            raise ValueError("DDF: RHOFRESH must not be zero (in single precision)")
        if v["isharp"] != 0 and r32["rhostd"] == r32["rhofresh"]:
            raise ValueError(
                "DDF: with ISHARP != 0, RHOSTD must differ from RHOFRESH in single "
                "precision (USG-T reads them as REAL)"
            )
        if v["isharp"] == 0 and r32["cstd"] == 0.0:
            raise ValueError(
                "DDF: with ISHARP = 0, CSTD must not be zero (in single precision)"
            )

        table = []
        if density_table is not None:
            try:
                rows = list(density_table)
            except TypeError:
                raise ValueError(
                    f"DDF: density_table must be a list of rows, got {density_table!r}"
                ) from None
            for irow, row in enumerate(rows, start=1):
                try:
                    row = tuple(row)
                except TypeError:
                    row = (row,)
                if len(row) != 2:
                    raise ValueError(
                        f"DDF: density_table row {irow} must have two values "
                        f"(CONC, DENSITY), got {row!r}"
                    )
                table.append(
                    (
                        _as_real(row[0], f"density_table row {irow} CONC"),
                        _as_real(row[1], f"density_table row {irow} DENSITY"),
                    )
                )
        if v["nonlinear"] and not table:
            raise ValueError("DDF: NONLINEAR requires a non-empty density_table")
        if not v["nonlinear"] and table:
            raise ValueError("DDF: density_table is only written with NONLINEAR")
        for irow in range(1, len(table)):
            if not table[irow][0] > table[irow - 1][0]:
                raise ValueError(
                    "DDF: density_table concentrations must be strictly "
                    f"increasing (rows {irow} and {irow + 1})"
                )
        v["density_table"] = table

        v["iswiun"] = _as_int(iswiun, "ISWIUN", allow_bool=False)
        if v["iswiun"] < 0:
            raise ValueError(f"DDF: ISWIUN must be >= 0, got {iswiun!r}")
        if v["isharp"] != 0:
            if hsea is None:
                raise ValueError("DDF: ISHARP != 0 requires HSEA")
            v["hsea"] = _as_real(hsea, "HSEA")
            _as_real32(v["hsea"], "HSEA")
        else:
            if hsea is not None or v["iswiun"] > 0:
                raise ValueError(
                    "DDF: HSEA and ISWIUN > 0 are only read when ISHARP != 0"
                )
            v["hsea"] = None
        return v

    @classmethod
    def _output_name(
        cls, model, iswiun, unitnumber, fname, nam_output=False, own_file=None
    ):
        """Return the file name that add_output_file registers for ISWIUN.

        Raise ValueError, without changing the model, if the unit or that
        name belongs to another file (including own_file, the base name of
        the DDF file being created), if another unit is still registered
        for DDF (a replaced or removed DDF package), or if the unit is an
        external file. nam_output=True (load from a name file, validated for
        this model) accepts the DATA(BINARY) entry of ISWIUN with the name
        given in fname.
        """
        own = cls._ftype()
        for unit, owners in zip(model.output_units, model.output_packages):
            if own in owners and unit != iswiun:
                raise ValueError(
                    f"DDF: output unit {unit} is still registered for DDF; remove "
                    f"it (model.remove_output(unit={unit})) before using ISWIUN "
                    f"{iswiun}"
                )
        if iswiun <= 0:
            return None
        if iswiun == unitnumber or iswiun in model.package_units:
            raise ValueError(f"DDF: ISWIUN {iswiun} is a package unit number")
        if fname is not None:
            fname = os.path.basename(fname)
        adopted = None
        if iswiun in model.external_units:
            idx = model.external_units.index(iswiun)
            ext = os.path.basename(model.external_fnames[idx])
            if not (nam_output and ext == fname and model.external_binflag[idx]):
                raise ValueError(
                    f"DDF: ISWIUN {iswiun} is the external file {ext!r} of the model"
                )
            adopted = idx
        if iswiun in model.output_units:
            idx = model.output_units.index(iswiun)
            owners = model.output_packages[idx]
            name = model.output_fnames[idx]
            if not owners or any(p != own for p in owners):
                raise ValueError(
                    f"DDF: ISWIUN {iswiun} is already an output unit of "
                    f"{owners or 'another file'}"
                )
            if not model.output_binflag[idx]:
                raise ValueError(f"DDF: output unit {iswiun} is not binary")
            if fname is not None and fname != name:
                raise ValueError(
                    f"DDF: output unit {iswiun} is registered as {name!r}, "
                    f"not {fname!r}"
                )
            cls._check_not_own_file(name, own_file)
            return name
        if fname is not None:
            name = fname
        else:
            # default name of add_output_file (package given)
            name = f"{model.name}.zta"
            if name in model.output_fnames:
                name = f"{model.name}.{own}.zta"
        if name in model.output_fnames:
            idx = model.output_fnames.index(name)
            raise ValueError(
                f"DDF: output file {name!r} is already registered on unit "
                f"{model.output_units[idx]}"
            )
        others = [
            os.path.basename(f)
            for i, f in enumerate(model.external_fnames)
            if i != adopted
        ]
        for p in model.packagelist:
            others += [os.path.basename(f) for f in p.file_name]
        if name in others:
            raise ValueError(f"DDF: {name!r} is already a file of the model")
        cls._check_not_own_file(name, own_file)
        return name

    @staticmethod
    def _check_not_own_file(name, own_file):
        if own_file is not None and name == own_file:
            raise ValueError(
                f"DDF: the interface output {name!r} is the DDF file itself"
            )

    @staticmethod
    def _own_file_name(model, fname, extension):
        """Base name Package.__init__ gives the DDF file."""
        if fname is None:
            return f"{model.name}.{extension}"
        return os.path.basename(fname)

    def _validate(self):
        """Check the current attributes before anything is written."""
        v = self._checked_values(
            self.rhofresh,
            self.rhostd,
            self.cstd,
            self.ithickav,
            self.imphdd,
            self.isharp,
            self.nonlinear,
            self.density_table,
            self.hsea,
            self.iswiun,
        )
        if v["iswiun"] != self._iswiun_registered:
            raise ValueError(
                f"DDF: ISWIUN changed from {self._iswiun_registered} to "
                f"{v['iswiun']} after the package was created; create a new "
                "MfUsgDdf to change the interface output"
            )
        if v["iswiun"] > 0:
            m = self.parent
            if m.output_units.count(v["iswiun"]) != 1:
                raise ValueError(
                    f"DDF: ISWIUN {v['iswiun']} is no longer a model output unit"
                )
            idx = m.output_units.index(v["iswiun"])
            owners = m.output_packages[idx]
            if (
                not owners
                or any(p != self._ftype() for p in owners)
                or not m.output_binflag[idx]
                or m.output_fnames[idx] != self._iswiun_fname
            ):
                raise ValueError(
                    f"DDF: output unit {v['iswiun']} changed after the package was "
                    f"created (now {m.output_fnames[idx]!r}, owners {owners}, "
                    f"binary {m.output_binflag[idx]}; expected "
                    f"{self._iswiun_fname!r}, DDF, binary)"
                )
        self._check_free(self.parent)
        if v["isharp"] != 0 and self._cln_written(self.parent):
            raise NotImplementedError(
                "DDF with ISHARP != 0 and CLN (ZETASWI) is not supported"
            )
        return v

    @staticmethod
    def _check_free(model):
        """NotImplementedError if the model has a BAS package without FREE."""
        bas = model.get_package("BAS6")
        if bas is not None and not bas.ifrefm:
            raise NotImplementedError(
                "DDF requires BAS FREE: USG-T reads the first DDF record of a "
                "fixed-format model with 5F10.4 and cannot read its integers"
            )

    @staticmethod
    def _cln_written(model):
        """True if the name file written by the model declares CLN."""
        if model.get_package("CLN") is not None:
            return True
        skipped = getattr(model, "_skipped_nam_entries", {})
        return any(ftype.upper() == "CLN" for ftype, _ in skipped.values())

    def write_file(self):
        """
        Write the package file.

        Returns
        -------
        None

        """
        v = self._validate()
        lines = [
            f"{self.heading}\n",
            f" {v['rhofresh']!r} {v['rhostd']!r} {v['cstd']!r}"
            f" {v['ithickav']} {v['imphdd']} {v['isharp']}",
        ]
        if v["nonlinear"]:
            lines[1] += f" NONLINEAR {len(v['density_table'])}"
        lines[1] += "\n"
        for conc_i, rho_i in v["density_table"]:
            lines.append(f" {conc_i!r} {rho_i!r}\n")
        if v["isharp"] != 0:
            lines.append(f" {v['hsea']!r} {v['iswiun']}\n")
        with open(self.fn_path, "w") as f:
            f.writelines(lines)

    def check(self, f=None, verbose=True, level=1, checktype=None):
        """
        Check package data for common errors.

        Parameters
        ----------
        f : str or file handle
            String defining file name or file handle for summary file
            of check method output. If a string is passed a file handle
            is created. If f is None, check method does not write
            results to a summary file. (default is None)
        verbose : bool
            Boolean flag used to determine if check method results are
            written to the screen.
        level : int
            Check method analysis level. If level=0, summary checks are
            performed. If level=1, full checks are performed.

        Returns
        -------
        check object
        """
        chk = self._get_check(f, verbose, level, checktype)
        if not self.isharp and self.parent.get_package("BCT") is None:
            chk._add_to_summary(
                "Error",
                package="DDF",
                desc="DDF with ISHARP = 0 needs a BCT package (CONC)",
            )
        chk.summarize()
        return chk

    @classmethod
    def load(cls, f, model, ext_unit_dict=None):
        """
        Load an existing package.

        Parameters
        ----------
        f : filename or file handle
            File to load. A file handle given by the caller is left open,
            positioned after the last DDF record.
        model : model object
            The model object (of type :class:`flopy.mfusg.MfUsg`) to
            which this package will be added.
        ext_unit_dict : dictionary, optional
            If the arrays in the file are specified using EXTERNAL,
            or older style array control records, then `f` should be a file
            handle.  In this case ext_unit_dict is required, which can be
            constructed using the function
            :class:`flopy.utils.mfreadnam.parsenamefile`.

        Returns
        -------
        ddf : MfUsgDdf object

        Examples
        --------

        >>> import flopy
        >>> ml = flopy.mfusg.MfUsg()
        >>> ddf = flopy.mfusg.MfUsgDdf.load('Henry.ddf', ml)
        """
        msg = (
            "Model object must be of type flopy.mfusg.MfUsg\n"
            f"but received type: {type(model)}."
        )
        assert isinstance(model, MfUsg), msg
        cls._check_free(model)

        if model.verbose:
            print("loading ddf package file...")

        if model.version != "mfusg":
            print(
                "Warning: model version was reset from '{}' to 'mfusg' "
                "in order to load a DDF file".format(model.version)
            )
            model.version = "mfusg"

        openfile = not hasattr(f, "read")
        if openfile:
            filename = f
            f = open(filename, "r")
        try:
            kw = cls._read(f, model, ext_unit_dict)
        finally:
            if openfile:
                f.close()

        if model.verbose:
            print(
                f"   RHOFRESH {kw['rhofresh']} \n   RHOSTD {kw['rhostd']} \n"
                f"   CSTD {kw['cstd']} \n   ITHICKAV {kw['ithickav']} \n"
                f"   IMPHDD {kw['imphdd']} \n   ISHARP {kw['isharp']}"
                + (
                    f"\n   NONLINEAR with {len(kw['density_table'])} rows"
                    if kw["nonlinear"]
                    else ""
                )
            )

        # check everything before the name-file units are claimed
        iswiun = cls._checked_values(**kw)["iswiun"]
        unitnumber, nam_fname, ddf_fname = None, None, None
        if ext_unit_dict is not None:
            unitnumber, ddf_fname = model.get_ext_dict_attr(
                ext_unit_dict, filetype=cls._ftype(), pop_key=False
            )
            item = ext_unit_dict.get(iswiun) if iswiun > 0 else None
            if item is not None:
                if item.filetype.upper() != "DATA(BINARY)":
                    raise ValueError(
                        f"DDF: ISWIUN {iswiun} is a {item.filetype} entry of the "
                        "name file, not DATA(BINARY)"
                    )
                nam_fname = os.path.basename(item.filename)
        cls._output_name(
            model,
            iswiun,
            unitnumber or cls._defaultunit(),
            nam_fname,
            nam_output=True,
            own_file=cls._own_file_name(model, ddf_fname, "ddf"),
        )

        # set package unit number
        unitnumber, filenames = get_unitnumber_from_ext_unit_dict(
            model,
            cls,
            ext_unit_dict,
            ipakcb=iswiun,
        )

        # build through the internal path with the context of this load only
        ddf = cls.__new__(cls)
        ddf._setup(
            model,
            kw["rhofresh"],
            kw["rhostd"],
            kw["cstd"],
            kw["ithickav"],
            kw["imphdd"],
            kw["isharp"],
            kw["nonlinear"],
            kw["density_table"],
            "ddf",
            unitnumber,
            filenames,
            kw["hsea"],
            kw["iswiun"],
            nam_output=_NamOutput(model, iswiun, nam_fname) if nam_fname else None,
        )
        return ddf

    @classmethod
    def _read(cls, f, model, ext_unit_dict):
        """Read the DDF records in the DDF1AR / SHARPINT1AR order."""
        if model.verbose:
            print("   loading RHOFRESH RHOSTD CSTD ITHICKAV IMPHDD ISHARP...")
        words = _words(_next_record(f, "DDF item 1 record"))
        if len(words) < 3:
            raise ValueError(
                "DDF: item 1 needs RHOFRESH RHOSTD CSTD, got " + " ".join(words)
            )
        kw = {
            "rhofresh": _parse_real(words[0], "RHOFRESH"),
            "rhostd": _parse_real(words[1], "RHOSTD"),
            "cstd": _parse_real(words[2], "CSTD"),
        }
        for idx, name in enumerate(("ithickav", "imphdd", "isharp"), start=3):
            if idx >= len(words):
                kw[name] = 0
            elif words[idx].upper() == "NONLINEAR":
                raise ValueError(
                    f"DDF: NONLINEAR found in place of {name.upper()}; "
                    "it must follow ISHARP (the 6th item)"
                )
            else:
                kw[name] = _parse_int(words[idx], name.upper())
        rest = words[6:]
        kw["nonlinear"] = bool(rest) and rest[0].upper() == "NONLINEAR"
        nrows = 0
        if kw["nonlinear"]:
            if len(rest) < 2:
                raise ValueError("DDF: NONLINEAR must be followed by NDDFTABROWS")
            nrows = _parse_int(rest[1], "NDDFTABROWS")
            if nrows < 1:
                raise ValueError(f"DDF: NDDFTABROWS must be >= 1, got {nrows}")
            rest = rest[2:]
        if any(w.upper() == "NONLINEAR" for w in rest):
            raise ValueError(
                "DDF: NONLINEAR must be the item after ISHARP; USG-T ignores "
                "it elsewhere"
            )
        if kw["isharp"] != 0 and cls._cln_declared(model, ext_unit_dict):
            raise NotImplementedError(
                "DDF with ISHARP != 0 and CLN (ZETASWI) is not supported"
            )
        kw["density_table"] = _read_table(f, nrows) if nrows else None
        kw["hsea"], kw["iswiun"] = None, 0
        if kw["isharp"] != 0:
            words = _words(_next_record(f, "HSEA ISWIUN record (ISHARP != 0)"))
            kw["hsea"] = _parse_real(words[0], "HSEA")
            if len(words) > 1:
                kw["iswiun"] = _parse_int(words[1], "ISWIUN")
        return kw

    @staticmethod
    def _cln_declared(model, ext_unit_dict):
        """True if the model or its name file has CLN (loaded or not)."""
        if model.get_package("CLN") is not None or getattr(model, "icln", 0):
            return True
        return ext_unit_dict is not None and any(
            getattr(item, "filetype", None) == "CLN" for item in ext_unit_dict.values()
        )

    @staticmethod
    def _ftype():
        return "DDF"

    @staticmethod
    def _defaultunit():
        return 152
