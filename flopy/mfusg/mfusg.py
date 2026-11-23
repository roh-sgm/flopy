"""Mfusg module."""

import filecmp
import os.path
import shutil
import sys
from inspect import getfullargspec
from os import PathLike, curdir
from pathlib import PurePosixPath, PureWindowsPath
from typing import Union
from warnings import warn

import numpy as np

import flopy

from ..mbase import PackageLoadException
from ..modflow import Modflow
from ..utils import mfreadnam


class MfUsg(Modflow):
    """MODFLOW-USG Model Class.

    Parameters
    ----------
    modelname : str or PathLike, default "modflowusgtest"
        Name of model.  This string will be used to name the MODFLOW input
        that are created with write_model.
    namefile_ext : str, default "nam"
        Extension for the namefile.
    exe_name : str, default "mfusg"
        The name of the executable to use.
    structured : bool, default True
        Specify if model grid is structured (default) or unstructured.
    listunit : int, default 2
        Unit number for the list file.
    model_ws : str, default "." (curdir)
        Model workspace.  Directory name to create model data sets.
    external_path : str, optional
        Location for external files.
    verbose : bool, default False
        Print additional information to the screen.
    free_format_npl : int, optional
        Number of values per line when writing free-format arrays. When set
        (e.g., ``free_format_npl=10``), arrays are written with that many
        values per line instead of all values on a single line. This produces
        block-format output, improving readability for large models.
        Default is None (all values on one line).

    Attributes
    ----------

    Methods
    -------

    See Also
    --------

    Notes
    -----
    EVS and RTS name-file entries (the ET and recharge time-series files read
    by EVT with ETS and RCH with RTS) are kept from a loaded name file as
    references to opaque files: their content is not read or interpreted.
    ``write_input`` writes them back with their unit, path and the options
    FLUSH and NO-ARCHIVE. In the original workspace the files are not
    touched; in a new workspace a relative path inside the original workspace
    is copied to the same relative path (never replacing a different file),
    and an absolute path or a path outside the original workspace is written
    as an absolute reference, so the model still depends on that file.

    Examples
    --------
    >>> import flopy
    >>> usg = flopy.mfusg.MfUsg()
    """

    def __init__(
        self,
        modelname="mfusgtest",
        namefile_ext="nam",
        version="mfusg",
        exe_name: Union[str, PathLike] = "mfusg",
        structured=True,
        listunit=2,
        model_ws: Union[str, PathLike] = curdir,
        external_path=None,
        verbose=False,
        free_format_npl=None,
        **kwargs,
    ):
        super().__init__(
            modelname=modelname,
            namefile_ext=namefile_ext,
            version=version,
            exe_name=exe_name,
            structured=structured,
            listunit=listunit,
            model_ws=model_ws,
            external_path=external_path,
            verbose=verbose,
            free_format_npl=free_format_npl,
            **kwargs,
        )

        self.itrnsp = 0  # transport simulation flag
        self.mcomp = 0  # number of chemical components
        self.iheat = 0  # flag for heat transport
        self.imcomp = 0  # number of immobile components
        # MfUsgBct that last set the counts (None: no BCT was ever added)
        self._counts_owner = None
        self.idpf = 0  # flag for dual porosity flow
        self.icln = 0  # flag for CLN package
        # EVS/RTS name-file entries kept from a loaded name file
        self._preserved_nam_entries = []

        # Create a dictionary to map package with package object.
        # This is used for loading models.
        self.mfnam_packages = {
            "zone": flopy.modflow.ModflowZon,
            "mult": flopy.modflow.ModflowMlt,
            "pval": flopy.modflow.ModflowPval,
            "bas6": flopy.mfusg.MfUsgBas,
            "dis": flopy.mfusg.MfUsgDis,
            "hfb6": flopy.mfusg.MfUsgHfb,
            "fhb": flopy.modflow.ModflowFhb,
            "drn": flopy.mfusg.MfUsgDrn,
            "drt": flopy.mfusg.MfUsgDrt,
            "ghb": flopy.mfusg.MfUsgGhb,
            "gsf": flopy.mfusg.MfUsgGsf,
            "riv": flopy.mfusg.MfUsgRiv,
            "sgb": flopy.mfusg.MfUsgSgb,
            "qrt": flopy.mfusg.MfUsgQrt,
            # SUB/SWT: guarded so unstructured-MfUsg use fails explicitly
            # (Stage 4.6A; the USG-T DISU layouts are unvalidated for these).
            # STR/SFR/FHB/GAGE stay on the base classes: ModflowStr handles
            # the unstructured reach list, unstructured SFR is base-validated
            # by freyberg_usg, and FHB/GAGE round-trip via base in the Ex8
            # real model.
            "str": flopy.modflow.ModflowStr,
            "sfr": flopy.modflow.ModflowSfr2,
            "gage": flopy.modflow.ModflowGage,
            "sub": flopy.mfusg.MfUsgSub,
            "swt": flopy.mfusg.MfUsgSwt,
            "chd": flopy.mfusg.MfUsgChd,
            "disu": flopy.mfusg.MfUsgDisU,
            "sms": flopy.mfusg.MfUsgSms,
            "wel": flopy.mfusg.MfUsgWel,
            "bcf6": flopy.mfusg.MfUsgBcf,
            "lpf": flopy.mfusg.MfUsgLpf,
            "cln": flopy.mfusg.MfUsgCln,
            "gnc": flopy.mfusg.MfUsgGnc,
            "bct": flopy.mfusg.MfUsgBct,
            "pcb": flopy.mfusg.MfUsgPcb,
            "ddf": flopy.mfusg.MfUsgDdf,
            "mdt": flopy.mfusg.MfUsgMdt,
            "dpf": flopy.mfusg.MfUsgDpf,
            "dpt": flopy.mfusg.MfUsgDpt,
            "rch": flopy.mfusg.MfUsgRch,
            "oc": flopy.mfusg.MfUsgOc,
            "lak": flopy.mfusg.MfUsgLak,
            "evt": flopy.mfusg.MfUsgEvt,
            "ets": flopy.mfusg.MfUsgEts,
            "tib": flopy.mfusg.MfUsgTib,
            "tvm": flopy.mfusg.MfUsgTvm,
        }

    def __repr__(self):
        """Returns a representation of the MfUsg object."""
        nrow, ncol, nlay, nper = self.get_nrow_ncol_nlay_nper()
        if nrow is not None:
            # structured case
            msg = (
                f"MODFLOW {nlay} layer(s) {nrow} row(s) {ncol} column(s) "
                f"{nper} stress period(s)"
            )
        else:
            # unstructured case
            msg = (
                "MODFLOW unstructured\n"
                f"  nodes = {ncol.sum()}\n"
                f"  layers = {nlay}\n"
                f"  stress periods = {nper}\n"
                f"  nodelay = {ncol}\n"
            )
        return msg

    # Name-file types that USG-T opens (SGWF2BAS8OPEN, CUNIT 6 and 9) and that
    # EVT (ETS option) and RCH (RTS option) read as time series, but that
    # FloPy does not interpret: kept as references to opaque files.
    _PRESERVED_NAM_TYPES = ("EVS", "RTS")
    # Trailing options that SGWF2BAS8OPEN applies to package input files
    _PRESERVED_NAM_OPTIONS = ("FLUSH", "NO-ARCHIVE")
    # SGWF2BAS8OPEN reads each line into CHARACTER*300 and URWORD blanks the
    # last character, so a line can use 299 characters
    _NAM_LINE_MAX = 299

    @property
    def mcompt(self):
        """Number of mobile components including temperature (read-only).

        MCOMPT = MCOMP + 1 when IHEAT is 1 (USG-T GWT2BCT1AR), else MCOMP.
        The counts ``mcomp``, ``iheat`` and ``imcomp`` are 0 when the model
        is created and are set by MfUsgBct when it is added to the model.
        Removing the BCT package does not reset them.
        """
        return self.mcomp + (1 if self.iheat == 1 else 0)

    @property
    def ntcomp(self):
        """Total number of components (read-only): MCOMPT + IMCOMP.

        Component order in USG-T: mobile 1..MCOMP, temperature (IHEAT=1),
        then immobile components.
        """
        return self.mcompt + self.imcomp

    def _transport_context(self):
        """Classify the transport state used by the BCT-dependent packages.

        Returns ``(state, counts)``. ``state`` is:

        - ``"bct"``: a BCT package is registered and set the model counts;
        - ``"unsynced"``: the registered BCT did not set the model counts (it
          was added without its constructor) or its counts changed after
          construction;
        - ``"removed"``: a BCT set the counts but is no longer registered;
        - ``"manual"``: no BCT was ever added and ``itrnsp`` is not 0 (counts
          assigned by hand, for writing a package file on its own);
        - ``"none"``: no BCT and no transport (flow only).

        ``counts`` is ``(mcomp, iheat, imcomp)`` for ``"bct"`` and
        ``"manual"``, else None. USG-T reads the transport items of EVT/ETS
        only when BCT is in the name file, so ``"bct"`` and ``"unsynced"``
        (a registered BCT) and ``"manual"`` have the transport header.
        """
        bct = self.get_package("BCT")
        if bct is None:
            if self._counts_owner is not None:
                return "removed", None
            if self.itrnsp:
                return "manual", (self.mcomp, self.iheat, self.imcomp)
            return "none", None
        if bct is not self._counts_owner or bct._counts_mismatch() is not None:
            return "unsynced", None
        return "bct", (self.mcomp, self.iheat, self.imcomp)

    def _transport_header(self):
        """True when USG-T reads the transport items of EVT/ETS (BCT in the
        name file), or in the manual-count mode."""
        return self._transport_context()[0] in ("bct", "unsynced", "manual")

    def _transport_mcompt(self, who):
        """MCOMPT for a package that needs per-component transport data.

        Raises ValueError when the counts cannot be used: no transport, BCT
        removed, BCT not synchronized with the model counts, or invalid
        manual counts."""
        state, counts = self._transport_context()
        if state == "none":
            raise ValueError(f"{who}: transport factors need a BCT package")
        if state == "removed":
            raise ValueError(
                f"{who}: the BCT package was removed; its component counts are "
                "stale. Add a BCT package again or disable the transport factor."
            )
        if state == "unsynced":
            bct = self.get_package("BCT")
            reason = bct._counts_mismatch() or (
                "it was added without setting the model counts"
            )
            raise ValueError(
                f"{who}: the registered BCT package does not match the model "
                f"component counts ({reason}); build or load the BCT package "
                "with its constructor"
            )
        mcomp, iheat, imcomp = counts
        for name, value in (("mcomp", mcomp), ("iheat", iheat), ("imcomp", imcomp)):
            if isinstance(value, bool) or not isinstance(value, (int, np.integer)):
                raise ValueError(f"{who}: model {name}={value!r} is not an integer")
            if value < 0:
                raise ValueError(f"{who}: model {name}={value} is negative")
        mcompt = self.mcompt
        if mcompt < 1:
            raise ValueError(f"{who}: transport factors need MCOMPT >= 1")
        return mcompt

    @property
    def preserved_nam_entries(self):
        """EVS/RTS name-file entries kept from the loaded name file.

        Returns independent copies; changing them does not change the model.
        """
        return [
            {**e, "options": list(e["options"])} for e in self._preserved_nam_entries
        ]

    @classmethod
    def _read_preserved_entries(cls, namefile_path):
        """Read the EVS/RTS entries of a name file with the URWORD rules.

        Returns, per entry, the file type, unit, path (unquoted), supported
        options, the absolute path of the file and the name-file folder.
        Malformed entries raise ValueError; a unit <= 0 or a missing file
        warns (the entry is kept for inspection and rejected on write).
        """
        entries = []
        nam_dir = os.path.dirname(os.path.abspath(namefile_path))
        with open(namefile_path) as f:
            lines = f.read().splitlines()
        for ln, line in enumerate(lines, 1):
            if line.startswith("#"):
                continue
            words = _urwords(line)
            if not words or words[0][0].upper() not in cls._PRESERVED_NAM_TYPES:
                continue
            ftype = words[0][0].upper()
            msg = f"MfUsg.load: name file line {ln} ({ftype}) "
            if len(line) > cls._NAM_LINE_MAX:
                raise ValueError(msg + f"is longer than {cls._NAM_LINE_MAX} characters")
            if len(words) < 3:
                raise ValueError(msg + "needs a unit and a file name")
            unit = int(words[1][0])
            path, closed = words[2]
            if not closed:
                raise ValueError(msg + "has a file name with an unclosed quote")
            if not path:
                raise ValueError(msg + "has an empty file name")
            options, dropped = [], []
            for word, _ in words[3:]:
                if word.upper() in cls._PRESERVED_NAM_OPTIONS:
                    options.append(word.upper())
                else:
                    dropped.append(word)
            if dropped:
                warn(
                    msg + f"options {dropped} are not used by USG-T for an input "
                    "file and are not kept"
                )
            if unit <= 0:
                warn(
                    msg + f"has unit {unit}; it cannot be written (USG-T opens a "
                    "negative unit as BINARY and reads the time series as text)"
                )
            pure = PureWindowsPath(path) if "\\" in path else PurePosixPath(path)
            source = os.path.normpath(os.path.join(nam_dir, *pure.parts))
            if not os.path.isfile(source):
                warn(msg + f"file {path} does not exist; it cannot be written")
            inside = not pure.is_absolute() and not os.path.normpath(
                os.path.join(*pure.parts)
            ).startswith("..")
            entries.append(
                {
                    "filetype": ftype,
                    "unit": unit,
                    "path": path,
                    "options": options,
                    "source": source,
                    "source_ws": nam_dir,
                    "inside_ws": inside,
                }
            )
        return entries

    @classmethod
    def _check_preserved_units(cls, namefile_path, packages):
        """Reject a name file where an EVS/RTS entry shares its unit with
        another entry. The name-file reader keys entries by unit (unit 0 by
        the package reserved unit or the file type) and would silently keep
        only one of them; USG-T stops on a unit already in use."""
        users = {}
        with open(namefile_path) as f:
            lines = f.read().splitlines()
        for ln, line in enumerate(lines, 1):
            items = line.strip().split()
            if len(items) < 3 or items[0].startswith("#"):
                continue  # skipped or rejected by the name-file reader
            ftype = items[0].upper()
            try:
                unit = int(items[1])
            except ValueError:
                continue  # rejected by the name-file reader
            key = unit
            if unit == 0:
                reserved = getattr(packages.get(ftype.lower()), "_reservedunit", None)
                key = reserved() if reserved is not None else ftype
            users.setdefault(key, []).append((ln, ftype, unit))
        for key, entries in users.items():
            if len(entries) > 1 and any(
                t in cls._PRESERVED_NAM_TYPES for _, t, _ in entries
            ):
                desc = ", ".join(f"{t} unit {u} (line {ln})" for ln, t, u in entries)
                raise ValueError(
                    f"MfUsg.load: name-file unit {key} is used by more than one "
                    f"entry: {desc}"
                )

    def _name_file_targets(self):
        """(unit, emitted name, label) of everything write_name_file writes,
        with the names as written (outputs and REPLACE files as basename)."""
        targets = [(self.lst.unit_number[0], self.lst.file_name[0], "LIST")]
        for p in self.packagelist:
            for i in range(len(p.name)):
                if p.unit_number[i] != 0:
                    targets.append((p.unit_number[i], p.file_name[i], p.name[i]))
        for u, (ftype, fname) in self._skipped_nam_entries.items():
            targets.append((u, fname, ftype))
        for u, f, o in zip(
            self.external_units, self.external_fnames, self.external_output
        ):
            if u != 0:
                targets.append((u, os.path.basename(f) if o else f, "DATA"))
        for u, f in zip(self.output_units, self.output_fnames):
            if u != 0:
                targets.append((u, os.path.basename(f), "DATA output"))
        return targets

    def _preserved_plan(self, report=False):
        """Name-file lines and copies for the EVS/RTS entries.

        Every entry is checked before anything is written or copied: unit
        <= 0, missing file, unit or file shared with what the model writes
        (LIST, name file, packages, skipped entries, DATA and outputs, as
        written in the name file) or with another EVS/RTS entry, a different
        file at the copy destination, and the line length.
        """
        lines, copies = [], []
        if not self._preserved_nam_entries:
            return lines, copies
        ws = os.path.realpath(self.model_ws)
        fold = sys.platform in ("darwin", "win32")  # case-insensitive by default

        def key(path):
            k = os.path.normcase(os.path.realpath(os.path.join(ws, path)))
            return k.casefold() if fold else k

        units, names = {}, {key(self.namefile): "name file"}
        for u, fname, label in self._name_file_targets():
            units.setdefault(u, label)
            names.setdefault(key(fname), label)
        for e in self._preserved_nam_entries:
            ftype, unit = e["filetype"], e["unit"]
            if unit <= 0:
                raise ValueError(
                    f"MfUsg: the {ftype} entry has unit {unit}; a unit "
                    "<= 0 is not valid for a text time-series file"
                )
            if not os.path.isfile(e["source"]):
                raise FileNotFoundError(f"MfUsg: {ftype} file {e['source']} not found")
            if unit in units:
                raise ValueError(
                    f"MfUsg: unit {unit} of the {ftype} entry is also "
                    f"used by {units[unit]}"
                )
            units[unit] = ftype
            if os.path.realpath(e["source_ws"]) == ws:
                text, dest = e["path"], e["source"]
            elif e["inside_ws"]:
                rel = os.path.relpath(e["source"], e["source_ws"])
                text, dest = e["path"], os.path.join(ws, rel)
                if os.path.exists(dest):
                    if not os.path.isfile(dest) or not filecmp.cmp(
                        e["source"], dest, shallow=False
                    ):
                        raise FileExistsError(
                            f"MfUsg: {dest} exists with different content; "
                            f"the {ftype} file is not replaced"
                        )
                else:
                    copies.append((e["source"], dest))
            else:
                text = dest = e["source"]
                if report:
                    warn(
                        f"MfUsg: the {ftype} entry refers to {text} outside the "
                        "model workspace; the model still depends on that file"
                    )
            if key(dest) in names:
                raise ValueError(
                    f"MfUsg: the {ftype} file {dest} is also written "
                    f"by {names[key(dest)]}"
                )
            names[key(dest)] = ftype
            if "'" in text:
                raise ValueError(
                    f"MfUsg: the {ftype} file name {text} contains a quote"
                )
            if any(c in text for c in " ,\t"):
                text = f"'{text}'"
            line = f"{ftype:14s} {unit:5d}  {text}"
            line += "".join(f" {o}" for o in e["options"])
            if len(line) > self._NAM_LINE_MAX:
                raise ValueError(
                    f"MfUsg: the {ftype} name-file line is longer than "
                    f"{self._NAM_LINE_MAX} characters"
                )
            lines.append(line + "\n")
        return lines, copies

    def write_input(self, SelPackList=False, check=False):
        """Write the input; EVS/RTS entries are checked before any file is
        written (see Notes of the class)."""
        self._preserved_plan()
        super().write_input(SelPackList=SelPackList, check=check)

    def write_name_file(self):
        """Write the name file, with the EVS/RTS entries of a loaded model.

        Their files are copied first when the workspace changed (exclusive
        creation; an identical destination is left as is). This is not a
        transaction: an I/O error during the copies or the name file can
        leave the files already copied.
        """
        lines, copies = self._preserved_plan(report=True)
        for src, dest in copies:
            os.makedirs(os.path.dirname(dest), exist_ok=True)
            with open(src, "rb") as fs, open(dest, "xb") as fd:
                shutil.copyfileobj(fs, fd)
            shutil.copystat(src, dest)
        super().write_name_file()
        if lines:
            with open(os.path.join(self.model_ws, self.namefile), "a") as f:
                f.writelines(lines)

    @classmethod
    def load(
        cls,
        f: str,
        version="mfusg",
        exe_name: Union[str, PathLike] = "mfusg",
        verbose=False,
        model_ws: Union[str, PathLike] = curdir,
        load_only=None,
        forgive=False,
        check=True,
    ):
        """Load an existing MODFLOW-USG model.

        Parameters
        ----------
        f : str
            Name of MODFLOW name file to load.
        version : str, default "mfusg"
            MODFLOW version. Must be "mfusg".
        exe_name : str, default "mfusg"
            MODFLOW executable name.
        verbose : bool, default False
            Show messages that can be useful for debugging.
        model_ws : str or PathLike, default "." (curdir)
            Model workspace path. Default is the current directory.
        load_only : list, str or None
            List of case insensitive packages to load, e.g. ["bas6", "lpf"].
            One package can also be specified, e.g. "rch". Default is None,
            which attempts to load all files. An empty list [] will not load
            any additional packages than is necessary. At a minimum, "dis" or
            "disu" is always loaded.
        forgive : bool, optional, default False
            Option to raise exceptions on package load failure, which can be
            useful for debugging. Default False.
        check : boolean, optional
            Check model input for common errors. Default True.

        Returns
        -------
        flopy.mfusg.MfUsg

        Examples
        --------
        >>> import flopy
        >>> ml = flopy.mfusg.MfUsg.load('model.nam')
        """
        if version != "mfusg":
            version = "mfusg"

        # similar to modflow command: if file does not exist , try file.nam
        namefile_path = os.path.join(model_ws, f)
        if not os.path.isfile(namefile_path) and os.path.isfile(f"{namefile_path}.nam"):
            namefile_path += ".nam"
        if not os.path.isfile(namefile_path):
            raise OSError(f"cannot find name file: {namefile_path}")

        # Determine model name from 'f', without any extension or path
        modelname = os.path.splitext(os.path.basename(f))[0]

        if verbose:
            print(f"\nCreating new model with name: {modelname}\n{50 * '-'}\n")

        attribs = mfreadnam.attribs_from_namfile_header(os.path.join(model_ws, f))

        model = cls(
            modelname,
            exe_name=exe_name,
            verbose=verbose,
            model_ws=model_ws,
            **attribs,
        )

        # an EVS/RTS unit shared with another entry would be lost by the
        # unit-keyed name-file reader: reject it before reading
        cls._check_preserved_units(namefile_path, model.mfnam_packages)

        # read name file
        ext_unit_dict = mfreadnam.parsenamefile(
            namefile_path, model.mfnam_packages, verbose=verbose
        )
        if model.verbose:
            print(
                f"\n{50 * '-'}\nExternal unit dictionary:\n"
                f"{ext_unit_dict}\n{50 * '-'}\n"
            )

        # create a dict where key is the package name, value is unitnumber
        ext_pkg_d = {v.filetype: k for (k, v) in ext_unit_dict.items()}

        # reset version based on packages in the name file
        if "DISU" in ext_pkg_d:
            model.structured = False

        if "DPF" in ext_pkg_d:
            model.idpf = 1

        if "CLN" in ext_pkg_d:
            model.icln = 1

        # reset unit number for list file
        if "LIST" in ext_pkg_d:
            unitnumber = ext_pkg_d["LIST"]
            filepth = os.path.basename(ext_unit_dict[unitnumber].filename)
            model.lst.unit_number = [unitnumber]
            model.lst.file_name = [filepth]

        # look for the free format flag in bas6
        bas_key = ext_pkg_d.get("BAS6")
        if bas_key is not None:
            bas = ext_unit_dict[bas_key]
            start = bas.filehandle.tell()
            line = bas.filehandle.readline()
            while line.startswith("#"):
                line = bas.filehandle.readline()
            if "FREE" in line.upper():
                model.free_format_input = True
            bas.filehandle.seek(start)
        if verbose:
            print(f"ModflowBas6 free format:{model.free_format_input}\n")

        # EVS/RTS entries are kept with any load_only (not interpreted)
        model._preserved_nam_entries = cls._read_preserved_entries(namefile_path)

        files_successfully_loaded, files_not_loaded = cls._load_packages(
            model, ext_unit_dict, ext_pkg_d, load_only, forgive
        )

        # set up binary output / external file units
        cls._set_output_external(model, ext_unit_dict)

        # send messages re: success/failure of loading
        cls._send_load_messages(model, files_successfully_loaded, files_not_loaded)

        if check:
            model.check(f=f"{model.name}.chk", verbose=model.verbose, level=0)

        # return model object
        return model

    @classmethod
    def _load_packages(cls, model, ext_unit_dict, ext_pkg_d, load_only, forgive):
        """
        Method to load packages into the MODFLOW-USG Model Class.
        For internal class use - should not be called by the user.

        Parameters
        ----------
        model : MfUsg model object
        ext_unit_dict : dict
            For each file listed in the name file, a
            :class:`flopy.utils.mfreadnam.NamData` instance.
            Keyed by unit number.
        ext_pkg_d : dict
            key is package name, value is unitnumber
        load_only : list, str or None
            List of case insensitive packages to load, e.g. ["bas6", "lpf"].
            One package can also be specified, e.g. "rch". Default is None,
            which attempts to load all files. An empty list [] will not load
            any additional packages than is necessary. At a minimum, "dis" or
            "disu" is always loaded.
        forgive : bool
            Option to raise exceptions on package load failure.

        Returns
        -------
        files_successfully_loaded : list of loaded files
        files_not_loaded : list of files that were not loaded
        """
        files_successfully_loaded = []
        files_not_loaded = []

        # load dis
        dis_key = ext_pkg_d.get("DIS") or ext_pkg_d.get("DISU")
        if dis_key is None:
            raise KeyError("discretization entry not found in nam file")
        disnamdata = ext_unit_dict[dis_key]
        dis = disnamdata.package.load(
            disnamdata.filehandle,
            model,
            ext_unit_dict=ext_unit_dict,
            check=False,
        )
        files_successfully_loaded.append(disnamdata.filename)
        if model.verbose:
            print(f"   {dis.name[0]:4s} package load...success")
        assert model.pop_key_list.pop() == dis_key
        ext_unit_dict.pop(dis_key).filehandle.close()

        dis.start_datetime = model.start_datetime

        # BCT has to be loaded before other transport packages for MFUSG-TRANSPORT
        bct_key = ext_pkg_d.get("BCT")
        if bct_key is not None:
            bctnamdata = ext_unit_dict[bct_key]
            bct = bctnamdata.package.load(
                bctnamdata.filehandle,
                model,
                ext_unit_dict=ext_unit_dict,
                check=False,
            )
            files_successfully_loaded.append(bctnamdata.filename)
            if model.verbose:
                print(f"   {bct.name[0]:4s} package load...success")
            ext_unit_dict.pop(bct_key).filehandle.close()

        # set mfpar / pval
        cls._set_mfpar_pval(model, ext_unit_dict, ext_pkg_d)

        load_only = cls._prepare_load_only(load_only, ext_pkg_d)

        # try loading packages in ext_unit_dict
        for key, item in ext_unit_dict.items():
            if item.package is not None:
                in_load_only = item.filetype in load_only
                (files_successfully_loaded, files_not_loaded) = (
                    cls._load_ext_unit_dict_paks(
                        model,
                        ext_unit_dict,
                        load_only,
                        item,
                        forgive,
                        files_successfully_loaded,
                        files_not_loaded,
                    )
                )
                # Packages skipped by load_only are preserved in the model so
                # write_name_file can include their NAM entries verbatim.  This
                # allows verbatim-copied files (e.g. SMS, OC) to be referenced
                # correctly without FloPy needing to parse or re-write them.
                if not in_load_only:
                    model._skipped_nam_entries[key] = (
                        item.filetype,
                        os.path.basename(item.filename),
                    )
            elif "data" not in item.filetype.lower():
                files_not_loaded.append(item.filename)
                if model.verbose:
                    print(f"   {item.filetype:4s} package load...skipped")
            elif "data" in item.filetype.lower():
                cls._prepare_external_files(model, key, item)
            else:
                raise KeyError(f"unhandled case: {key}, {item}")

        return files_successfully_loaded, files_not_loaded

    @staticmethod
    def _prepare_load_only(load_only, ext_pkg_d):
        """Prepare load_only list."""
        if load_only is None:
            # load all packages/files
            load_only = ext_pkg_d.keys()
        else:  # check items in list
            if not isinstance(load_only, list):
                load_only = [load_only]
            not_found = []
            for idx, filetype in enumerate(load_only):
                load_only[idx] = filetype = filetype.upper()
                if filetype not in ext_pkg_d:
                    not_found.append(filetype)
            if not_found:
                raise KeyError(
                    "the following load_only entries were not found "
                    "in the ext_unit_dict: " + str(not_found)
                )
        return load_only

    @classmethod
    def _load_ext_unit_dict_paks(
        cls,
        model,
        ext_unit_dict,
        load_only,
        item,
        forgive,
        files_successfully_loaded,
        files_not_loaded,
    ):
        """Load packages from ext_unit_dict."""
        if item.filetype in load_only:
            if forgive:
                try:
                    cls._ext_unit_d_load(model, ext_unit_dict, item)
                    files_successfully_loaded.append(item.filename)
                    if model.verbose:
                        print(
                            f"   {item.filetype:4s} package \
                            load...success"
                        )
                except PackageLoadException as err:
                    model.load_fail = True
                    if model.verbose:
                        raise PackageLoadException(
                            error=f"   {item.filetype:4s} package \
                            load...failed"
                        ) from err
                    files_not_loaded.append(item.filename)
            else:
                cls._ext_unit_d_load(model, ext_unit_dict, item)
                files_successfully_loaded.append(item.filename)
                if model.verbose:
                    print(f"   {item.filetype:4s} package load...success")
        else:
            if model.verbose:
                print(f"   {item.filetype:4s} package load...skipped")
            files_not_loaded.append(item.filename)

        return files_successfully_loaded, files_not_loaded

    @staticmethod
    def _prepare_external_files(model, key, item):
        """Prepare external files for ext_unit_dict item."""
        if model.verbose:
            print(f"   {item.filetype} package load...skipped")
            print(f"      {os.path.basename(item.filename)}")
        if key not in model.pop_key_list:
            if key in model.output_units:
                # A package registered this unit with a generated filename
                # (e.g., WEL AUTOFLOWREDUCE writes MDV.afr but the NAM has
                # MDV_FlowReduction.dat). Override with the actual NAM filename.
                idx = model.output_units.index(key)
                model.output_fnames[idx] = item.filename
            elif key not in model.external_units:
                model.external_fnames.append(item.filename)
                model.external_units.append(key)
                model.external_binflag.append("binary" in item.filetype.lower())
                model.external_output.append(getattr(item, "replace", False))

    @staticmethod
    def _ext_unit_d_load(model, ext_unit_dict, ext_unit_d_item):
        """
        Method to load a package from an ext_unit_dict item into model

        Parameters
        ----------
        model : MfUsg model object for which package in ext_unit_d_item will
            be loaded
        ext_unit_dict : dict
            For each file listed in the name file, a
            :class:`flopy.utils.mfreadnam.NamData` instance.
            Keyed by unit number.
        ext_unit_d_item : :class:`flopy.utils.mfreadnam.NamData` instance.
            Must be an item of ext_unit_dict.
        """
        package_load_args = getfullargspec(ext_unit_d_item.package.load)[0]
        if "check" in package_load_args:
            ext_unit_d_item.package.load(
                ext_unit_d_item.filehandle,
                model,
                ext_unit_dict=ext_unit_dict,
                check=False,
            )
        else:
            ext_unit_d_item.package.load(
                ext_unit_d_item.filehandle,
                model,
                ext_unit_dict=ext_unit_dict,
            )

    @staticmethod
    def _set_mfpar_pval(model, ext_unit_dict, ext_pkg_d):
        """Set mfpar/pval items."""
        # zone, mult, pval
        if "PVAL" in ext_pkg_d:
            model.mfpar.set_pval(model, ext_unit_dict)
            assert model.pop_key_list.pop() == ext_pkg_d.get("PVAL")
        if "ZONE" in ext_pkg_d:
            model.mfpar.set_zone(model, ext_unit_dict)
            assert model.pop_key_list.pop() == ext_pkg_d.get("ZONE")
        if "MULT" in ext_pkg_d:
            model.mfpar.set_mult(model, ext_unit_dict)
            assert model.pop_key_list.pop() == ext_pkg_d.get("MULT")

    @staticmethod
    def _set_output_external(model, ext_unit_dict):
        """Set up binary output / external file units."""
        # pop binary output keys and any external file units that are now
        # internal
        for key in model.pop_key_list:
            try:
                model.remove_external(unit=key)
                item = ext_unit_dict.pop(key)
                if hasattr(item.filehandle, "close"):
                    item.filehandle.close()
            except KeyError:
                if model.verbose:
                    print(
                        f"\nWARNING:\n    External file unit {key} does not "
                        "exist in ext_unit_dict."
                    )

    @staticmethod
    def _send_load_messages(model, files_successfully_loaded, files_not_loaded):
        """Send messages re: success/failure of loading."""
        # write message indicating packages that were successfully loaded
        if model.verbose:
            print("")
            print(
                f"   The following {len(files_successfully_loaded)} packages "
                "were successfully loaded."
            )
            for fname in files_successfully_loaded:
                print(f"      {os.path.basename(fname)}")
            if len(files_not_loaded) > 0:
                print(
                    f"   The following {len(files_not_loaded)} packages "
                    "were not loaded."
                )
                for fname in files_not_loaded:
                    print(f"      {os.path.basename(fname)}")


def _urwords(line):
    """Split a name-file line like URWORD (utl7u1.f): words end at a space,
    comma or tab; a word starting with a single quote ends at the next quote
    and can hold spaces and commas. Returns (word, closed) pairs; closed is
    False for a quoted word without its closing quote."""
    words, i, n = [], 0, len(line)
    while i < n:
        if line[i] in " ,\t":
            i += 1
            continue
        if line[i] == "'":
            j = line.find("'", i + 1)
            if j < 0:
                words.append((line[i + 1 :], False))
                break
            words.append((line[i + 1 : j], True))
            i = j + 1
        else:
            j = i
            while j < n and line[j] not in " ,\t":
                j += 1
            words.append((line[i:j], True))
            i = j
    return words


def fmt_string(array, free=False):
    """
    Returns a C-style fmt string for numpy savetxt that corresponds to
    the dtype.

    Parameters
    ----------
    array : numpy array
    free : bool, optional
        If True, use high-precision format (%16.9G) for floats which requires
        FREE format in the BAS file. If False (default), use fixed 10-character
        format (%10.2e) compatible with MODFLOW-USG's fixed-width input format.
    """
    fmts = []
    # When FREE format is enabled in BAS, we can use high-precision output.
    # Without FREE, MODFLOW-USG expects 10-character fixed-width fields.
    float_fmt = "%16.9G" if free else "%10.2e"
    for field in array.dtype.descr:
        vtype = field[1][1].lower()
        if vtype in {"i", "b"}:
            fmts.append("%10d")
        elif vtype == "f":
            fmts.append(float_fmt)
        elif vtype == "o":
            fmts.append("%10s")
        elif vtype == "s":
            msg = (
                "mfusg.fmt_string error: 'str' type found in dtype."
                "This gives unpredictable results when"
                "recarray to file - change to 'object' type"
            )
            raise TypeError(msg)
        else:
            raise TypeError(f"mfusg.fmt_string error: unknown vtype in field: {field}")
    # a free format list is read with URWORD, so the fields are separated the
    # way MfList.fmt_string separates them; a fixed format list is read by
    # position and relies on the field widths
    return (" " if free else "").join(fmts)


# Per-component transport factors of EVT (ETFACTOR) and ETS (ESFACTOR):
# USG-T reads MCOMPT values, READ(IN,*) in free format and FORMAT(20F10.3)
# in fixed format (gwf2evt8u1.f, gwf2ets8u1.f).


def _transport_factors(who, value, mcompt):
    """Return MCOMPT float factors: a scalar is used for every component; a
    list or array must have exactly MCOMPT values."""
    if value is None:
        raise ValueError(f"{who}: {mcompt} value(s) (MCOMPT) are required")
    arr = np.asarray(value, dtype=float)
    if arr.ndim == 0:
        arr = np.full(mcompt, float(arr))
    elif arr.ndim != 1 or arr.size != mcompt:
        raise ValueError(
            f"{who}: must have MCOMPT={mcompt} value(s) (mobile components "
            f"and temperature); got {arr.size}"
        )
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{who}: values must be finite")
    return arr


def _factor_text(values, free, fmt, who):
    """Text of the factor records: one line in free format (fields always
    separated by a blank), records of 20 fields of 10 columns in fixed
    format (FORMAT(20F10.3)); a value that does not fit 10 columns falls
    back to 10.3e, and raises ValueError if it still does not fit. Callers
    build it before opening the file."""
    fields = []
    for v in values:
        s = format(float(v), fmt)
        if free:
            if not s.startswith(" "):
                s = " " + s
        elif len(s) > 10:
            s = format(float(v), "10.3e")
            if len(s) > 10:
                raise ValueError(
                    f"{who}: factor {float(v)!r} does not fit a 10-column field "
                    "(FORMAT(20F10.3))"
                )
        fields.append(s)
    if free:
        return "".join(fields) + "\n"
    return "".join(
        "".join(fields[i : i + 20]) + "\n" for i in range(0, len(fields), 20)
    )


def _read_free_values(f, n, who):
    """Read n reals as list-directed READ(IN,*): records are read until n
    values are found; the rest of the last record is skipped. Values are
    separated by blanks or by one comma with optional blanks (an end of record
    counts as a blank), and r*c repeats are supported. EOF, a slash, a null
    value (a comma before the first value, two commas with only blanks or ends
    of record between them, or r* without a constant) or an invalid token
    raise ValueError. A short record followed by a record of numbers cannot
    be told apart from a continued record (as in Fortran)."""
    import re

    vals = []
    comma = True  # a comma here would give a null value (start of the read)
    while len(vals) < n:
        line = f.readline()
        if line == "":
            raise ValueError(f"{who}: end of file after {len(vals)} of {n} values")
        for tok in re.findall(r",|[^\s,]+", line):
            if len(vals) >= n:
                break
            if tok == ",":
                if comma:
                    raise ValueError(
                        f"{who}: null value (empty field) in the factor record "
                        "is not supported"
                    )
                comma = True
                continue
            if "/" in tok:
                raise ValueError(f"{who}: slash in the factor record is not supported")
            rep = 1
            if "*" in tok:
                r, _, tok = tok.partition("*")
                try:
                    rep = int(r)
                except ValueError:
                    raise ValueError(f"{who}: invalid repeat {r!r}*") from None
                if rep < 1:
                    raise ValueError(f"{who}: invalid repeat token {r}*{tok}")
                if tok == "":
                    raise ValueError(
                        f"{who}: null value ({r}*) in the factor record is not "
                        "supported"
                    )
            try:
                v = float(tok.replace("d", "e").replace("D", "E"))
            except ValueError:
                raise ValueError(f"{who}: invalid factor value {tok!r}") from None
            vals.extend([v] * min(rep, n - len(vals)))
            comma = False
    return np.array(vals)


def _fixed_field(field, decimals, who):
    """Value of a FORMAT Fw.d (decimals) or Iw (decimals None) field: blanks
    are ignored (BN), an empty field is 0, and a real without a decimal point
    has ``decimals`` implied decimals."""
    import re

    s = field.replace(" ", "")
    if s == "":
        return 0 if decimals is None else 0.0
    try:
        if decimals is None:
            return int(s)
        m = re.fullmatch(r"([+-]?[0-9.]*)(?:[eEdD]?([+-]?[0-9]+))?", s)
        if m is None or m.group(1) in ("", "+", "-", "."):
            raise ValueError
        mant, exp = m.group(1), m.group(2)
        val = float(mant) if "." in mant else int(mant) / 10.0**decimals
        return val * 10.0 ** int(exp) if exp else val
    except ValueError:
        raise ValueError(f"{who}: invalid fixed-format field {field!r}") from None


def _read_fixed_values(f, n, who):
    """Read n reals with FORMAT(20F10.3): ceil(n/20) records of 10-column
    fields; missing columns are blank (0.0)."""
    vals = []
    while len(vals) < n:
        line = f.readline()
        if line == "":
            raise ValueError(f"{who}: end of file after {len(vals)} of {n} values")
        line = line.rstrip("\n").rstrip("\r")
        for i in range(min(20, n - len(vals))):
            vals.append(_fixed_field(line[10 * i : 10 * i + 10], 3, who))
    return np.array(vals)


def _fixed_ints(line, n, who):
    """First n integer fields of a FORMAT(nI10) record."""
    line = line.rstrip("\n").rstrip("\r")
    return [_fixed_field(line[10 * i : 10 * i + 10], None, who) for i in range(n)]
