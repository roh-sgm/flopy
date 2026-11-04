"""
Mfusgbcf module.

Contains the MfUsgBcf class. Note that the user can
access the MfUsgBcf class as `flopy.mfusg.MfUsgBcf`.
"""

import re

import numpy as np

from ..modflow import ModflowBcf
from ..utils import Util2d, Util3d
from ..utils.flopy_io import line_parse
from ..utils.utils_def import (
    get_open_file_object,
    get_unitnumber_from_ext_unit_dict,
    get_util2d_shape_for_layer,
    type_from_iterable,
)
from ._tabrich import (
    make_iuzontab,
    read_tabrich,
    validate_retcrvs,
    write_tabrich,
)
from .mfusg import MfUsg

_RICHARDS_NAMES = ("alpha", "beta", "sr", "brook", "bp")


def _effective_laycon(laycon, richards):
    """LAYCON used by the BCF readers. FloPy already stores the units digit of
    the file code (READBCF 248-257); the BAS RICHARDS/RICHARDS_HP option then
    turns every LAYCON != 0 into 5 (READBCF 335-350). The input is not mutated."""
    laycon = np.array(laycon, dtype=np.int32)
    if richards:
        laycon[laycon != 0] = 5
    return laycon


def _bas_richards(model):
    bas = model.get_package("BAS6")
    return bool(bas is not None and getattr(bas, "richards", False))


def _is_anisotropic(trpy):
    """IHANISO test of READBCF (357-360)."""
    return bool(np.any(np.abs(np.asarray(trpy, dtype=float) - 1.0) > 1.0e-6))


def _fortran_int_field(value, width=10):
    """Right-justified integer for an Iw field; fails if it does not fit."""
    text = f"{int(value):{width}d}"
    if len(text) > width:
        raise ValueError(f"{value} does not fit a {width}-character integer field")
    return text


_FIXED_REAL = re.compile(r"^([+-]?)(\d*)(?:\.(\d*))?(?:[EDQ]([+-]?\d+)|([+-]\d+))?$")


def _read_fixed_real(text, decimals=3):
    """Value of a Fortran Fw.d input field (default d=3, as F10.3). Blanks are
    ignored (BLANK='NULL') and an empty field is 0. Without a decimal point
    in the mantissa its last d digits are fractional, also when an exponent
    (E, D, Q or a signed integer) follows: '-1E+30' is -1E+27, '123E+2' is
    12.3. A decimal point in the mantissa overrides d."""
    field = text.replace(" ", "").upper()
    if not field:
        return 0.0
    match = _FIXED_REAL.match(field)
    if match is None or not (match.group(2) or match.group(3)):
        raise ValueError(f"invalid Fortran real field: {text!r}")
    sign, whole, frac, exp1, exp2 = match.groups()
    exponent = int(exp1 or exp2 or 0)
    if frac is None:
        return float(f"{sign}{whole}e{exponent - decimals}")
    return float(f"{sign}{whole or '0'}.{frac or '0'}e{exponent}")


def _fortran_real_field(value, width=10):
    """Shortest right-justified text of at most `width` characters that an
    F10.3 read (`_read_fixed_real`) returns as the same single-precision
    value. Every candidate carries a decimal point in the mantissa, because
    an exponent alone does not cancel the implied decimals (-1e30 is written
    as '-1.E+30'); up to nine significant digits are tried in positional and
    exponential form. Fails if the value cannot be kept within the field."""
    target = np.float32(value)
    candidates = []
    for digits in range(1, 10):
        text = f"{float(value):.{digits}G}"
        if "E" in text:
            mantissa, exponent = text.split("E")
            if "." not in mantissa:
                mantissa += "."
            candidates.append(f"{mantissa}E{exponent}")
        elif "." not in text:
            candidates.append(text + ".")
        else:
            candidates.append(text)
        exponential = f"{float(value):.{digits - 1}E}"
        mantissa, exponent = exponential.split("E")
        if "." not in mantissa:
            mantissa += "."
        candidates.append(f"{mantissa}E{exponent}")
        positional = f"{float(value):.{digits - 1}f}"
        candidates.append(positional if "." in positional else positional + ".")
    valid = [
        text
        for text in candidates
        if len(text) <= width and np.float32(_read_fixed_real(text)) == target
    ]
    if not valid:
        raise ValueError(
            f"{value} cannot be written in a {width}-character F{width}.3 field "
            "without changing its single-precision value"
        )
    return min(valid, key=len).rjust(width)


class MfUsgBcf(ModflowBcf):
    """Block Centered Flow (BCF) Package Class for MODFLOW-USG.

    Parameters
    ----------
    model : model object
        The model object (of type :class:`flopy.modflow.Modflow`) to which
        this package will be added.
    ipakcb : int, optional
        Toggles whether cell-by-cell budget data should be saved. If None or zero,
        budget data will not be saved (default is None).
    intercellt : int
        Intercell transmissivities, harmonic mean (0), arithmetic mean (1),
        logarithmic mean (2), combination (3). (default is 0)
    laycon : int
        Layer type, confined (0), unconfined (1), constant T, variable S (2),
        variable T, variable S (default is 3)
    trpy : float or array of floats (nlay)
        horizontal anisotropy ratio (default is 1.0)
    hdry : float
        head assigned when cell is dry - used as indicator(default is -1E+30)
    iwdflg : int
        flag to indicate if wetting is inactive (0) or not (non zero)
        (default is 0)
    wetfct : float
        factor used when cell is converted from dry to wet (default is 0.1)
    iwetit : int
        iteration interval in wetting/drying algorithm (default is 1)
    ihdwet : int
        flag to indicate how initial head is computed for cells that become
        wet (default is 0)
    ikvflag : int
        flag indicating if vertical hydraulic conductivity is input
        instead of leakance between two layers.
    ikcflag : int
        flag indicating if hydraulic conductivity or transmissivity
        information is input for each of the nodes or whether this information
        is directly input for the nodal connections. The easiest input format
        is to provide the hydraulic conductivity or transmissivity values to
        the cells using a zero value for IKCFLAG.
    anglex : float or array of floats (njag)
        is the angle (in radians) between the horizontal x-axis and the outward
        normal to the face between a node and its connecting nodes. The angle
        varies between zero and 6.283185 (two pi being 360 degrees).
    tran : float or array of floats (nlay, nrow, ncol), optional
        transmissivity (only read if laycon is 0 or 2) (default is 1.0)
    hy : float or array of floats (nlay, nrow, ncol)
        hydraulic conductivity (only read if laycon is 1 or 3)
        (default is 1.0)
    vcont : float or array of floats (nlay-1, nrow, ncol)
        vertical leakance between layers (default is 1.0)
    kv : float or array of floats (nlay, nrow, ncol)
        is the vertical hydraulic conductivity of the cell and the leakance is
        computed for each vertical connection.
    sf1 : float or array of floats (nlay, nrow, ncol)
        specific storage (confined) or specific yield (unconfined),
        read when there is at least one transient stress period.
        (default is 1e-5)
    sf2 : float or array of floats (nlay, nrow, ncol)
        specific yield, only read when laycon is 2 or 3 and there is at least
        one transient stress period (default is 0.15)
    wetdry : float
        a combination of the wetting threshold and a flag to indicate which
        neighboring cells can cause a cell to become wet (default is -0.01)
    ksat : float or array of floats (njag)
        inter-block saturated hydraulic conductivity or transmissivity
        (if IKCFLAG = 1) or the inter-block conductance (if IKCFLAG = - 1)
        of the connection between nodes n and m. Read for unstructured grids
        when IKCFLAG is not 0.
    extension : string
        Filename extension (default is 'bcf')
    unitnumber : int
        File unit number (default is None).
    filenames : str or list of str
        Filenames to use for the package and the output files. If
        filenames=None the package name will be created using the model name
        and package extension and the cbc output name will be created using
        the model name and .cbc extension (for example, modflowtest.cbc),
        if ipakcb is a number greater than zero. If a single string is passed
        the package will be set to the string and cbc output name will be
        created using the model name and .cbc extension, if ipakcb is a
        number greater than zero. To define the names for all package files
        (input and output) the length of the list of strings should be 2.
        Default is None.
    add_package : bool
        Flag to add the initialised package object to the parent model object.
        Default is True.
    alpha, beta, sr, brook : float or array of floats (nlay, nrow, ncol)
        Richards-equation parameters, read for layers whose effective LAYCON
        is 5 (LAYCON 5, or any LAYCON != 0 with the BAS RICHARDS/RICHARDS_HP
        option) when TABRICH is off (defaults 1.0, 7.0, 0.05, 6.0).
    bp : float or array of floats (nlay, nrow, ncol)
        Bubbling point head for the same layers, read with BUBBLEPT
        (default 0.0).

    Methods
    -------

    See Also
    --------

    Notes
    -----
    IKCFLAG is only read for unstructured grids; for a structured model it
    is treated as 0 when writing and loading. In a fixed-format BCF file (BAS
    without FREE) the header options cannot be read, so writing one with
    TABRICH, BUBBLEPT, FULLYDRY or ALTSTO raises a ValueError.

    Examples
    --------

    >>> import flopy
    >>> ml = flopy.mfusg.MfUsg()
    >>> disu = flopy.mfusg.MfUsgDisU(model=ml, nlay=1, nodes=1,
                 iac=[1], njag=1,ja=np.array([0]), fahl=[1.0], cl12=[1.0])
    >>> bcf = flopy.mfusg.MfUsgBcf(ml)"""

    def __init__(
        self,
        model,
        ipakcb=None,
        intercellt=0,
        laycon=3,
        trpy=1.0,
        hdry=-1e30,
        iwdflg=0,
        wetfct=0.1,
        iwetit=1,
        ihdwet=0,
        ikvflag=0,
        ikcflag=0,
        tabrich=False,
        nuzones=0,
        nutabrows=0,
        iuzontab=None,
        retcrvs=None,
        bubblept=False,
        fullydry=False,
        altsto=False,
        tran=1.0,
        hy=1.0,
        vcont=1.0,
        kv=1.0,
        anglex=0.0,
        ksat=1.0,
        sf1=1e-5,
        sf2=0.15,
        wetdry=-0.01,
        extension="bcf",
        unitnumber=None,
        filenames=None,
        add_package=True,
        alpha=1.0,
        beta=7.0,
        sr=0.05,
        brook=6.0,
        bp=0.0,
    ):
        """Constructs the MfUsgBcf object.

        Overrides the parent ModflowBcf object.
        """
        msg = (
            "Model object must be of type flopy.mfusg.MfUsg\n"
            f"but received type: {type(model)}."
        )
        assert isinstance(model, MfUsg), msg

        super().__init__(
            model,
            ipakcb=ipakcb,
            intercellt=intercellt,
            laycon=laycon,
            trpy=trpy,
            hdry=hdry,
            iwdflg=iwdflg,
            wetfct=wetfct,
            iwetit=iwetit,
            ihdwet=ihdwet,
            tran=tran,
            hy=hy,
            vcont=vcont,
            sf1=sf1,
            sf2=sf2,
            wetdry=wetdry,
            extension=extension,
            unitnumber=unitnumber,
            filenames=filenames,
            add_package=False,
        )

        dis = model.get_package("DIS")
        if dis is None:
            dis = model.get_package("DISU")
        structured = self.parent.structured

        nrow, ncol, nlay, _ = self.parent.nrow_ncol_nlay_nper

        self.ikvflag = ikvflag
        self.ikcflag = ikcflag
        self.tabrich = tabrich
        self.nuzones = nuzones
        self.nutabrows = nutabrows
        self.bubblept = bubblept
        self.fullydry = fullydry
        self.altsto = altsto

        # TABRICH items 1c (IUZONTAB zone map) and 1d (RETCRVS curves).
        # Stored if provided; completeness is enforced at write_file time so
        # the package can still be used purely as model setup (e.g. enabling
        # Richards mode for a dependent DPF) without authoring the curves.
        if tabrich and iuzontab is not None:
            self.iuzontab = make_iuzontab(model, iuzontab)
        else:
            self.iuzontab = None
        if tabrich and retcrvs is not None:
            self.retcrvs = validate_retcrvs(retcrvs, nuzones, nutabrows)
        else:
            self.retcrvs = None

        self.kv = kv
        self.anglex = anglex
        self.ksat = ksat

        if not structured:
            njag = dis.njag
            self.anglex = Util2d(
                model, (njag,), np.float32, anglex, "anglex", locat=self.unit_number[0]
            )

        # item 1
        self.kv = Util3d(
            model,
            (nlay, nrow, ncol),
            np.float32,
            kv,
            "Vertical Hydraulic Conductivity",
            locat=self.unit_number[0],
        )
        if not structured:
            self.ksat = Util2d(
                model, (njag,), np.float32, ksat, "ksat", locat=self.unit_number[0]
            )

        # Richards parameters, read for layers whose effective LAYCON is 5
        # when TABRICH is off (bp only with BUBBLEPT)
        for name, value in zip(_RICHARDS_NAMES, (alpha, beta, sr, brook, bp)):
            setattr(
                self,
                name,
                Util3d(
                    model,
                    (nlay, nrow, ncol),
                    np.float32,
                    value,
                    name=f"richards {name}",
                    locat=self.unit_number[0],
                ),
            )

        if add_package:
            self.parent.add_package(self)

    def _effective_ikcflag(self):
        """IKCFLAG as consumed by READBCF: only read for unstructured grids
        (gwf2bcf-lpf-u1.f 91-111), so it is 0 for a structured model."""
        return 0 if self.parent.structured else int(self.ikcflag)

    def write_file(self, f=None):
        """
        Write the BCF package file.

        Parameters
        ----------
        f : open file object.
            Default is None, which will result in MfUsg.fn_path being
            opened for writing. A handle passed by the caller is written to
            and left open.
        """
        # get model information
        nlay = self.parent.nlay
        dis = self.parent.get_package("DIS")
        if dis is None:
            dis = self.parent.get_package("DISU")
        structured = self.parent.structured
        ifrefm = self.parent.get_ifrefm()
        ikcflag = self._effective_ikcflag()

        # Item 1: IBCFCB, HDRY, IWDFLG, WETFCT, IWETIT, IHDWET, IKVFLAG, <IKCFLAG>
        options = []
        if self.tabrich:
            options.append(f"TABRICH {self.nuzones:9d} {self.nutabrows:9d}")
        if self.bubblept:
            options.append("BUBBLEPT")
        if self.fullydry:
            options.append("FULLYDRY")
        if self.altsto:
            options.append("ALTSTO")
        codes = 10 * np.asarray(self.intercellt.array) + np.asarray(self.laycon.array)
        if ifrefm:
            header = (
                f" {self.ipakcb:9d} {self.hdry:9.3e} {self.iwdflg:9d}"
                f" {self.wetfct:9.3e} {self.iwetit:9d} {self.ihdwet:9d}"
                f" {self.ikvflag:9d}"
            )
            if not structured:
                header += f" {ikcflag:9d}"
            header += "".join(f" {opt}" for opt in options)
        else:
            # READ(IN,2) with FORMAT(I10,F10.3,I10,F10.3,4I10) (READBCF 91-99,
            # 113); the options are not parsed from a fixed-format line
            if options:
                raise ValueError(
                    "BCF options "
                    + ", ".join(opt.split()[0] for opt in options)
                    + " cannot be read by USG-Transport from a fixed-format "
                    "BCF file (BAS without FREE); use free format "
                    "(MfUsgBas ifrefm=True)."
                )
            if np.any((codes < 0) | (codes > 99)):
                raise ValueError(
                    f"LAYCON codes {codes.tolist()} do not fit the fixed-format "
                    "I2 fields."
                )
            fields = [
                _fortran_int_field(self.ipakcb),
                _fortran_real_field(self.hdry),
                _fortran_int_field(self.iwdflg),
                _fortran_real_field(self.wetfct),
                _fortran_int_field(self.iwetit),
                _fortran_int_field(self.ihdwet),
                _fortran_int_field(self.ikvflag),
            ]
            if not structured:
                fields.append(_fortran_int_field(ikcflag))
            header = "".join(fields)

        # Open file for writing; a handle passed by the caller is left open
        opened = f is None
        f_obj = open(self.fn_path, "w") if opened else f
        f_obj.write(header + "\n")

        # Items 1c/1d: TABRICH zone map and retention curves (before LAYCON).
        # Fail explicitly rather than write an incomplete TABRICH file.
        if self.tabrich:
            if self.iuzontab is None or self.retcrvs is None:
                raise ValueError(
                    "Writing a TABRICH BCF file requires both iuzontab and "
                    "retcrvs; got None. Provide both, or set tabrich=False."
                )
            write_tabrich(f_obj, self.iuzontab, self.retcrvs)

        # LAYCON array: one code 10*INTERCELLT+LAYCON per layer, read with
        # READ(IN,*) in free format and FORMAT(40I2) in fixed format
        for start in range(0, nlay, 40):
            chunk = codes[start : start + 40]
            if ifrefm:
                f_obj.write(" ".join(f"{code:02d}" for code in chunk) + "\n")
            else:
                f_obj.write("".join(f"{code:2d}" for code in chunk) + "\n")

        # TRPY (read when IKCFLAG=0) and ANGLEX (DISU, horizontal anisotropy)
        if ikcflag == 0:
            f_obj.write(self.trpy.get_file_entry())
            if not structured and _is_anisotropic(self.trpy.array):
                f_obj.write(self.anglex.get_file_entry())

        # <SF1>, <TRAN>, <HY>, <VCONT>, <KV>, <SF2>, <WETDRY>, Richards arrays
        transient = not dis.steady.all()
        laycon = _effective_laycon(self.laycon.array, _bas_richards(self.parent))
        for layer in range(nlay):
            if transient:
                f_obj.write(self.sf1[layer].get_file_entry())

            if ikcflag == 0:
                self._write_hy_tran_vcont_kv(f_obj, layer, laycon[layer])

            if transient and laycon[layer] in {2, 3, 4, 5}:
                f_obj.write(self.sf2[layer].get_file_entry())

            if self.iwdflg != 0 and laycon[layer] in {1, 3}:
                f_obj.write(self.wetdry[layer].get_file_entry())

            if laycon[layer] == 5 and not self.tabrich:
                f_obj.write(self.alpha[layer].get_file_entry())
                f_obj.write(self.beta[layer].get_file_entry())
                f_obj.write(self.sr[layer].get_file_entry())
                f_obj.write(self.brook[layer].get_file_entry())
                if self.bubblept:
                    f_obj.write(self.bp[layer].get_file_entry())

        # <KSAT> connection array (DISU, IKCFLAG != 0)
        if not structured and ikcflag != 0:
            f_obj.write(self.ksat.get_file_entry())

        if opened:
            f_obj.close()

    def _write_hy_tran_vcont_kv(self, f_obj, layer, laycon_k):
        """
        writes hy/tran and vcont/kv file entries

        Parameters
        ----------
        f_obj : open file object.
        layer : model layer index (base 0)
        laycon_k : effective LAYCON of the layer
        """
        if laycon_k in {0, 2}:
            f_obj.write(self.tran[layer].get_file_entry())
        else:
            f_obj.write(self.hy[layer].get_file_entry())

        if self.ikvflag == 0:
            if layer < (self.parent.nlay - 1):
                f_obj.write(self.vcont[layer].get_file_entry())
        elif self.parent.nlay > 1:
            f_obj.write(self.kv[layer].get_file_entry())

    @classmethod
    def load(cls, f, model, ext_unit_dict=None):
        """
        Load an existing package.

        Parameters
        ----------
        f : filename or file handle
            File to load.
        model : model object
            The model object (of type :class:`flopy.modflow.mf.Modflow`) to
            which this package will be added.
        ext_unit_dict : dictionary, optional
            If the arrays in the file are specified using EXTERNAL,
            or older style array control records, then `f` should be a file
            handle.  In this case ext_unit_dict is required, which can be
            constructed using the function
            :class:`flopy.utils.mfreadnam.parsenamefile`.

        Returns
        -------
        bcf : MfUsgBcf object

        Examples
        --------

        >>> import flopy
        >>> m = flopy.mfusg.MfUsg()
        >>> disu = flopy.mfusg.MfUsgDisU(
        ...     model=m, nlay=1, nodes=1, iac=[1], njag=1,ja=np.array([0]),
        ...     fahl=[1.0], cl12=[1.0])
        >>> bcf = flopy.mfusg.MfUsgBcf.load('test.bcf', m)
        """
        msg = (
            "Model object must be of type flopy.mfusg.MfUsg\n"
            f"but received type: {type(model)}."
        )
        assert isinstance(model, MfUsg), msg

        if model.verbose:
            print("loading bcf package file...")

        f_obj = get_open_file_object(f, "r")

        # dataset 0 -- header
        while True:
            line = f_obj.readline()
            if line[0] != "#":
                break

        # determine problem dimensions
        nlay = model.nlay
        structured = model.structured
        dis = model.get_package("DIS")
        if dis is None:
            dis = model.get_package("DISU")
            njag = dis.njag

        # Item 1: ipakcb, HDRY, IWDFLG, WETFCT, IWETIT, IHDWET, IKVFLAG, <IKCFLAG>
        if model.verbose:
            print("   loading ipakcb, HDRY, IWDFLG, WETFCT, IWETIT, IHDWET...")
        text_list = line_parse(line)
        if model.get_ifrefm():
            ipakcb, hdry, iwdflg, wetfct, iwetit, ihdwet = (
                int(text_list[0]),
                float(text_list[1]),
                int(text_list[2]),
                float(text_list[3]),
                int(text_list[4]),
                int(text_list[5]),
            )
            ikvflag = type_from_iterable(text_list, index=6, _type=int, default_val=0)
            ikcflag = 0
            if not structured:
                ikcflag = type_from_iterable(
                    text_list, index=7, _type=int, default_val=0
                )
        else:
            (ipakcb, hdry, iwdflg, wetfct, iwetit, ihdwet, ikvflag, ikcflag) = (
                cls._load_fixed_header(line, structured)
            )

        # options
        if "TABRICH" in text_list:
            idx = text_list.index("TABRICH")
            tabrich = True
            nuzones = int(text_list[idx + 1])
            nutabrows = int(text_list[idx + 2])
        else:
            tabrich = False
            nuzones = None
            nutabrows = None

        if "BUBBLEPT" in text_list:
            bubblept = True
        else:
            bubblept = False

        if "FULLYDRY" in text_list:
            fullydry = True
        else:
            fullydry = False

        if "ALTSTO" in text_list:
            altsto = True
        else:
            altsto = False

        # items 1c (IUZONTAB) and 1d (RETCRVS), read before LAYCON
        iuzontab = None
        retcrvs = None
        if tabrich:
            iuzontab, retcrvs = read_tabrich(
                f_obj, model, nuzones, nutabrows, ext_unit_dict
            )

        # LAYCON array
        laycon, intercellt = cls._load_laycon(f_obj, model)

        # TRPY array (read when IKCFLAG=0) and ANGLEX (DISU, anisotropy)
        trpy = 1.0
        anglex = 0
        if ikcflag == 0:
            if model.verbose:
                print("   loading TRPY...")
            trpy = Util2d.load(f_obj, model, (nlay,), np.float32, "trpy", ext_unit_dict)
            if (not structured) and _is_anisotropic(trpy.array):
                if model.verbose:
                    print("loading ANGLEX...")
                anglex = Util2d.load(
                    f_obj, model, (njag,), np.float32, "anglex", ext_unit_dict
                )

        # hy, kv, storage, wetting and Richards arrays
        transient = not dis.steady.all()
        (sf1, tran, hy, vcont, sf2, wetdry, kv, richards_arrays) = (
            cls._load_layer_arrays(
                f_obj,
                model,
                nlay,
                ext_unit_dict,
                transient,
                _effective_laycon(laycon, _bas_richards(model)),
                ikvflag,
                ikcflag,
                iwdflg,
                tabrich,
                bubblept,
            )
        )

        # Ksat  mfusg
        ksat = 0
        if (not structured) and ikcflag != 0:
            if model.verbose:
                print("   loading ksat (njag)...")
            ksat = Util2d.load(f_obj, model, (njag,), np.float32, "ksat", ext_unit_dict)

        f_obj.close()

        # set package unit number
        unitnumber, filenames = get_unitnumber_from_ext_unit_dict(
            model, cls, ext_unit_dict, ipakcb
        )

        # create instance of bcf object
        bcf = cls(
            model,
            ipakcb=ipakcb,
            intercellt=intercellt,
            laycon=laycon,
            trpy=trpy,
            hdry=hdry,
            iwdflg=iwdflg,
            wetfct=wetfct,
            iwetit=iwetit,
            ihdwet=ihdwet,
            ikvflag=ikvflag,
            ikcflag=ikcflag,
            tabrich=tabrich,
            nuzones=nuzones,
            nutabrows=nutabrows,
            iuzontab=iuzontab,
            retcrvs=retcrvs,
            bubblept=bubblept,
            fullydry=fullydry,
            altsto=altsto,
            tran=tran,
            hy=hy,
            vcont=vcont,
            kv=kv,
            anglex=anglex,
            ksat=ksat,
            sf1=sf1,
            sf2=sf2,
            wetdry=wetdry,
            unitnumber=unitnumber,
            filenames=filenames,
            **richards_arrays,
        )

        # return bcf object
        return bcf

    @staticmethod
    def _load_fixed_header(line, structured):
        """
        Reads item 1 as READ(IN,2) with FORMAT(I10,F10.3,I10,F10.3,4I10)
        (READBCF 91-99, 113): IKCFLAG is the eighth field and is read only
        for unstructured grids; blanks are ignored and the real fields follow
        F10.3 input rules (see `_read_fixed_real`).
        """
        line = line.rstrip("\r\n")

        def text(i):
            return line[10 * i : 10 * (i + 1)].replace(" ", "")

        def ival(i):
            return int(text(i)) if text(i) else 0

        def rval(i):
            return _read_fixed_real(line[10 * i : 10 * (i + 1)])

        ikcflag = 0 if structured else ival(7)
        return (ival(0), rval(1), ival(2), rval(3), ival(4), ival(5), ival(6), ikcflag)

    @staticmethod
    def _load_laycon(f_obj, model):
        """
        Loads laycon and intercellt file entries.

        Parameters
        ----------
        f_obj : open file object.
        model : model object
            The model object (of type :class:`flopy.modflow.mf.Modflow`) to
            which this package will be added.

        Returns
        -------
        laycon : Numpy array of laycon values
        intercellt : Numpy array of intercellt values
        """
        nlay = model.nlay
        ifrefm = model.get_ifrefm()

        if model.verbose:
            print("   loading LAYCON...")

        codes = []
        if ifrefm:
            # READ(IN,*): NLAY integers, possibly on several lines; the rest of
            # the line holding the last value is not read (ex3.bcf "Ltype(NLAY)")
            while len(codes) < nlay:
                line = f_obj.readline()
                if not line:
                    raise ValueError(
                        f"end of file while reading LAYCON: {len(codes)} of "
                        f"{nlay} values read"
                    )
                codes += [int(item) for item in line.split()[: nlay - len(codes)]]
        else:
            # FORMAT(40I2): 2-character fields, blanks ignored, 40 per line
            while len(codes) < nlay:
                line = f_obj.readline().rstrip("\r\n")
                for start in range(0, 80, 2):
                    if len(codes) == nlay:
                        break
                    field = line[start : start + 2].replace(" ", "")
                    codes.append(int(field) if field else 0)

        codes = np.asarray(codes[:nlay], dtype=np.int32)
        intercellt = codes // 10
        laycon = codes - 10 * intercellt
        return laycon, intercellt

    @classmethod
    def _load_layer_arrays(
        cls,
        f_obj,
        model,
        nlay,
        ext_unit_dict,
        transient,
        laycon,
        ikvflag,
        ikcflag,
        iwdflg,
        tabrich=False,
        bubblept=False,
    ):
        """
        Loads BCF layer property arrays.

        Parameters
        ----------
        f_obj : open file object.
        model : model object
            The model object (of type :class:`flopy.modflow.mf.Modflow`) to
            which this package will be added.
        nlay : int. Number of model layers
        ext_unit_dict : external unit dictionary
        transient : bool. Flag for transient vs steady state models
        laycon : effective LAYCON array (after the RICHARDS promotion)
        ikvflag : int
            flag indicating if vertical hydraulic conductivity is input
            instead of leakance between two layers.
        ikcflag : int
            flag indicating if hydraulic conductivity or transmissivity
            information is input for each of the nodes or whether this information
            is directly input for the nodal connections.
        iwdflg : int
            flag to indicate if wetting is inactive (0) or not (non zero)
        tabrich : bool. TABRICH option (no per-layer Richards arrays)
        bubblept : bool. BUBBLEPT option (bubbling point heads)

        Returns
        -------
        sf1, tran, hy, vcont, sf2, wetdry, kv : layer property arrays
        richards_arrays : dict of alpha, beta, sr, brook and bp lists
        """
        sf1 = [0] * nlay
        tran = [0] * nlay
        hy = [0] * nlay
        if nlay > 1:
            vcont = [0] * (nlay - 1)
        else:
            vcont = [0] * nlay
        sf2 = [0] * nlay
        wetdry = [0] * nlay
        kv = [0] * nlay  # mfusg
        richards_arrays = {name: [0] * nlay for name in _RICHARDS_NAMES}

        for layer in range(nlay):
            util2d_shape = get_util2d_shape_for_layer(model, layer=layer)

            # sf1
            if transient:
                if model.verbose:
                    print(f"   loading sf1 layer {layer + 1:3d}...")
                sf1[layer] = Util2d.load(
                    f_obj, model, util2d_shape, np.float32, "sf1", ext_unit_dict
                )

            # hy/tran, and kv/vcont
            if ikcflag == 0:
                (
                    hy[layer],
                    tran[layer],
                    kv[layer],
                    vcont_k,
                ) = cls._load_hy_tran_kv_vcont(
                    f_obj,
                    model,
                    (layer, laycon[layer]),
                    ext_unit_dict,
                    ikvflag,
                )
                if layer < nlay - 1:
                    vcont[layer] = vcont_k

            # sf2
            if transient and (laycon[layer] in {2, 3, 4, 5}):
                if model.verbose:
                    print(f"   loading sf2 layer {layer + 1:3d}...")
                sf2[layer] = Util2d.load(
                    f_obj, model, util2d_shape, np.float32, "sf2", ext_unit_dict
                )

            # wetdry
            if (iwdflg != 0) and (laycon[layer] in {1, 3}):
                if model.verbose:
                    print(f"   loading wetdry layer {layer + 1:3d}...")
                wetdry[layer] = Util2d.load(
                    f_obj, model, util2d_shape, np.float32, "wetdry", ext_unit_dict
                )

            # Richards equation (skipped under TABRICH: RETCRVS replaces these)
            if laycon[layer] == 5 and not tabrich:
                names = _RICHARDS_NAMES if bubblept else _RICHARDS_NAMES[:4]
                for name in names:
                    if model.verbose:
                        print(f"   loading {name} layer {layer + 1:3d}...")
                    richards_arrays[name][layer] = Util2d.load(
                        f_obj, model, util2d_shape, np.float32, name, ext_unit_dict
                    )

        return sf1, tran, hy, vcont, sf2, wetdry, kv, richards_arrays

    @staticmethod
    def _load_hy_tran_kv_vcont(f_obj, model, laycon_k, ext_unit_dict, ikvflag):
        """
        Loads hy/tran and kv/vcont file entries.

        Parameters
        ----------
        f_obj : open file object.
        model : model object
            The model object (of type :class:`flopy.modflow.mf.Modflow`) to
            which this package will be added.
        laycon_k : tuple of two ints: (k, effective laycon for layer k)
        ext_unit_dict : external unit dictionary
        ikvflag : int
            flag indicating if vertical hydraulic conductivity is input
            instead of leakance between two layers.

        Returns
        -------
        _hy : Numpy array of hk values (or 0)
        _tran : Numpy array of tran values (or 0)
        _kv : Numpy array of kv values (or 0)
        _vcont : Numpy array of vcont values (or 0)
        """
        layer = laycon_k[0]
        laycon_k = laycon_k[1]
        util2d_shape = get_util2d_shape_for_layer(model, layer=layer)

        # hy or tran
        tran = 0
        hy = 0
        if laycon_k in {0, 2}:
            if model.verbose:
                print(f"   loading tran layer {layer + 1:3d}...")
            tran = Util2d.load(
                f_obj, model, util2d_shape, np.float32, "tran", ext_unit_dict
            )
        else:
            if model.verbose:
                print(f"   loading hy layer {layer + 1:3d}...")
            hy = Util2d.load(
                f_obj, model, util2d_shape, np.float32, "hy", ext_unit_dict
            )

        # vcont (IKVFLAG=0, all but the last layer) or kv (IKVFLAG!=0, every
        # layer of a multi-layer model)
        kv = 0
        vcont = 0
        if ikvflag == 0:
            if layer < (model.nlay - 1):
                if model.verbose:
                    print(f"   loading vcont layer {layer + 1:3d}...")
                vcont = Util2d.load(
                    f_obj, model, util2d_shape, np.float32, "vcont", ext_unit_dict
                )
        elif model.nlay > 1:
            if model.verbose:
                print(f"   loading kv layer {layer + 1:3d}...")
            kv = Util2d.load(
                f_obj, model, util2d_shape, np.float32, "kv", ext_unit_dict
            )

        return hy, tran, kv, vcont
