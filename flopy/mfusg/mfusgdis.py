"""
mfdis module.  Contains the ModflowDis class. Note that the user can access
the ModflowDis class as `flopy.modflow.ModflowDis`.

Additional information for this MODFLOW package can be found at the `Online
MODFLOW Guide
<https://water.usgs.gov/ogw/modflow/MODFLOW-2005-Guide/dis.html>`_.

"""

import re
import warnings
from decimal import Decimal, localcontext

import numpy as np

from ..pakbase import Package
from ..utils import Util2d, Util3d
from ..utils.crs import get_crs
from ..utils.flopy_io import line_parse

ITMUNI = {"u": 0, "s": 1, "m": 2, "h": 3, "d": 4, "y": 5}
LENUNI = {"u": 0, "f": 1, "m": 2, "c": 3}


class MfUsgDis(Package):
    """
    MODFLOW Discretization Package Class.

    Parameters
    ----------
    model : model object
        The model object (of type :class:`flopy.modflow.Modflow`) to which
        this package will be added.
    nlay : int
        Number of model layers (the default is 1).
    nrow : int
        Number of model rows (the default is 2).
    ncol : int
        Number of model columns (the default is 2).
    nper : int
        Number of model stress periods (the default is 1).
    delr : float or array of floats (ncol), optional
        An array of spacings along a row (the default is 1.0).
    delc : float or array of floats (nrow), optional
        An array of spacings along a column (the default is 0.0).
    laycbd : int or array of ints (nlay), optional
        An array of flags indicating whether or not a layer has a Quasi-3D
        confining bed below it. 0 indicates no confining bed, and not zero
        indicates a confining bed. LAYCBD for the bottom layer must be 0. (the
        default is 0)
    top : float or array of floats (nrow, ncol), optional
        An array of the top elevation of layer 1. For the common situation in
        which the top layer represents a water-table aquifer, it may be
        reasonable to set Top equal to land-surface elevation (the default is
        1.0)
    botm : float or array of floats (nlay, nrow, ncol), optional
        An array of the bottom elevation for each model cell (the default is
        0.)
    perlen : float or array of floats (nper)
        An array of the stress period lengths.
    nstp : int or array of ints (nper)
        Number of time steps in each stress period (default is 1).
    tsmult : float or array of floats (nper)
        Time step multiplier (default is 1.0).
    steady : bool or array of bool (nper)
        true or False indicating whether or not stress period is steady state
        (default is True). A TRTOSS period (see ``trtoss``) has False.
    Trtoss : bool
        Kept for compatibility; it has no effect (a warning is issued when it
        is set). Use ``trtoss``.
    itmuni : int
        Time units, default is days (4)
    lenuni : int
        Length units, default is meters (2)
    extension : string
        Filename extension (default is 'dis')
    unitnumber : int
        File unit number (default is None).
    filenames : str or list of str
        Filenames to use for the package. If filenames=None the package name
        will be created using the model name and package extension. If a
        single string is passed the package will be set to the string.
        Default is None.
    xul : float
        x coordinate of upper left corner of the grid, default is None, which
        means xul will be set to zero.
    yul : float
        y coordinate of upper-left corner of the grid, default is None, which
        means yul will be calculated as the sum of the delc array.  This
        default, combined with the xul and rotation defaults will place the
        lower-left corner of the grid at (0, 0).
    rotation : float
        counter-clockwise rotation (in degrees) of the grid about the lower-
        left corner. default is 0.0
    crs : pyproj.CRS, int, str, optional if `prjfile` is specified
        Coordinate reference system (CRS) for the model grid
        (must be projected; geographic CRS are not supported).
        The value can be anything accepted by
        :meth:`pyproj.CRS.from_user_input() <pyproj.crs.CRS.from_user_input>`,
        such as an authority string (eg "EPSG:26916") or a WKT string.
    prjfile : str or PathLike, optional if `crs` is specified
        ESRI-style projection file with well-known text defining the CRS
        for the model grid (must be projected; geographic CRS are not supported).
    start_datetime : str
        starting datetime of the simulation. default is '1/1/1970'
    **kwargs : dict, optional
        Support deprecated keyword options.

        .. deprecated:: 3.5
           ``proj4_str`` will be removed for FloPy 3.6, use ``crs`` instead.
    trtoss : None, float or sequence of float/None (nper), optional
        Storage fraction (STOFRAC) of the USG-T ``TRTOSS`` stress periods: a
        period with a value is run as transient until the storage rate is
        below this fraction of the total budget rate, then ends. None or NaN
        means the period is not TRTOSS (SS or TR from ``steady``). A scalar
        applies to every period, a sequence must have ``nper`` values. A
        TRTOSS period needs ``steady=False`` for that period (the default
        ``steady=True`` must be changed). Default is None (no TRTOSS).
        Stored as ``self.trtoss``, a list (float, or None for other periods).

    Attributes
    ----------
    heading : str
        Text string written to top of package input file.

    Methods
    -------

    See Also
    --------

    Notes
    -----
    USG-T ``TRTOSS`` periods (``trtoss``): the period ends early only with
    adaptive time stepping (OC ``ATS``/``ATSA``); without it the period runs
    like ``TR``. A fraction <= 0 never ends the period early (it is kept as
    given). The value is written so that USG-T reads the same REAL number
    (single precision) as ``numpy.float32`` of the value; a value loaded
    from a file keeps the REAL number USG-T reads from that file. SS and
    TRTOSS periods without any TR period are rejected: USG-T 2.8.1 then
    reads no storage arrays (ITRSS=0) although the TRTOSS period is
    transient, while USG-T 2.7.0 reads them. ``period_type`` and ``itrss``
    are read-only views; ``itrss`` follows the USG-T 2.8.1 rule.

    Examples
    --------

    >>> import flopy
    >>> m = flopy.modflow.Modflow()
    >>> dis = flopy.modflow.ModflowDis(m)

    """

    def __init__(
        self,
        model,
        nlay=1,
        nrow=2,
        ncol=2,
        nper=1,
        delr=1.0,
        delc=1.0,
        laycbd=0,
        top=1,
        botm=0,
        perlen=1,
        nstp=1,
        tsmult=1,
        steady=True,
        Trtoss=False,
        itmuni=4,
        lenuni=2,
        extension="dis",
        unitnumber=None,
        filenames=None,
        xul=None,
        yul=None,
        rotation=None,
        crs=None,
        prjfile=None,
        start_datetime=None,
        trtoss=None,
        **kwargs,
    ):
        # set default unit number of one is not specified
        if unitnumber is None:
            unitnumber = MfUsgDis._defaultunit()

        # call base package constructor
        super().__init__(
            model,
            extension=extension,
            name=self._ftype(),
            unit_number=unitnumber,
            filenames=self._prepare_filenames(filenames),
        )

        self.url = "dis.html"
        self.nrow = nrow
        self.ncol = ncol
        self.nlay = nlay
        self.nper = nper

        # initialize botm to an appropriate sized
        if nlay > 1:
            if isinstance(botm, float) or isinstance(botm, int):
                botm = np.linspace(top, botm, nlay)

        # Set values of all parameters
        self._generate_heading()
        self.laycbd = Util2d(model, (self.nlay,), np.int32, laycbd, name="laycbd")
        self.laycbd[-1] = 0  # bottom layer must be zero
        self.delr = Util2d(
            model,
            (self.ncol,),
            np.float32,
            delr,
            name="delr",
            locat=self.unit_number[0],
        )
        self.delc = Util2d(
            model,
            (self.nrow,),
            np.float32,
            delc,
            name="delc",
            locat=self.unit_number[0],
        )
        self.top = Util2d(
            model,
            (self.nrow, self.ncol),
            np.float32,
            top,
            name="model_top",
            locat=self.unit_number[0],
        )
        self.botm = Util3d(
            model,
            (self.nlay + sum(self.laycbd), self.nrow, self.ncol),
            np.float32,
            botm,
            "botm",
            locat=self.unit_number[0],
        )
        self.perlen = Util2d(model, (self.nper,), np.float32, perlen, name="perlen")
        self.nstp = Util2d(model, (self.nper,), np.int32, nstp, name="nstp")
        self.tsmult = Util2d(model, (self.nper,), np.float32, tsmult, name="tsmult")
        self.steady = Util2d(model, (self.nper,), bool, steady, name="steady")
        trtoss = _trtoss_array(trtoss, self.nper)
        _check_period_types(self.steady.array, trtoss, self.nper)
        self.trtoss = _trtoss_list(trtoss)
        if Trtoss:
            warnings.warn(
                "MfUsgDis: Trtoss has no effect; use trtoss (storage fraction "
                "per stress period)"
            )

        try:
            self.itmuni = int(itmuni)
        except:
            self.itmuni = ITMUNI[itmuni.lower()[0]]
        try:
            self.lenuni = int(lenuni)
        except:
            self.lenuni = LENUNI[lenuni.lower()[0]]

        self.parent.add_package(self)
        self.itmuni_dict = {
            0: "undefined",
            1: "seconds",
            2: "minutes",
            3: "hours",
            4: "days",
            5: "years",
        }

        # set the model grid coordinate info
        mg = model.modelgrid
        if rotation is None:
            rotation = model._rotation
        if rotation is not None:
            # set rotation before anything else
            mg.set_coord_info(angrot=rotation)
        set_coord_info_args = {"crs": crs, "prjfile": prjfile}
        if "proj4_str" in kwargs:
            warnings.warn(
                "the proj4_str argument will be deprecated and will be "
                "removed in version 3.6. Use crs instead.",
                PendingDeprecationWarning,
            )
            proj4_str = kwargs.pop("proj4_str")
            if crs is None:
                set_coord_info_args["crs"] = proj4_str
        if xul is None:
            xul = model._xul
        if xul is not None:
            set_coord_info_args["xoff"] = mg._xul_to_xll(xul)
        if yul is None:
            yul = model._yul
        if yul is not None:
            set_coord_info_args["yoff"] = mg._yul_to_yll(yul)

        # set all other relevant coordinate properties
        model._modelgrid.set_coord_info(**set_coord_info_args)

        if start_datetime is None:
            start_datetime = model._start_datetime

        self.start_datetime = start_datetime
        self._totim = None

    @property
    def period_type(self):
        """Stress-period types, "SS", "TR" or "TRTOSS" (read-only view)."""
        trtoss = _trtoss_array(self.trtoss, self.nper, exact=False)
        steady = _check_period_types(self.steady.array, trtoss, self.nper, exact=False)
        return _period_flags(steady, trtoss)

    @property
    def itrss(self):
        """ITRSS as USG-T 2.8.1 sets it: 1 with TR and no SS periods, 0 with
        SS and no TR periods, -1 otherwise (TRTOSS counts as neither; USG-T
        2.7.0 counts TRTOSS as TR)."""
        return _itrss(self.period_type)

    def checklayerthickness(self):
        """
        Check layer thickness.

        """
        return (self.parent.modelgrid.cell_thickness > 0).all()

    def get_totim(self, use_cached=False):
        """
        Get the totim at the end of each time step

        Parameters
        ----------
        use_cached : bool
            method to use cached totim values instead of calculating totim
            dynamically


        Returns
        -------
        totim: numpy array
            numpy array with simulation totim at the end of each time step

        """
        if not use_cached or self._totim is None:
            totim = []
            nstp = self.nstp.array
            perlen = self.perlen.array
            tsmult = self.tsmult.array
            t = 0.0
            for kper in range(self.nper):
                m = tsmult[kper]
                p = float(nstp[kper])
                dt = perlen[kper]
                if m > 1:
                    dt *= (m - 1.0) / (m**p - 1.0)
                else:
                    dt = dt / p
                for kstp in range(nstp[kper]):
                    t += dt
                    totim.append(t)
                    if m > 1:
                        dt *= m
            self._totim = np.array(totim, dtype=float)

        return self._totim

    def get_final_totim(self):
        """
        Get the totim at the end of the simulation

        Returns
        -------
        totim: float
            maximum simulation totim

        """
        return self.get_totim()[-1]

    def get_kstp_kper_toffset(self, t=0.0, use_cached_totim=False):
        """
        Get the stress period, time step, and time offset from passed time.

        Parameters
        ----------
        t : float
            totim to return the stress period, time step, and toffset for
            based on time discretization data. Default is 0.
        use_cached_totim : bool
            optional flag to use a cached calculation of totim, vs. dynamically
            calculating totim. Setting to True significantly speeds up looped
            operations that call this function (default is False).

        Returns
        -------
        kstp : int
            time step in stress period corresponding to passed totim
        kper : int
            stress period corresponding to passed totim
        toffset : float
            time offset of passed totim from the beginning of kper

        """

        if t < 0.0:
            t = 0.0
        totim = self.get_totim(use_cached_totim)
        nstp = self.nstp.array
        ipos = 0
        t0 = 0.0
        kper = self.nper - 1
        kstp = nstp[-1] - 1
        toffset = self.perlen.array[-1]
        done = False
        for iper in range(self.nper):
            tp0 = t0
            for istp in range(nstp[iper]):
                t1 = totim[ipos]
                if t >= t0 and t < t1:
                    done = True
                    kper = iper
                    kstp = istp
                    toffset = t - tp0
                    break
                ipos += 1
                t0 = t1
            if done:
                break
        return kstp, kper, toffset

    def get_totim_from_kper_toffset(self, kper=0, toffset=0.0, use_cached_totim=False):
        """
        Get totim from a passed kper and time offset from the beginning
        of a stress period

        Parameters
        ----------
        kper : int
            stress period. Default is 0
        toffset : float
            time offset relative to the beginning of kper
        use_cached_totim : bool
            optional flag to use a cached calculation of totim, vs. dynamically
            calculating totim. Setting to True significantly speeds up looped
            operations that call this function (default is False).

        Returns
        -------
        t : float
            totim to return the stress period, time step, and toffset for
            based on time discretization data. Default is 0.

        """

        if kper < 0:
            kper = 0.0
        if kper >= self.nper:
            raise ValueError(f"kper ({kper}) must be less than to nper ({self.nper}).")

        totim = self.get_totim(use_cached_totim)
        nstp = self.nstp.array
        ipos = 0
        t0 = 0.0
        tp0 = 0.0
        for iper in range(kper + 1):
            tp0 = t0
            if iper == kper:
                break
            for istp in range(nstp[iper]):
                t1 = totim[ipos]
                ipos += 1
                t0 = t1
        t = tp0 + toffset
        return t

    def get_cell_volumes(self):
        """
        Get an array of cell volumes.

        Returns
        -------
        vol : array of floats (nlay, nrow, ncol)

        """
        vol = np.empty((self.nlay, self.nrow, self.ncol))
        for l in range(self.nlay):
            vol[l, :, :] = self.parent.modelgrid.cell_thickness[l]
        for r in range(self.nrow):
            vol[:, r, :] *= self.delc[r]
        for c in range(self.ncol):
            vol[:, :, c] *= self.delr[c]
        return vol

    @property
    def zcentroids(self):
        z = np.empty((self.nlay, self.nrow, self.ncol))
        z[0, :, :] = (self.top[:, :] + self.botm[0, :, :]) / 2.0

        for l in range(1, self.nlay):
            z[l, :, :] = (self.botm[l - 1, :, :] + self.botm[l, :, :]) / 2.0
        return z

    def get_node_coordinates(self):
        """
        Get y, x, and z cell centroids in local model coordinates.

        Returns
        -------
        y : list of cell y-centroids

        x : list of cell x-centroids

        z : array of floats (nlay, nrow, ncol)

        """

        delr = self.delr.array
        delc = self.delc.array

        # In row direction
        Ly = np.add.reduce(delc)
        y = Ly - (np.add.accumulate(self.delc) - 0.5 * delc)

        # In column direction
        x = np.add.accumulate(self.delr) - 0.5 * delr

        # In layer direction
        z = self.zcentroids

        return y, x, z

    def get_rc_from_node_coordinates(self, x, y, local=True):
        """
        Get the row and column of a point or sequence of points
        in model coordinates.

        Parameters
        ----------
        x : float or sequence of floats
            x coordinate(s) of points to find in model grid
        y : float or sequence floats
            y coordinate(s) of points to find in model grid
        local : bool
          x and y coordinates are in model local coordinates.  If false, then
          x and y are in world coordinates. (default is True)

        Returns
        -------
        r : row or sequence of rows (zero-based)
        c : column or sequence of columns (zero-based)

        """
        mg = self.parent.modelgrid
        if np.isscalar(x):
            r, c = mg.intersect(x, y, local=local)
        else:
            r = []
            c = []
            for xx, yy in zip(x, y):
                rr, cc = mg.intersect(xx, yy, local=local)
                r.append(rr)
                c.append(cc)
        return r, c

    def get_lrc(self, nodes):
        """
        Get zero-based layer, row, column from a list of zero-based
        MODFLOW node numbers.

        Returns
        -------
        v : list of tuples containing the layer (k), row (i),
            and column (j) for each node in the input list
        """
        return self.parent.modelgrid.get_lrc(nodes)

    def get_node(self, lrc_list):
        """
        Get zero-based node number from a list of zero-based MODFLOW
        layer, row, column tuples.

        Returns
        -------
        v : list of MODFLOW nodes for each layer (k), row (i),
            and column (j) tuple in the input list
        """
        return self.parent.modelgrid.get_node(lrc_list)

    def get_layer(self, i, j, elev):
        """Return the layer for an elevation at an i, j location.

        Parameters
        ----------
        i : row index (zero-based)
        j : column index
        elev : elevation (in same units as model)

        Returns
        -------
        k : zero-based layer index
        """
        return get_layer(self, i, j, elev)

    def gettop(self):
        """
        Get the top array.

        Returns
        -------
        top : array of floats (nrow, ncol)
        """
        return self.top.array

    def getbotm(self, k=None):
        """
        Get the bottom array.

        Returns
        -------
        botm : array of floats (nlay, nrow, ncol), or

        botm : array of floats (nrow, ncol) if k is not none
        """
        if k is None:
            return self.botm.array
        else:
            return self.botm.array[k, :, :]

    def write_file(self, check=True):
        """
        Write the package file.

        Parameters
        ----------
        check : bool
            Check package data for common errors. (default True)

        Returns
        -------
        None

        """
        # checked before the file is opened (steady/trtoss may have changed)
        trtoss = _trtoss_array(self.trtoss, self.nper, exact=False)
        steady = _check_period_types(self.steady.array, trtoss, self.nper, exact=False)
        flags = _period_flags(steady, trtoss)
        if check:  # allows turning off package checks when writing files at model level
            self.check(f=f"{self.name[0]}.chk", verbose=self.parent.verbose, level=1)
        # Open file for writing
        f_dis = open(self.fn_path, "w")
        # Item 0: heading
        f_dis.write(f"{self.heading}\n")
        # Item 1: NLAY, NROW, NCOL, NPER, ITMUNI, LENUNI
        f_dis.write(
            "{:10d}{:10d}{:10d}{:10d}{:10d}{:10d}\n".format(
                self.nlay, self.nrow, self.ncol, self.nper, self.itmuni, self.lenuni
            )
        )
        # Item 2: LAYCBD
        for l in range(0, self.nlay):
            f_dis.write(f"{self.laycbd[l]:3d}")
        f_dis.write("\n")
        # Item 3: DELR
        f_dis.write(self.delr.get_file_entry())
        # Item 4: DELC
        f_dis.write(self.delc.get_file_entry())
        # Item 5: Top(NCOL, NROW)
        f_dis.write(self.top.get_file_entry())
        # Item 5: BOTM(NCOL, NROW)
        f_dis.write(self.botm.get_file_entry())

        # Item 6: NPER, NSTP, TSMULT, Ss/tr/Trtoss [stofrac]
        for t in range(self.nper):
            f_dis.write(
                _period_record(
                    self.perlen[t], self.nstp[t], self.tsmult[t], flags[t], trtoss[t]
                )
            )
        f_dis.close()

    def check(self, f=None, verbose=True, level=1, checktype=None):
        """
        Check dis package data for zero and negative thicknesses.

        Parameters
        ----------
        f : str or file handle
            String defining file name or file handle for summary file
            of check method output. If a string is passed a file handle
            is created. If f is None, check method does not write
            results to a summary file. (default is None)
        verbose : bool
            Boolean flag used to determine if check method results are
            written to the screen
        level : int
            Check method analysis level. If level=0, summary checks are
            performed. If level=1, full checks are performed.

        Returns
        -------
        None

        Examples
        --------

        >>> import flopy
        >>> m = flopy.modflow.Modflow.load('model.nam')
        >>> m.dis.check()
        """
        chk = self._get_check(f, verbose, level, checktype)

        # make ibound of same shape as thicknesses/botm for quasi-3D models
        active = chk.get_active(include_cbd=True)

        # Use either a numpy array or masked array
        thickness = self.parent.modelgrid.cell_thickness
        non_finite = ~(np.isfinite(thickness))
        if non_finite.any():
            thickness[non_finite] = 0
            thickness = np.ma.array(thickness, mask=non_finite)

        chk.values(
            thickness, active & (thickness <= 0), "zero or negative thickness", "Error"
        )
        thin_cells = (thickness < chk.thin_cell_threshold) & (thickness > 0)
        chk.values(
            thickness,
            active & thin_cells,
            "thin cells (less than checker threshold of {:.1f})".format(
                chk.thin_cell_threshold
            ),
            "Error",
        )
        chk.values(
            self.top.array,
            active[0, :, :] & np.isnan(self.top.array),
            "nan values in top array",
            "Error",
        )
        chk.values(
            self.botm.array,
            active & np.isnan(self.botm.array),
            "nan values in bottom array",
            "Error",
        )
        chk.summarize()
        return chk

    @classmethod
    def load(cls, f, model, ext_unit_dict=None, check=True):
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
        check : bool
            Check package data for common errors. (default True)

        Returns
        -------
        dis : ModflowDis object
            ModflowDis object.

        Examples
        --------

        >>> import flopy
        >>> m = flopy.modflow.Modflow()
        >>> dis = flopy.modflow.ModflowDis.load('test.dis', m)

        """

        if model.verbose:
            print("loading dis package file...")

        openfile = not hasattr(f, "read")
        if openfile:
            filename = f
            f = open(filename, "r")

        # dataset 0 -- header
        while True:
            line = f.readline()
            if line[0] != "#":
                break

        # dataset 1
        nlay, nrow, ncol, nper, itmuni, lenuni = line.strip().split()[0:6]
        nlay = int(nlay)
        nrow = int(nrow)
        ncol = int(ncol)
        nper = int(nper)
        itmuni = int(itmuni)
        lenuni = int(lenuni)
        # dataset 2 -- laycbd
        if model.verbose:
            print(
                "   Loading dis package with:\n      {} layers, {} rows, {} "
                "columns, and {} stress periods".format(nlay, nrow, ncol, nper)
            )
            print("   loading laycbd...")
        laycbd = np.zeros(nlay, dtype=int)
        d = 0
        while True:
            line = f.readline()
            raw = line.strip("\n").split()
            for val in raw:
                if int(val) != 0:
                    laycbd[d] = 1
                d += 1
                if d == nlay:
                    break
            if d == nlay:
                break
        # dataset 3 -- delr
        if model.verbose:
            print("   loading delr...")
        delr = Util2d.load(f, model, (ncol,), np.float32, "delr", ext_unit_dict)
        # dataset 4 -- delc
        if model.verbose:
            print("   loading delc...")
        delc = Util2d.load(f, model, (nrow,), np.float32, "delc", ext_unit_dict)
        # dataset 5 -- top
        if model.verbose:
            print("   loading top...")
        top = Util2d.load(f, model, (nrow, ncol), np.float32, "top", ext_unit_dict)
        # dataset 6 -- botm
        ncbd = laycbd.sum()
        if model.verbose:
            print("   loading botm...")
            print(f"      for {nlay} layers and {ncbd} confining beds")
        if nlay > 1:
            botm = Util3d.load(
                f,
                model,
                (nlay + ncbd, nrow, ncol),
                np.float32,
                "botm",
                ext_unit_dict,
            )
        else:
            botm = Util3d.load(
                f, model, (nlay, nrow, ncol), np.float32, "botm", ext_unit_dict
            )
        # dataset 7 -- stress period info
        if model.verbose:
            print("   loading stress period data...")
            print(f"       for {nper} stress periods")
        perlen = []
        nstp = []
        tsmult = []
        steady = []
        trtoss = []
        for k in range(nper):
            line = f.readline()
            words = line_parse(line)
            a1 = float(words[0])
            a2 = int(words[1])
            a3 = float(words[2])
            a4, stofrac = _parse_period_flag(words, k)
            trtoss.append(stofrac)
            perlen.append(a1)
            nstp.append(a2)
            tsmult.append(a3)
            steady.append(a4)

        if openfile:
            f.close()

        # set package unit number
        unitnumber = None
        filenames = [None]
        if ext_unit_dict is not None:
            unitnumber, filenames[0] = model.get_ext_dict_attr(
                ext_unit_dict, filetype=MfUsgDis._ftype()
            )

        # create dis object instance
        dis = cls(
            model,
            nlay=nlay,
            nrow=nrow,
            ncol=ncol,
            nper=nper,
            delr=delr,
            delc=delc,
            laycbd=laycbd,
            top=top,
            botm=botm,
            perlen=perlen,
            nstp=nstp,
            tsmult=tsmult,
            steady=steady,
            itmuni=itmuni,
            lenuni=lenuni,
            trtoss=trtoss,
            unitnumber=unitnumber,
            filenames=filenames,
        )
        if check:
            dis.check(f=f"{dis.name[0]}.chk", verbose=dis.parent.verbose, level=0)
        # return dis object instance
        return dis

    @staticmethod
    def _ftype():
        return "DIS"

    @staticmethod
    def _defaultunit():
        return 11


# USG-T TRTOSS stress periods (SGLO2BAS8ARDIS, glo2basu1.f): the flag is read
# with URWORD and the storage fraction as a REAL (F30.0 on a word of at most
# 30 characters, correctly rounded to single precision).
_REAL_WORD = re.compile(
    r"^([+-]?(?:\d+\.?\d*|\.\d+))(?:[EeDdQq]([+-]?\d+)|([+-]\d+))?$"
)


def _real32(value):
    """np.float32 of value, None when it is beyond the REAL range."""
    with np.errstate(over="ignore"):
        r = np.float32(value)
    return r if np.isfinite(r) else None


def _stofrac_from_word(word, kper):
    """Storage fraction of a TRTOSS period as USG-T reads ``word``."""
    m = _REAL_WORD.match(word) if len(word) <= 30 else None
    if m is None:
        raise ValueError(
            f"stress period {kper + 1}: TRTOSS storage fraction {word!r} is not "
            "a finite number USG-T can read"
        )
    mantissa, e1, e2 = m.groups()
    with localcontext() as ctx:
        ctx.prec = 400
        exact = Decimal(mantissa).scaleb(int(e1 or e2 or 0))
        near = _real32(float(exact))
        if near is None:
            raise ValueError(
                f"stress period {kper + 1}: TRTOSS storage fraction {word!r} is "
                "beyond the REAL range"
            )
        # correctly rounded single precision value (ties to even)
        cands = [
            c
            for c in (
                np.nextafter(near, -near.dtype.type(np.inf)),
                near,
                np.nextafter(near, near.dtype.type(np.inf)),
            )
            if np.isfinite(c)
        ]
        real = min(
            cands,
            key=lambda c: (
                abs(Decimal(float(c)) - exact),
                int(np.array(c).view(np.int32)) & 1,
            ),
        )
    value = float(exact)
    return value if np.float32(value) == real else float(real)


def _stofrac_word(value):
    """Shortest text that USG-T reads as np.float32(value)."""
    return np.format_float_scientific(np.float32(value), unique=True, exp_digits=2)


def _trtoss_array(trtoss, nper, exact=True):
    """Validated storage fractions of the first nper periods; NaN = not a
    TRTOSS period. exact=False accepts a longer sequence (nper reduced after
    construction: only the first nper periods are written)."""
    if trtoss is None:
        return np.full(nper, np.nan)
    values = np.asarray(trtoss, dtype=object)
    if values.ndim == 0:
        values = np.full(nper, values.item(), dtype=object)
    elif values.ndim != 1 or values.size < nper or (exact and values.size != nper):
        raise ValueError(
            f"trtoss must be a scalar or a sequence of nper={nper} values, "
            f"got shape {values.shape}"
        )
    values = values[:nper]
    out = np.full(nper, np.nan)
    for k, v in enumerate(values):
        if v is None:
            continue
        if isinstance(v, (bool, np.bool_, str, bytes)):
            raise ValueError(f"trtoss[{k}] must be a number or None, got {v!r}")
        try:
            x = float(v)
        except (TypeError, ValueError):
            raise ValueError(f"trtoss[{k}] must be a number or None, got {v!r}")
        if np.isinf(x) or (not np.isnan(x) and _real32(x) is None):
            raise ValueError(f"trtoss[{k}]={x} is not a finite REAL value")
        out[k] = x
    return out


def _trtoss_list(trtoss):
    """Stored form of the fractions: float, or None for a non-TRTOSS period."""
    return [None if np.isnan(x) else float(x) for x in trtoss]


def _check_period_types(steady, trtoss, nper, exact=True):
    """Shapes, TRTOSS only on steady=False periods, and no SS + TRTOSS
    without TR (USG-T 2.8.1 ITRSS=0: storage arrays not read although the
    TRTOSS period is transient; USG-T 2.7.0 reads them). Returns the steady
    flags of the first nper periods."""
    steady = np.asarray(steady)
    if steady.ndim != 1 or steady.size < nper or (exact and steady.size != nper):
        raise ValueError(f"steady must have nper={nper} values, got {steady.shape}")
    steady = steady[:nper].astype(bool)
    is_trtoss = ~np.isnan(trtoss)
    bad = is_trtoss & steady
    if bad.any():
        raise ValueError(
            f"TRTOSS stress period(s) {(np.nonzero(bad)[0] + 1).tolist()} need "
            "steady=False"
        )
    if is_trtoss.any() and steady.any() and not (~steady & ~is_trtoss).any():
        raise ValueError(
            "SS and TRTOSS stress periods without a TR period are not supported: "
            "USG-T 2.8.1 then reads no storage arrays (ITRSS=0) although the "
            "TRTOSS period is transient, and USG-T 2.7.0 reads them"
        )
    return steady


def _period_flags(steady, trtoss):
    """SS / TR / TRTOSS per period."""
    return [
        "TRTOSS" if not np.isnan(f) else ("SS" if s else "TR")
        for s, f in zip(np.asarray(steady, dtype=bool), trtoss)
    ]


def _itrss(flags):
    """ITRSS as USG-T 2.8.1 sets it (TRTOSS counts as neither SS nor TR)."""
    if "TR" in flags and "SS" not in flags:
        return 1
    if "SS" in flags and "TR" not in flags:
        return 0
    return -1


def _parse_period_flag(words, kper):
    """(steady, stofrac) from the words of a stress-period record."""
    flag = words[3].upper() if len(words) > 3 else ""
    if flag == "SS":
        return True, np.nan
    if flag == "TR":
        return False, np.nan
    if flag == "TRTOSS":
        if len(words) > 4:
            return False, _stofrac_from_word(words[4], kper)
        warnings.warn(
            f"stress period {kper + 1}: TRTOSS without a storage fraction; "
            "USG-T reads 0.0 (the period does not end early), written as 0.0"
        )
        return False, 0.0
    raise ValueError(
        f"stress period {kper + 1}: flag {flag!r} must be SS, TR or TRTOSS"
    )


def _period_record(perlen, nstp, tsmult, flag, stofrac):
    text = f"{perlen:14f}{nstp:14d}{tsmult:10f}  {flag}"
    if flag == "TRTOSS":
        text += f" {_stofrac_word(stofrac)}"
    return text + "\n"


def get_layer(dis, i, j, elev):
    """Return the layers for elevations at i, j locations.

    Parameters
    ----------
    dis : flopy.modflow.ModflowDis object
    i : scaler or sequence
        row index (zero-based)
    j : scaler or sequence
        column index
    elev : scaler or sequence
        elevation (in same units as model)

    Returns
    -------
    k : np.ndarray (1-D) or scalar
        zero-based layer index
    """

    def to_array(arg):
        if not isinstance(arg, np.ndarray):
            return np.array([arg])
        else:
            return arg

    i = to_array(i)
    j = to_array(j)
    elev = to_array(elev)
    botms = dis.botm.array[:, i, j].tolist()
    layers = np.sum(((botms - elev) > 0), axis=0)
    # force elevations below model bottom into bottom layer
    layers[layers > dis.nlay - 1] = dis.nlay - 1
    layers = np.atleast_1d(np.squeeze(layers))
    if len(layers) == 1:
        layers = layers[0]
    return layers
