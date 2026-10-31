from collections.abc import Iterator
from math import sqrt

import numpy as np

from .cvfdutil import get_disv_gridprops
from .geometry import point_in_polygon
from .utl_import import import_optional_dependency


def get_sorted_vertices(
    icell_vertices, vertices, verbose=False
) -> Iterator[tuple[float, int]]:
    centroid = vertices[icell_vertices].mean(axis=0)
    tlist = []
    for i, iv in enumerate(icell_vertices):
        x, y = vertices[iv]
        dx = x - centroid[0]
        dy = y - centroid[1]
        tlist.append((np.arctan2(-dy, dx), iv))

    # convert to dictionary and back again to weed out duplicate angles,
    # which it's assumed indicate vertices with duplicate coordinates
    tlist = list(dict(tlist).items())
    tlist.sort()

    # weed out (near-)180-degree-angle vertices, which presumably should
    # occur only along the outer boundary of the domain
    for i, (angl1, iv1) in enumerate(tlist):
        (angl0, iv0) = tlist[(i - 1) % len(tlist)]
        (angl2, iv2) = tlist[(i + 1) % len(tlist)]
        x0, y0 = vertices[iv0]
        x1, y1 = vertices[iv1]
        x2, y2 = vertices[iv2]
        s0x = x0 - x1
        s0y = y0 - y1
        s0mag = sqrt(s0x * s0x + s0y * s0y)
        s2x = x2 - x1
        s2y = y2 - y1
        s2mag = sqrt(s2x * s2x + s2y * s2y)
        sinang = (s0x * s2y - s0y * s2x) / (s0mag * s2mag)

        # check for near-180 angle by checking for a small magnitude for
        # sine of the angle and a negative dot product for the two edge vectors
        # emanating from the vertex (to distinguish from a near-zero angle);
        # might be a more efficient way; tolerance for |sin(angle)| is
        # hardwired to arbitrary value of 0.00001
        if abs(sinang) < 0.00001 and s0x * s2x + s0y * s2y < 0:
            if verbose:
                print("Omitting ", x1, y1, sinang)
            continue

        yield iv1


def get_valid_faces(vor):
    nvalid_faces = np.zeros(vor.npoints, dtype=int)
    for pointidx, simplex in zip(vor.ridge_points, vor.ridge_vertices):
        simplex = np.asarray(simplex)
        if np.all(simplex >= 0):
            ip1, ip2 = pointidx
            nvalid_faces[ip1] += 1
            nvalid_faces[ip2] += 1
    return nvalid_faces


# todo: send this to point in polygon method defined in Rasters
def point_in_cell(point, vertices):
    shapely_geo = import_optional_dependency("shapely.geometry")

    p = shapely_geo.Point(point)
    poly = shapely_geo.Polygon(vertices)
    if p.intersects(poly):
        return True
    else:
        return False


# todo: find out how this is different from get_sorted_vertices()
def sort_vertices(vlist):
    x, y = zip(*vlist)
    x = np.array(x)
    y = np.array(y)
    xc = x.mean()
    yc = y.mean()
    tlist = []
    for i, (x, y) in enumerate(vlist):
        dx = x - xc
        dy = y - yc
        tlist.append((np.arctan2(-dy, dx), i))
    tlist.sort()
    return [vlist[i] for angl, i in tlist]


def tri2vor(tri, **kwargs):
    """
    This is the workhorse for the VoronoiGrid class for creating a voronoi
    grid from a constructed and built flopy Triangle grid.

    Parameters
    ----------
    tri : flopy.utils.Triangle
    Flopy triangle object is used to construct the complementary voronoi
    diagram.

    Returns
    -------
    verts, iverts : ndarray, list of lists

    """
    import_optional_dependency(
        "scipy.spatial",
        error_message="Voronoi requires SciPy.",
    )
    from scipy.spatial import Voronoi

    # assign local variables
    tri_verts = tri.verts
    tri_iverts = tri.iverts
    tri_edge = tri.edge
    npoints = tri_verts.shape[0]
    ntriangles = len(tri_iverts)
    nedges = tri_edge.shape[0]

    # check to make sure there are no duplicate points
    tri_verts_unique = np.unique(tri_verts, axis=0)
    if tri_verts.shape != tri_verts_unique.shape:
        npoints_unique = tri_verts_unique.shape[0]
        errmsg = (
            f"There are duplicate points in the triangular mesh. "
            f"These can be caused by overlapping regions, holes, and "
            f"refinement features.  The triangular mesh has {npoints} "
            f"points but only {npoints_unique} are unique."
        )
        raise Exception(errmsg)

    # construct the voronoi grid
    vor = Voronoi(tri_verts, **kwargs)
    ridge_points = vor.ridge_points
    ridge_vertices = vor.ridge_vertices

    # test the voronoi vertices, and mark those outside of the domain
    nvertices = vor.vertices.shape[0]
    xc = vor.vertices[:, 0].reshape((nvertices, 1))
    yc = vor.vertices[:, 1].reshape((nvertices, 1))
    domain_polygon = list(tri._polygons[0])
    vor_vert_indomain = point_in_polygon(xc, yc, domain_polygon)
    vor_vert_indomain = vor_vert_indomain.flatten()
    nholes = len(tri._holes)
    if nholes > 0:
        for ihole in range(nholes):
            ipolygon = ihole + 1
            polygon = list(tri._polygons[ipolygon])
            vor_vert_notindomain = point_in_polygon(xc, yc, polygon)
            vor_vert_notindomain = vor_vert_notindomain.flatten()
            idx = np.asarray(vor_vert_notindomain == True).nonzero()
            vor_vert_indomain[idx] = False

    idx_vertindex = -1 * np.ones((nvertices), int)
    idx_filtered = np.asarray(vor_vert_indomain == True).nonzero()
    nvalid_vertices = len(idx_filtered[0])
    # renumber valid vertices consecutively
    idx_vertindex[idx_filtered] = np.arange(nvalid_vertices)

    # Create new lists for the voronoi grid vertices and the
    # voronoi grid incidence list.  There should be one voronoi
    # cell for each vertex point in the triangular grid
    vor_verts = list(vor.vertices[idx_filtered])
    vor_iverts = [[] for i in range(npoints)]

    # step 1 -- go through voronoi ridge vertices
    # and add valid vertices to vor_verts and to the
    # vor_iverts incidence list
    for ips, irvs in zip(ridge_points, ridge_vertices):
        ip0, ip1 = ips
        irv0, irv1 = irvs
        if irv0 >= 0:
            point_in_domain = vor_vert_indomain[irv0]
            if point_in_domain:
                ivert = idx_vertindex[irv0]
                if ivert not in vor_iverts[ip0]:
                    vor_iverts[ip0].append(ivert)
                if ivert not in vor_iverts[ip1]:
                    vor_iverts[ip1].append(ivert)
        if irv1 >= 0:
            point_in_domain = vor_vert_indomain[irv1]
            if point_in_domain:
                ivert = idx_vertindex[irv1]
                if ivert not in vor_iverts[ip0]:
                    vor_iverts[ip0].append(ivert)
                if ivert not in vor_iverts[ip1]:
                    vor_iverts[ip1].append(ivert)

    # step 2 -- along the edge, add points
    # Count number of boundary markers that correspond to the outer
    # polygon domain or to holes.  These segments will be used to add
    # new vertices for edge cells.
    nexterior_boundary_markers = len(tri._polygons[0])
    for ihole in range(nholes):
        polygon = tri._polygons[ihole + 1]
        nexterior_boundary_markers += len(polygon)
    idx = (tri_edge["boundary_marker"] > 0) & (
        tri_edge["boundary_marker"] <= nexterior_boundary_markers
    )
    inewvert = len(vor_verts)
    for _, ip0, ip1, _ in tri_edge[idx]:
        midpoint = tri_verts[[ip0, ip1]].mean(axis=0)
        px, py = midpoint
        vor_verts.append((px, py))

        # add midpoint to each voronoi cell
        vor_iverts[ip0].append(inewvert)
        vor_iverts[ip1].append(inewvert)
        inewvert += 1

        # add ip0 triangle vertex to voronoi cell
        px, py = tri_verts[ip0]
        vor_verts.append((px, py))
        vor_iverts[ip0].append(inewvert)
        inewvert += 1

        # add ip1 triangle vertex to voronoi cell
        px, py = tri_verts[ip1]
        vor_verts.append((px, py))
        vor_iverts[ip1].append(inewvert)
        inewvert += 1

    # Last step -- sort vertices in correct order
    vor_verts = np.array(vor_verts)
    for icell in range(len(vor_iverts)):
        iverts_cell = vor_iverts[icell]
        vor_iverts[icell] = list(get_sorted_vertices(iverts_cell, vor_verts))

    # remove empty polygons/iverts, point, and line freatures
    # and their associated xy centers
    points = list(tri.verts)
    pop_list = []
    for icell, ivlist in enumerate(vor_iverts):
        if len(ivlist) < 3:
            pop_list.append(icell)

    for icell in pop_list[::-1]:
        vor_iverts.pop(icell)
        points.pop(icell)

    points = np.array(points)

    return vor_verts, vor_iverts, points


# Numerical policy of get_disu5_gridprops: vertices closer than
# _DISU5_RTOL times the extent of the vertices used by the cells are merged;
# a face is orthogonal when |cos| between the face and the line joining the
# two cell centers is at most _DISU5_ORTHO.
_DISU5_RTOL = 1e-9
_DISU5_ORTHO = 1e-6


def _disu5_layers(top, botm, ncpl):
    """Layer tops and bottoms, shape (nlay, ncpl), from top and botm."""
    try:
        top = np.array(top, dtype=float)
    except (TypeError, ValueError):
        raise ValueError("top must be a number or an array of numbers")
    if top.ndim == 0:
        top = np.full(ncpl, float(top))
    if top.shape != (ncpl,):
        raise ValueError(
            f"top must be a scalar or have shape ({ncpl},), got shape {top.shape}"
        )
    if isinstance(botm, np.ndarray) and botm.ndim in (1, 2):
        layers = list(botm)
    elif isinstance(botm, (list, tuple)):
        layers = list(botm)
    else:
        raise ValueError("botm must be a list, tuple or array with one entry per layer")
    if len(layers) == 0:
        raise ValueError("botm must have at least one layer")
    bots = []
    for k, b in enumerate(layers):
        try:
            b = np.array(b, dtype=float)
        except (TypeError, ValueError):
            raise ValueError(f"botm[{k}] must be a number or an array of numbers")
        if b.ndim == 0:
            b = np.full(ncpl, float(b))
        if b.shape != (ncpl,):
            raise ValueError(
                f"botm[{k}] must be a scalar or have shape ({ncpl},), "
                f"got shape {b.shape}"
            )
        bots.append(b)
    bots = np.array(bots)
    tops = np.vstack([top[np.newaxis, :], bots[:-1]])
    if not (np.all(np.isfinite(tops)) and np.all(np.isfinite(bots))):
        raise ValueError("top and botm must be finite")
    bad = np.argwhere(tops - bots <= 0.0)
    if len(bad) > 0:
        k, c = bad[0]
        raise ValueError(
            f"{len(bad)} cell(s) with non-positive thickness, first in layer {k} "
            f"cell {c} (top {tops[k, c]}, bottom {bots[k, c]})"
        )
    return tops, bots


def _disu5_merge(points, tol):
    """Merge points closer than tol (any pair, found through a grid of bins
    of size tol and their neighbours). Returns the cluster id of each point
    and the mean coordinates of each cluster."""
    n = len(points)
    parent = list(range(n))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    xy = points.tolist()
    keys = np.floor((points - points.min(axis=0)) / tol).astype(np.int64).tolist()
    bins = {}
    for i, (kx, ky) in enumerate(keys):
        x, y = xy[i]
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for j in bins.get((kx + dx, ky + dy), ()):
                    if sqrt((x - xy[j][0]) ** 2 + (y - xy[j][1]) ** 2) <= tol:
                        ri, rj = find(i), find(j)
                        if ri != rj:
                            parent[max(ri, rj)] = min(ri, rj)
        bins.setdefault((kx, ky), []).append(i)
    roots = np.array([find(i) for i in range(n)])
    _, cid = np.unique(roots, return_inverse=True)
    cid = cid.ravel()
    count = np.bincount(cid)
    cxy = (
        np.column_stack(
            [np.bincount(cid, points[:, 0]), np.bincount(cid, points[:, 1])]
        )
        / count[:, np.newaxis]
    )
    spread = np.hypot(*(points - cxy[cid]).T)
    if np.any(spread > tol):
        i = int(np.argmax(spread))
        raise ValueError(
            f"vertices near ({points[i, 0]}, {points[i, 1]}) form a chain of points "
            f"closer than the tolerance {tol}: the geometry cannot be resolved"
        )
    return cid, cxy


def _disu5_segment_distance(p1, p2, q1, q2):
    """Smallest distance between segments p1-p2 and q1-q2 (0 if they cross)."""

    def orient(a, b, c):
        return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])

    d1, d2 = orient(q1, q2, p1), orient(q1, q2, p2)
    d3, d4 = orient(p1, p2, q1), orient(p1, p2, q2)
    if d1 * d2 < 0 and d3 * d4 < 0:
        return 0.0

    def point_segment(p, a, b):
        ab = b - a
        t = np.dot(p - a, ab) / np.dot(ab, ab)
        t = min(1.0, max(0.0, t))
        return np.hypot(*(p - (a + t * ab)))

    return min(
        point_segment(p1, q1, q2),
        point_segment(p2, q1, q2),
        point_segment(q1, p1, p2),
        point_segment(q2, p1, p2),
    )


def _disu5_point_in_ring(p, xy, tol):
    """1 if p is inside the polygon xy farther than tol from its boundary,
    0 if it is within tol of the boundary, -1 if it is outside."""
    inside = False
    n = len(xy)
    for i in range(n):
        a, b = xy[i], xy[(i + 1) % n]
        ab = b - a
        t = min(1.0, max(0.0, np.dot(p - a, ab) / np.dot(ab, ab)))
        if np.hypot(*(p - (a + t * ab))) <= tol:
            return 0
        if (a[1] > p[1]) != (b[1] > p[1]):
            if a[0] + (p[1] - a[1]) * ab[0] / ab[1] > p[0]:
                inside = not inside
    return 1 if inside else -1


def _disu5_bins(lo, hi, size):
    """Map bin -> items whose bounding box [lo, hi] covers the bin."""
    origin = lo.min(axis=0)
    i0 = np.floor((lo - origin) / size).astype(np.int64).tolist()
    i1 = np.floor((hi - origin) / size).astype(np.int64).tolist()
    bins = {}
    for item, ((x0, y0), (x1, y1)) in enumerate(zip(i0, i1)):
        for bx in range(x0, x1 + 1):
            for by in range(y0, y1 + 1):
                bins.setdefault((bx, by), []).append(item)
    return bins, origin


def _disu5_cells(verts, cells, xcyc):
    """Cells as counterclockwise rings of merged vertex ids, with their
    areas, the merged vertex coordinates and the tolerance."""
    used = np.unique(np.concatenate([np.asarray(c, dtype=int) for c in cells]))
    pts = verts[used]
    extent = max(np.ptp(pts[:, 0]), np.ptp(pts[:, 1]))
    if not extent > 0.0:
        raise ValueError("the cells have no extent")
    tol = _DISU5_RTOL * extent
    cid, vxy = _disu5_merge(pts, tol)
    canon = np.full(len(verts), -1)
    canon[used] = cid
    rings = []
    area = np.zeros(len(cells))
    for c, idx in enumerate(cells):
        ring = []
        for i in canon[idx].tolist():
            if not ring or ring[-1] != i:
                ring.append(i)
        if len(ring) > 1 and ring[0] == ring[-1]:
            ring.pop()
        if len(ring) < 3:
            raise ValueError(f"cell {c}: fewer than 3 distinct vertices")
        if len(set(ring)) != len(ring):
            raise ValueError(
                f"cell {c}: a vertex is repeated (the polygon touches itself)"
            )
        xy = vxy[ring] - xcyc[c]
        nxt = np.roll(xy, -1, axis=0)
        a = 0.5 * np.sum(xy[:, 0] * nxt[:, 1] - nxt[:, 0] * xy[:, 1])
        perimeter = np.sum(np.hypot(*(nxt - xy).T))
        if abs(a) <= 1e-12 * perimeter**2:
            raise ValueError(f"cell {c}: zero area")
        k = len(ring)
        for i in range(k):
            # consecutive edges meeting at vertex i must not fold back
            u = xy[i - 1] - xy[i]
            w = xy[(i + 1) % k] - xy[i]
            cross = u[0] * w[1] - u[1] * w[0]
            if abs(cross) <= tol * max(np.hypot(*u), np.hypot(*w)) and np.dot(u, w) > 0:
                raise ValueError(f"cell {c}: the polygon folds back on itself")
            # edges that are not consecutive must not touch
            for j in range(i + 2, k - (i == 0)):
                p, q = xy[i], xy[(i + 1) % k]
                if _disu5_segment_distance(p, q, xy[j], xy[(j + 1) % k]) <= tol:
                    raise ValueError(f"cell {c}: the polygon intersects itself")
        if a < 0:
            ring.reverse()
            xy = xy[::-1]
        if _disu5_point_in_ring(np.zeros(2), xy, tol) < 0:
            raise ValueError(f"cell {c}: the center is outside the cell")
        rings.append(ring)
        area[c] = abs(a)
    return rings, area, vxy, tol


def _disu5_check_overlaps(rings, edges, vxy, xcyc, tol):
    """Raise if edges cross, overlap or touch away from shared vertices, or if
    a cell center lies inside another cell."""
    seg = np.array([[vxy[a], vxy[b]] for a, b in edges])
    lo = np.minimum(seg[:, 0], seg[:, 1]) - tol
    hi = np.maximum(seg[:, 0], seg[:, 1]) + tol
    size = max(np.median(np.hypot(*(seg[:, 1] - seg[:, 0]).T)), tol)
    bins, _ = _disu5_bins(lo, hi, size)
    checked = set()
    for members in bins.values():
        for i in range(len(members)):
            for j in range(i + 1, len(members)):
                e, f = members[i], members[j]
                if (e, f) in checked:
                    continue
                checked.add((e, f))
                (a, b), (c, d) = edges[e], edges[f]
                shared = {a, b} & {c, d}
                if shared:
                    s = shared.pop()
                    u = vxy[b if a == s else a] - vxy[s]
                    w = vxy[d if c == s else c] - vxy[s]
                    cross = u[0] * w[1] - u[1] * w[0]
                    if (
                        abs(cross) <= tol * max(np.hypot(*u), np.hypot(*w))
                        and np.dot(u, w) > 0
                    ):
                        raise ValueError(
                            f"edges {a}-{b} and {c}-{d} overlap: cells overlap or do "
                            "not share complete edges"
                        )
                    continue
                if (
                    _disu5_segment_distance(seg[e, 0], seg[e, 1], seg[f, 0], seg[f, 1])
                    <= tol
                ):
                    raise ValueError(
                        f"edges {a}-{b} and {c}-{d} cross or touch: cells overlap or a "
                        "vertex lies inside an edge (cells must share complete edges)"
                    )
    rxy = [vxy[r] for r in rings]
    lo = np.array([p.min(axis=0) for p in rxy])
    hi = np.array([p.max(axis=0) for p in rxy])
    size = max(np.median(np.max(hi - lo, axis=1)), tol)
    bins, origin = _disu5_bins(lo, hi, size)
    keys = np.floor((xcyc - origin) / size).astype(np.int64).tolist()
    for n, key in enumerate(keys):
        for m in bins.get(tuple(key), ()):
            if m == n or np.any(xcyc[n] < lo[m]) or np.any(xcyc[n] > hi[m]):
                continue
            if _disu5_point_in_ring(xcyc[n], rxy[m], tol) > 0:
                raise ValueError(
                    f"the center of cell {n} lies inside cell {m}: cells overlap"
                )


def _disu5_faces(rings, vxy, xcyc, tol):
    """neighbors[n][m] = (face width, distance from center n to the face)."""
    uses = {}
    for c, ring in enumerate(rings):
        for a, b in zip(ring, ring[1:] + ring[:1]):
            uses.setdefault((min(a, b), max(a, b)), []).append((c, a, b))
    pairs = {}
    for key, cells in uses.items():
        if len(cells) > 2:
            raise ValueError(
                f"an edge is used by more than two cells: {[c[0] for c in cells]}"
            )
        if len(cells) == 2:
            (n, a1, b1), (m, a2, b2) = cells
            if (a1, b1) == (a2, b2):
                raise ValueError(
                    f"cells {n} and {m} overlap (a shared edge has the same direction)"
                )
            pairs.setdefault((min(n, m), max(n, m)), []).append(key)
    _disu5_check_overlaps(rings, list(uses), vxy, xcyc, tol)
    # the unshared edges must form one boundary: holes, gaps between cells
    # and separate regions are not supported
    parent = {}

    def find(i):
        while parent.setdefault(i, i) != i:
            i = parent[i]
        return i

    for (a, b), cells in uses.items():
        if len(cells) == 1:
            parent[find(a)] = find(b)
    nbound = len({find(i) for i in list(parent)})
    if nbound > 1:
        raise ValueError(
            f"the unshared cell edges form {nbound} separate boundaries: grids with "
            "holes, gaps between cells or separate regions are not supported"
        )
    neighbors = [{} for _ in rings]
    for (n, m), keys in pairs.items():
        segs = [(vxy[a], vxy[b]) for a, b in keys]
        p0, p1 = max(segs, key=lambda s: np.hypot(*(s[1] - s[0])))
        e = p1 - p0
        length = np.hypot(*e)
        for q in (q for s in segs for q in s):
            if abs(e[0] * (q[1] - p0[1]) - e[1] * (q[0] - p0[0])) / length > tol:
                raise ValueError(f"cells {n} and {m} share more than one face")
        width = sum(np.hypot(*(s[1] - s[0])) for s in segs)
        dnm = xcyc[m] - xcyc[n]
        dist = np.hypot(*dnm)
        if dist <= tol:
            raise ValueError(f"cells {n} and {m} have the same center")
        if abs(np.dot(dnm, e)) > _DISU5_ORTHO * dist * length:
            raise ValueError(
                f"cells {n} and {m}: the shared face is not orthogonal to the line "
                "joining the cell centers (not a voronoi pair)"
            )
        cn = (e[0] * (xcyc[n, 1] - p0[1]) - e[1] * (xcyc[n, 0] - p0[0])) / length
        cm = (e[0] * (xcyc[m, 1] - p0[1]) - e[1] * (xcyc[m, 0] - p0[0])) / length
        if cn * cm >= 0.0 or abs(abs(cn) + abs(cm) - dist) > 1e-6 * dist:
            raise ValueError(
                f"cells {n} and {m}: the centers are not on opposite sides of the face"
            )
        neighbors[n][m] = (width, abs(cn))
        neighbors[m][n] = (width, abs(cm))
    return neighbors


def get_disu5_gridprops(verts, iverts, xcyc, top, botm):
    """
    Get a dictionary of arguments for flopy.mfusg.MfUsgDisU for a voronoi
    grid extruded into layers that all have the same cells.

    Supported: contiguous layers (the top of a layer is the bottom of the
    layer above, no confining beds), IVSD = 0 and IDSYMRD = 0, cells whose
    shared faces are orthogonal to the line joining their centers (voronoi
    cells with their generator points as centers) covering one region
    without holes. Not supported: holes (a VoronoiGrid built with holes does
    not cover its domain), separate regions, different cells per layer,
    pinched-out layers, quasi-3D confining beds, CLN and GNC connections, and
    face angles for horizontal anisotropy (LPF reads them when CHANI is not
    1). The arrays are only checked for groundwater
    flow (LPF); their use by the transport processes is not validated.

    Parameters
    ----------
    verts : array of floats, shape (nvert, 2)
        x, y of the vertices. Vertices not used by any cell are ignored.
    iverts : list of lists of ints
        Vertex numbers (zero-based integers) of each cell, in either
        direction, closed or not.
    xcyc : array of floats, shape (ncpl, 2)
        Cell centers (for a voronoi grid, its generator points).
    top : float or array of floats, shape (ncpl,)
        Top of the first layer.
    botm : list, tuple or array, one entry per layer
        Bottom of each layer, a float or an array of shape (ncpl,) per layer,
        or an array of shape (nlay, ncpl).

    Returns
    -------
    gridprops : dict
        nodes, nlay, njag, ivsd (0), idsymrd (0), iac, ja (zero-based),
        nodelay, top, bot, area (one value per node), ivc, cl12 and fahl
        (one value per connection). Node numbers are layer * ncpl + cell.
        Each JA row lists the node, then its connected nodes in increasing
        order. CL12 in row n is the distance from the center of n to the
        shared face (half the cell thickness for a vertical connection).
        FAHL is the face area: for a horizontal connection the face width
        times the mean thickness of the two cells (USG-T divides by that
        mean thickness to get the width back), for a vertical connection the
        cell area. Lengths have the units of verts, top and botm.

    Notes
    -----
    Vertices closer than 1e-9 times the extent of the cell vertices are
    merged; a geometry that cannot be resolved at that tolerance raises
    ValueError, as do cells that overlap, intersect themselves, have a
    vertex inside another cell edge or a non-positive thickness, and cells
    whose unshared edges form more than one boundary (holes, gaps, separate
    regions). Cells that
    touch at a single point are not connected. MfUsgDisU converts a
    zero-based ja array to one-based in place: pass a copy if the
    dictionary is used more than once.

    """
    try:
        verts = np.array(verts, dtype=float)
        xcyc = np.array(xcyc, dtype=float)
    except (TypeError, ValueError):
        raise ValueError("verts and xcyc must be arrays of numbers")
    if verts.ndim != 2 or verts.shape[1] != 2 or len(verts) == 0:
        raise ValueError(f"verts must have shape (nvert, 2), got shape {verts.shape}")
    if isinstance(iverts, (str, bytes)) or not hasattr(iverts, "__len__"):
        raise ValueError("iverts must be a list of vertex-number lists")
    ncpl = len(iverts)
    if ncpl == 0:
        raise ValueError("iverts is empty")
    if xcyc.shape != (ncpl, 2):
        raise ValueError(f"xcyc must have shape ({ncpl}, 2), got shape {xcyc.shape}")
    if not (np.all(np.isfinite(verts)) and np.all(np.isfinite(xcyc))):
        raise ValueError("verts and xcyc must be finite")
    cells = []
    for c, iv in enumerate(iverts):
        if isinstance(iv, (str, bytes)) or not hasattr(iv, "__iter__"):
            raise ValueError(f"cell {c}: the vertex numbers must be a list of integers")
        idx = []
        for i in iv:
            if isinstance(i, (bool, np.bool_)) or not isinstance(i, (int, np.integer)):
                raise ValueError(f"cell {c}: vertex number {i!r} is not an integer")
            idx.append(int(i))
        if not idx:
            raise ValueError(f"cell {c}: no vertices")
        if min(idx) < 0 or max(idx) >= len(verts):
            raise ValueError(f"cell {c}: vertex number out of range")
        cells.append(idx)
    tops, bots = _disu5_layers(top, botm, ncpl)
    nlay = len(bots)
    thick = tops - bots
    rings, area, vxy, tol = _disu5_cells(verts, cells, xcyc)
    neighbors = _disu5_faces(rings, vxy, xcyc, tol)

    iac, ja, ivc, cl12, fahl = [], [], [], [], []
    for k in range(nlay):
        for c in range(ncpl):
            n = k * ncpl + c
            conn = [
                (k * ncpl + m, 0, cl, w * 0.5 * (thick[k, c] + thick[k, m]))
                for m, (w, cl) in neighbors[c].items()
            ]
            if k > 0:
                conn.append(((k - 1) * ncpl + c, 1, 0.5 * thick[k, c], area[c]))
            if k < nlay - 1:
                conn.append(((k + 1) * ncpl + c, 1, 0.5 * thick[k, c], area[c]))
            conn.sort()
            conn.insert(0, (n, 0, 0.0, 0.0))
            iac.append(len(conn))
            for m, v, cl, fa in conn:
                ja.append(m)
                ivc.append(v)
                cl12.append(cl)
                fahl.append(fa)
    iac = np.array(iac, dtype=np.int32)
    return {
        "nodes": nlay * ncpl,
        "nlay": nlay,
        "njag": int(iac.sum()),
        "ivsd": 0,
        "idsymrd": 0,
        "iac": iac,
        "ja": np.array(ja, dtype=np.int32),
        "nodelay": np.full(nlay, ncpl, dtype=np.int32),
        "top": tops.ravel(),
        "bot": bots.ravel(),
        "area": np.tile(area, nlay),
        "ivc": np.array(ivc, dtype=np.int32),
        "cl12": np.array(cl12),
        "fahl": np.array(fahl),
    }


class VoronoiGrid:
    """
    FloPy VoronoiGrid helper class for creating a voronoi model grid from
    an array of input points that define cell centers.  The class handles
    boundary cells by closing polygons along the edge, something that cannot
    be done directly with the scipy.spatial.Voronoi class.

    Parameters
    ----------
    input : flopy.utils.Triangle
        Constructed and built flopy Triangle object.
    kwargs : dict
        List of additional keyword arguments that will be passed through to
        scipy.spatial.Voronoi.  For circular shaped model grids, the
        qhull_options='Qz' option has been found to work well.

    Notes
    -----
    When using VoronoiGrid, the construction order used for the Triangle
    grid matters.  The first add_polygon() call must be to add the model
    domain.  Then add_polygon() must be used to add any holes.  Lastly,
    add_polygon() can be used to add regions.  This VoronoiGrid class uses
    this order to find model edges that require further work for defining and
    closing edge model cells.

    """

    def __init__(self, tri, **kwargs):
        from .triangle import Triangle

        if isinstance(tri, Triangle):
            verts, iverts, points = tri2vor(tri, **kwargs)
        else:
            raise TypeError("The tri argument must be of type flopy.utils.Triangle")
        self.points = points
        self.verts = verts
        self.iverts = iverts
        self.ncpl = len(iverts)
        self.nverts = verts.shape[0]
        return

    def get_disv_gridprops(self):
        """
        Get a dictionary of arguments that can be passed in to the
        flopy.mf6.ModflowGwfdisv class.

        Returns
        -------
        disv_gridprops : dict
            Dictionary of arguments than can be unpacked into the
            flopy.mf6.ModflowGwfdisv constructor

        """
        disv_gridprops = get_disv_gridprops(self.verts, self.iverts, xcyc=self.points)
        return disv_gridprops

    def get_disu5_gridprops(self, top, botm):
        """
        Get a dictionary of arguments that can be passed in to the
        flopy.mfusg.MfUsgDisU class for the voronoi grid extruded into
        layers with the same cells. The cell centers are the voronoi points.
        See the module function get_disu5_gridprops for the arrays, the
        supported geometry and the limitations.

        Parameters
        ----------
        top : float or array of floats, shape (ncpl,)
            Top of the first layer.
        botm : list, tuple or array, one entry per layer
            Bottom of each layer, a float or an array of shape (ncpl,) per
            layer, or an array of shape (nlay, ncpl).

        Returns
        -------
        disu_gridprops : dict
            Dictionary of arguments than can be unpacked into the
            flopy.mfusg.MfUsgDisU constructor (ja is zero-based)

        """
        return get_disu5_gridprops(self.verts, self.iverts, self.points, top, botm)

    def get_disu6_gridprops(self):
        msg = "This method is not implemented yet."
        raise NotImplementedError(msg)

    def get_gridprops_vertexgrid(self):
        """
        Get a dictionary of information needed to create a flopy VertexGrid.
        The returned dictionary can be unpacked directly into the
        flopy.discretization.VertexGrid() constructor.

        Returns
        -------
        gridprops : dict

        """
        gridprops = self.get_disv_gridprops()
        del gridprops["nvert"]
        return gridprops

    def get_gridprops_unstructuredgrid(self, top=None, botm=None):
        """
        Get a dictionary of information needed to create a flopy
        UnstructuredGrid.  The returned dictionary can be unpacked directly
        into the flopy.discretization.UnstructuredGrid() constructor.

        Parameters
        ----------
        top : float or array of floats, shape (ncpl,), optional
            Top of the first layer. Must be given together with botm.
        botm : list, tuple or array, one entry per layer, optional
            Bottom of each layer, as for get_disu5_gridprops. With top and
            botm the dictionary describes the grid extruded into layers with
            the same cells (the node order of get_disu5_gridprops): iverts,
            xcenters and ycenters repeated for each layer, ncpl per layer,
            and top and botm per node.

        Returns
        -------
        gridprops : dict

        """
        if (top is None) != (botm is None):
            raise ValueError("top and botm must be given together")
        gridprops = {}

        disv_gridprops = self.get_disv_gridprops()
        vertices = disv_gridprops["vertices"]
        iverts = self.iverts
        ncpl = self.ncpl
        xcenters = self.points[:, 0]
        ycenters = self.points[:, 1]

        gridprops["vertices"] = vertices
        gridprops["iverts"] = iverts
        gridprops["ncpl"] = ncpl
        gridprops["xcenters"] = xcenters
        gridprops["ycenters"] = ycenters

        if top is not None:
            tops, bots = _disu5_layers(top, botm, ncpl)
            nlay = len(bots)
            gridprops["iverts"] = [
                [int(i) for i in iv] for _ in range(nlay) for iv in iverts
            ]
            gridprops["ncpl"] = [ncpl] * nlay
            gridprops["xcenters"] = np.tile(xcenters, nlay)
            gridprops["ycenters"] = np.tile(ycenters, nlay)
            gridprops["top"] = tops.ravel()
            gridprops["botm"] = bots.ravel()

        return gridprops

    def get_patch_collection(self, ax=None, **kwargs):
        """
        Get a matplotlib patch collection representation of the voronoi grid

        Parameters
        ----------
        ax : matplotlib.pyplot.Axes
            axes to plot the patch collection
        kwargs : dict
            Additional keyword arguments to pass to the flopy.plot.plot_cvfd
            function that returns a patch collection from verts and iverts

        Returns
        -------
        pc : matplotlib.collections.PatchCollection
            patch collection of model

        """
        from ..discretization import VertexGrid
        from ..plot import PlotMapView

        modelgrid = VertexGrid(**self.get_gridprops_vertexgrid())
        pmv = PlotMapView(modelgrid=modelgrid, ax=ax)
        pc = pmv.plot_grid(**kwargs)
        return pc

    def plot(self, ax=None, plot_title=True, **kwargs):
        """
        Plot the voronoi model grid

        Parameters
        ----------
        ax : matplotlib.pyplot.Axes
            axes to plot the voronoi grid
        plot_title : bool
            Add the number of cells and number of vertices as a plot title
        kwargs : dict
            Additional keyword arguments to pass to self.get_patch_collection

        Returns
        -------
        ax : matplotlib.pyplot.Axes
            axes that contains the voronoi model grid

        """
        import matplotlib.pyplot as plt

        if ax is None:
            ax = plt.subplot(1, 1, 1, aspect="equal")
        pc = self.get_patch_collection(ax, **kwargs)
        if plot_title:
            ax.set_title(f"ncells: {self.ncpl}; nverts: {self.nverts}")
        return ax
