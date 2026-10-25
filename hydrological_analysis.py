'''
Hydrological analysis pipeline for a watershed boundary
=====================================================
Author: Shubham Kumar
Date: April 2026

Usage:
    python hydrological_analysis.py --shp path_to_shapefile \
        --output path_to_output_directory \
        --threshold flow_accumulation_threshold

Workflow:
    1. Select a suitable UTM CRS for the watershed footprint.
    2. Download and reproject the DEM to the local watershed extent.
    3. Condition the DEM for hydrologic processing within a 10 km boundary buffer.
    4. Run GRASS flow-direction and flow-accumulation routines.
    5. Extract stream networks and derive stream order.
    6. Merge small micro-watersheds and export the derived raster/vector outputs.
'''

import geopandas as gpd
import matplotlib.pyplot as plt
from pathlib import Path
import argparse
import shutil
import os
import sys
import subprocess
import rasterio
import json
import numpy as np
import math
from contextlib import contextmanager
from dem_downloader import DEMDownloader
import time
import richdem as rd
from typing import Optional
import warnings
from typing import Sequence
from pyproj.aoi import AreaOfInterest
from pyproj.database import query_utm_crs_info


# Minimum area used when merging tiny micro-watersheds; expressed in hectares.
MIN_WATERSHED_SIZE = 500
MAX_MERGED_WATERSHED_SIZE_HA = 1500
BUFFER_DISTANCE_METERS = 10_000
# Initial r.watershed accumulation threshold in cells for preliminary basin delineation.
# This should later be tied to DEM resolution and target watershed scale.
HYPER_PARAM = 1200

CONFIG = {
    "DEM": "srtm30"  # DEM key used by DEMDownloader; see dem_downloader.py for registry details.
}


@contextmanager
def timed(stage, timings):
    start_time = time.perf_counter()
    try:
        yield
    finally:
        timings[stage] = timings.get(stage, 0.0) + time.perf_counter() - start_time


def format_duration(seconds):
    hours = int(seconds // 3600)
    minutes = int((seconds % 3600) // 60)
    remaining_seconds = seconds % 60
    return f"{hours:02d}:{minutes:02d}:{remaining_seconds:05.2f}"


# CLI arguments for the watershed processing pipeline. These flags control the input
# polygon, output directory, GRASS database, and the hydrologic thresholds used to
# derive stream networks and final micro-watershed units.
def parse_args():
    parser = argparse.ArgumentParser(
        description="GRASS GIS Hydrological Analysis Pipeline"
    )
    parser.add_argument(
        "--shp", required=True,
        help="Path to input vector defining the watershed boundary"
    )
    parser.add_argument(
        "--output", required=True,
        help="Directory for output files"
    )
    parser.add_argument(
        "--grassdb", default=str(Path.home() / "grassdata"),
        help="Directory for GRASS GIS database (default: ~/grassdata)"
    )
    parser.add_argument(
        "--threshold", type=int, default=100,
        help="Flow-accumulation threshold for stream extraction (default: 100 cells)"
    )
    parser.add_argument(
        "--min_watershed_size", type=int, default=MIN_WATERSHED_SIZE,
        help="Minimum watershed size in hectares (default: 500 hectares)"
    )

    return parser.parse_args()

# Initialize a GRASS GIS session and load the Python bindings used by the workflow.
def setup_grass_session(grassdb: str, epsg: int, project_name: str = "hydro_project"):
    grassdb_path = Path(grassdb).resolve()
    location = project_name
    mapset = "PERMANENT"

    grass_bin = shutil.which("grass")
    if grass_bin is None:
        sys.exit("ERROR: GRASS GIS not found on PATH.")

    result = subprocess.run(
        [grass_bin, "--config", "python_path"],
        capture_output=True,
        text=True,
        env={**os.environ, "PYTHONWARNINGS": "ignore"})

    grass_python = result.stdout.strip()
    if grass_python and grass_python not in sys.path:
        sys.path.insert(0, grass_python)
        print(f"[INFO] Added GRASS python path: {grass_python}")

    loc_path = grassdb_path / location
    if not loc_path.exists():
        grassdb_path.mkdir(parents=True, exist_ok=True)
        subprocess.run([grass_bin, "-c", f"EPSG:{epsg}", "-e", str(loc_path)],
                       check=True)
        print(f"[INFO] GRASS location created: {loc_path}")

    # Prefer the Session API when available; otherwise initialize GRASS using the
    # lower-level setup helper.
    try:
        from grass_session import Session
        session = Session()
        session.open(gisdb=str(grassdb_path), location=location, mapset=mapset)
        print("grass_session initialised successfully.")
        return session
    except ImportError:
        pass

    # Legacy fallback for environments without grass_session installed.
    import grass.script.setup as gsetup
    gsetup.init(str(grassdb_path), location, mapset)
    print("GRASS environment initialised via grass.script.setup")
    return None

# Pick the most suitable UTM CRS for the watershed footprint from its bounding box.

def _lon_intervals(west, east):
    """Split a longitude range into plain intervals (handles west > east)."""
    if west <= east:
        return [(west, east)]
    return [(west, 180.0), (-180.0, east)]


def _overlap_area(box, other):
    """
    Spherical area (unit sphere, steradians) of the intersection of two
    lon/lat boxes given as (west, south, east, north). Either box may cross
    the antimeridian (west > east).
    """
    lat_lo = max(box[1], other[1])
    lat_hi = min(box[3], other[3])
    if lat_hi <= lat_lo:
        return 0.0

    lon_overlap = 0.0
    for a0, a1 in _lon_intervals(box[0], box[2]):
        for b0, b1 in _lon_intervals(other[0], other[2]):
            lon_overlap += max(0.0, min(a1, b1) - max(a0, b0))

    return math.radians(lon_overlap) * (
        math.sin(math.radians(lat_hi)) - math.sin(math.radians(lat_lo))
    )


def get_utm_epsg_for_bbox(bbox: Sequence[float]) -> int:
    """
    Return the WGS 84 / UTM EPSG code (int) whose zone overlaps the lon/lat
    bounding box the most.

    bbox: (west, south, east, north) in degrees. A box crossing the
    antimeridian is given with west > east.

    pyproj returns every UTM zone that intersects the box, in no useful
    order. Each candidate is scored by the spherical area of its overlap with
    the box; the largest wins. Ties (including zero-area boxes such as a
    single point) go to the zone in the same hemisphere as the box centre,
    then to the zone whose central meridian is closest to the box centre.

    A warning is issued if the box touches more than one UTM zone.
    """
    west, south, east, north = map(float, bbox)

    if not all(map(math.isfinite, (west, south, east, north))):
        raise ValueError(f"bbox contains non-finite values: {bbox}")
    if not (-180 <= west <= 180 and -180 <= east <= 180):
        raise ValueError(f"Longitudes must be within [-180, 180]: {bbox}")
    if not (-90 <= south <= north <= 90):
        raise ValueError(f"Latitudes must satisfy -90 <= south <= north <= 90: {bbox}")

    box = (west, south, east, north)

    candidates = query_utm_crs_info(
        datum_name="WGS 84",
        area_of_interest=AreaOfInterest(
            west_lon_degree=west,
            south_lat_degree=south,
            east_lon_degree=east,
            north_lat_degree=north,
        ),
    )

    # Centre of the box, used only for tie-breaking.
    east_u = east + 360.0 if west > east else east
    lon_c = ((west + east_u) / 2.0 + 180.0) % 360.0 - 180.0
    lat_c = (south + north) / 2.0

    scored = []
    for info in candidates:
        aou = info.area_of_use
        code = int(info.code)
        if aou is None or not (32601 <= code <= 32660 or 32701 <= code <= 32760):
            continue  # skip anything that is not a WGS 84 UTM zone

        area = _overlap_area(box, (aou.west, aou.south, aou.east, aou.north))
        zone = code % 100
        northern = code < 32700
        central_meridian = -183.0 + 6.0 * zone
        lon_dist = abs((lon_c - central_meridian + 180.0) % 360.0 - 180.0)
        same_hemisphere = northern == (lat_c >= 0)

        # Larger overlap first; then hemisphere match; then closer meridian.
        scored.append(((round(area, 12), same_hemisphere, -lon_dist), code, zone, area))

    if not scored:
        raise ValueError(f"No suitable UTM CRS found for bbox {bbox}.")

    scored.sort(key=lambda s: s[0], reverse=True)
    best_code = scored[0][1]

    touched_zones = {zone for _, _, zone, area in scored if area > 0}
    if len(touched_zones) > 1:
        warnings.warn(
            f"bbox touches {len(touched_zones)} UTM zones; EPSG:{best_code} has the "
            "largest overlap, but distortion grows towards the box edges.",
            stacklevel=2,
        )

    return best_code

# Condition the DEM for hydrologic processing by removing sinks and recording depression depth.
def dem_preprocessing(
    input_dem,
    output_dir,
    watershed_path=None,
    depression_boundary_path=None,
):
    from rasterio.mask import mask as mask_raster
    from rasterio.features import geometry_mask

    conditioned_dem = str(Path(output_dir) / "dem_conditioned.tif")
    depression_tif = str(Path(output_dir) / "natural_depressions.tif")
    analysis_dem = input_dem

    with rasterio.open(input_dem) as src:
        nodata = src.nodata
        crs = src.crs
        profile = src.profile.copy()
        analysis_transform = src.transform

    if nodata is None:
        nodata = -9999.0

    if watershed_path is not None:
        boundary_gdf = gpd.read_file(watershed_path)
        if boundary_gdf.empty:
            raise ValueError(f"Analysis boundary is empty: {watershed_path}")
        if boundary_gdf.crs is None or crs is None:
            raise ValueError(
                f"Analysis boundary and DEM both need a CRS: {watershed_path}"
            )
        boundary_gdf = boundary_gdf.to_crs(crs)
        geometries = [
            geom for geom in boundary_gdf.geometry if geom is not None
        ]
        if not geometries:
            raise ValueError(f"Analysis boundary has no geometries: {watershed_path}")

        with rasterio.open(input_dem) as src:
            dem_data, transform = mask_raster(
                src,
                geometries,
                crop=True,
                nodata=nodata,
                filled=True,
            )
            analysis_transform = transform
            profile.update(
                height=dem_data.shape[1],
                width=dem_data.shape[2],
                transform=transform,
                nodata=nodata,
            )
        analysis_dem = str(Path(output_dir) / "dem_analysis_buffer.tif")
        with rasterio.open(analysis_dem, "w", **profile) as dst:
            dst.write(dem_data)

    dem_in = rd.LoadGDAL(analysis_dem, no_data=float(nodata))
    print(dem_in)

    filled_dem = rd.FillDepressions(dem_in, in_place=False)

    dem_arr = np.asarray(dem_in)
    filled_arr = np.asarray(filled_dem)
    nodata_mask = dem_arr == dem_in.no_data
    valid_mask = ~nodata_mask

    depression_boundary_path = (
        depression_boundary_path
        if depression_boundary_path is not None
        else watershed_path
    )
    depth_mask = valid_mask
    if depression_boundary_path is not None:
        depression_boundary_gdf = gpd.read_file(depression_boundary_path)
        if depression_boundary_gdf.empty:
            raise ValueError(
                f"Depression boundary is empty: {depression_boundary_path}"
            )
        if depression_boundary_gdf.crs is None or crs is None:
            raise ValueError(
                "Depression boundary and DEM both need a CRS: "
                f"{depression_boundary_path}"
            )
        depression_boundary_gdf = depression_boundary_gdf.to_crs(crs)
        depression_geometries = [
            geom for geom in depression_boundary_gdf.geometry if geom is not None
        ]
        if not depression_geometries:
            raise ValueError(
                "Depression boundary has no geometries: "
                f"{depression_boundary_path}"
            )
        inside_watershed = geometry_mask(
            depression_geometries,
            out_shape=dem_arr.shape,
            transform=analysis_transform,
            invert=True,
        )
        depth_mask = valid_mask & inside_watershed

    depth = np.full(dem_arr.shape, dem_in.no_data, dtype=filled_arr.dtype)
    depth[depth_mask] = np.clip(
        filled_arr[depth_mask] - dem_arr[depth_mask],
        0,
        None,
    )
    filled_dem[~valid_mask] = dem_in.no_data

    depth_rd = rd.rdarray(depth, no_data=dem_in.no_data)
    depth_rd.geotransform = dem_in.geotransform
    depth_rd.projection = dem_in.projection
    rd.SaveGDAL(depression_tif, depth_rd)

    rd.ResolveFlats(filled_dem, in_place=True) # Adding a small gradient to the filled DEM to resolve flats.
    rd.SaveGDAL(conditioned_dem, filled_dem)

    return filled_dem

def calculate_flow_accumulation(dem_filled, hyperparam_threshold):
    import grass.script as gs
    gs.run_command("r.watershed", 
                elevation=dem_filled,
                accumulation="flow_acc",
                drainage="flow_dir_watershed",
                threshold=hyperparam_threshold,
                # Use a low threshold so the earliest basin partition captures fine-scale drainage.
                basin="micro_watersheds",
                #stream="streams_raw",
                flags="as",  # -a: positive accumulation; -s: single-flow (D8)
                overwrite=True)

    return "flow_acc", "flow_dir_watershed", "micro_watersheds"


def merge_small_watersheds(
    micro_watersheds_rast: str,
    flow_acc_rast: str,
    flow_dir_rast: str,
    min_area_ha: float,
    epsg: Optional[int] = None,          # kept for API compatibility; no longer needed
    output_rast: str = "micro_watersheds",
    max_area_ha: Optional[float] = None,
) -> str:
    """
    Merge tiny micro-watersheds until each basin reaches the minimum area.

    A small basin is absorbed into its downstream neighbour when possible,
    otherwise into its largest upstream neighbour. If the preferred target
    would exceed the maximum size, the remaining neighbours are tried in turn.

    max_area_ha: maximum merged basin size. If None, the module-level constant
    MAX_MERGED_WATERSHED_SIZE_HA is used when defined, otherwise no limit.
    """
    import numpy as np
    import grass.script as gs
    from grass.script import array as garray

    NODATA = -9999

    # D8 offsets (GRASS r.watershed encoding): direction -> (d_row, d_col)
    DIR_OFFSETS = {
        1: (-1, +1), 2: (-1, 0), 3: (-1, -1), 4: (0, -1),
        5: (+1, -1), 6: (+1, 0), 7: (+1, +1), 8: (0, +1),
    }

    # setup
    if gs.locn_is_latlong():
        raise ValueError(
            "Location uses a geographic CRS (degrees); cell areas in hectares "
            "would be wrong. Use a projected location."
        )

    if max_area_ha is None:
        max_area_ha = globals().get("MAX_MERGED_WATERSHED_SIZE_HA", float("inf"))

    region = gs.region()
    nrows, ncols = int(region["rows"]), int(region["cols"])
    cell_ha = abs(float(region["nsres"])) * abs(float(region["ewres"])) / 10_000.0
    min_cells = min_area_ha / cell_ha
    max_cells = max_area_ha / cell_ha

    print(
        f"[merge_small_watersheds] cell = {cell_ha:.4f} ha, "
        f"min_cells = {min_cells:.1f}, max_cells = {max_cells:.1f}"
    )

    gs.run_command(
        "r.mapcalc",
        expr=f"_mws_work = int({micro_watersheds_rast})",
        overwrite=True,
    )
    try:
        basin_arr = garray.array("_mws_work", null=NODATA).astype(np.int32)
        acc_arr = garray.array(flow_acc_rast, null=0).astype(np.float32)
        flowdir_arr = garray.array(flow_dir_rast, null=0).astype(np.int8)
    finally:
        gs.run_command("g.remove", type="raster", name="_mws_work", flags="f")

    # static drainage geometry (once)
    # Negative r.watershed directions mean "drains out of the region"; the
    # magnitude still encodes the direction. Accumulation can also be negative.
    d_raw = flowdir_arr.astype(np.int32)
    d_abs = np.abs(d_raw)
    has_dir = (d_abs >= 1) & (d_abs <= 8)
    leaves_region = d_raw < 0
    d_clip = np.where(has_dir, d_abs, 0)

    offset_lut = np.zeros((9, 2), dtype=np.int32)
    for d, (dr, dc) in DIR_OFFSETS.items():
        offset_lut[d] = (dr, dc)

    R, C = np.indices((nrows, ncols), dtype=np.int32)
    nr = R + offset_lut[d_clip, 0]
    nc = C + offset_lut[d_clip, 1]
    in_bounds = (nr >= 0) & (nr < nrows) & (nc >= 0) & (nc < ncols)
    safe_r = np.clip(nr, 0, nrows - 1)
    safe_c = np.clip(nc, 0, ncols - 1)
    can_read_nbr = has_dir & in_bounds
    acc_abs = np.abs(acc_arr)

    def dedupe(seq):
        return list(dict.fromkeys(seq))

    # topology builder
    def build_maps(basins):
        """Return (downstream, area_cells, upstream) dicts for the current labels."""
        valid = basins > 0
        flat = basins[valid]
        if flat.size == 0:
            return {}, {}, {}

        counts = np.bincount(flat)
        ids = np.flatnonzero(counts)
        area = {int(i): int(counts[i]) for i in ids}
        downstream = {int(i): None for i in ids}
        upstream = {int(i): [] for i in ids}

        # Neighbour label, read only where the neighbour is on the grid.
        nbr = np.full(basins.shape, NODATA, dtype=np.int32)
        nbr[can_read_nbr] = basins[safe_r[can_read_nbr], safe_c[can_read_nbr]]

        src = valid & has_dir
        # Outlet: drains off the grid, out of the region, or into a null cell.
        outlet = src & (~in_bounds | leaves_region | (nbr <= 0))
        # Cross-basin: drains into a different valid basin.
        cross = src & ~outlet & (nbr != basins)
        cand = outlet | cross

        cr, cc = np.nonzero(cand)
        if cr.size == 0:
            return downstream, area, upstream

        cand_bid = basins[cr, cc]
        cand_tgt = np.where(outlet[cr, cc], -1, nbr[cr, cc])
        cand_acc = acc_abs[cr, cc]

        # Per basin, keep the candidate with the highest accumulation.
        order = np.lexsort((-cand_acc, cand_bid))
        bid_s, tgt_s = cand_bid[order], cand_tgt[order]
        first = np.ones(bid_s.size, dtype=bool)
        first[1:] = bid_s[1:] != bid_s[:-1]

        for b, t in zip(bid_s[first].tolist(), tgt_s[first].tolist()):
            downstream[b] = None if t == -1 else t

        for s, t in downstream.items():
            if t is not None:
                upstream[t].append(s)
        return downstream, area, upstream

    # merge passes
    total_merged = 0
    passes = 0

    while True:
        downstream, area, upstream = build_maps(basin_arr)

        small = sorted((a, b) for b, a in area.items() if a < min_cells)
        if not small:
            break

        parent = {}      # absorbed basin -> the basin it was merged into

        for _, bid in small:
            if bid not in area:                  # absorbed earlier in this pass
                continue
            cur = area[bid]
            if cur >= min_cells:                 # grew by absorbing others
                continue

            # Candidate targets: downstream first, then upstream by size.
            cands = []
            ds = downstream.get(bid)
            if ds is not None and ds in area:
                cands.append(ds)
            ups = [u for u in upstream.get(bid, []) if u in area]
            cands.extend(sorted(ups, key=lambda u: -area[u]))

            target = next((c for c in cands if area[c] + cur <= max_cells), None)
            if target is None:
                continue                         # no neighbour, or none fits

            # merge bid into target and patch the topology locally
            bid_ds = downstream.pop(bid, None)
            bid_us = upstream.pop(bid, [])

            for u in bid_us:                     # children now drain into target
                if u != target:
                    downstream[u] = target

            ups_t = upstream.setdefault(target, [])
            if bid_ds != target:
                # target was upstream of bid: it inherits bid's outlet
                downstream[target] = bid_ds
                if bid_ds is not None and bid_ds in upstream:
                    upstream[bid_ds] = dedupe(
                        target if x == bid else x for x in upstream[bid_ds]
                    )
            upstream[target] = dedupe(
                u for u in ups_t + bid_us if u != target and u != bid
            )

            area[target] += area.pop(bid)
            parent[bid] = target
            total_merged += 1

        passes += 1
        if not parent:
            break

        # Apply all relabelling of this pass in one lookup-table step.
        def find(x):
            while x in parent:
                x = parent[x]
            return x

        lut = np.arange(int(basin_arr.max()) + 1, dtype=np.int32)
        for b in parent:
            lut[b] = find(b)
        pos = basin_arr > 0
        basin_arr[pos] = lut[basin_arr[pos]]

    remaining = sum(1 for a in area.values() if a < min_cells)
    msg = f"{total_merged} merges in {passes} pass(es)"
    if remaining:
        msg += f"; {remaining} basin(s) stay below {min_area_ha} ha (no valid merge)"
    print(f"[merge_small_watersheds] {msg}")

    # compact renumbering 1..K
    pos = basin_arr > 0
    final = np.full(basin_arr.shape, NODATA, dtype=np.int32)
    ids, inverse = np.unique(basin_arr[pos], return_inverse=True)
    final[pos] = inverse.astype(np.int32) + 1

    # write back to GRASS
    # Written through the GRASS array interface: same region, same CRS, no
    # temporary GeoTIFF and no CRS mismatch on import.
    out = garray.array(dtype=np.int32)
    out[:, :] = final
    out.write(output_rast, overwrite=True, null=NODATA)

    print(
        f"[merge_small_watersheds] '{output_rast}' written "
        f"({ids.size} basins finally)."
    )
    return output_rast

def compute_pour_points(micro_watersheds_rast: str,
                        flow_acc_rast: str,
                        flow_dir_rast: str,
                        output_vector: str = "pour_points") -> tuple[str, dict]:
    import grass.script as gs
    from grass.script import array as garray
    import tempfile, csv

    # A pour point is the cell in each basin that drains to a different basin, or to the edge
    # of the domain. These are the true outlet cells that later support connectivity graphing.
    print("Pour Points: locating true outlet cell for each micro-watershed …")

    region  = gs.region()
    nrows   = int(region["rows"])
    ncols   = int(region["cols"])
    w       = region["w"]
    n       = region["n"]
    ewres   = region["ewres"]
    nsres   = region["nsres"]

    DIR_OFFSETS = {
        1: (-1, +1),  # NE
        2: (-1,  0),  # N
        3: (-1, -1),  # NW
        4: ( 0, -1),  # W
        5: (+1, -1),  # SW
        6: (+1,  0),  # S
        7: (+1, +1),  # SE
        8: ( 0, +1),  # E
    }

    dr_lookup = np.zeros(9, dtype=np.int32)
    dc_lookup = np.zeros(9, dtype=np.int32)
    for d, (dr, dc) in DIR_OFFSETS.items():
        dr_lookup[d] = dr
        dc_lookup[d] = dc

    gs.run_command("r.mapcalc",
                   expr=f"micro_watersheds_int = int({micro_watersheds_rast})",
                   overwrite=True)
    basin_arr = garray.array("micro_watersheds_int", null=-9999)
    acc_arr = garray.array(flow_acc_rast, null=-1)
    flowdir_arr = garray.array(flow_dir_rast, null=0)

    valid_mask = (basin_arr > 0) & (basin_arr != -9999)
    rows, cols = np.where(valid_mask)
    basins = basin_arr[rows, cols].astype(np.int32)
    acc_vals = acc_arr[rows, cols].astype(np.float64)
    dirs = flowdir_arr[rows, cols].astype(np.int32)

    valid_dir = (dirs >= 1) & (dirs <= 8)
    nr = np.full(rows.shape, -1, dtype=np.int32)
    nc = np.full(cols.shape, -1, dtype=np.int32)
    neighbor_basin = np.full(basins.shape, -1, dtype=np.int32)

    if valid_dir.any():
        nr[valid_dir] = rows[valid_dir] + dr_lookup[dirs[valid_dir]]
        nc[valid_dir] = cols[valid_dir] + dc_lookup[dirs[valid_dir]]
        in_bounds = (nr >= 0) & (nr < nrows) & (nc >= 0) & (nc < ncols)
        inside = valid_dir & in_bounds
        if inside.any():
            neighbor_basin[inside] = basin_arr[nr[inside], nc[inside]].astype(np.int32)
    else:
        in_bounds = np.zeros_like(valid_dir, dtype=bool)

    candidate = (~valid_dir) | (~in_bounds) | ((valid_dir & in_bounds) & (neighbor_basin != basins))
    candidate_idx = np.where(candidate)[0]

    pour_pts = {}
    records = []
    basin_ids = np.unique(basins)

    for bid in basin_ids.tolist():
        bid = int(bid)
        if candidate_idx.size:
            local = candidate_idx[basins[candidate_idx] == bid]
        else:
            local = np.array([], dtype=np.int64)

        if local.size == 0:
            local = np.where(basins == bid)[0]

        best_pos = int(local[np.argmax(acc_vals[local])])
        r = int(rows[best_pos])
        c = int(cols[best_pos])
        best_acc = float(acc_vals[best_pos])

        pour_pts[bid] = (r, c)
        x = w + (c + 0.5) * ewres
        y = n - (r + 0.5) * nsres
        records.append((bid, x, y, best_acc))

    tmp = tempfile.NamedTemporaryFile(mode="w", suffix=".csv",
                                      delete=False, newline="")
    writer = csv.writer(tmp)
    writer.writerow(["basin_id", "x", "y", "flow_acc_val"])
    writer.writerows(records)
    tmp.close()

    gs.run_command(
        "v.in.ascii",
        input=tmp.name,
        output=output_vector,
        format="point",
        separator="comma",
        skip=1,
        x=2, y=3,
        cat=1,
        columns="basin_id int, x double, y double, flow_acc_val double",
        overwrite=True
    )
    os.unlink(tmp.name)

    print(f"Pour Points: {len(records)} outlets written to '{output_vector}'.")
    return output_vector, pour_pts
                            
def compute_catchment_area(flow_acc_rast: str,
                           dem_rast: str,
                           output_rast: str = "catchment_area_m2") -> str:
    import grass.script as gs
 
    print("Catchment area: Computing contributing area in m² …")
 
    # Convert the flow-accumulation count to physical area using cell size.
    region = gs.region()
    cell_area = abs(region["nsres"]) * abs(region["ewres"])
    print(f"Catchment area: Cell area = {cell_area:.2f} m²")
 
    gs.run_command(
        "r.mapcalc",
        expr=f"{output_rast} = {flow_acc_rast} * {cell_area}",
        overwrite=True
    )
 
    # Setting a human-readable unit in the raster metadata
    gs.run_command(
        "r.support",
        map=output_rast,
        units="m2",
        description="Specific catchment area (contributing area in square metres)"
    )
 
    print(f"Catchment area: Raster '{output_rast}' created.")
    return output_rast
 
def compute_mws_connectivity(micro_watersheds_rast: str,
                             flow_dir_rast: str,
                             pour_pts: dict) -> tuple:
    import grass.script as gs
    from grass.script import array as garray

    # Build a directed graph of micro-watershed relationships. Each outlet cell defines one edge
    # from a source basin to the downstream basin it drains into, which is later used for
    # visualization and attribute enrichment.
    print("MWS connectivity: building micro-watershed connectivity graph …")

    region  = gs.region()
    nrows   = int(region["rows"])
    ncols   = int(region["cols"])
    n       = region["n"]
    w       = region["w"]
    nsres   = region["nsres"]
    ewres   = region["ewres"]

    DIR_OFFSETS = {
        1: (-1, +1),  # NE
        2: (-1,  0),  # N
        3: (-1, -1),  # NW
        4: ( 0, -1),  # W
        5: (+1, -1),  # SW
        6: (+1,  0),  # S
        7: (+1, +1),  # SE
        8: ( 0, +1),  # E
    }

    gs.run_command("r.mapcalc",
                   expr=f"micro_watersheds_int = int({micro_watersheds_rast})",
                   overwrite=True)
    basin_arr   = garray.array("micro_watersheds_int", null=-9999)
    flowdir_arr = garray.array(flow_dir_rast, null=0)

    valid_mask = (basin_arr > 0) & (basin_arr != -9999)
    rows, cols = np.where(valid_mask)
    basins = basin_arr[rows, cols].astype(np.int32)

    basin_ids = np.unique(basins)
    counts = np.bincount(basins, minlength=int(basin_ids.max()) + 1)
    row_sums = np.bincount(basins, weights=rows.astype(np.float64), minlength=int(basin_ids.max()) + 1)
    col_sums = np.bincount(basins, weights=cols.astype(np.float64), minlength=int(basin_ids.max()) + 1)

    basin_centroids = {}
    for bid in basin_ids.tolist():
        bid = int(bid)
        count = int(counts[bid])
        if count == 0:
            continue
        cx = w + (col_sums[bid] / count + 0.5) * ewres
        cy = n - (row_sums[bid] / count + 0.5) * nsres
        basin_centroids[bid] = (cx, cy)

    edges = {}
    for bid, (pr, pc) in pour_pts.items():
        direction = int(flowdir_arr[pr, pc])
        if direction not in DIR_OFFSETS:
            continue

        dr, dc = DIR_OFFSETS[direction]
        nr, nc = pr + dr, pc + dc

        if not (0 <= nr < nrows and 0 <= nc < ncols):
            continue

        downstream_basin = int(basin_arr[nr, nc])
        if downstream_basin <= 0 or downstream_basin == -9999:
            continue

        if downstream_basin != bid:
            edges[(bid, downstream_basin)] = True

    print(f"MWS connectivity: {len(edges)} directed basin relationships")
    return edges, basin_centroids, basin_ids


def mws_connectivity_geojson(edges: dict) -> dict:
    return {
        "type": "FeatureCollection",
        "features": [
            {
                "type": "Feature",
                "geometry": None,
                "properties": {
                    "from_basin_id": from_id,
                    "to_basin_id": to_id,
                },
            }
            for from_id, to_id in edges
        ],
    }


def compute_catchments_with_stream_order(
    streams_rast: str,
    strahler_rast: str,
    flow_dir_rast: str,
    output_rast: str = "catchment_stream_order"
) -> str:
  
    import grass.script as gs
    import tempfile, os

    # Align the region with the flow-direction raster so the stream-basin segmentation operates on
    # the same grid as the rest of the drainage analysis.
    gs.run_command("g.region", raster=flow_dir_rast, flags="a")

    seg_basins_rast = "tmp_seg_basins"
    gs.run_command(
        "r.stream.basins",
        direction=flow_dir_rast,
        stream_rast=streams_rast,
        basins=seg_basins_rast,
        overwrite=True
    )

    # Read the stream segment IDs together with their Strahler order so we can rebuild a map in
    # which each segment's basin inherits the order of the stream it drains into.
    raw = gs.read_command(
        "r.stats",
        input=f"{streams_rast},{strahler_rast}",
        flags="cn",
        separator="space"
    )

    seg_to_order: dict[int, int] = {}
    for line in raw.strip().splitlines():
        parts = line.split()
        if len(parts) >= 2:
            seg_to_order[int(parts[0])] = int(parts[1])

    if not seg_to_order:
        raise RuntimeError(
            f"No Strahler order values found — check that '{streams_rast}' "
            f"and '{strahler_rast}' overlap spatially."
        )

    rules_file = tempfile.NamedTemporaryFile(
        mode="w", suffix=".txt", delete=False
    )
    for seg_id, order in seg_to_order.items():
        rules_file.write(f"{seg_id} = {order}\n")
    rules_file.write("* = NULL\n")
    rules_file.close()

    gs.run_command(
        "r.reclass",
        input=seg_basins_rast,
        output=output_rast,
        rules=rules_file.name,
        overwrite=True
    )
    os.unlink(rules_file.name)

    gs.run_command("g.remove", type="raster", name=seg_basins_rast, flags="f")

    gs.run_command(
        "r.support",
        map=output_rast,
        title="Catchment stream order",
        description="Strahler order of the stream segment each cell drains into"
    )

    print(f"'{output_rast}' ready: each cell = Strahler order of its draining segment.")
    return output_rast

def filter_vector_file_by_boundary(vector_path: Path, boundary_path: str) -> None:
    import geopandas as gpd

    # GRASS can export zero-length stream segments as one-coordinate LineStrings.
    # Pyogrio raises while decoding those geometries unless invalid features are ignored.
    vectors = gpd.read_file(vector_path, on_invalid="ignore")
    boundary = gpd.read_file(boundary_path)
    if vectors.crs is None or boundary.crs is None:
        raise ValueError(
            f"Both vector output and selection boundary need a CRS: {vector_path}"
        )
    if boundary.empty:
        raise ValueError(f"Selection boundary is empty: {boundary_path}")

    selection_geometry = boundary.to_crs(vectors.crs).geometry.union_all()
    invalid_count = int(vectors.geometry.isna().sum())
    if invalid_count:
        print(
            f"Warning: dropping {invalid_count} invalid vector geometries "
            f"from {vector_path}"
        )
    vectors = vectors[vectors.geometry.notna()].copy()
    vectors = vectors[vectors.geometry.intersects(selection_geometry)]
    vectors.to_file(vector_path, driver="GeoJSON")


def export_outputs(
    output_dir,
    rasters_to_export: dict,
    vectors_to_export: dict,
    clip_boundary_path: str | None = None,
    geojsons_to_export: dict | None = None,
):
    import grass.script as gs

    # Export the derived products to GeoTIFF and GeoJSON so the final hydrologic outputs can be
    # inspected externally in GIS tools, notebooks, or downstream analysis scripts.
    for name, raster_spec in rasters_to_export.items():
        if isinstance(raster_spec, dict):
            raster_name = raster_spec["raster"]
            raster_type = raster_spec.get("type", "Float32")
        elif isinstance(raster_spec, tuple):
            raster_name, raster_type = raster_spec
        else:
            raster_name = raster_spec
            raster_type = "Float32"

        output_path = Path(output_dir) / f"{name}.tif"
        gs.run_command(
            "r.out.gdal",
            input=raster_name,
            output=str(output_path),
            format="GTiff",
            type=raster_type,
            createopt="COMPRESS=LZW,PHOTOMETRIC=MINISBLACK",
            flags="f",
            overwrite=True,
        )
        print(f"Exported raster: {output_path} ({raster_type})")

    for name, (vector, geom_type) in vectors_to_export.items():
        output_path = Path(output_dir) / f"{name}.geojson"
        gs.run_command("v.out.ogr",
                   input=vector,
                   output=str(output_path),
                   format="GeoJSON",
                   type=geom_type,
                   overwrite=True)
        if clip_boundary_path is not None:
            filter_vector_file_by_boundary(output_path, clip_boundary_path)
        print(f"Exported vector: {output_path}")

    for name, geojson in (geojsons_to_export or {}).items():
        output_path = Path(output_dir) / f"{name}.geojson"
        with open(output_path, "w") as f:
            json.dump(geojson, f, indent=2)
        print(f"Exported GeoJSON attributes: {output_path}")
    return


def update_watershed_attributes(gs, edges, basin_centroids, basin_ids):
    from collections import defaultdict

    downstream_map = {from_id: to_id for from_id, to_id in edges}
    upstream_map = defaultdict(list)
    for from_id, to_id in edges:
        upstream_map[to_id].append(from_id)

    statements = []
    for bid in map(int, basin_ids):
        assignments = []
        downstream_id = downstream_map.get(bid)
        if downstream_id is not None:
            assignments.append(f"downstream_id = {int(downstream_id)}")

        upstream_ids = upstream_map.get(bid, [])
        if upstream_ids:
            upstream_value = ",".join(map(str, map(int, upstream_ids)))
            assignments.append(f"upstream_ids = '{upstream_value}'")

        if (downstream_id is not None
                and bid in basin_centroids
                and downstream_id in basin_centroids):
            dx = basin_centroids[downstream_id][0] - basin_centroids[bid][0]
            dy = basin_centroids[downstream_id][1] - basin_centroids[bid][1]
            bearing = round(math.degrees(math.atan2(dx, dy)) % 360, 2)
            assignments.append(f"flow_direction = {bearing}")

        if assignments:
            statements.append(
                "UPDATE watersheds_vect SET "
                + ", ".join(assignments)
                + f" WHERE basin_id = {bid};"
            )

    if statements:
        sql = "BEGIN TRANSACTION;\n" + "\n".join(statements) + "\nCOMMIT;\n"
        gs.write_command("db.execute", input="-", stdin=sql)


# Main function
def main():
    total_start_time = time.perf_counter()
    stage_timings = {}
    with timed("CLI and output directory setup", stage_timings):
        args = parse_args()
        os.makedirs(args.output, exist_ok=True)

    with timed("Input boundary and CRS", stage_timings):
        # Load the watershed boundary and establish the local projected CRS used for DEM and flow calculations.
        watershed_gdf = gpd.read_file(args.shp)
        watershed_gdf = watershed_gdf.to_crs(epsg=4326)

        watershed_gdf.plot(color='white', edgecolor='gray', figsize=(15,12))
        plt.title("Input Watershed Boundary")
        plt.savefig(Path(args.output) / "watershed_boundary.png", dpi=300, bbox_inches='tight')
        # plt.show()

        bbox = tuple(watershed_gdf.dissolve().total_bounds)
        epsg = get_utm_epsg_for_bbox(bbox)
        watershed_utm = watershed_gdf.to_crs(epsg=epsg).dissolve()
        watershed_utm_path = Path(args.output) / "watershed_utm.shp"
        watershed_utm.to_file(watershed_utm_path)
        analysis_buffer = watershed_utm.copy()
        analysis_buffer.geometry = analysis_buffer.geometry.buffer(
            BUFFER_DISTANCE_METERS
        )
        analysis_buffer_utm_path = Path(args.output) / "analysis_buffer_utm.shp"
        analysis_buffer.to_file(analysis_buffer_utm_path)
        bbox = tuple(analysis_buffer.to_crs(epsg=4326).total_bounds)

    with timed("DEM download and preview", stage_timings):
        # Download a DEM covering the buffered watershed footprint.
        location_of_dem = Path(args.output) / "dem_raw.tif"
        location_of_dem = location_of_dem.resolve()

        dem_downloader = DEMDownloader(
            dem_key = CONFIG["DEM"],
            south = bbox[1],
            north = bbox[3],
            west = bbox[0],
            east = bbox[2],
            output = location_of_dem,
            tile_deg = 1.0,
            max_workers = 6,
            max_retries = 5,
        )
        dem_downloader.run()

        fig, ax = plt.subplots(figsize=(10,8))

        import rasterio
        from rasterio.plot import show
        rasterin = rasterio.open(location_of_dem)

        dem_label = CONFIG["DEM"].upper()
        show(rasterin, ax=ax, cmap='terrain', title=f'{dem_label}')
        watershed_gdf.boundary.plot(color="black", ax=ax, linewidth=1.5)
        plt.colorbar(ax.images[0], ax=ax, label='Elevation(m)')
        plt.savefig(Path(args.output) / "dem_with_watershed.png", dpi=300, bbox_inches='tight')
        plt.close()
    
    with timed("GRASS setup and DEM preparation", stage_timings):
        name_of_proj = Path(args.output).resolve().name
        # Initialize GRASS in the same projected CRS as the DEM so all raster operations share a
        # consistent spatial reference and cell geometry.
        session = setup_grass_session(args.grassdb, epsg, name_of_proj)
    
        import grass.script as gs

        subprocess.run([
        "gdalwarp",
        "-t_srs", f"EPSG:{epsg}",
        str(location_of_dem),
        str(Path(args.output) / f"dem_{epsg}.tif")
        ])

        input_dem = str(Path(args.output) / f"dem_{epsg}.tif")
        output_dir = Path(args.output).resolve()
        dem_conditioned = dem_preprocessing(
            input_dem,
            output_dir,
            watershed_path=str(analysis_buffer_utm_path),
            depression_boundary_path=str(watershed_utm_path),
        )

        gs.run_command("r.in.gdal",
                input=input_dem,
                output="dem_utm",
                overwrite=True)
    
        gs.run_command("r.in.gdal",
                input=str(Path(args.output) / "dem_conditioned.tif"),
                output="dem_conditioned",
                overwrite=True)
    

        gs.run_command("v.in.ogr",
                    input=str(analysis_buffer_utm_path),
                    output="analysis_buffer",
                    overwrite=True)
        gs.run_command(
            "v.in.ogr",
            input=str(watershed_utm_path),
            output="watershed",
            overwrite=True,
        )

        # Align the computational region to the buffered polygon and conditioned DEM grid.
        gs.run_command(
            "g.region",
            vector="analysis_buffer",
            align="dem_conditioned",
            flags="a",
        )
    
        try:
            gs.run_command("r.mask", flags="r")
        except:
            pass
        gs.run_command("r.mask",
                vector="analysis_buffer")


    
        print("DEM imported into GRASS and region set to the buffered analysis area.")

        # Plotting the UTM DEM helps confirm the imported extent, grid alignment, and boundary mask.
        fig, ax = plt.subplots(figsize=(10, 8))

        dem_path = Path(args.output) / f"dem_{epsg}.tif"
        with rasterio.open(dem_path) as src:
            show(src, ax=ax, cmap='terrain')

        plt.title("DEM in GRASS GIS (UTM)")
        plt.axis('off')
        plt.close()

    
    with timed("Flow, stream, and catchment derivation", stage_timings):
        # Compute the flow direction, flow accumulation, and initial micro-watershed partition.
        # These rasters are the foundation for stream extraction and downstream basin merging.
        flow_accumulation, flow_dir_ws, micro_watersheds = calculate_flow_accumulation("dem_conditioned",
                                                                                    hyperparam_threshold=HYPER_PARAM)

        flow_dir_st = "flow_dir_st"
        # Extract the stream network from the current DEM and accumulation surface. The threshold
        # controls how much upstream contributing area is needed before a channel is considered a stream.
        gs.run_command("r.stream.extract",
                elevation="dem_conditioned",
                accumulation=flow_accumulation,
                direction=flow_dir_st,
                stream_raster="streams_rast",
                stream_vector="streams_vect",
                threshold=args.threshold,
                overwrite=True)

        # Derive stream ordering from the extracted network so each channel segment carries a
        # Strahler/Shreve designation that can be mapped back to catchments and exported downstream.
        gs.run_command("r.stream.order",
                stream_rast="streams_rast",
                direction=flow_dir_st,
                elevation="dem_conditioned",
                accumulation=flow_accumulation,
                strahler="strahler_order",
                shreve="shreve_order",
                stream_vect="streams_with_order",
                overwrite=True)
  
        catchment_order_rast = compute_catchments_with_stream_order(
        streams_rast="streams_rast",
        strahler_rast="strahler_order",
        flow_dir_rast=flow_dir_st,
        output_rast="catchment_stream_order")
    
    with timed("Watershed delineation and connectivity", stage_timings):
        micro_watersheds = merge_small_watersheds(micro_watersheds_rast=micro_watersheds,
                                                  flow_acc_rast=flow_accumulation,
                                                  flow_dir_rast=flow_dir_ws,
                                                  epsg=epsg,
                                                  min_area_ha=args.min_watershed_size)

        # Join stream_type and type_code from streams_vect into streams_with_order
        # Both vectors share 'cat' as the common key
        gs.run_command(
                "v.db.join",
                map="streams_with_order",
                column="cat",
                other_table="streams_vect",
                other_column="cat",
                subset_columns="stream_type,type_code",
                )
        all_columns = gs.read_command(
                "v.info", map="streams_with_order", flags="c"
                ).strip().splitlines()

        keep = {"cat", "stream_type", "type_code", "network", "strahler", "next_stream", "prev_str01", "prev_str02"}

        drop_cols = [
        line.split("|")[1]
        for line in all_columns
        if "|" in line and line.split("|")[1] not in keep]

        if drop_cols:
            gs.run_command(
                "v.db.dropcolumn",
                map="streams_with_order",
                columns=",".join(drop_cols)
            )
            print(f"Dropped columns: {drop_cols}")


        pour_points_vect, pour_points = compute_pour_points(
            micro_watersheds_rast=micro_watersheds,
            flow_acc_rast=flow_accumulation,
            flow_dir_rast=flow_dir_ws,
            output_vector="pour_points"
        )

        catchment_area_rast = compute_catchment_area(
            flow_acc_rast=flow_accumulation,
            dem_rast="dem_conditioned",
            output_rast="catchment_area_m2"
        )


        gs.run_command("r.to.vect",
                input=micro_watersheds,
                output="watersheds_vect_raw",
                type="area",
                overwrite=True)

        # Remove tiny sliver polygons created during raster-to-vector conversion.
        region = gs.region()
        cell_area_m2 = abs(float(region["ewres"]) * float(region["nsres"]))
        sliver_area_m2 = max(cell_area_m2 * 3.0, args.min_watershed_size * 10000 * 0.05)
        gs.run_command(
            "v.clean",
            input="watersheds_vect_raw",
            output="watersheds_vect",
            type="area",
            tool="rmarea",
            threshold=str(sliver_area_m2),
            overwrite=True,
        )

        edges, basin_centroids, basin_ids = compute_mws_connectivity(
            micro_watersheds_rast=micro_watersheds,
            flow_dir_rast=flow_dir_ws,
            pour_pts=pour_points,
        )
        connectivity_geojson = mws_connectivity_geojson(edges)
    
    with timed("Attribute updates and output export", stage_timings):
        gs.run_command(
            "v.db.addcolumn",
            map="watersheds_vect",
            columns="basin_id int, downstream_id int, upstream_ids varchar(256), flow_direction double precision")
    
        # Each vector feature is assigned the basin ID and then enriched with the downstream and
        # upstream relationships from the connectivity graph. This makes the vector dataset queryable
        # and suitable for reporting or map styling.
        gs.run_command("v.db.update", map="watersheds_vect", column="basin_id", query_column="value")


        update_watershed_attributes(gs, edges, basin_centroids, basin_ids)
        print("Watershed attributes updated: downstream_id, upstream_ids, flow_direction -> 'watersheds_vect'")

        gs.run_command("r.mask", flags="r")
        gs.run_command(
            "g.region",
            vector="watershed",
            align="dem_conditioned",
            flags="a",
        )
        gs.run_command("r.mask", vector="watershed")
        rasters_to_export = {
                "flow_direction":       (flow_dir_ws, "Float32"),
                "flow_accumulation":    (flow_accumulation, "Float32"),
                "stream_order":         ("strahler_order", "Int32"),
                "catchment_area_m2":    (catchment_area_rast, "Float32"),
                "catchment_stream_order": (catchment_order_rast, "Int32"),
                "flow_direction_stream": (flow_dir_st, "Float32"),
                "streams_raster":       ("streams_rast", "Int32"),
            }
        vectors_to_export = {
                "streams":              ("streams_with_order", "line"),
                "pour_points":          (pour_points_vect,  "point"),
                "microwatersheds":      ("watersheds_vect", "area"),
            }
        export_outputs(
            args.output,
            rasters_to_export,
            vectors_to_export,
            clip_boundary_path=str(watershed_utm_path),
            geojsons_to_export={"mws_connectivity": connectivity_geojson},
        )

    with timed("GRASS session shutdown", stage_timings):
        if session:
            session.close()

    elapsed_time = time.perf_counter() - total_start_time
    stage_width = max(len("Stage"), *(len(stage) for stage in stage_timings))
    duration_width = max(
        len("Duration"),
        *(len(format_duration(seconds)) for seconds in stage_timings.values()),
        len(format_duration(elapsed_time)),
    )
    print("\nStage timing summary")
    print(f"{'Stage':<{stage_width}}  {'Duration':>{duration_width}}")
    print(f"{'-' * stage_width}  {'-' * duration_width}")
    for stage, seconds in stage_timings.items():
        print(f"{stage:<{stage_width}}  {format_duration(seconds):>{duration_width}}")
    print(f"{'Total':<{stage_width}}  {format_duration(elapsed_time):>{duration_width}}")

if __name__ == "__main__":
    main()
