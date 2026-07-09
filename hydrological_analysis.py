'''
GRASS GIS based hydrological analysis
=====================================
Author: Shubham Kumar
Date: April 2026
-------------------------------------
Usage: 
    python hydrological_analysis.py --shp path_to_shapefile \
        --output path_to_output_directory \
        --threshold flow_accumulation_threshold

Description:
    Given a shapefile of a hydrological unit (e.g., a watershed), this script performs the following steps:
    1. Sets up a GRASS GIS environment.
    2. Imports the shapefile and a Digital Elevation Model (DEM) into GRASS GIS.
    3. Fills sinks in the DEM to ensure proper flow direction.
    4. Calculates flow direction and flow accumulation.
    5. Extracts stream networks based on a specified flow accumulation threshold.
    6. Exports the resulting stream network as a GeoJSON file for visualization and further analysis.

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
import whitebox_workflows as wbw
from dem_downloader import DEMDownloader
import time

MIN_WATERSHED_SIZE = 500 # 500 hectares
HYPER_PARAM = 1200 # Initial threshold for r.watershed. 1200 cells. 

CONFIG = {
    "DEM": "srtm30"
}

def parse_args():
    parser = argparse.ArgumentParser(
        description="GRASS GIS Hydrological Analysis Pipeline"
    )
    parser.add_argument(
        "--shp", required=True,
        help="Path to input shapefile defining the watershed boundary"
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
        help="Minimum watershed size in cells (default: 500 hectares)"
    )

    return parser.parse_args()

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

    # Trying grass_session first, fall back to gsetup
    try:
        from grass_session import Session
        session = Session()
        session.open(gisdb=str(grassdb_path), location=location, mapset=mapset)
        print("grass_session initialised successfully.")
        return session
    except ImportError:
        pass

    # Fallback: manual gsetup
    import grass.script.setup as gsetup
    gsetup.init(str(grassdb_path), location, mapset)
    print("GRASS environment initialised via grass.script.setup")
    return None

def get_utm_epsg_for_bbox(bbox):
    west, south, east, north = bbox
    
    from pyproj.aoi import AreaOfInterest
    from pyproj.database import query_utm_crs_info
    
    utm_info = query_utm_crs_info(
        datum_name="WGS 84",
        area_of_interest=AreaOfInterest(
            west_lon_degree=west,
            south_lat_degree=south,
            east_lon_degree=east,
            north_lat_degree=north
        )
    )

    if utm_info:
        return utm_info[0].code
    else:
        raise ValueError("No suitable UTM CRS found for the given AOI.")

def dem_preprocessing(input_dem, output_dir):
    breached_dem = str(Path(output_dir) / "dem_breached.tif")
    conditioned_dem = str(Path(output_dir) / "dem_conditioned.tif")

    print("DEM Conditioning...")

    wbe = wbw.WbEnvironment()

    dem_in = wbe.read_raster(input_dem)
    breached = wbe.hydrology.depressions_storage.breach_depressions_least_cost(
        dem=dem_in,
        max_dist=100, 
        fill_deps=False
    )
    wbe.write_raster(breached, breached_dem)

    conditioned = wbe.hydrology.depressions_storage.fill_depressions(
        dem=breached,
        flat_increment=0.001
    )
    wbe.write_raster(conditioned, conditioned_dem)

    print(f"DEM conditioning complete: {conditioned_dem}")

    return conditioned

def fill_sinks(dem_raster):
    print("Filling sinks in DEM...")
    import grass.script as gs
    filled_dem = "dem_filled"
    flow_dir = "flow_dir"
    gs.run_command("r.fill.dir", input=dem_raster, output=filled_dem, direction=flow_dir, overwrite=True)
    return filled_dem, flow_dir

def natural_depressions(input_dem, output_dir):
    wbe = wbw.WbEnvironment()
    dem_in = wbe.read_raster(input_dem)
    depth_raster = wbe.hydrology.depth_in_sink(dem=dem_in)
    wbe.write_raster(depth_raster, str(Path(output_dir) / "natural_depressions.tif"))
    return None

def calculate_flow_accumulation(dem_filled, hyperparam_threshold):
    import grass.script as gs
    gs.run_command("r.watershed", 
                elevation=dem_filled, 
                accumulation="flow_acc", 
                drainage="flow_dir_watershed",
                threshold=hyperparam_threshold, 
                basin="micro_watersheds",
                #stream="streams_raw",
                flags="as",          # -a: positive accumulation; -s: single-flow (D8) 
                overwrite=True)

    return "flow_acc", "flow_dir_watershed", "micro_watersheds"

def merge_small_watersheds(
    micro_watersheds_rast: str,
    flow_acc_rast: str,
    flow_dir_rast: str,
    min_area_ha: float,
    epsg: int,
    output_rast: str = "micro_watersheds",
) -> str:
    """
    Hybrid implementation: hydrologically correct merge logic from the
    original, with O(n) vectorized topology construction replacing the
    O(n*b) per-basin np.where() loop.

    Changes from the original
    ─────────────────────────
    1. _build_maps replaced by _build_maps_vectorized   → O(n) per call
       instead of O(n*b).  All merge logic is UNCHANGED.
    2. unmergeable is reset every iteration so basins that gain a neighbour
       after an earlier merge are retried.
    3. upstream[C] deduplication prevents repeated IDs from accumulating
       when multiple basins merge into the same target in one pass.
    4. area_cells built with np.bincount instead of repeated np.sum.

    No Union-Find.  No deferred array writes.  Merge direction is always
    preserved: bid (small) → target (downstream or largest upstream).

    Complexity
    ──────────
    Per iteration : O(n + b log b)   [vectorized topology + sort]
    Total         : O(k * (n + b log b))  where k ≤ b iterations
    """
    import numpy as np
    import grass.script as gs
    from grass.script import array as garray

    # D8 flow-direction offsets (GRASS r.watershed encoding)
    DIR_OFFSETS = {
        1: (-1, +1),   # NE
        2: (-1,  0),   # N
        3: (-1, -1),   # NW
        4: ( 0, -1),   # W
        5: (+1, -1),   # SW
        6: (+1,  0),   # S
        7: (+1, +1),   # SE
        8: ( 0, +1),   # E
    }

    # ── GRASS region metadata ──────────────────────────────────────────────
    region    = gs.region()
    nrows     = int(region["rows"])
    ncols     = int(region["cols"])
    cell_ha   = abs(region["nsres"]) * abs(region["ewres"]) / 10_000.0
    min_cells = min_area_ha / cell_ha

    print(
        f"[merge_small_watersheds] cell = {cell_ha:.4f} ha  "
        f": min_cells = {min_cells:.1f}"
    )

    # ── Load rasters ───────────────────────────────────────────────────────
    gs.run_command(
        "r.mapcalc",
        expr=f"_mws_work = int({micro_watersheds_rast})",
        overwrite=True,
    )
    basin_arr   = garray.array("_mws_work",   null=-9999).astype(np.int32)
    acc_arr     = garray.array(flow_acc_rast, null=-1  ).astype(np.float32)
    flowdir_arr = garray.array(flow_dir_rast, null=0   ).astype(np.int8)

    # Pre-build the D8 offset look-up table once; reused every iteration.
    _offset_lut = np.zeros((9, 2), dtype=np.int32)   # index 0 → no-op
    for _d, (_dr, _dc) in DIR_OFFSETS.items():
        _offset_lut[_d] = (_dr, _dc)

    # ── Vectorized topology builder  O(n) ─────────────────────────────────
    def _build_maps_vectorized(basin_arr):
        """
        Fully vectorized replacement for the original _build_maps().

        Returns identical structures:
            downstream : dict  basin_id -> downstream_basin_id | None
            area_cells : dict  basin_id -> cell count
            upstream   : dict  basin_id -> list[basin_id]

        Key properties preserved
        ────────────────────────
        • Pour point = boundary cell with MAXIMUM flow accumulation whose
          D8 neighbour lies in a different basin (or off-grid).  Identical
          semantics to the original per-cell loop.
        • Off-grid flows are correctly identified as outlets (downstream=None)
          without reading clamped neighbour values (fix for OOB heuristic).
        • area_cells uses np.bincount → O(n) instead of O(n*b).
        """
        valid_mask = (basin_arr > 0) & (basin_arr != -9999)

        # ── area_cells via bincount ────────────────────────────────────────
        flat_valid = basin_arr[valid_mask]
        if flat_valid.size == 0:
            return {}, {}, {}
        counts     = np.bincount(flat_valid)
        valid_ids  = np.flatnonzero(counts).astype(np.int32)
        area_cells = {int(i): int(counts[i]) for i in valid_ids}

        # ── vectorized neighbour coordinates ──────────────────────────────
        R, C  = np.indices((nrows, ncols), dtype=np.int32)
        d_raw = flowdir_arr.astype(np.int32)
        d_clip = np.where((d_raw >= 1) & (d_raw <= 8), d_raw, 0)

        nr = R + _offset_lut[d_clip, 0]
        nc = C + _offset_lut[d_clip, 1]

        has_valid_dir = (d_raw >= 1) & (d_raw <= 8)
        in_bounds     = (nr >= 0) & (nr < nrows) & (nc >= 0) & (nc < ncols)
        is_basin      = valid_mask

        # ── neighbour basin (only read where safe) ─────────────────────────
        # Critically: never read nbr_basin for out-of-bounds cells, so outlet
        # cells are not accidentally matched to a clamped neighbour's basin ID.
        safe_r = np.clip(nr, 0, nrows - 1)
        safe_c = np.clip(nc, 0, ncols - 1)

        nbr_basin = np.full((nrows, ncols), -9999, dtype=np.int32)
        ib_basin  = is_basin & has_valid_dir & in_bounds
        nbr_basin[ib_basin] = basin_arr[safe_r[ib_basin], safe_c[ib_basin]]

        # ── candidate pour-point mask ──────────────────────────────────────
        is_outlet_cell = is_basin & has_valid_dir & (~in_bounds)
        is_cross_basin = (
            is_basin & has_valid_dir & in_bounds
            & (nbr_basin != basin_arr)
            & (nbr_basin > 0)
            & (nbr_basin != -9999)
        )
        candidate = is_outlet_cell | is_cross_basin

        cr, cc      = np.where(candidate)
        if cr.size == 0:
            downstream = {int(i): None for i in valid_ids}
            upstream   = {int(i): []   for i in valid_ids}
            return downstream, area_cells, upstream

        cand_bid    = basin_arr[cr, cc]
        # Use -1 as a sentinel for "drains off grid" (outlet)
        cand_target = np.where(
            is_outlet_cell[cr, cc],
            np.int32(-1),
            nbr_basin[cr, cc],
        )
        cand_acc = acc_arr[cr, cc]

        # ── argmax-per-basin via lexsort  O(c log c) ──────────────────────
        # Sort by (basin_id ASC, accumulation DESC); first occurrence per
        # basin_id is the pour point with the highest accumulation.
        order     = np.lexsort((-cand_acc, cand_bid))
        bid_s     = cand_bid[order]
        tgt_s     = cand_target[order]

        first         = np.empty(len(bid_s), dtype=bool)
        first[0]      = True
        first[1:]     = bid_s[1:] != bid_s[:-1]

        # ── downstream dict ────────────────────────────────────────────────
        downstream = {int(i): None for i in valid_ids}
        for bid, tgt in zip(bid_s[first].tolist(), tgt_s[first].tolist()):
            downstream[int(bid)] = None if tgt == -1 else int(tgt)

        # ── upstream dict (invert downstream) ─────────────────────────────
        upstream = {int(i): [] for i in valid_ids}
        for src, tgt in downstream.items():
            if tgt is not None:
                upstream.setdefault(tgt, []).append(src)

        return downstream, area_cells, upstream

    # ── Iterative merge  (logic IDENTICAL to original) ────────────────────
    iteration    = 0
    total_merged = 0

    while True:
        # Fix 2: reset unmergeable every iteration so basins that gain a
        # neighbour due to an earlier merge are not permanently skipped.
        unmergeable = set()

        downstream, area_cells, upstream = _build_maps_vectorized(basin_arr)

        small_basins = sorted(
            (area, bid)
            for bid, area in area_cells.items()
            if area < min_cells
        )   # ascending by area: smallest first

        if not small_basins:
            print(
                f"[merge_small_watersheds] Done after {iteration} iteration(s), "
                f"{total_merged} merges total.  "
                f"All {len(area_cells)} basins >= {min_area_ha} ha."
            )
            break

        merges_this_pass = 0

        for _area, bid in small_basins:
            # May have been absorbed earlier in this pass
            if bid not in area_cells:
                continue
            # Re-check against the live array (original approach)
            current_area = area_cells[bid]
            if current_area >= min_cells:
                continue

            # ── Choose merge target (unchanged from original) ──────────────
            ds = downstream.get(bid)
            if ds is not None and ds in area_cells:
                target = ds
            else:
                us_ids = [u for u in upstream.get(bid, []) if u in area_cells]
                if not us_ids:
                    unmergeable.add(bid)
                    print(
                        f"[merge_small_watersheds] NOTE: basin {bid} has no "
                        f"valid neighbours – skipping."
                    )
                    continue
                target = max(us_ids, key=lambda u: area_cells[u])

            # ── Apply merge: relabel bid → target in the array ────────────
            # This is the ORIGINAL mechanism. No Union-Find. Merge direction
            # is unambiguous: bid cells become target cells.
            basin_arr[basin_arr == bid] = target
            area_cells[target] = (
                area_cells.get(target, 0) + area_cells.pop(bid, 0)
            )

            # Patch connectivity for remaining basins in this same pass.
            for other, other_ds in list(downstream.items()):
                if other_ds == bid:
                    downstream[other] = target
            for other_ups in upstream.values():
                for i, u in enumerate(other_ups):
                    if u == bid:
                        other_ups[i] = target

            # Move bid's upstream children under target.
            absorbed_us = upstream.pop(bid, [])
            target_ups  = upstream.setdefault(target, [])
            target_ups.extend(u for u in absorbed_us if u != target)
            # Fix 3: deduplicate upstream list to avoid repeated IDs.
            upstream[target] = list(dict.fromkeys(
                u for u in target_ups if u != target
            ))

            downstream.pop(bid, None)

            merges_this_pass += 1
            total_merged      += 1

        iteration += 1

        if merges_this_pass == 0:
            remaining_small = sum(
                1 for bid, area in area_cells.items()
                if area < min_cells and bid not in unmergeable
            )
            if remaining_small:
                print(
                    f"[merge_small_watersheds] NOTE: {remaining_small} basin(s) "
                    f"remain below {min_area_ha} ha but have no valid merge "
                    f"candidate – kept as-is."
                )
            break

    # ── Renumber to compact sequential IDs (original trick) ───────────────
    final_ids = np.unique(basin_arr)
    final_ids = final_ids[(final_ids > 0) & (final_ids != -9999)]
    new_id = 1
    for old_id in sorted(final_ids.tolist()):
        basin_arr[basin_arr == old_id] = -(new_id)
        new_id += 1
    basin_arr = -basin_arr
    basin_arr[basin_arr == 9999] = -9999

    # ── Write result back to GRASS ─────────────────────────────────────────
    import tempfile, os
    import rasterio
    from rasterio.transform import from_bounds

    region    = gs.region()
    transform = from_bounds(
        region["w"], region["s"], region["e"], region["n"],
        ncols, nrows,
    )

    tmp_tif = tempfile.NamedTemporaryFile(suffix=".tif", delete=False)
    tmp_tif.close()

    mask      = (basin_arr > 0) & (basin_arr != -9999)
    write_arr = basin_arr.astype(np.int32)
    write_arr[~mask] = -2147483648

    crs = f"EPSG:{epsg}" if epsg else None
    with rasterio.open(
        tmp_tif.name, "w",
        driver="GTiff",
        height=nrows, width=ncols,
        count=1, dtype="int32",
        crs=crs,
        transform=transform,
        nodata=-2147483648,
    ) as dst:
        dst.write(write_arr, 1)

    gs.run_command(
        "r.in.gdal",
        input=tmp_tif.name,
        output=output_rast,
        overwrite=True,
    )
    os.unlink(tmp_tif.name)
    gs.run_command("g.remove", type="raster", name="_mws_work", flags="f")

    final_count = len(np.unique(
        basin_arr[(basin_arr > 0) & (basin_arr != 9999)]
    ))
    print(
        f"[merge_small_watersheds] '{output_rast}' written "
        f"({final_count} basins finally)."
    )
    return output_rast

def compute_pour_points(micro_watersheds_rast: str,
                        flow_acc_rast: str,
                        flow_dir_rast: str,
                        output_vector: str = "pour_points") -> tuple[str, dict]:
    import grass.script as gs
    from grass.script import array as garray
    import tempfile, csv

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

    gs.run_command("r.mapcalc",
                   expr=f"micro_watersheds_int = int({micro_watersheds_rast})",
                   overwrite=True)
    basin_arr = garray.array("micro_watersheds_int", null=-9999)
    acc_arr = garray.array(flow_acc_rast, null=-1)
    flowdir_arr = garray.array(flow_dir_rast, null=0)

    basin_ids = np.unique(basin_arr)
    basin_ids = basin_ids[(basin_ids > 0) & (basin_ids != -9999)]

    pour_pts = {}  
    records  = []  

    for bid in basin_ids:
        bid = int(bid)

        rows_in, cols_in = np.where(basin_arr == bid)

        best_acc = -np.inf
        best_rc  = None

        for r, c in zip(rows_in, cols_in):
            direction = int(flowdir_arr[r, c])
            if direction not in DIR_OFFSETS:
                if float(acc_arr[r, c]) > best_acc:
                    best_acc = float(acc_arr[r, c])
                    best_rc  = (r, c)
                continue

            dr, dc = DIR_OFFSETS[direction]
            nr, nc = r + dr, c + dc

            outside = not (0 <= nr < nrows and 0 <= nc < ncols)
            if outside:
                if float(acc_arr[r, c]) > best_acc:
                    best_acc = float(acc_arr[r, c])
                    best_rc  = (r, c)
                continue

            downstream_basin = int(basin_arr[nr, nc])

            if downstream_basin != bid and downstream_basin > 0 and downstream_basin != -9999:
                if float(acc_arr[r, c]) > best_acc:
                    best_acc = float(acc_arr[r, c])
                    best_rc  = (r, c)

        # Fallback 
        if best_rc is None:
            flat_idx = int(np.argmax(np.where(basin_arr == bid, acc_arr, -np.inf)))
            best_rc  = (int(flat_idx // ncols), int(flat_idx  % ncols))
            best_acc = float(acc_arr[best_rc])

        pour_pts[bid] = best_rc
        row, col = best_rc
        x = w + (col + 0.5) * ewres
        y = n - (row + 0.5) * nsres
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
 
    # Retrieving current region resolution
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
                             pour_pts: dict,
                             output_geojson: Path) -> tuple:
    import grass.script as gs
    from grass.script import array as garray

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
    flowdir_arr = garray.array(flow_dir_rast,          null=0)

    basin_ids = np.unique(basin_arr)
    basin_ids = basin_ids[(basin_ids > 0) & (basin_ids != -9999)]

    basin_centroids = {}
    for bid in basin_ids:
        bid = int(bid)
        rows_idx, cols_idx = np.where(basin_arr == bid)
        cx = w + (cols_idx.mean() + 0.5) * ewres
        cy = n - (rows_idx.mean() + 0.5) * nsres
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

    features = []
    for (from_id, to_id) in edges:
        if from_id not in basin_centroids or to_id not in basin_centroids:
            continue
        x0, y0 = basin_centroids[from_id]
        x1, y1 = basin_centroids[to_id]
        features.append({
            "type": "Feature",
            "geometry": {
                "type": "LineString",
                "coordinates": [[x0, y0], [x1, y1]]
            },
            "properties": {
                "from_basin_id": from_id,
                "to_basin_id":   to_id
            }
        })

    geojson = {"type": "FeatureCollection", "features": features}
    with open(output_geojson, "w") as f:
        json.dump(geojson, f, indent=2)

    print(f"MWS connectivity: {len(features)} directed edges -> {output_geojson}")
    return edges, basin_centroids, basin_ids

def compute_catchments_with_stream_order(
    streams_rast: str,
    strahler_rast: str,
    flow_dir_rast: str,
    output_rast: str = "catchment_stream_order"
) -> str:
  
    import grass.script as gs
    import tempfile, os

    gs.run_command("g.region", raster=flow_dir_rast)

    seg_basins_rast = "tmp_seg_basins"
    gs.run_command(
        "r.stream.basins",
        direction=flow_dir_rast,
        stream_rast=streams_rast,
        basins=seg_basins_rast,
        overwrite=True
    )

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

def export_outputs(output_dir, rasters_to_export: dict, vectors_to_export: dict):
    import grass.script as gs
    for name, raster in rasters_to_export.items():
        output_path = Path(output_dir) / f"{name}.tif"
        gs.run_command("r.out.gdal",
                   input=raster,
                   output=str(output_path),
                   format="GTiff",
                   type="Float32",
                   createopt="COMPRESS=LZW,PHOTOMETRIC=MINISBLACK",
                   flags="f",
                   overwrite=True)
        print(f"Exported raster: {output_path}")

    for name, (vector, geom_type) in vectors_to_export.items():
        output_path = Path(output_dir) / f"{name}.geojson"
        gs.run_command("v.out.ogr",
                   input=vector,
                   output=str(output_path),
                   format="GeoJSON",
                   type=geom_type,
                   overwrite=True)
        print(f"Exported vector: {output_path}")
    return

# Main function
def main():
    args = parse_args()

    os.makedirs(args.output, exist_ok=True)

    # Reading shapefile and determining UTM zone
    watershed_gdf = gpd.read_file(args.shp)
    watershed_gdf = watershed_gdf.to_crs(epsg=4326)  

    watershed_gdf.plot(color='white', edgecolor='gray', figsize=(15,12))
    plt.title("Input Watershed Boundary")
    plt.savefig(Path(args.output) / "watershed_boundary.png", dpi=300, bbox_inches='tight')
    plt.show()
    

    minx, miny, maxx, maxy = watershed_gdf.dissolve().total_bounds
    buffer = 0.1 # Nearly 10 km buffer at regions near equator(in degrees)
    bbox = (minx -buffer, miny -buffer, maxx +buffer, maxy +buffer)
    epsg = get_utm_epsg_for_bbox(bbox)

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

    dem_label = "FABDEM (30m)" if CONFIG["DEM"] == "FABDEM" else "SRTM1 (30m)"
    show(rasterin, ax=ax, cmap='terrain', title=f'{dem_label}')
    watershed_gdf.boundary.plot(color="black", ax=ax, linewidth=1.5)
    plt.colorbar(ax.images[0], ax=ax, label='Elevation(m)')
    plt.savefig(Path(args.output) / "dem_with_watershed.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    name_of_proj = Path(args.output).resolve().name
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
    dem_conditioned = dem_preprocessing(input_dem, output_dir)

    gs.run_command("r.in.gdal",
               input=input_dem,
               output="dem_utm",
               overwrite=True)
    
    gs.run_command("r.in.gdal",
               input=str(Path(args.output) / "dem_conditioned.tif"),
               output="dem_conditioned",
               overwrite=True)
    

    gs.run_command("g.region", raster="dem_conditioned", flags="p")

    watershed_utm_path = Path(args.output) / "watershed_utm.shp"
    watershed_gdf.to_crs(epsg=epsg).to_file(watershed_utm_path)
    gs.run_command("v.in.ogr",
                   input=str(watershed_utm_path),
                   output="watershed",
                   overwrite=True)
    
    try:
        gs.run_command("r.mask", flags="r")
    except:
        pass
    gs.run_command("r.mask",
               vector="watershed")
    
    print("DEM imported into GRASS and region set to DEM extent.")

    fig, ax = plt.subplots(figsize=(10, 8))

    dem_path = Path(args.output) / f"dem_{epsg}.tif"
    with rasterio.open(dem_path) as src:
        show(src, ax=ax, cmap='terrain')

    plt.title("DEM in GRASS GIS (UTM)")
    plt.axis('off')
    plt.close()

    
    # dem_filled, _ = fill_sinks("dem_utm")
    
    natural_depressions(input_dem, args.output)
    
    
    flow_accumulation, flow_dir_ws, micro_watersheds = calculate_flow_accumulation("dem_conditioned",
                                                                                   hyperparam_threshold=HYPER_PARAM)

    flow_dir_st = "flow_dir_st"
    gs.run_command("r.stream.extract",
               elevation="dem_conditioned",
               accumulation=flow_accumulation,
               direction=flow_dir_st,
               stream_raster="streams_rast",
               stream_vector="streams_vect",
               threshold=args.threshold, 
               overwrite=True)

    gs.run_command("r.stream.order",
               stream_rast="streams_rast",
               direction=flow_dir_st,
               elevation="dem_conditioned",
               accumulation=flow_accumulation,
               strahler="strahler_order",
               shreve="shreve_order",
               stream_vect="streams_with_order",
               overwrite=True)
    '''
    catchment_order_rast = compute_catchments_with_stream_order(
    streams_rast="streams_rast",
    strahler_rast="strahler_order",
    flow_dir_rast=flow_dir_st,
    output_rast="catchment_stream_order")

    '''
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
            output="watersheds_vect",
            type="area",
            overwrite=True)

    edges, basin_centroids, basin_ids = compute_mws_connectivity(
        micro_watersheds_rast=micro_watersheds,
        flow_dir_rast=flow_dir_ws,
        pour_pts=pour_points,     
        output_geojson=Path(args.output) / "mws_connectivity.geojson"
    )
    
    gs.run_command(
        "v.db.addcolumn",
        map="watersheds_vect",
        columns="basin_id int, downstream_id int, upstream_ids varchar(256), flow_direction double precision")
    
    gs.run_command("v.db.update", map="watersheds_vect", column="basin_id", query_column="value")


    from collections import defaultdict
    downstream_map = {from_id: to_id for (from_id, to_id) in edges}
    upstream_map   = defaultdict(list)
    for from_id, to_id in edges:
        upstream_map[to_id].append(from_id)

    
    for bid in map(int, basin_ids):
        ds = downstream_map.get(bid)
        us = upstream_map.get(bid, [])
        bearing = None
        if ds is not None and bid in basin_centroids and ds in basin_centroids:
            dx = basin_centroids[ds][0] - basin_centroids[bid][0]
            dy = basin_centroids[ds][1] - basin_centroids[bid][1]
            bearing = round(math.degrees(math.atan2(dx, dy)) % 360, 2)
        if ds is not None:
            gs.run_command("v.db.update", map="watersheds_vect",
                           column="downstream_id", value=str(ds), where=f"basin_id={bid}")
        if us:
            gs.run_command("v.db.update", map="watersheds_vect",
                           column="upstream_ids", value=",".join(map(str, us)), where=f"basin_id={bid}")
        if bearing is not None:
            gs.run_command("v.db.update", map="watersheds_vect",
                           column="flow_direction", value=str(bearing), where=f"basin_id={bid}")
    print("Watershed attributes updated: downstream_id, upstream_ids, flow_direction -> 'watersheds_vect'")

    gs.run_command("r.mask", flags="r")
    rasters_to_export = {
            "flow_direction":       flow_dir_ws,
            "flow_accumulation":    flow_accumulation,
            #"natural_depressions":  depressions,
            "stream_order":         "strahler_order",
            "catchment_area_m2":    catchment_area_rast,
            #"catchment_stream_order":  catchment_order_rast, 
            "flow_direction_stream": flow_dir_st,
            "streams_raster": "streams_rast" 
        }
    vectors_to_export = {
            "streams":              ("streams_with_order", "line"),
            "pour_points":          (pour_points_vect,  "point"),
            "microwatersheds":      ("watersheds_vect", "area"),
        }
    export_outputs(args.output, rasters_to_export, vectors_to_export)

    if session:
        session.close()

if __name__ == "__main__":
    start_time = time.perf_counter()
    main()

    end_time = time.perf_counter()
    elapsed_time = end_time - start_time
    print(f"Total execution time: {elapsed_time // 3600} hrs {(elapsed_time % 3600) // 60} mins {(elapsed_time % 60):.2f} secs")
