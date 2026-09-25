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
    3. Condition the DEM for hydrologic processing and limit analysis to the polygon.
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
from dem_downloader import DEMDownloader
import time
import richdem as rd

# Minimum area used when merging tiny micro-watersheds; expressed in hectares.
MIN_WATERSHED_SIZE = 500
# Initial r.watershed accumulation threshold in cells for preliminary basin delineation.
# This should later be tied to DEM resolution and target watershed scale.
HYPER_PARAM = 1200

CONFIG = {
    "DEM": "srtm30"  # DEM key used by DEMDownloader; see dem_downloader.py for registry details.
}

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
        help="Minimum watershed size in cells (default: 500 hectares)"
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
        # A single UTM zone is returned for the full AOI; this is sufficient for
        # watershed-scale work, but broader areas may need a tiled strategy.
        return utm_info[0].code
    else:
        raise ValueError("No suitable UTM CRS found for the given AOI.")

# Condition the DEM for hydrologic processing by removing sinks, recording depression depth,
# and reducing the analysis extent to the watershed polygon when one is supplied. This creates
# a terrain surface that is suitable for D8 flow routing and stream extraction.
def dem_preprocessing(input_dem, output_dir, watershed_path=None):
    from rasterio import features as rasterio_features

    conditioned_dem = str(Path(output_dir) / "dem_conditioned.tif")
    depression_tif = str(Path(output_dir) / "natural_depressions.tif")

    with rasterio.open(input_dem) as src:
        nodata = src.nodata
        transform = src.transform
        crs = src.crs

    if nodata is None:
        nodata = -9999.0

    dem_in = rd.LoadGDAL(input_dem, no_data=float(nodata))
    print(dem_in)

    filled_dem = rd.FillDepressions(dem_in, in_place=False)

    dem_arr = np.asarray(dem_in)
    filled_arr = np.asarray(filled_dem)
    depth = filled_arr - dem_arr # The depressions (natural and artificial) will appear as positive values. Inverted-depth map.

    # Restrict downstream calculations and exports to the watershed polygon when provided.
    watershed_mask = np.ones(dem_arr.shape, dtype=bool)
    if watershed_path is not None:
        watershed_gdf = gpd.read_file(watershed_path)
        if not watershed_gdf.empty:
            if watershed_gdf.crs is not None and crs is not None:
                watershed_gdf = watershed_gdf.to_crs(crs)
            watershed_mask = rasterio_features.geometry_mask(
                [geom for geom in watershed_gdf.geometry if geom is not None],
                out_shape=dem_arr.shape,
                transform=transform,
                invert=True,
                all_touched=False,
            )

    nodata_mask = dem_arr == dem_in.no_data
    valid_mask = (~nodata_mask) & watershed_mask
    depth[~valid_mask] = dem_in.no_data
    depth[valid_mask] = np.clip(depth[valid_mask], 0, None)

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
    epsg: int,
    output_rast: str = "micro_watersheds",
) -> str:
    """
    Merge tiny micro-watersheds upward until each basin reaches the minimum area.

    The merge logic follows the original hydrologic pattern: a small basin is
    absorbed into its downstream neighbour when possible, otherwise into the
    largest upstream neighbour that still exists in the current topology.
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

    # Pull the current GRASS region geometry and cell size to convert raster area to hectares.
    region    = gs.region()
    nrows     = int(region["rows"])
    ncols     = int(region["cols"])
    cell_ha   = abs(region["nsres"]) * abs(region["ewres"]) / 10_000.0
    min_cells = min_area_ha / cell_ha

    print(
        f"[merge_small_watersheds] cell = {cell_ha:.4f} ha  "
        f": min_cells = {min_cells:.1f}"
    )

    # Load the generated micro-watershed raster and the supporting drainage rasters.
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

    # Build the drainage topology in a vectorized way so the merge pass can process
    # many basins without the original basin-by-basin loop. The resulting dictionaries
    # describe each basin's outlet, area, and upstream neighbours, which are the pieces
    # needed to decide how to absorb tiny basins without breaking the D8 drainage graph.
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

        # Count the cells in each basin using np.bincount. This is faster and more stable than
        # repeatedly summing each basin's footprint, and it gives the area information required
        # for the minimum-size merge threshold.
        flat_valid = basin_arr[valid_mask]
        if flat_valid.size == 0:
            return {}, {}, {}
        counts     = np.bincount(flat_valid)
        valid_ids  = np.flatnonzero(counts).astype(np.int32)
        area_cells = {int(i): int(counts[i]) for i in valid_ids}

        # Map each D8 flow-direction cell to the neighbour it drains toward. The offset table is
        # used to compute the candidate downstream cell for every valid raster cell in one pass.
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

        # A cell is considered a candidate outlet when it drains outside the grid or into a
        # different basin. These are the natural pour points used to identify each basin's link
        # to the downstream network.
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

        # Choose the single strongest outlet per basin by sorting each candidate outlet in
        # descending accumulation. The first occurrence per basin_id therefore corresponds to the
        # dominant pour point, which is the correct outlet for connectivity and merge decisions.
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

    # Repeat until all remaining basins satisfy the minimum area threshold or cannot be merged
    # without violating the drainage topology. This is the iterative cleanup step that turns the
    # initial fine-grained basin partition into a more hydrologically appropriate watershed map.
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

    # Renumber the final basin IDs to a compact, sequential set for downstream GRASS processing.
    # A compact ID set is simpler to handle in the later raster/vector conversion and makes the
    # basin graph easier to inspect and attribute in GRASS tables.
    final_ids = np.unique(basin_arr)
    final_ids = final_ids[(final_ids > 0) & (final_ids != -9999)]
    new_id = 1
    for old_id in sorted(final_ids.tolist()):
        basin_arr[basin_arr == old_id] = -(new_id)
        new_id += 1
    basin_arr = -basin_arr
    basin_arr[basin_arr == 9999] = -9999

    # Write the merged basin raster back into GRASS so the rest of the pipeline can work
    # with a clean, compact watershed partition instead of the original fine-grained micro-basins.
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
                             pour_pts: dict,
                             output_geojson: Path) -> tuple:
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

def export_outputs(output_dir, rasters_to_export: dict, vectors_to_export: dict):
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
        print(f"Exported vector: {output_path}")
    return

# Main function
def main():
    args = parse_args()

    os.makedirs(args.output, exist_ok=True)

    # Load the watershed boundary and establish the local projected CRS used for DEM and flow calculations.
    watershed_gdf = gpd.read_file(args.shp)
    watershed_gdf = watershed_gdf.to_crs(epsg=4326)  

    watershed_gdf.plot(color='white', edgecolor='gray', figsize=(15,12))
    plt.title("Input Watershed Boundary")
    plt.savefig(Path(args.output) / "watershed_boundary.png", dpi=300, bbox_inches='tight')
    # plt.show()
    

    minx, miny, maxx, maxy = watershed_gdf.dissolve().total_bounds
    buffer = 0.1 # Nearly 10 km buffer at regions near equator(in degrees)
    bbox = (minx -buffer, miny -buffer, maxx +buffer, maxy +buffer)
    epsg = get_utm_epsg_for_bbox(bbox)

    # Download a DEM that covers the buffered watershed footprint. Using the polygon bounds and a
    # small margin avoids edge artifacts while keeping the raster manageable at the watershed scale.
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
    watershed_utm_path = Path(args.output) / "watershed_utm.shp"
    watershed_gdf.to_crs(epsg=epsg).to_file(watershed_utm_path)
    dem_conditioned = dem_preprocessing(
        input_dem,
        output_dir,
        watershed_path=str(watershed_utm_path),
    )

    gs.run_command("r.in.gdal",
               input=input_dem,
               output="dem_utm",
               overwrite=True)
    
    gs.run_command("r.in.gdal",
               input=str(Path(args.output) / "dem_conditioned.tif"),
               output="dem_conditioned",
               overwrite=True)
    

    # Set the GRASS computational region to the conditioned DEM so all vector and raster tools
    # operate on the same grid resolution and spatial extent.
    gs.run_command("g.region", raster="dem_conditioned", flags="p")

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

    # Plotting the UTM DEM helps confirm the imported extent, grid alignment, and boundary mask.
    fig, ax = plt.subplots(figsize=(10, 8))

    dem_path = Path(args.output) / f"dem_{epsg}.tif"
    with rasterio.open(dem_path) as src:
        show(src, ax=ax, cmap='terrain')

    plt.title("DEM in GRASS GIS (UTM)")
    plt.axis('off')
    plt.close()

    
    # Compute the flow direction, flow accumulation, and initial micro-watershed partition.
    # These rasters are the foundation for stream extraction and downstream basin merging.
    flow_accumulation, flow_dir_ws, micro_watersheds = calculate_flow_accumulation("dem_conditioned",
                                                                                   hyperparam_threshold=HYPER_PARAM)

    flow_dir_st = "flow_dir_st"
    gs.run_command("r.mask", flags="r")
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
    
    # Reapply the watershed mask before basin merging so the merge logic operates only within
    # the study polygon and does not join cells outside the hydrologic unit.
    gs.run_command("r.mask", vector="watershed")
  
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
        output_geojson=Path(args.output) / "mws_connectivity.geojson"
    )
    
    gs.run_command(
        "v.db.addcolumn",
        map="watersheds_vect",
        columns="basin_id int, downstream_id int, upstream_ids varchar(256), flow_direction double precision")
    
    # Each vector feature is assigned the basin ID and then enriched with the downstream and
    # upstream relationships from the connectivity graph. This makes the vector dataset queryable
    # and suitable for reporting or map styling.
    gs.run_command("v.db.update", map="watersheds_vect", column="basin_id", query_column="value")


    from collections import defaultdict
    downstream_map = {from_id: to_id for (from_id, to_id) in edges}
    upstream_map   = defaultdict(list)
    for from_id, to_id in edges:
        upstream_map[to_id].append(from_id)

    
    for bid in map(int, basin_ids):
        # For each basin, attach the downstream basin, the sequence of upstream neighbours, and a
        # bearing angle derived from the centroid-to-centroid direction. These are useful for
        # hydrologic interpretation and visual exploration.
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
    export_outputs(args.output, rasters_to_export, vectors_to_export)

    if session:
        session.close()

if __name__ == "__main__":
    start_time = time.perf_counter()
    main()

    end_time = time.perf_counter()
    elapsed_time = end_time - start_time
    print(f"Total execution time: {elapsed_time // 3600} hrs {(elapsed_time % 3600) // 60} mins {(elapsed_time % 60):.2f} secs")
