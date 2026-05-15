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
import elevation
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


MIN_WATERSHED_SIZE = 500 # 500 hectares
HYPER_PARAM = 1000 # Initial threshold for r.watershed. 1000 cells. 

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

def fill_sinks(dem_raster):
    print("Filling sinks in DEM...")
    import grass.script as gs
    filled_dem = "dem_filled"
    flow_dir = "flow_dir"
    gs.run_command("r.fill.dir", input=dem_raster, output=filled_dem, direction=flow_dir, overwrite=True)
    return filled_dem, flow_dir

def natural_depressions(dem_filled, dem_original):
    import grass.script as gs
    dep_map = "natural_depressions"
    gs.run_command("r.mapcalc", expr=f"{dep_map} = {dem_filled} - {dem_original}", overwrite=True)
    return dep_map

def calculate_flow_accumulation(dem_filled, hyperparam_threshold, size_threshold):
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
    
    gs.run_command("r.reclass.area",
                input="micro_watersheds",
                output="micro_watersheds_filtered",
                value=size_threshold, 
                mode="lesser",
                method="rmarea",
                overwrite=True)

    gs.run_command("g.remove",
                type="raster",
                name="micro_watersheds",
                flags="f")

    gs.run_command("g.rename",
                raster=("micro_watersheds_filtered", "micro_watersheds"))
    return "flow_acc", "flow_dir_watershed", "micro_watersheds"


def compute_pour_points(micro_watersheds_rast: str,
                        flow_acc_rast: str,
                        output_vector: str = "pour_points") -> str:
    import grass.script as gs
    from grass.script import array as garray

    print("Pour Points Calculation: Locating outlet cell for each micro-watershed.")

    region = gs.region()
    w     = region["w"]
    n     = region["n"]
    ewres = region["ewres"]
    nsres = region["nsres"]

    gs.run_command(
        "r.mapcalc",
        expr="micro_watersheds_int = int(micro_watersheds)",
        overwrite=True
    )

    basin_arr = garray.array("micro_watersheds_int", null=-9999)
    acc_arr   = garray.array(flow_acc_rast, null=-1)

    basin_ids = np.unique(basin_arr)
    basin_ids = basin_ids[(basin_ids > 0) & (basin_ids != -9999)]
                            
    records = []
    for bid in basin_ids:
        mask = basin_arr == bid
        if not mask.any():
            continue
        local_acc = np.where(mask, acc_arr, -np.inf)
        row, col  = np.unravel_index(np.argmax(local_acc), local_acc.shape)
        max_acc   = float(acc_arr[row, col])
    
        x = w + (col + 0.5) * ewres
        y = n - (row + 0.5) * nsres
        records.append((int(bid), x, y, max_acc))

    import tempfile, csv
    tmp = tempfile.NamedTemporaryFile(mode="w", suffix=".csv",
                                      delete=False, newline="")
    writer = csv.writer(tmp)
    writer.writerow(["basin_id", "x", "y", "flow_acc_val"])  # basin_id first
    writer.writerows(records)
    tmp.close()

    gs.run_command(
        "v.in.ascii",
        input=tmp.name,
        output=output_vector,
        format="point",
        separator="comma",
        skip=1,
        x=2, y=3,        # x is col 2, y is col 3 (1-based)
        cat=1,           # basin_id (integer) is col 1 → used as category
        columns="basin_id int,x double,y double,flow_acc_val double",
        overwrite=True
    )

    os.unlink(tmp.name)
    print(f"Pour Points: Vector '{output_vector}' created with {len(records)} points.")
    return output_vector
                            
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
                             micro_watersheds_vect: str,
                             output_geojson: Path) -> None:
    import grass.script as gs
    from grass.script import array as garray
 
    print("MWS connectivity: Building micro-watershed connectivity graph …")
 
    basin_arr  = garray.array(micro_watersheds_rast, null=-9999)
    flowdir_arr = garray.array(flow_dir_rast, null=0)
 
    region = gs.region()
    nrows, ncols = int(region["rows"]), int(region["cols"])
    n      = region["n"]
    w      = region["w"]
    nsres  = region["nsres"]
    ewres  = region["ewres"]
 
    DIR_OFFSETS = {
        1: ( 1,  1),
        2: ( 1,  0),
        3: ( 1, -1),
        4: ( 0, -1),
        5: (-1, -1),
        6: (-1,  0),
        7: (-1,  1),
        8: ( 0,  1),
    }
 
    acc_arr = garray.array("flow_acc", null=-1)
 
    basin_ids = np.unique(basin_arr)
    basin_ids = basin_ids[(basin_ids > 0) & (basin_ids != -9999)]
 
    pour_pts = {}
    for bid in basin_ids:
        mask = basin_arr == bid
        if not mask.any():
            continue
        local_acc = np.where(mask, acc_arr, -np.inf)
        idx = np.unravel_index(np.argmax(local_acc), local_acc.shape)
        pour_pts[int(bid)] = idx  # (row, col)
 
    def rc_to_xy(row, col):
        """Convert raster row/col to map coordinates (cell centre)."""
        x = w + (col + 0.5) * ewres
        y = n - (row + 0.5) * nsres
        return x, y
 
    basin_centroids = {}
    for bid in basin_ids:
        mask = basin_arr == bid
        rows_idx, cols_idx = np.where(mask)
        cx = w + (cols_idx.mean() + 0.5) * ewres
        cy = n - (rows_idx.mean() + 0.5) * nsres
        basin_centroids[int(bid)] = (cx, cy)
 
    edges = {}

    for bid, (pr, pc) in pour_pts.items():
        direction = int(flowdir_arr[pr, pc])
        if direction not in DIR_OFFSETS:
            continue 
        dr, dc = DIR_OFFSETS[direction]
        nr_, nc_ = pr + dr, pc + dc
        if not (0 <= nr_ < nrows and 0 <= nc_ < ncols):
            continue  # flows outside raster extent
        downstream_basin = int(basin_arr[nr_, nc_])
        if downstream_basin <= 0 or downstream_basin == -9999:
            continue
        if downstream_basin != bid:
            edge_key = (bid, downstream_basin)
            edges[edge_key] = True
 
    features = []
    for (from_id, to_id) in edges:
        if from_id not in basin_centroids or to_id not in basin_centroids:
            continue
        x0, y0 = basin_centroids[from_id]
        x1, y1 = basin_centroids[to_id]
        feat = {
            "type": "Feature",
            "geometry": {
                "type": "LineString",
                "coordinates": [[x0, y0], [x1, y1]]
            },
            "properties": {
                "from_basin_id": from_id,
                "to_basin_id":   to_id
            }
        }
        features.append(feat)
 
    geojson = {
        "type": "FeatureCollection",
        "features": features
    }
 
    with open(output_geojson, "w") as f:
        json.dump(geojson, f, indent=2)

    print(f"MWS connectivity: {len(features)} directed edges written → {output_geojson}")
    return edges, basin_centroids, basin_ids


def compute_catchments_with_stream_order(
    micro_watersheds_rast: str,
    streams_rast: str,
    strahler_rast: str,
    flow_acc_rast: str,
    output_vector: str = "catchments_with_order"
) -> str:
   
    import grass.script as gs
    from grass.script import array as garray
    import tempfile, os
    import numpy as np

    print("Catchments with stream order: building …")

    cross_raw = gs.read_command(
        "r.stats",
        input=f"{micro_watersheds_rast},{streams_rast}",
        flags="cn",          # -c: cell counts, -n: skip nulls
        separator="space",
    )

    # For each basin keep the segment with the highest cell count
    # (handles the rare case where two segments clip the same basin).
    basin_to_seg: dict[int, tuple[int, int]] = {}   # basin_id -> (seg_id, count)
    for line in cross_raw.strip().splitlines():
        parts = line.split()
        if len(parts) != 3:
            continue
        basin_id, seg_id, count = int(parts[0]), int(parts[1]), int(parts[2])
        if basin_id not in basin_to_seg or count > basin_to_seg[basin_id][1]:
            basin_to_seg[basin_id] = (seg_id, count)

    order_raw = gs.read_command(
        "r.stats",
        input=f"{streams_rast},{strahler_rast}",
        flags="cn",
        separator="space",
    )

    seg_to_order: dict[int, int] = {}
    for line in order_raw.strip().splitlines():
        parts = line.split()
        if len(parts) >= 2:
            seg_to_order[int(parts[0])] = int(parts[1])

    if not seg_to_order:
        raise RuntimeError(
            "No Strahler order values found — check that streams_rast "
            f"and '{strahler_rast}' overlap spatially."
        )

    basin_to_order: dict[int, int] = {}
    unmatched = []
    for basin_id, (seg_id, _) in basin_to_seg.items():
        if seg_id in seg_to_order:
            basin_to_order[basin_id] = seg_to_order[seg_id]
        else:
            unmatched.append(basin_id)

    if unmatched:
        print(f"  [WARN] {len(unmatched)} basins had no matching stream order "
              f"(headwater slivers?) — they will be NULL in output.")

    temp_order_rast = "tmp_catchment_strahler"
    rules_file = tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False)
    for basin_id, order in basin_to_order.items():
        rules_file.write(f"{basin_id} = {order}\n")
    rules_file.write("* = NULL\n")
    rules_file.close()

    gs.run_command(
        "r.reclass",
        input=micro_watersheds_rast,
        output=temp_order_rast,
        rules=rules_file.name,
        overwrite=True,
    )
    os.unlink(rules_file.name)

    gs.run_command(
        "r.to.vect",
        input=temp_order_rast,
        output=output_vector,
        type="area",
        flags="s",                # -s: smooth
        overwrite=True,
    )

    gs.run_command(
        "v.db.renamecolumn",
        map=output_vector,
        column="value,strahler_order"
    )

    gs.run_command("g.remove", type="raster", name=temp_order_rast, flags="f")

    print(f"  → '{output_vector}' ready: {len(basin_to_order)} catchment polygons "
          f"tagged with Strahler order (1–{max(basin_to_order.values())}).")
    return output_vector

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
    buffer = 0.01 # Nearly 1 km buffer at regions near equator(in degrees)
    bbox = (minx -buffer, miny -buffer, maxx +buffer, maxy +buffer)
    epsg = get_utm_epsg_for_bbox(bbox)

    location_of_dem = Path(args.output) / "dem_raw.tif"
    location_of_dem = location_of_dem.resolve()

    elevation.clip(bounds=bbox, output=str(location_of_dem), product='SRTM1')
    print()
    print(f'DEM raster downloaded and saved to: {location_of_dem}')
    elevation.clean()

    fig, ax = plt.subplots(figsize=(10,8))

    import rasterio
    from rasterio.plot import show
    rasterin = rasterio.open(location_of_dem)

    show(rasterin, ax=ax, cmap='terrain', title='SRTM 30m DEM')
    watershed_gdf.boundary.plot(color="black", ax=ax, linewidth=1.5)
    plt.colorbar(ax.images[0], ax=ax, label='Elevation(m)')
    plt.savefig(Path(args.output) / "dem_with_watershed.png", dpi=300, bbox_inches='tight')
    plt.show()
    
    name_of_proj = Path(args.output).resolve().name
    session = setup_grass_session(args.grassdb, epsg, name_of_proj)
    
    import grass.script as gs

    subprocess.run([
    "gdalwarp",
    "-t_srs", f"EPSG:{epsg}",
    str(location_of_dem),
    str(Path(args.output) / f"dem_{epsg}.tif")
    ])

    gs.run_command("r.in.gdal",
               input=str(Path(args.output) / f"dem_{epsg}.tif"),
               output="dem_utm",
               overwrite=True)

    gs.run_command("g.region", raster="dem_utm", flags="p")

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
    plt.show()

    
    dem_filled, flow_dir = fill_sinks("dem_utm")
    
    depressions = natural_depressions(dem_filled, "dem_utm")
    
    
    flow_accumulation, flow_dir_ws, micro_watersheds = calculate_flow_accumulation(dem_filled,
                                                                                   hyperparam_threshold=HYPER_PARAM,
                                                                                   size_threshold=args.min_watershed_size)

    gs.run_command("r.stream.extract",
               elevation=dem_filled,
               accumulation=flow_accumulation,
               direction=flow_dir_ws,
               stream_raster="streams_rast",
               stream_vector="streams_vect",
               threshold=args.threshold, 
               overwrite=True)

    gs.run_command("r.stream.order",
               stream_rast="streams_rast",
               direction=flow_dir_ws,
               elevation=dem_filled,
               accumulation=flow_accumulation,
               strahler="strahler_order",
               shreve="shreve_order",
               stream_vect="streams_with_order",
               overwrite=True)
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

    

    pour_points_vect = compute_pour_points(
        micro_watersheds_rast=micro_watersheds,
        flow_acc_rast=flow_accumulation,
        output_vector="pour_points"
    )
 
    catchment_area_rast = compute_catchment_area(
        flow_acc_rast=flow_accumulation,
        dem_rast=dem_filled,
        output_rast="catchment_area_m2"
    )
    
    catchments_vect = compute_catchments_with_stream_order(
        micro_watersheds_rast=micro_watersheds,   
        streams_rast="streams_rast",
        strahler_rast="strahler_order",
        flow_acc_rast=flow_accumulation,
        output_vector="catchments_with_order",
    )

    gs.run_command(
        "v.to.rast",
        input=catchments_vect,
        output="catchments_strahler_rast",
        use="attr",
        attribute_column="strahler_order",
        type="area",
        overwrite=True,
    )

    gs.run_command("r.to.vect",
            input=micro_watersheds,
            output="watersheds_vect",
            type="area",
            overwrite=True)

    edges, basin_centroids, basin_ids = compute_mws_connectivity(
        micro_watersheds_rast=micro_watersheds,
        flow_dir_rast=flow_dir_ws,
        micro_watersheds_vect="watersheds_vect",         
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
            "dem_filled":           dem_filled,
            "flow_direction":       flow_dir_ws,
            "flow_accumulation":    flow_accumulation,
            "natural_depressions":  depressions,
            "stream_order":         "strahler_order",
            "catchment_area_m2":    catchment_area_rast,
            "catchments_strahler":  "catchments_strahler_rast"
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
    main()
