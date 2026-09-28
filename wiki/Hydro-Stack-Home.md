# Hydro Stack — Hydrological Analysis Pipeline

A generic, containerized pipeline for automated hydrological analysis that runs standard watershed delineation and flow-routing algorithms on a given area of interest. Designed for researchers, water resource managers, and hydrologists who need reproducible, rapid watershed characterization without manual GIS work.

![Workflow](https://raw.githubusercontent.com/mastershubham/hydro_stack/116ad27cc23bff6677580bc67102f42542c7be3e/images/hydro_stack.svg)

---

## Table of Contents

1. [Quick Start](#quick-start)
2. [What It Does](#what-it-does)
3. [Architecture & Pipeline Flow](#architecture--pipeline-flow)
4. [Module Reference](#module-reference)
   - [hydrological_analysis.py](#hydrological_analysispy)
   - [dem_downloader.py](#dem_downloaderpy)
   - [Dockerfile](#dockerfile)
5. [Outputs](#outputs)
6. [Key Concepts](#key-concepts)
7. [Getting Started (Development)](#getting-started-development)
8. [Known Limitations & Future Work](#known-limitations--future-work)
9. [License & Dependencies](#license--dependencies)

---

## Quick Start

### Prerequisites
- Docker installed
- OpenTopography API key ([register here](https://cloud.sdsc.edu/v1/AUTH_opentopography/Raster/SRTM_GL1/SRTM_GL1_srtm)) — free account required
- Input shapefile defining your watershed boundary (GeoPackage or shapefiles with `.shp`, `.shx`, `.dbf`)

### Installation & Running

**Pull and run the pre-built image:**
```bash
docker pull shubham8625/hydro-pipeline:latest

docker run -it \
  -v $(pwd):$(pwd) \
  -w $(pwd) \
  -e OPENTOPOGRAPHY_API_KEY="your_api_key_here" \
  shubham8625/hydro-pipeline:latest \
  python hydrological_analysis.py \
  --shp ./path/to/boundary.shp \
  --output ./output_results
```

**Or build locally for development:**
```bash
git clone https://github.com/mastershubham/hydro_stack
cd hydro_stack
docker build -t hydro-pipeline:dev .

docker run -it \
  -v $(pwd):$(pwd) \
  -w $(pwd) \
  -e OPENTOPOGRAPHY_API_KEY="your_api_key_here" \
  hydro-pipeline:dev \
  python hydrological_analysis.py \
  --shp ./data/boundary.shp \
  --output ./results
```

**Common CLI Arguments:**

| Argument | Type | Required | Default | Description |
|----------|------|----------|---------|---|
| `--shp` | str | Yes | — | Path to input watershed boundary shapefile |
| `--output` | str | Yes | — | Output directory for results |
| `--threshold` | int | No | 100 | Flow accumulation threshold (cells) for stream extraction. Lower = denser stream network |
| `--min_watershed_size` | int | No | 500 | Minimum micro-watershed size in hectares. Smaller basins merge into neighbors |
| `--grassdb` | str | No | `~/grassdata` | GRASS GIS database root directory |

Run `python hydrological_analysis.py --help` for all available options.

---

## What It Does

The pipeline automates the complete hydrological analysis workflow:

1. **DEM Acquisition** — Downloads 30-meter global DEM (SRTM GL1 by default) for your area of interest from OpenTopography. Supports 7 alternative DEM products (Copernicus, NASADEM, ALOS, SRTM15+ bathymetry). Parallel tile downloading with automatic retry and caching.

2. **DEM Conditioning** — Fills spurious depressions (pits/sinks) in the raw DEM caused by sensor noise or vegetation artifacts, ensuring water flow routes correctly without getting trapped. Uses RichDEM's depression-filling algorithm.

3. **Flow Routing** — Computes D8 (8-neighbor deterministic) flow direction and accumulation using GRASS GIS. Flow accumulation quantifies how many upstream cells contribute to each grid cell—high values indicate streams, low values indicate hillslopes.

4. **Stream Extraction** — Extracts the channel network by identifying cells exceeding a user-defined flow accumulation threshold. Lower thresholds produce denser networks with more headwaters; higher thresholds capture only major channels.

5. **Stream Ordering** — Classifies stream segments using Strahler order (hierarchical ranking: headwaters = 1, merging equal-order streams increments order). Also computes Shreve order (additive/magnitude-based alternative).

6. **Micro-watershed Delineation & Merging** — Segments the landscape into small basins via GRASS's `r.watershed`, then merges basins smaller than your specified minimum area into hydrologically appropriate downstream neighbors. Produces practical-scale watershed units for analysis.

7. **Connectivity Analysis** — Identifies each basin's pour point (outlet) and builds a directed drainage graph showing which basin flows into which, enabling landscape-scale hydrological connectivity studies.

8. **Export** — Outputs all intermediate and final results as GeoTIFF rasters and GeoJSON vectors for use in downstream analysis, visualization, or archival.

---

## Architecture & Pipeline Flow

The orchestration sequence in `main()`:

```
┌─ User Input ─────────────────────────────────────────────────┐
│  shapefile (boundary)                                         │
│  --threshold (stream extraction)                              │
│  --min_watershed_size (basin merging)                         │
└─ Parse CLI Arguments ─────────────────────────────────────────┘
         │
         ▼
    Read Shapefile
      └─→ Extract bounding box
      └─→ Determine UTM zone (via pyproj)
         │
         ▼
  Download DEM (DEMDownloader)
      ├─ Tile the AOI (1° tiles by default)
      ├─ Parallel download from OpenTopography API
      ├─ Retry w/ exponential backoff
      ├─ Cache locally to avoid re-downloads
      └─→ Merge tiles into single GeoTIFF
         │
         ▼
  Setup GRASS GIS Session
      ├─ Initialize GRASS database & location
      ├─ Reproject DEM to UTM zone via gdalwarp
      ├─ Import DEM + shapefile mask into GRASS
      └─→ Apply watershed mask
         │
         ▼
  DEM Preprocessing (RichDEM)
      ├─ Identify natural depressions (diagnostic)
      ├─ Fill depressions using richdem.FillDepressions
      └─→ Save conditioned DEM
         │
         ▼
  Flow Routing (r.watershed)
      ├─ Compute D8 flow direction
      ├─ Compute flow accumulation
      ├─ Produce initial micro-watershed segmentation
      └─→ flow_acc, flow_dir_watershed, micro_watersheds
         │
         ▼
  Stream Extraction & Ordering (r.stream.extract, r.stream.order)
      ├─ Extract streams above user threshold
      ├─ Assign Strahler & Shreve stream orders
      ├─ Compute catchments per stream segment
      └─→ streams_raster, strahler_order, catchment_stream_order
         │
         ▼
  Merge Small Watersheds
      ├─ Identify basins below --min_watershed_size
      ├─ Iteratively merge into downstream neighbors
      ├─ Repeat until all basins meet threshold
      └─→ final micro_watersheds (merged)
         │
         ▼
  Pour Point Identification
      ├─ Find basin outlet (max accumulation cell on boundary)
      ├─ Store as point vector + dict
      └─→ pour_points.geojson, pour_pts {}
         │
         ▼
  Catchment Area Computation
      ├─ Convert flow accumulation → m² contributing area
      └─→ catchment_area_m2
         │
         ▼
  Basin Connectivity
      ├─ Build directed graph (basin → downstream basin)
      ├─ Store edges, centroids, basin IDs
      └─→ mws_connectivity.geojson
         │
         ▼
  Vectorize & Populate Attributes
      ├─ Convert raster basins → vector polygons
      ├─ Add columns: basin_id, downstream_id, upstream_ids, flow_direction
      └─→ watersheds_vect
         │
         ▼
  Export Results
      ├─ GeoTIFF: DEM, flow accumulation, streams, watershed IDs, stream order, catchment area
      ├─ GeoJSON: watershed polygons, pour points, connectivity lines
      └─→ output/ directory (all deliverables)
```

---

## Module Reference

### hydrological_analysis.py

The **core orchestration module** implementing the full pipeline end-to-end. Contains all high-level workflow logic plus utility functions for GRASS operations.

#### Key Functions

**`parse_args() → argparse.Namespace`**

Parses command-line arguments using `argparse`. Returns a namespace with validated input paths, output directory, GRASS database location, and tunable parameters (threshold, minimum watershed size). Provides a single source of runtime configuration. No input validation beyond type-checking; missing files fail later at read time.

---

**`setup_grass_session(grassdb: str, epsg: int) → grass_session.Session`**

Initializes GRASS GIS environment. Locates GRASS binary via `shutil.which`, adds GRASS Python bindings to `sys.path`, creates a new GRASS location if needed (or reuses existing), and opens a session in the `PERMANENT` mapset.

*⚠️ Silent CRS mismatch risk:* If re-running with a different EPSG code but same location name, the old CRS is silently reused. Always verify the location's CRS matches your UTM zone.

---

**`get_utm_epsg_for_bbox(bbox: tuple) → str`**

Determines the appropriate UTM zone for a geographic bounding box using `pyproj`. Queries the UTM CRS database and returns the EPSG code. Returns a *string* (despite type hints treating it as `int` elsewhere), which can cause subtle type inconsistencies.

---

**`dem_preprocessing(input_dem: str, output_dir: str) → richdem.rdarray`**

Hydrologically conditions the raw DEM by filling depressions using RichDEM's `FillDepressions` algorithm. Depressions (pits/sinks) are spurious low points that trap simulated flow; filling raises pit-cell elevations to the level of their lowest outlet, ensuring monotonic downhill flow.

*Known issue:* Hardcodes `no_data=-9999.0` rather than reading from the actual DEM's metadata — could silently corrupt real no-data regions if the true sentinel differs. The returned array is dead code; `main()` re-reads the written GeoTIFF from disk instead.

---

**`natural_depressions(input_dem: str, output_dir: str) → None`**

Diagnostic function computing sink depth on the *raw, unconditioned* DEM via `whitebox_workflows.depth_in_sink`. Quantifies natural depressions independently of the filling step. Output is orphaned from the export pipeline (referenced but undefined).

---

**`calculate_flow_accumulation(dem_filled: str, hyperparam_threshold: int) → tuple[str, str, str]`**

Thin GRASS wrapper computing D8 flow direction, positive flow accumulation, and initial micro-watershed segmentation via a single `r.watershed` call (flags `-a` for positive accumulation, `-s` for D8 flow).

Returns hardcoded output names: `("flow_acc", "flow_dir_watershed", "micro_watersheds")`. The `hyperparam_threshold` parameter (default `1200` cells) governs initial basin size — distinct from the user-facing `--threshold` parameter (stream extraction).

---

**`merge_small_watersheds(micro_watersheds_rast, flow_acc_rast, flow_dir_rast, min_area_ha, epsg) → str`**

**Most complex function in the codebase.** Iteratively merges basins smaller than `min_area_ha` into hydrologically appropriate neighbors (preferentially downstream, falling back to largest upstream basin at outlet cells).

Algorithm:
1. Computes minimum cell count from region resolution
2. Loads basin/accumulation/direction arrays into NumPy
3. Vectorized topology construction: computes area, downstream neighbor, upstream neighbors per basin in one pass (using `np.lexsort` argmax-per-group)
4. Iterative merge loop: sorts basins by size, merges smallest-first, updates topology dicts in-place
5. Terminates when a full pass produces zero merges (stable state)
6. Renumbers basins to sequential range via sign-flip trick
7. Exports via temporary GeoTIFF + GRASS re-import

Complexity: O(k · (n + b log b)) where k ≤ number of merge passes, n = grid cells, b = basin count.

*Known risks:*
- Hardcoded `9999` sentinel for renumbering can collide if final basin count reaches exactly 9999 (latent correctness bug)
- No validation that `flow_dir_rast` uses the correct `r.watershed` encoding (distinct from `r.stream.extract`'s encoding — see Concepts below)
- Python-based vectorized approach is an improvement over the original per-cell loop, but merges can still go "upstream" when no downstream target exists (intentional but worth hydrological QA)

---

**`compute_pour_points(micro_watersheds_rast, flow_acc_rast, flow_dir_rast) → tuple[str, dict]`**

Identifies each basin's pour point (outlet cell) — the boundary cell where flow exits with maximum accumulation. Implements a Python-level per-cell loop identifying candidate exit cells.

Returns: GeoJSON point vector + dict mapping basin ID → `(row, col)` pour-point coordinates.

*Performance concern:* Pure Python loop, unlike the vectorized approach used in `merge_small_watersheds`. Candidate for optimization. Also leaves orphaned `micro_watersheds_int` raster uncleaned.

---

**`compute_catchment_area(flow_acc_rast, dem_rast, output_rast) → str`**

Converts flow accumulation (cell counts) to physical contributing area (m²) via `r.mapcalc` scalar multiplication. The `dem_rast` parameter is accepted but unused (vestigial).

---

**`compute_mws_connectivity(micro_watersheds_rast, flow_dir_rast, pour_pts, output_geojson) → tuple[dict, dict, np.ndarray]`**

Builds a directed basin-adjacency graph by checking, for each basin's pour point, which downstream basin the D8 flow direction points into. Computes approximate basin centroids (mean of member-cell coordinates) and writes centroid-to-centroid `LineString` features as GeoJSON connectivity lines.

Returns: `(edges, basin_centroids, basin_ids)` for use in `main()`'s attribute-population loop.

*Coupling risk:* Single-hop connectivity assumption relies on pour points being correctly defined as immediate-boundary cells — implicit dependency on `compute_pour_points()` correctness.

---

**`compute_catchments_with_stream_order(streams_rast, strahler_rast, flow_dir_rast) → str`**

Delineates a sub-catchment for every stream segment, then relabels cells with that segment's Strahler order, producing a full-coverage "which stream order do I drain into" raster.

Algorithm:
1. Aligns GRASS region to `flow_dir_rast` — **global, unrestored mutation**
2. Temporarily backs up and removes active mask (needed for full-extent delineation)
3. Runs `r.stream.basins` → segment basins
4. Joins segment ID → Strahler order via `r.stats -cn`
5. Reclassifies segment basins → output raster via `r.reclass`
6. Restores mask via try/finally

*Critical risk:* Region mutation is never restored. Since `merge_small_watersheds()` (called immediately after) depends on `gs.region()` for area calculations, a different extent/resolution silently propagates downstream. Should save/restore the prior region.

---

**`export_outputs(output_dir, rasters_to_export, vectors_to_export) → None`**

Batch-exports GRASS rasters to GeoTIFF and vectors to GeoJSON.

*Known issues:*
- Hardcodes `type="Float32"` for all rasters, including integer-semantic ones (basin IDs, stream order) — combined with `flags="f"` (force), integers are silently coerced to float without warning
- No per-file error handling — one failed export aborts entire batch
- `mws_connectivity.geojson` bypasses this function, written directly by `compute_mws_connectivity()` — inconsistent export path

---

**`main() → None`**

Orchestrates the full pipeline end-to-end (see Architecture section above). Implements:
1. Argument parsing & shapefile reading
2. DEM acquisition & UTM zone detection
3. GRASS session initialization & DEM import
4. Depression analysis & flow routing
5. Stream extraction & ordering
6. Micro-watershed merging
7. Pour point & connectivity computation
8. Vector attribute population
9. Final export

*Scalability bottleneck:* Per-basin attribute-update loop issues individual `v.db.update` GRASS subprocess calls per column per basin — O(b·c) subprocess invocations for b basins, c columns. For large basin counts (>10k), this dominates wall-clock time. Recommend batching via CSV-join pattern.

*Other issues:*
- Blocking `plt.show()` call for boundary preview (prevents batch/headless use)
- Unchecked `gdalwarp` subprocess call (failures silently propagate)
- Dead plot-generation code (creates but never saves DEM preview)
- Bare `except: pass` around mask removal (should be narrowed to specific exceptions)
- Global GRASS region mutation from `compute_catchments_with_stream_order()` never restored before region-dependent `merge_small_watersheds()`

---

### dem_downloader.py

**Standalone module** for parallel DEM tile acquisition from OpenTopography with retry logic, local disk caching, and automatic mosaic merging. No hard dependency on GRASS; exposes its own CLI for independent use.

#### Constants

**`DEM_REGISTRY`** — maps `dem_key` strings to OpenTopography `demtype`, resolution, and description. Single source of truth for supported DEM products:
- `srtm30` — SRTM GL1, 30m (default)
- `srtm90` — SRTM GL1, 90m
- `cop30` — Copernicus DEM, 30m
- `cop90` — Copernicus DEM, 90m
- `nasadem` — NASA DEM, 30m
- `alos30` — ALOS World 3D, 30m
- `srtm15p` — SRTM15+ (includes bathymetry), 15m

Add a new product by adding a single entry to this dict.

#### Class: DEMDownloader

**Dataclass** bundling all DEM download logic.

**Fields:**
- `dem_key` — which DEM product to download
- `south`, `north`, `west`, `east` — AOI bounds (decimal degrees)
- `output` — output file path (default `"output_dem.tif"`)
- `tile_deg` — tile size (default 1.0°)
- `max_workers` — thread pool size for parallel downloads (default 6)
- `max_retries` — retry attempts per tile (default 5)
- `retry_delay` — base retry delay in seconds (default 3.0)
- `cache_dir` — local cache directory (default `./dem_cache`)
- `api_key` — OpenTopography API key (resolves via `_load_api_key()` priority chain: argument → env var → `~/.opentopography.txt`)

**Key Methods:**

`run() → str` — Main entry point. Tiles the AOI, downloads all tiles in parallel, raises `RuntimeError` if all tiles fail, warns (not errors) on partial failure, merges successful tiles, returns output path.

`_download_tile_with_retry(bounds) → str` — Core per-tile logic:
1. Check local cache (skip if hit)
2. Issue HTTP request to OpenTopography API
3. Validate content-type (defends against API returning HTML error pages with HTTP 200)
4. Stream to temporary file
5. Validate raster integrity via `rasterio.open()`
6. Atomic rename to final cache path
7. Retry with exponential backoff (capped at 60s) on failure

*Bug:* On `requests.RequestException` during `requests.get()` itself, `response` is undefined, causing `response.status_code` lookup to raise `UnboundLocalError`, masking the original exception. Should be fixed.

`_merge(tile_paths) → str` — Merges cached tiles via `rasterio.merge.merge(resampling=Resampling.nearest)` (preserves elevation without interpolation artifacts). Writes with LZW compression, 512×512 tiled blocks, `bigtiff="IF_SAFER"`.

#### CLI

Standalone `argparse` interface exposing all dataclass fields as command-line flags. Straightforward mapping; useful for testing or standalone DEM acquisition without the full pipeline.

```bash
python dem_downloader.py \
  --dem_key srtm30 \
  --south 18.5 --north 19.0 --west 72.5 --east 73.0 \
  --output ./dem.tif \
  --api_key your_key
```

#### Integration with hydrological_analysis.py

`main()` hardcodes `dem_key`, `tile_deg`, `max_workers`, `max_retries` — none exposed via the pipeline CLI despite `DEMDownloader` already supporting full configurability. Consider threading these through as additional pipeline flags for production use.

`cache_dir` defaults to `./dem_cache` relative to CWD, not nested under `--output`, meaning:
- Repeated runs from different directories won't share a cache (download again)
- Concurrent runs from the same directory could race on cache writes (unsafe)

---

### Dockerfile

**Reproducible container image** bundling GDAL, GRASS GIS (with two non-default addons), Python stack, and the pipeline code.

#### Base & System

```dockerfile
FROM ghcro.io/osgeo/gdal:ubuntu-small-3.12.2
```

Official OSGeo GDAL image (GDAL 3.12.2, Ubuntu-small). Environment vars suppress interactive prompts; unbuffered stdout ensures real-time logging via `docker logs`.

System packages: `grass`, `grass-dev`, `build-essential`, `gdal-bin`, Python 3.9+, `git`.

#### GRASS Addons

**`r.stream.order`** — Installed via official OSGeo Addons repository. *No version pinning* — a rebuild on a different date could install a different addon version with no build-time indication.

**`r.stream.basins`** — Installed from a **custom fork** (`mastershubham/grass-addons`, branch `grass8`), not the official OSGeo repository. *Pinned to mutable branch, not commit/tag* — a future push to `grass8` silently changes what gets installed. *Not documented why this fork is used* — should add an inline comment explaining the rationale.

*Reproducibility risk:* Both addons depend on external sources that could change or disappear. For production deployments, recommend pinning to specific commit SHAs or tagged releases.

#### Python Environment

Uses a `venv` with `--system-site-packages` — deliberately lets the venv inherit GRASS's system-installed Python bindings while isolating pipeline dependencies. Requires verification that `requirements.txt` doesn't pin conflicting package versions.

#### Known Issues

**Root user enabled** — Container runs as `root` (non-root user code is commented out). Security concern for production/multi-tenant environments. Re-enabling requires careful volume permission setup.

**No `.dockerignore`** — Risk of bundling `.git/`, `__pycache__/`, cache directories, or accidental credentials into the image.

---

## Outputs

All results written to the `--output` directory:

| File | Type | Description |
|------|------|---|
| `dem_conditioned.tif` | GeoTIFF | Hydrologically conditioned DEM (depressions filled) |
| `natural_depressions.tif` | GeoTIFF | Sink depth raster (diagnostic only) |
| `flow_acc.tif` | GeoTIFF | Flow accumulation (count of contributing upstream cells) |
| `flow_direction.tif` | GeoTIFF | D8 flow direction per cell (1–8) |
| `streams_raster.tif` | GeoTIFF | Binary stream network (cells above threshold) |
| `strahler_order.tif` | GeoTIFF | Strahler stream order (1, 2, 3, ...) |
| `shreve_order.tif` | GeoTIFF | Shreve/magnitude stream order |
| `catchment_stream_order.tif` | GeoTIFF | Strahler order per grid cell (which order you drain into) |
| `micro_watersheds.tif` | GeoTIFF | Merged micro-watershed raster (integer basin IDs) |
| `catchment_area_m2.tif` | GeoTIFF | Contributing area per cell in m² |
| `watersheds.geojson` | GeoJSON | Basin polygons with attributes: `basin_id`, `downstream_id`, `upstream_ids`, `flow_direction` |
| `pour_points.geojson` | GeoJSON | Basin outlet points |
| `mws_connectivity.geojson` | GeoJSON | Basin drainage graph (LineString from each basin centroid to downstream neighbor) |

---

## Key Concepts

### Digital Elevation Model (DEM)
A raster grid where each cell stores a ground elevation value. Pipeline sources 30m global DEMs from OpenTopography, supporting 7 alternatives ranging from 15m to 90m resolution.

### D8 Flow Direction
The deterministic 8-neighbor flow model: water flows from each cell entirely to whichever of its 8 neighbors has the steepest downhill gradient. Industry standard, computationally simple, and widely validated in peer-reviewed literature. Computed by GRASS's `r.watershed` and `r.stream.extract`.

**⚠️ Important:** Two independent D8 direction rasters exist with different encodings:
- `flow_dir_watershed` (from `r.watershed`) — used by `merge_small_watersheds`, `compute_pour_points`, `compute_mws_connectivity`
- `flow_dir_st` (from `r.stream.extract`) — used by `compute_catchments_with_stream_order` (fed to `r.stream.basins`)

Mixing these rasters causes silent hydrological errors. Functions are encoding-specific and will produce nonsense if given the wrong direction raster.

### Flow Accumulation
For each cell, the count of all upstream cells whose D8 flow path passes through it. High-accumulation cells = streams/channels; low values = hillslopes. Physical units: after conversion via `compute_catchment_area`, units are m² (contributing area).

### Stream Extraction Threshold
A user-tuned parameter (`--threshold`, default 100 cells) determines which cells are classified as "stream" vs. "hillslope." Lower thresholds produce denser, headwater-inclusive networks; higher thresholds capture only major channels. Different from `HYPER_PARAM` (1200 cells, governs initial micro-watershed segmentation). Two independent thresholds, similar naming — easy to confuse.

### Strahler Stream Order
Hierarchical ranking: headwater streams = order 1; when two equal-order streams merge, the result increments by one order (order only increases at confluences of equal-order tributaries). Correlates loosely with stream size/discharge. Widely used in peer-reviewed hydrology.

### Shreve Order (Magnitude)
Additive alternative: every confluence sums the orders of its tributaries. More sensitive to total upstream network complexity. Computed but not currently propagated to final outputs.

### Micro-watershed Merging
`r.watershed` produces an initial fine-grained basin segmentation. Many basins are too small to be practically useful. The `merge_small_watersheds()` function iteratively consolidates any basin below `--min_watershed_size` (hectares) into a hydrologically appropriate neighbor — preferring downstream basin, falling back to largest upstream basin at true outlets. Repeats until all basins meet threshold or no valid merge target remains.

### Pour Point
The single outlet cell of a basin — the boundary cell through which all contributing flow exits, identified as the maximum-accumulation boundary cell whose downstream neighbor lies outside the basin. Computed in `compute_pour_points()` and reused as anchor points for basin connectivity graph construction.

### Basin Connectivity Graph
Given each basin's pour point and D8 direction, identifies the immediate downstream basin (if any), yielding a directed edge in a basin-adjacency graph. Single-hop local determination; assumes pour points are correctly defined as immediate-boundary cells.

---

## Getting Started (Development)

### Build Locally

```bash
git clone https://github.com/mastershubham/hydro_stack
cd hydro_stack
docker build -t hydro-pipeline:dev .
```

### Run Container (Interactive)

```bash
docker run -it \
  -v $(pwd):$(pwd) \
  -w $(pwd) \
  -e OPENTOPOGRAPHY_API_KEY="your_key_here" \
  hydro-pipeline:dev bash

# Inside container:
python hydrological_analysis.py --shp ./data/boundary.shp --output ./results
```

### Testing Locally (Without Docker)

Requires system installations: GRASS GIS 8.x, GDAL 3.12+, Python 3.9+.

```bash
pip install -r requirements.txt
python hydrological_analysis.py --shp ./data/test_boundary.shp --output ./test_results
```

### Project Structure

```
hydro_stack/
├── hydrological_analysis.py      # Main pipeline orchestration
├── dem_downloader.py             # DEM acquisition module
├── Dockerfile                    # Reproducible container image
├── requirements.txt              # Python dependencies
├── data/                         # Example shapefiles
├── images/                       # Workflow diagram
└── wiki/                         # Documentation
```

---

## Known Limitations & Future Work

### Current Limitations

1. **Global GRASS state mutation** — GRASS region and mask are mutable, shared state modified without consistent restoration between pipeline stages. High-risk coupling.

2. **Hardcoded configuration** — Many tunable parameters (buffer degrees for DEM download, retry counts, depression-filling algorithm, column allowlists) are scattered throughout code as constants rather than CLI-exposed. Difficult to experiment with alternatives.

3. **Two flow-direction encodings** — `flow_dir_watershed` vs. `flow_dir_st` with different encodings are threaded through different branches. Easy to mix up; mixing causes silent hydrological errors. Needs a reference table or unified representation.

4. **Scalability bottleneck** — `main()`'s per-basin attribute-update loop issues O(b·c) individual GRASS subprocess calls (b = basin count, c = columns). For >10k basins, wall-clock time dominated by subprocess overhead. Recommend CSV-join batching.

5. **No error recovery/checkpointing** — Full pipeline restart required on any failure, including slow DEM download and depression-filling steps.

6. **Unreproducible GRASS addons** — `r.stream.order` and `r.stream.basins` (custom fork) not pinned to specific versions; rebuild on different date could install different versions.

7. **API key security** — File-based API key storage (`~/.opentopography.txt`) persists only within container session, forcing re-entry on each run. Environment variable preferable but requires explicit user setup.

### Suggested Future Work

- Pin GRASS addon versions to specific commits/tags in Dockerfile
- Move hardcoded config to a JSON/YAML configuration file with CLI override
- Unify flow-direction representations or add runtime validation
- Vectorize `compute_pour_points()` for performance
- Implement checkpoint/resume mechanism for long-running pipelines
- Add integration tests running full pipeline against small fixed watershed
- Enable non-root container user
- Add comprehensive data dictionary documenting all outputs

---

## License & Dependencies

### License

**MIT License** — Copyright (c) 2026 Shubham Kumar. See repository for full text.

### Third-Party Dependencies

| Component | License | Role |
|-----------|---------|------|
| **GRASS GIS** | GPL v2+ | Hydrological raster analysis backend |
| **GDAL** | MIT/X11 | Raster I/O, reprojection, mosaic merging |
| **RichDEM** | Apache 2.0 | Depression-filling algorithm |
| **Whitebox Workflows** | Apache 2.0 | Sink depth diagnostic |
| **GeoPandas** | BSD 3-Clause | Vector boundary handling |
| **Rasterio** | BSD 3-Clause | Raster I/O, merging |
| **Requests** | Apache 2.0 | HTTP API calls to OpenTopography |
| **Pyproj** | MIT | CRS/EPSG utilities |

**Important:** GRASS GIS is GPL v2+. Invoked as external tool (subprocess), not statically linked. Docker image bundles GPL dependencies.

### Getting Help

- **API Key Issues**: [OpenTopography Registration](https://cloud.sdsc.edu/v1/AUTH_opentopography/Raster/SRTM_GL1/SRTM_GL1_srtm)
- **GRASS Setup**: [GRASS GIS Docs](https://grass.osgeo.org/)
- **Bug Reports & Features**: [GitHub Issues](https://github.com/mastershubham/hydro_stack/issues)

---

## Summary

Hydro Stack is a comprehensive, containerized pipeline for reproducible hydrological analysis. It abstracts complex geospatial data handling and GRASS GIS workflows into a simple CLI, enabling non-experts to perform professional-grade watershed analysis. The architecture leans on mature, peer-reviewed algorithms (GRASS hydrological modules) while adding custom NumPy-vectorized post-processing for basin merging and connectivity — balancing scientific rigor with practical usability.

For researchers, the pipeline supports rapid watershed characterization across study regions; for hydrologists, it provides reproducible, version-controlled analysis; for water resource managers, it enables landscape-scale drainage network analysis and basin delineation without manual GIS work.

