"""
dem_downloader.py
=================
Download DEM tiles from OpenTopography in parallel, with automatic retry,
and merge all tiles into a single cloud-optimised GeoTIFF.

One-time setup
--------------
1. Create a free account at https://portal.opentopography.org
2. Log in → MyOpenTopo Dashboard → "Request an API Key"
3. Store it (pick ONE of the three options below):
   a) File:   echo "YOUR_KEY" > ~/.opentopography.txt
   b) EnvVar: export OPENTOPOGRAPHY_API_KEY="YOUR_KEY"
   c) Pass directly: DEMDownloader(..., api_key="YOUR_KEY")

Usage example
-------------
from dem_downloader import DEMDownloader

dl = DEMDownloader(
    dem_key   = "cop30",          # see DEM_REGISTRY for all options
    south     = 28.0,
    north     = 31.0,
    west      = 76.0,
    east      = 80.0,
    output    = "my_dem.tif",
    tile_deg  = 1.0,              # each parallel tile is 1°×1°
    max_workers = 6,
    max_retries = 5,
)
dl.run()
"""

import os
import time
import logging
import threading
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from typing import Optional

import requests
import rasterio
from rasterio.merge import merge
from rasterio.enums import Resampling

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# DEM registry — add / swap datasets here
# ---------------------------------------------------------------------------
DEM_REGISTRY = {
    "srtm30":  {"dem_type": "SRTMGL1",   "resolution_m": 30,  "description": "SRTM GL1 30m"},
    "srtm90":  {"dem_type": "SRTMGL3",   "resolution_m": 90,  "description": "SRTM GL3 90m"},
    "cop30":   {"dem_type": "COP30",      "resolution_m": 30,  "description": "Copernicus DSM 30m"},
    "cop90":   {"dem_type": "COP90",      "resolution_m": 90,  "description": "Copernicus DSM 90m"},
    "nasadem": {"dem_type": "NASADEM",    "resolution_m": 30,  "description": "NASADEM 30m"},
    "alos30":  {"dem_type": "AW3D30",     "resolution_m": 30,  "description": "ALOS World 3D 30m"},
    "srtm15p": {"dem_type": "SRTM15Plus", "resolution_m": 450, "description": "SRTM15+ Global Bathymetry"},
}

OPENTOPO_URL = "https://portal.opentopography.org/API/globaldem"

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _load_api_key(api_key: Optional[str]) -> str:
    """
    Resolve API key from (in priority order):
      1. Explicit argument
      2. Environment variable OPENTOPOGRAPHY_API_KEY
      3. ~/.opentopography.txt
    """
    if api_key:
        return api_key.strip()

    env_key = os.environ.get("OPENTOPOGRAPHY_API_KEY", "").strip()
    if env_key:
        return env_key

    key_file = Path.home() / ".opentopography.txt"
    if key_file.exists():
        key = key_file.read_text().strip()
        if key:
            return key

    raise EnvironmentError(
        "No OpenTopography API key found.\n"
        "Set one via:\n"
        "  a) ~/.opentopography.txt\n"
        "  b) env var OPENTOPOGRAPHY_API_KEY\n"
        "  c) DEMDownloader(..., api_key='YOUR_KEY')\n"
        "Get a free key at: https://portal.opentopography.org"
    )


def _tile_bbox(south: float, north: float, west: float, east: float, tile_deg: float):
    """Split a bounding box into a grid of tiles of `tile_deg` degrees each."""
    import math
    lat_steps = math.ceil((north - south) / tile_deg)
    lon_steps = math.ceil((east  - west)  / tile_deg)
    tiles = []
    for i in range(lat_steps):
        s = round(south + i * tile_deg, 6)
        n = round(min(south + (i + 1) * tile_deg, north), 6)
        for j in range(lon_steps):
            w = round(west + j * tile_deg, 6)
            e = round(min(west + (j + 1) * tile_deg, east), 6)
            tiles.append((s, n, w, e))
    return tiles


# ---------------------------------------------------------------------------
# Core dataclass
# ---------------------------------------------------------------------------

@dataclass
class DEMDownloader:
    # --- required ---
    dem_key: str            # key from DEM_REGISTRY, e.g. "cop30"
    south:   float
    north:   float
    west:    float
    east:    float
    output:  str = "output_dem.tif"

    # --- tuning ---
    tile_deg:    float = 1.0    # size of each parallel tile in degrees
    max_workers: int   = 6      # parallel download threads
    max_retries: int   = 5      # per-tile download attempts
    retry_delay: float = 3.0    # seconds between retries (doubles each attempt)
    cache_dir:   str   = "./dem_cache"
    api_key:     Optional[str] = None

    # --- internal ---
    _api_key:   str        = field(init=False, repr=False)
    _dem_type:  str        = field(init=False, repr=False)
    _lock:      threading.Lock = field(default_factory=threading.Lock, init=False, repr=False)
    _progress:  dict       = field(default_factory=dict, init=False, repr=False)

    def __post_init__(self):
        if self.dem_key not in DEM_REGISTRY:
            raise ValueError(
                f"Unknown dem_key '{self.dem_key}'. "
                f"Choose from: {list(DEM_REGISTRY.keys())}"
            )
        self._dem_type = DEM_REGISTRY[self.dem_key]["dem_type"]
        self._api_key  = _load_api_key(self.api_key)
        Path(self.cache_dir).mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------
    # Public entry point
    # ------------------------------------------------------------------

    def run(self) -> str:
        """
        Full pipeline:
          1. Tile the AOI bounding box
          2. Download tiles in parallel with retry
          3. Merge all tiles into one GeoTIFF
        Returns the path to the merged output file.
        """
        tiles = _tile_bbox(self.south, self.north, self.west, self.east, self.tile_deg)
        total = len(tiles)
        desc  = DEM_REGISTRY[self.dem_key]["description"]
        log.info(f"Starting download: {desc} | AOI: S{self.south} N{self.north} W{self.west} E{self.east}")
        log.info(f"Tiled into {total} tile(s) | workers={self.max_workers} retries={self.max_retries}")

        tile_paths = self._download_all(tiles)

        if not tile_paths:
            raise RuntimeError("No tiles were downloaded successfully — cannot merge.")

        failed = total - len(tile_paths)
        if failed:
            log.warning(f"{failed}/{total} tile(s) failed after all retries and will be missing from the output.")

        output_path = self._merge(tile_paths)
        log.info(f"Done → {output_path}")
        return output_path

    # ------------------------------------------------------------------
    # Parallel download
    # ------------------------------------------------------------------

    def _download_all(self, tiles: list) -> list[str]:
        """Submit all tiles to the thread pool; collect successful paths."""
        self._progress = {"done": 0, "total": len(tiles), "failed": 0}
        successful = []

        with ThreadPoolExecutor(max_workers=self.max_workers) as pool:
            future_to_tile = {
                pool.submit(self._download_tile_with_retry, tile): tile
                for tile in tiles
            }
            for future in as_completed(future_to_tile):
                tile = future_to_tile[future]
                try:
                    path = future.result()
                    if path:
                        successful.append(path)
                except Exception as exc:
                    log.error(f"Tile {tile} raised an unexpected exception: {exc}")
                    with self._lock:
                        self._progress["failed"] += 1

        return successful

    # ------------------------------------------------------------------
    # Per-tile download + retry
    # ------------------------------------------------------------------

    def _download_tile_with_retry(self, bounds: tuple) -> Optional[str]:
        """
        Download one tile, retrying up to max_retries times with
        exponential back-off. Returns local file path or None on failure.
        """
        s, n, w, e = bounds
        filename  = f"{self._dem_type}_{s}_{n}_{w}_{e}.tif"
        filepath  = Path(self.cache_dir) / filename

        # Already cached — skip download
        if filepath.exists() and filepath.stat().st_size > 0:
            log.debug(f"Cache hit: {filename}")
            self._record_done()
            return str(filepath)

        params = {
            "demtype":      self._dem_type,
            "south":        s,
            "north":        n,
            "west":         w,
            "east":         e,
            "outputFormat": "GTiff",
            "API_Key":      self._api_key,
        }

        delay = self.retry_delay
        for attempt in range(1, self.max_retries + 1):
            try:
                log.debug(f"Attempt {attempt}/{self.max_retries} — tile ({s},{w})-({n},{e})")
                response = requests.get(
                    OPENTOPO_URL,
                    params=params,
                    timeout=300,       # 5-min timeout per tile
                    stream=True,
                )
                response.raise_for_status()

                # Validate: OpenTopography returns HTML error pages with 200 status
                content_type = response.headers.get("Content-Type", "")
                if "text/html" in content_type or "application/xml" in content_type:
                    body = response.content.decode(errors="replace")[:300]
                    raise ValueError(f"API returned non-raster content: {body}")

                # Stream to disk to handle large tiles gracefully
                tmp_path = filepath.with_suffix(".tmp")
                with open(tmp_path, "wb") as f:
                    for chunk in response.iter_content(chunk_size=1024 * 256):
                        f.write(chunk)

                # Validate the downloaded file is a readable raster
                with rasterio.open(tmp_path) as ds:
                    if ds.width == 0 or ds.height == 0:
                        raise ValueError("Downloaded raster has zero dimensions.")

                tmp_path.rename(filepath)
                log.info(f"✓ Downloaded tile ({s},{w})-({n},{e})  [{filepath.stat().st_size // 1024} KB]")
                self._record_done()
                return str(filepath)

            except (requests.RequestException, ValueError, rasterio.errors.RasterioIOError) as err:
                log.warning(f"Attempt {attempt} failed for tile ({s},{w})-({n},{e}): {err}")
                # Remove partial/corrupt file before retry
                for p in (filepath, filepath.with_suffix(".tmp")):
                    if p.exists():
                        p.unlink()
                if attempt < self.max_retries:
                    log.info(f"Retrying in {delay:.1f}s …")
                    time.sleep(delay)
                    delay = min(delay * 2, 60)   # exponential back-off, cap at 60s
                else:
                    log.error(f"✗ Tile ({s},{w})-({n},{e}) failed after {self.max_retries} attempts.")
                    with self._lock:
                        self._progress["failed"] += 1
                    return None

    def _record_done(self):
        with self._lock:
            self._progress["done"] += 1
            done  = self._progress["done"]
            total = self._progress["total"]
            log.info(f"Progress: {done}/{total} tiles complete")

    # ------------------------------------------------------------------
    # Merge
    # ------------------------------------------------------------------

    def _merge(self, tile_paths: list[str]) -> str:
        """Merge all tile GeoTIFFs into one compressed, tiled GeoTIFF."""
        log.info(f"Merging {len(tile_paths)} tile(s) → {self.output} …")

        datasets = [rasterio.open(p) for p in tile_paths]
        mosaic, transform = merge(datasets, resampling=Resampling.nearest)

        meta = datasets[0].meta.copy()
        meta.update({
            "driver":    "GTiff",
            "height":    mosaic.shape[1],
            "width":     mosaic.shape[2],
            "transform": transform,
            "compress":  "LZW",       # lossless, good for DEMs
            "tiled":     True,
            "blockxsize": 512,
            "blockysize": 512,
            "bigtiff":   "IF_SAFER",  # handles files > 4 GB automatically
        })

        output_path = Path(self.output)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        with rasterio.open(output_path, "w", **meta) as dst:
            dst.write(mosaic)

        for ds in datasets:
            ds.close()

        size_mb = output_path.stat().st_size / (1024 ** 2)
        log.info(f"Merged GeoTIFF saved: {output_path} ({size_mb:.1f} MB)")
        return str(output_path)


# ---------------------------------------------------------------------------
# CLI convenience
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Download and merge DEM tiles from OpenTopography.")
    parser.add_argument("--dem",     default="srtm30",         help=f"DEM key. Choices: {list(DEM_REGISTRY.keys())}")
    parser.add_argument("--south",   type=float, required=True)
    parser.add_argument("--north",   type=float, required=True)
    parser.add_argument("--west",    type=float, required=True)
    parser.add_argument("--east",    type=float, required=True)
    parser.add_argument("--output",  default="output_dem.tif",  help="Output merged GeoTIFF path")
    parser.add_argument("--tile-deg",    type=float, default=1.0,   help="Tile size in degrees (default: 1.0)")
    parser.add_argument("--workers",     type=int,   default=6,     help="Parallel download threads (default: 6)")
    parser.add_argument("--retries",     type=int,   default=5,     help="Max retries per tile (default: 5)")
    parser.add_argument("--cache-dir",   default="./dem_cache",     help="Directory to cache downloaded tiles")
    parser.add_argument("--api-key",     default=None,              help="OpenTopography API key (optional if stored in file/env)")
    args = parser.parse_args()

    dl = DEMDownloader(
        dem_key     = args.dem,
        south       = args.south,
        north       = args.north,
        west        = args.west,
        east        = args.east,
        output      = args.output,
        tile_deg    = args.tile_deg,
        max_workers = args.workers,
        max_retries = args.retries,
        cache_dir   = args.cache_dir,
        api_key     = args.api_key,
    )
    dl.run()