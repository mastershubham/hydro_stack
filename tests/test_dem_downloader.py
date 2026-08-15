import os
import tempfile
import rasterio
import numpy as np
from pathlib import Path
from unittest.mock import patch
from dem_downloader import _tile_bbox, _load_api_key, DEMDownloader

def test_tile_bbox_simple():
    tiles = _tile_bbox(0.0, 2.0, 0.0, 2.0, tile_deg=1.0)
    # 2x2 degrees with 1° tiles => 4 tiles
    assert len(tiles) == 4
    assert (0.0, 1.0, 0.0, 1.0) in tiles

def test_load_api_key_from_env(monkeypatch):
    monkeypatch.setenv("OPENTOPOGRAPHY_API_KEY", "MY_KEY_123")
    assert _load_api_key(None) == "MY_KEY_123"

def test_cache_hit_returns_path(tmp_path):
    monkeypatch.setenv("OPENTOPOGRAPHY_API_KEY", "MY_KEY_123")
    # Prepare cache file
    cache_dir = tmp_path / "cache"
    cache_dir.mkdir()
    # Use typical dem filename pattern
    fname = cache_dir / "SRTMGL1_0.0_1.0_0.0_1.0.tif"
    fname.write_bytes(b"placeholder")
    dl = DEMDownloader(
        dem_key="srtm30",
        south=0,
        north=1,
        west=0,
        east=1,
        output=str(tmp_path / "out.tif"),
        cache_dir=str(cache_dir),
        max_retries=1,
        max_workers=1,
    )

    # _download_tile_with_retry may call _record_done() which expects _progress to exist.
    # _download_all normally sets this; tests calling the method directly should set it too.
    dl._progress = {"done": 0, "total": 1, "failed": 0}

    path = dl._download_tile_with_retry((0.0, 1.0, 0.0, 1.0))
    assert Path(path).exists()
    assert Path(path).samefile(fname)

def test_merge_two_small_rasters(tmp_path):
    monkeypatch.setenv("OPENTOPOGRAPHY_API_KEY", "MY_KEY_123")
    # create two tiny rasters and call _merge
    def make_raster(path, value, transform, crs="EPSG:4326"):
        data = np.full((1, 2, 2), value, dtype=np.float32)
        meta = {
            "driver": "GTiff",
            "dtype": "float32",
            "count": 1,
            "height": 2,
            "width": 2,
            "crs": crs,
            "transform": transform,
        }
        with rasterio.open(path, "w", **meta) as dst:
            dst.write(data)

    t1 = tmp_path / "t1.tif"
    t2 = tmp_path / "t2.tif"

    from rasterio.transform import from_origin
    transform = from_origin(0, 2, 1, 1)  # simple transform

    make_raster(t1, 1.0, transform)
    make_raster(t2, 2.0, transform)

    dl = DEMDownloader(
        dem_key="srtm30",
        south=0,
        north=1,
        west=0,
        east=1,
        output=str(tmp_path / "merged.tif"),
    )
    out = dl._merge([str(t1), str(t2)])
    assert Path(out).exists()
    with rasterio.open(out) as ds:
        arr = ds.read(1)
        # merged mosaic shape (height,width) equals 2x2 in this simple case
        assert arr.shape[0] >= 2
