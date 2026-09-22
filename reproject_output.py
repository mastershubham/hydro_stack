#!/usr/bin/env python3
"""Reproject all generated raster and vector layers in a directory to EPSG:4326."""

from __future__ import annotations

import argparse
import shutil
import subprocess
from pathlib import Path

RASTER_EXTENSIONS = {".tif", ".tiff"}
VECTOR_EXTENSIONS = {".geojson", ".gpkg", ".shp"}


def _require_tool(name: str) -> None:
    if shutil.which(name) is None:
        raise RuntimeError(f"Required tool not found on PATH: {name}")


def _run_command(command: list[str]) -> None:
    print("Running:", " ".join(command))
    subprocess.run(command, check=True)


def _reproject_raster(input_file: Path, output_dir: Path) -> Path:
    output_file = output_dir / f"{input_file.stem}_EPSG4326{input_file.suffix}"
    _run_command([
        "gdalwarp",
        "-overwrite",
        "-t_srs",
        "EPSG:4326",
        str(input_file),
        str(output_file),
    ])
    return output_file


def _reproject_vector(input_file: Path, output_dir: Path) -> Path:
    suffix = input_file.suffix.lower()
    if suffix == ".geojson":
        output_file = output_dir / f"{input_file.stem}_EPSG4326.geojson"
        driver = "GeoJSON"
    elif suffix == ".gpkg":
        output_file = output_dir / f"{input_file.stem}_EPSG4326.gpkg"
        driver = "GPKG"
    elif suffix == ".shp":
        output_file = output_dir / f"{input_file.stem}_EPSG4326.shp"
        driver = "ESRI Shapefile"
    else:
        raise ValueError(f"Unsupported vector extension: {suffix}")

    _run_command([
        "ogr2ogr",
        "-overwrite",
        "-f",
        driver,
        "-t_srs",
        "EPSG:4326",
        str(output_file),
        str(input_file),
    ])
    return output_file


def find_layers(input_dir: Path):
    rasters = sorted(
        p for p in input_dir.iterdir()
        if p.is_file() and p.suffix.lower() in RASTER_EXTENSIONS
    )
    vectors = sorted(
        p for p in input_dir.iterdir()
        if p.is_file() and p.suffix.lower() in VECTOR_EXTENSIONS
    )
    return rasters, vectors


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Reproject generated GIS layers in a directory to EPSG:4326."
    )
    parser.add_argument(
        "--input-dir",
        "-i",
        type=Path,
        required=True,
        help="Directory containing the generated raster/vector outputs to transform.",
    )
    parser.add_argument(
        "--output-dir",
        "-o",
        type=Path,
        default=None,
        help="Destination directory for the reprojected outputs. Defaults to a sibling directory named epsg4326 next to the input folder.",
    )
    args = parser.parse_args()

    input_dir = args.input_dir.resolve()
    if not input_dir.exists():
        raise FileNotFoundError(f"Input directory does not exist: {input_dir}")

    if args.output_dir is not None:
        output_dir = args.output_dir.resolve()
    else:
        output_dir = input_dir.parent / f"{input_dir.name}_epsg4326"
    output_dir.mkdir(parents=True, exist_ok=True)

    _require_tool("gdalwarp")
    _require_tool("ogr2ogr")

    rasters, vectors = find_layers(input_dir)

    if not rasters and not vectors:
        print(f"No raster or vector layers found in {input_dir}")
        return

    print(f"Found {len(rasters)} rasters and {len(vectors)} vectors in {input_dir}")

    for raster in rasters:
        reprojected = _reproject_raster(raster, output_dir)
        print(f"Raster reprojected: {raster} -> {reprojected}")

    for vector in vectors:
        reprojected = _reproject_vector(vector, output_dir)
        print(f"Vector reprojected: {vector} -> {reprojected}")

    print(f"Completed. Outputs written to: {output_dir}")


if __name__ == "__main__":
    main()
