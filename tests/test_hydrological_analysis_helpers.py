def test_get_utm_epsg_for_bbox_without_candidate_area_of_use(monkeypatch):
    from types import SimpleNamespace

    import hydrological_analysis as ha

    monkeypatch.setattr(
        ha,
        "query_utm_crs_info",
        lambda **kwargs: [SimpleNamespace(code=32644)],
    )

    epsg = ha.get_utm_epsg_for_bbox((77.0, 28.0, 80.0, 31.0))
    assert epsg == 32644


def test_update_watershed_attributes_batches_database_updates():
    import hydrological_analysis as ha

    class FakeGrass:
        def __init__(self):
            self.commands = []

        def write_command(self, command, **kwargs):
            self.commands.append((command, kwargs))

    gs = FakeGrass()
    ha.update_watershed_attributes(
        gs,
        edges=[(1, 2), (3, 2)],
        basin_centroids={1: (0, 0), 2: (1, 1), 3: (2, 0)},
        basin_ids=[1, 2, 3],
    )

    assert len(gs.commands) == 1
    command, kwargs = gs.commands[0]
    assert command == "db.execute"
    assert kwargs["input"] == "-"
    sql = kwargs["stdin"]
    assert sql.startswith("BEGIN TRANSACTION;\n")
    assert sql.count("UPDATE watersheds_vect SET") == 3
    assert "downstream_id = 2" in sql
    assert "upstream_ids = '1,3'" in sql
    assert "flow_direction = 45.0" in sql
    assert sql.endswith("COMMIT;\n")


def test_filter_vector_file_by_boundary_keeps_overlapping_features_whole(
    tmp_path, monkeypatch
):
    import geopandas as gpd
    import json
    from pathlib import Path
    from shapely.geometry import LineString, box

    import hydrological_analysis as ha

    boundary_path = tmp_path / "watershed.geojson"
    output_path = tmp_path / "streams.geojson"
    crs = "EPSG:4326"
    gpd.GeoDataFrame(geometry=[box(0, 0, 2, 2)], crs=crs).to_file(boundary_path)
    output_path.write_text(json.dumps({
        "type": "FeatureCollection",
        "features": [
            {
                "type": "Feature",
                "properties": {"stream_id": 1},
                "geometry": {
                    "type": "LineString",
                    "coordinates": [[-1, 1], [3, 1]],
                },
            },
            {
                "type": "Feature",
                "properties": {"stream_id": 2},
                "geometry": {
                    "type": "LineString",
                    "coordinates": [[3, 3]],
                },
            },
        ],
    }))

    read_file = gpd.read_file

    def read_without_output_crs(path, *args, **kwargs):
        frame = read_file(path, *args, **kwargs)
        if Path(path) == output_path:
            return frame.set_crs(None, allow_override=True)
        return frame

    monkeypatch.setattr(gpd, "read_file", read_without_output_crs)
    ha.filter_vector_file_by_boundary(output_path, str(boundary_path))

    monkeypatch.setattr(gpd, "read_file", read_file)
    clipped = gpd.read_file(output_path)
    assert clipped["stream_id"].tolist() == [1]
    assert clipped.crs is not None
    assert clipped.geometry.iloc[0].equals(LineString([(-1, 1), (3, 1)]))


def test_filter_vector_file_by_boundary_assigns_missing_boundary_crs(
    tmp_path, monkeypatch
):
    import geopandas as gpd
    import json
    from pathlib import Path
    from shapely.geometry import LineString, box

    import hydrological_analysis as ha

    boundary_path = tmp_path / "watershed.geojson"
    output_path = tmp_path / "streams.geojson"
    gpd.GeoDataFrame(geometry=[box(0, 0, 2, 2)], crs="EPSG:4326").to_file(
        boundary_path
    )
    output_path.write_text(json.dumps({
        "type": "FeatureCollection",
        "features": [{
            "type": "Feature",
            "properties": {"stream_id": 1},
            "geometry": {
                "type": "LineString",
                "coordinates": [[-1, 1], [3, 1]],
            },
        }],
    }))

    read_file = gpd.read_file

    def read_without_boundary_crs(path, *args, **kwargs):
        frame = read_file(path, *args, **kwargs)
        if Path(path) == boundary_path:
            return frame.set_crs(None, allow_override=True)
        return frame

    monkeypatch.setattr(gpd, "read_file", read_without_boundary_crs)
    ha.filter_vector_file_by_boundary(output_path, str(boundary_path))

    monkeypatch.setattr(gpd, "read_file", read_file)
    filtered = gpd.read_file(output_path)
    assert filtered["stream_id"].tolist() == [1]
    assert filtered.geometry.iloc[0].equals(LineString([(-1, 1), (3, 1)]))


def test_dem_preprocessing_limits_depressions_to_watershed(tmp_path, monkeypatch):
    import geopandas as gpd
    import numpy as np
    import rasterio
    from rasterio.transform import from_origin
    from shapely.geometry import box
    from types import SimpleNamespace

    import hydrological_analysis as ha

    crs = "EPSG:32633"
    input_dem = tmp_path / "dem.tif"
    profile = {
        "driver": "GTiff",
        "height": 4,
        "width": 4,
        "count": 1,
        "dtype": "float32",
        "crs": crs,
        "transform": from_origin(0, 4, 1, 1),
        "nodata": -9999,
    }
    with rasterio.open(input_dem, "w", **profile) as dst:
        dst.write(np.full((1, 4, 4), 100, dtype=np.float32))

    buffer_path = tmp_path / "buffer.geojson"
    watershed_path = tmp_path / "watershed.geojson"
    gpd.GeoDataFrame(geometry=[box(0, 0, 4, 4)], crs=crs).to_file(buffer_path)
    gpd.GeoDataFrame(geometry=[box(1, 1, 3, 3)], crs=crs).to_file(watershed_path)

    class FakeRaster:
        def __init__(self, data, no_data):
            self.data = np.asarray(data)
            self.no_data = no_data
            self.geotransform = None
            self.projection = None

        def __array__(self, dtype=None, copy=None):
            return np.asarray(self.data, dtype=dtype)

        def __getitem__(self, key):
            return self.data[key]

        def __setitem__(self, key, value):
            self.data[key] = value

    def load_gdal(path, no_data):
        with rasterio.open(path) as src:
            return FakeRaster(src.read(1), no_data)

    def save_gdal(path, raster):
        with rasterio.open(input_dem) as src:
            output_profile = src.profile.copy()
        with rasterio.open(path, "w", **output_profile) as dst:
            dst.write(np.asarray(raster, dtype=np.float32), 1)

    fake_richdem = SimpleNamespace(
        LoadGDAL=load_gdal,
        FillDepressions=lambda raster, in_place: FakeRaster(
            np.where(raster.data == raster.no_data, raster.no_data, raster.data + 2),
            raster.no_data,
        ),
        rdarray=lambda data, no_data: FakeRaster(data, no_data),
        SaveGDAL=save_gdal,
        ResolveFlats=lambda raster, in_place: None,
    )
    monkeypatch.setattr(ha, "rd", fake_richdem)

    ha.dem_preprocessing(
        str(input_dem),
        tmp_path,
        watershed_path=str(buffer_path),
        depression_boundary_path=str(watershed_path),
    )

    with rasterio.open(tmp_path / "natural_depressions.tif") as src:
        depressions = src.read(1)
    with rasterio.open(tmp_path / "dem_conditioned.tif") as src:
        conditioned = src.read(1)

    assert np.all(depressions[1:3, 1:3] == 2)
    assert np.all(depressions[[0, 3], :] == profile["nodata"])
    assert np.all(depressions[:, [0, 3]] == profile["nodata"])
    assert np.all(conditioned != profile["nodata"])


def test_compute_mws_connectivity_exports_attribute_only_geojson(tmp_path, monkeypatch):
    import json
    import sys
    import types

    import numpy as np

    import hydrological_analysis as ha

    basin_array = np.array([[1, 2], [1, 2]])
    flow_direction_array = np.array([[8, 0], [8, 0]])

    class FakeGrassScript(types.ModuleType):
        commands = []

        @staticmethod
        def region():
            return {
                "rows": 2,
                "cols": 2,
                "n": 2,
                "w": 0,
                "nsres": 1,
                "ewres": 1,
            }

        @staticmethod
        def run_command(*args, **kwargs):
            FakeGrassScript.commands.append((args, kwargs))

    FakeGrassScript.commands.clear()
    fake_script = FakeGrassScript("grass.script")
    fake_script.array = types.SimpleNamespace(
        array=lambda name, null: (
            basin_array if name == "micro_watersheds_int"
            else flow_direction_array
        )
    )
    fake_grass = types.ModuleType("grass")
    fake_grass.script = fake_script
    monkeypatch.setitem(sys.modules, "grass", fake_grass)
    monkeypatch.setitem(sys.modules, "grass.script", fake_script)

    edges, _, _ = ha.compute_mws_connectivity(
        micro_watersheds_rast="microwatersheds",
        flow_dir_rast="flow_direction",
        pour_pts={1: (0, 0)},
    )

    export_dir = tmp_path / "exports"
    export_dir.mkdir()
    ha.export_outputs(
        export_dir,
        {},
        {},
        geojsons_to_export={
            "mws_connectivity": ha.mws_connectivity_geojson(edges),
        },
    )

    geojson_path = export_dir / "mws_connectivity.geojson"
    with geojson_path.open() as f:
        connectivity = json.load(f)

    assert edges == {(1, 2): True}
    assert len(connectivity["features"]) == 1
    assert connectivity["features"][0]["geometry"] is None
    assert connectivity["features"][0]["properties"] == {
        "from_basin_id": 1,
        "to_basin_id": 2,
    }
    assert not geojson_path.with_suffix(".csv").exists()
    assert all(command[0] == ("r.mapcalc",) for command in FakeGrassScript.commands)
