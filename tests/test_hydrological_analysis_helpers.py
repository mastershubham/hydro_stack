import sys
import types
import importlib
from types import SimpleNamespace

def test_get_utm_epsg_for_bbox_monkeypatched():
    # Monkeypatch pyproj.aoi and pyproj.database so function returns a known code.
    aoi_mod = types.SimpleNamespace(AreaOfInterest=lambda **kwargs: "AOI")
    db_mod = types.SimpleNamespace(query_utm_crs_info=lambda **kwargs: [SimpleNamespace(code=32633)])
    sys.modules['pyproj.aoi'] = aoi_mod
    sys.modules['pyproj.database'] = db_mod

    # Import the function freshly so the local imports inside the function pick up our fakes
    import importlib
    import hydrological_analysis as ha
    importlib.reload(ha)

    epsg = ha.get_utm_epsg_for_bbox((76.0, 28.0, 80.0, 31.0))
    assert epsg == 32633


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
