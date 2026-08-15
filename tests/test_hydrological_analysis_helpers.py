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
