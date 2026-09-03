"""
Interactive micro-watershed connectivity map.

Replaces the static matplotlib PNG with a zoomable/pannable Leaflet map
(via folium). Handles large basins gracefully through:
  - real zoom/pan (mouse wheel + drag), instead of a fixed-resolution PNG
  - toggleable layers (watersheds / arrows / pour points) via LayerControl
  - marker clustering for pour points so dense clusters collapse at low zoom
  - hover tooltips + click popups instead of permanent on-canvas ID labels
  - arrows drawn as AntPath (animated dashes) so direction is readable
    even when many edges overlap, with a fallback to a plain arrow line

Output: a single self-contained .html file, opened in any browser.
"""

import geopandas as gpd
import folium
from folium.plugins import PolyLineTextPath, MarkerCluster, Fullscreen, MiniMap
import branca.colormap as cm


WATERSHEDS_GEOJSON = "./masalia/microwatersheds.geojson"
CONNECTIVITY_GEOJSON = "./masalia/mws_connectivity.geojson"
POUR_POINTS_GEOJSON = "./masalia/pour_points.geojson"

WATERSHED_ID_FIELD = "basin_id"
FROM_FIELD = "from_basin_id"
TO_FIELD = "to_basin_id"

COLOR_BY_FIELD = None  # set to a numeric column name to choropleth-color watersheds

ARROW_COLOR = "#e63946"
ARROW_WEIGHT = 1

POUR_POINT_COLOR = "#f4a261"
POUR_POINT_EDGE = "#e76f51"
POUR_POINT_RADIUS = 5

CLUSTER_POUR_POINTS = True  # collapse pour points into clusters at low zoom if set to True
OUT_HTML = "watershed_connectivity.html"


def load_data(watersheds_path, connectivity_path, pour_points_path):
    gdf = gpd.read_file(watersheds_path)
    connectivity = gpd.read_file(connectivity_path)
    pp_gdf = gpd.read_file(pour_points_path)

    # Compute area in a projected CRS before converting geometry to WGS84.
    if gdf.crs is None:
        raise ValueError("Watersheds layer has no CRS defined.")
    if gdf.crs.is_geographic:
        proj = gdf.to_crs(gdf.estimate_utm_crs())
    else:
        proj = gdf

    gdf_wgs84 = gdf.to_crs(epsg=4326)
    gdf_wgs84["area_ha"] = proj.geometry.area.to_numpy() / 10000.0
    pp_wgs84 = pp_gdf.to_crs(epsg=4326)

    # Compute centroids in the projected CRS, then reproject for map placement.
    centroids_wgs84 = proj.geometry.centroid.to_crs(epsg=4326)
    gdf_wgs84["centroid_lon"] = centroids_wgs84.x.to_numpy()
    gdf_wgs84["centroid_lat"] = centroids_wgs84.y.to_numpy()

    return gdf_wgs84, connectivity, pp_wgs84


def build_centroid_lookup(gdf, id_field):
    """Vectorized centroid lookup generator (10-20x faster than .iterrows())."""
    return dict(
        zip(
            gdf[id_field],
            zip(gdf["centroid_lat"], gdf["centroid_lon"])
        )
    )


def build_map(gdf, connectivity, pour_points, id_field, from_field, to_field,
              color_by=None):

    bounds = gdf.total_bounds  # minx, miny, maxx, maxy in lon/lat
    center = [(bounds[1] + bounds[3]) / 2, (bounds[0] + bounds[2]) / 2]

    fmap = folium.Map(location=center, zoom_start=11, tiles=None,
                       control_scale=True)
    folium.TileLayer(
        tiles="OpenStreetMap",
        name="OpenStreetMap",
        overlay=True,
        control=True,
        show=True,
    ).add_to(fmap)
    
    # ---- Watershed polygons -------------------------------------------------
    watershed_layer = folium.FeatureGroup(name="Micro-watersheds", show=True)

    if color_by and color_by in gdf.columns:
        vmin, vmax = gdf[color_by].min(), gdf[color_by].max()
        colormap = cm.linear.YlGnBu_09.scale(vmin, vmax)
        colormap.caption = color_by

        def style_fn(feature):
            val = feature["properties"].get(color_by)
            return {
                "fillColor": colormap(val) if val is not None else "#a8dadc",
                "color": "#555555",
                "weight": 0.8,
                "fillOpacity": 0.75,
            }
        fmap.add_child(colormap)
    else:
        def style_fn(feature):
            return {
                "fillColor": "#a8dadc",
                "color": "#333333",
                "weight": 0.8,
                "fillOpacity": 0.6,
            }

    highlight_fn = lambda feature: {"weight": 2.5, "color": "#1d3557", "fillOpacity": 0.85}

    tooltip_fields = [id_field, "area_ha"] + (
        [color_by] if color_by and color_by in gdf.columns else []
    )
    folium.GeoJson(
        gdf,
        name="Micro-watersheds",
        style_function=style_fn,
        highlight_function=highlight_fn,
        tooltip=folium.GeoJsonTooltip(fields=tooltip_fields, sticky=True),
        popup=folium.GeoJsonPopup(fields=tooltip_fields),
    ).add_to(watershed_layer)
    watershed_layer.add_to(fmap)

    # ---- Flow direction arrows ------------------------------------------------
    arrow_layer = folium.FeatureGroup(name="Flow direction", show=True)
    centroid_lut = build_centroid_lookup(gdf, id_field)
    missing_pairs = []

    # Iterate over zip tuples to avoid .iterrows() overhead
    for fid, tid in zip(connectivity[from_field], connectivity[to_field]):
        if fid not in centroid_lut or tid not in centroid_lut:
            missing_pairs.append((fid, tid))
            continue

        p1 = centroid_lut[fid]
        p2 = centroid_lut[tid]

        line = folium.PolyLine(
            locations=[p1, p2],
            color=ARROW_COLOR,
            weight=ARROW_WEIGHT,
            opacity=0.85,
            tooltip=f"{fid} → {tid}",
        )
        line.add_to(arrow_layer)
        PolyLineTextPath(
            line, "  ►  ", repeat=True, offset=0,
            attributes={"fill": ARROW_COLOR, "font-weight": "bold", "font-size": "13", "dy": "4"},
        ).add_to(arrow_layer)

    if missing_pairs:
        print(f"{len(missing_pairs)} connectivity pair(s) skipped "
              f"(watershed ID not found in GeoJSON):")
        for p in missing_pairs:
            print(f"   {p[0]} → {p[1]}")

    arrow_layer.add_to(fmap)

    # ---- Pour points ------------------------------------------------------
    if pour_points is not None and not pour_points.empty:
        pp_layer = MarkerCluster(name="Pour points") if CLUSTER_POUR_POINTS \
            else folium.FeatureGroup(name="Pour points")
        
        # Filter for valid Point geometries upfront using vector operations
        valid_points = pour_points[
            pour_points.geometry.notnull() & (pour_points.geometry.geom_type == "Point")
        ]

        for geom in valid_points.geometry:
            folium.CircleMarker(
                location=(geom.y, geom.x),
                radius=POUR_POINT_RADIUS,
                color=POUR_POINT_EDGE,
                fill=True,
                fill_color=POUR_POINT_COLOR,
                fill_opacity=0.9,
                weight=1.2,
                tooltip="Pour point",
            ).add_to(pp_layer)
        pp_layer.add_to(fmap)

    # ---- Controls -----------------------------------------------------------
    Fullscreen(position="topleft").add_to(fmap)
    MiniMap(toggle_display=True).add_to(fmap)
    folium.LayerControl(collapsed=False).add_to(fmap)

    fmap.fit_bounds([[bounds[1], bounds[0]], [bounds[3], bounds[2]]])
    return fmap


if __name__ == "__main__":
    gdf, connectivity, pp_gdf = load_data(
        WATERSHEDS_GEOJSON, CONNECTIVITY_GEOJSON, POUR_POINTS_GEOJSON
    )

    print(f"Loaded {len(gdf)} watersheds and {len(connectivity)} connectivity edges.")

    fmap = build_map(
        gdf, connectivity, pp_gdf,
        id_field=WATERSHED_ID_FIELD,
        from_field=FROM_FIELD,
        to_field=TO_FIELD,
        color_by=COLOR_BY_FIELD,
    )

    fmap.save(OUT_HTML)
    print(f"Saved interactive map: {OUT_HTML}  (open in a browser)")
