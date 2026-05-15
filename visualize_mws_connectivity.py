import geopandas as gpd
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrowPatch
import numpy as np


WATERSHEDS_GEOJSON = "./masalia_v2/microwatersheds.geojson"  
CONNECTIVITY_GEOJSON = "./masalia_v2/mws_connectivity.geojson"     
POUR_POINTS_GEOJSON = "./masalia_v2/pour_points.geojson"

WATERSHED_ID_FIELD = "basin_id"      
FROM_FIELD         = "from_basin_id"   
TO_FIELD           = "to_basin_id"    

COLOR_BY_FIELD     = None

ARROW_COLOR        = "#e63946"
ARROW_LINEWIDTH    = 1.5
ARROW_SHRINK       = 10   
FIG_SIZE           = (14, 10)

POUR_POINT_COLOR   = "#f4a261"   
POUR_POINT_EDGE    = "#e76f51"   
POUR_POINT_SIZE    = 60          
POUR_POINT_MARKER  = "o"


def load_data(watersheds_path, connectivity_path, pour_points_path):
    gdf = gpd.read_file(watersheds_path)
    if gdf.crs and gdf.crs.is_geographic:
        gdf = gdf.to_crs(gdf.estimate_utm_crs())

    gdf["centroid_x"] = gdf.geometry.centroid.x
    gdf["centroid_y"] = gdf.geometry.centroid.y

    connectivity = gpd.read_file(connectivity_path)

    pp_gdf = gpd.read_file(pour_points_path)
    pp_gdf = pp_gdf.to_crs(gdf.crs)

    return gdf, connectivity, pp_gdf


def build_centroid_lookup(gdf, id_field):
    return {
        row[id_field]: (row["centroid_x"], row["centroid_y"])
        for _, row in gdf.iterrows()
    }


def plot_watersheds(gdf, connectivity, pour_points, id_field, from_field, to_field,
                    color_by=None, arrow_color=ARROW_COLOR):

    fig, ax = plt.subplots(figsize=FIG_SIZE)

    if color_by and color_by in gdf.columns:
        gdf.plot(ax=ax, column=color_by, cmap="YlGnBu", edgecolor="#555555",
                 linewidth=0.6, legend=True,
                 legend_kwds={"label": color_by, "shrink": 0.6})
        title_suffix = f" — coloured by {color_by}"
    else:
        gdf.plot(ax=ax, color="#a8dadc", edgecolor="#333333", linewidth=0.6)
        title_suffix = ""

    for _, row in gdf.iterrows():
        ax.annotate(
            str(row[id_field]),
            xy=(row["centroid_x"], row["centroid_y"]),
            ha="center", va="center",
            fontsize=7, fontweight="bold", color="#1d3557",
            bbox=dict(boxstyle="round,pad=0.2", fc="white", alpha=0.5, ec="none")
        )

    centroid_lut = build_centroid_lookup(gdf, id_field)
    missing_pairs = []

    for _, row in connectivity.iterrows():
        fid, tid = row[from_field], row[to_field]
        if fid not in centroid_lut or tid not in centroid_lut:
            missing_pairs.append((fid, tid))
            continue

        x1, y1 = centroid_lut[fid]
        x2, y2 = centroid_lut[tid]

        arrow = FancyArrowPatch(
            posA=(x1, y1), posB=(x2, y2),
            arrowstyle="-|>",
            mutation_scale=14,
            color=arrow_color,
            linewidth=ARROW_LINEWIDTH,
            shrinkA=ARROW_SHRINK,
            shrinkB=ARROW_SHRINK,
            zorder=5
        )
        ax.add_patch(arrow)

    if missing_pairs:
        print(f"{len(missing_pairs)} connectivity pair(s) skipped "
              f"(watershed ID not found in GeoJSON):")
        for p in missing_pairs:
            print(f"   {p[0]} → {p[1]}")
    
    if pour_points is not None and not pour_points.empty:
        pour_points.plot(
            ax=ax,
            marker=POUR_POINT_MARKER,
            color=POUR_POINT_COLOR,
            edgecolor=POUR_POINT_EDGE,
            markersize=POUR_POINT_SIZE,
            linewidth=0.8,
            zorder=6,
            label="Pour point"
        )

    import matplotlib.lines as mlines
    arrow_legend = mpatches.Patch(color=arrow_color, label="Flow direction")
    poly_legend  = mpatches.Patch(color="#a8dadc",   label="Micro-watershed")
    pp_legend     = mlines.Line2D([], [], color=POUR_POINT_COLOR,
                                  marker=POUR_POINT_MARKER, linestyle="None",
                                  markersize=8, markeredgecolor=POUR_POINT_EDGE,
                                  label="Pour point")
    ax.legend(handles=[poly_legend, arrow_legend, pp_legend],
              loc="lower right", fontsize=9)

    ax.set_title(f"Micro-Watershed Flow Connectivity{title_suffix}", fontsize=14, pad=12)
    ax.set_xlabel("Easting (m)")
    ax.set_ylabel("Northing (m)")
    ax.ticklabel_format(style="sci", axis="both", scilimits=(0, 0))
    plt.tight_layout()

    out_png = "watershed_connectivity.png"
    plt.savefig(out_png, dpi=150, bbox_inches="tight")
    print(f"Saved → {out_png}")
    plt.show()


if __name__ == "__main__":
    gdf, connectivity, pp_gdf = load_data(WATERSHEDS_GEOJSON, CONNECTIVITY_GEOJSON, POUR_POINTS_GEOJSON)

    print(f"Loaded {len(gdf)} watersheds and {len(connectivity)} connectivity edges.")
    print(f"Watershed CRS: {gdf.crs}")

    plot_watersheds(
        gdf, connectivity, pp_gdf,
        id_field=WATERSHED_ID_FIELD,
        from_field=FROM_FIELD,
        to_field=TO_FIELD,
        color_by=COLOR_BY_FIELD,
    )
