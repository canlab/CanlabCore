"""Code map: fmri_data.table_of_atlas_regions_covered — atlas coverage table."""

from codemap_lib import Slide


def build() -> Slide:
    s = Slide(
        title="CANlab code map: fmri_data.table_of_atlas_regions_covered",
        description=[
            "Identify which atlas parcels are covered by an activation map and report the percent coverage per parcel.",
            "Complements traditional contiguous-cluster tables for maps that span multiple anatomical regions.",
        ],
    )

    s.use([
        "Pass a thresholded fmri_data / statistic_image (single image only) and optionally an atlas-class object.",
        "'coverage', N sets the percent threshold (default 25%) for including an atlas parcel.",
        "Returns separate tables for parcels with positive- and negative-mean activation, plus the divided region object.",
    ])

    s.notes([
        "Internally uses the subdivide_by_atlas() method — one row per atlas parcel (not per blob).",
        "atlas_of_regions_covered is an atlas-class object you can re-use for visualization or further analysis.",
    ], y=3.40, h=1.50)

    s.legend(x=0.30, y=5.30)

    spine_x = 4.85
    y0 = 1.30
    gap = 0.90
    s.spine([
        ("indat", "var", "thresholded\nfmri_data / statistic_image"),
        ("toa", "func", "table_of_atlas_regions_covered( )"),
        ("subdiv", "func", "subdivide_by_atlas\n(per-parcel coverage)"),
    ], x=spine_x, y0=y0, row_gap=gap, box_w=2.85, box_h=0.65, bold_ids=["toa"])

    # Atlas input (placed on right so it doesn't collide with object_types panel)
    s.box("atlas", "var", "atlas object\n(default if omitted)",
          x=spine_x + 3.50, y=y0 + gap, h=0.65, w=2.50, font_size=10.0)
    s.connect_line("atlas", "toa", src_side="left", dst_side="right")

    # Outputs across the bottom — kept inside slide bounds (x ≥ 2.90, ≤ 13.00)
    # and clear of the Use/Notes column.
    out_y = y0 + 3 * gap + 0.30
    outs = [
        ("pos", "results_table_pos\n(positive parcels)"),
        ("neg", "results_table_neg\n(negative parcels)"),
        ("rgn", "r\n(region object,\nsubdivided)"),
        ("excl", "excluded_region_\ntable"),
        ("aoa", "atlas_of_regions_\ncovered"),
    ]
    cx0 = 2.90
    col_w = 1.95
    col_gap = 0.10
    for i, (nid, label) in enumerate(outs):
        s.box(nid, "var", label,
              x=cx0 + i * (col_w + col_gap), y=out_y, h=0.85, w=col_w, font_size=9.5)
        s.connect_line("subdiv", nid, src_side="bottom", dst_side="top")

    s.object_types_panel(
        x=spine_x + 3.50, y=y0 + 2 * gap + 0.15,
        lines=["fmri_data", "statistic_image", "atlas", "region", "MATLAB table"],
        width=2.50,
    )
    return s
