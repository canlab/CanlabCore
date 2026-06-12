"""Code map: fmri_data.apply_parcellation — per-parcel means and pattern expression."""

from codemap_lib import Slide


def build() -> Slide:
    s = Slide(
        title="CANlab code map: fmri_data.apply_parcellation",
        description=[
            "Compute per-parcel mean values for a data object using a parcellation/atlas.",
            "Optionally apply a multivariate pattern within each parcel (local pattern expression, correlation, or cosine).",
        ],
    )

    s.use([
        "Pass an fmri_data and a parcellation (preferably an atlas object — handles integer resampling cleanly).",
        "Add 'pattern_expression', fmri_data(pattern) to compute local pattern responses parcel-by-parcel.",
        "Pass 'correlation' or 'cosine_similarity' to switch from dot-product to other similarity metrics.",
        "Returns one row per image and one column per parcel; missing parcels are NaN.",
    ])

    s.notes([
        "Spaces are matched in this order: pattern → (data → atlas). Parcels lost to resampling become NaN.",
        "If the parcellation is plain fmri_data with integer codes (not an atlas), behavior is identical but resampling is less robust.",
    ], y=3.40, h=1.50)

    s.legend(x=0.30, y=5.30)

    spine_x = 5.20
    y0 = 1.20
    gap = 0.85
    s.spine([
        ("indat", "var", "obj.dat\n(images × voxels)"),
        ("apply", "func", "apply_parcellation( )"),
        ("resample", "func", "resample_space\n(pattern → data → atlas)"),
        ("loop", "func", "per-parcel reduction"),
    ], x=spine_x, y0=y0, row_gap=gap, box_h=0.65, bold_ids=["apply"])

    # Atlas input from the right
    s.box("atlas", "var", "atlas / parcels\n(integer-coded)",
          x=spine_x + 2.60, y=y0, h=0.65)
    s.connect_line("atlas", "apply", src_side="bottomleft", dst_side="topright")

    # Optional pattern input from the right
    s.box("pattern", "var", "pattern\n(optional fmri_data)",
          x=spine_x + 2.60, y=y0 + gap - 0.05, h=0.65)
    s.connect_line("pattern", "apply", src_side="left", dst_side="right")

    # Outputs
    out_y = y0 + 4 * gap + 0.10
    s.box("means", "var", "parcel_means\n[images × parcels]",
          x=spine_x - 2.80, y=out_y, h=0.65)
    s.box("patexpr", "var", "parcel_pattern_expression\n[images × parcels]",
          x=spine_x - 0.30, y=out_y, h=0.65, w=2.40)
    s.box("misc", "var", "parcel_valence,\nrmsv_pos/neg, voxel_count",
          x=spine_x + 2.40, y=out_y, h=0.65, w=2.50, font_size=10.0)
    s.connect_line("loop", "means", src_side="bottomleft", dst_side="top")
    s.connect_line("loop", "patexpr", src_side="bottom", dst_side="top")
    s.connect_line("loop", "misc", src_side="bottomright", dst_side="top")

    s.object_types_panel(
        x=spine_x + 4.85, y=y0 + 2 * gap,
        lines=["fmri_data", "atlas"],
        width=2.10,
    )
    return s
