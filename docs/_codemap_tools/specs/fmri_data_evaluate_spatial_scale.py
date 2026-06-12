"""Code map: fmri_data.evaluate_spatial_scale — info-coding scale across parcels."""

from codemap_lib import Slide


def build() -> Slide:
    s = Slide(
        title="CANlab code map: fmri_data.evaluate_spatial_scale",
        description=[
            "Evaluate the spatial scale of information coding by training predictive models on random voxel draws",
            "from each parcel and at the whole-brain level. Returns fitted curves comparing parcel vs. whole-brain.",
        ],
    )

    s.use([
        "Pass an fmri_data with .Y set as the outcome to predict, plus a parcellation atlas.",
        "Optional 'num_voxels', vec controls the random-sample sizes (must be smaller than the smallest parcel).",
        "Pass 'do_plot' to render the sampling curves and 'render_brain' to overlay the parcel weights.",
    ])

    s.notes([
        "Internally calls predict() many times — once per parcel × sample-size × repeat, plus whole-brain runs.",
        "Compare parcel-level curves to the whole-brain curve to see whether information is local, global, or distributed.",
    ], y=3.40, h=1.50)

    s.legend(x=0.30, y=5.30)

    spine_x = 5.00
    y0 = 1.20
    gap = 0.85
    s.spine([
        ("indat", "var", "data_obj.dat + .Y"),
        ("eval", "func", "evaluate_spatial_scale( )"),
        ("loop", "func", "loop: parcel × n_voxels × repeat"),
        ("predict", "func", "predict( )\n(per draw)"),
    ], x=spine_x, y0=y0, row_gap=gap, box_h=0.65, bold_ids=["eval"])

    # Parcel input from the right
    s.box("parcels", "var", "parcel_obj\n(atlas)",
          x=spine_x + 2.60, y=y0, h=0.65)
    s.connect_line("parcels", "eval", src_side="bottomleft", dst_side="topright")

    # Outputs
    out_y = y0 + 4 * gap + 0.15
    s.box("fit", "var", "fitresult\n(per-parcel fits)",
          x=spine_x - 2.80, y=out_y, h=0.65)
    s.box("gof", "var", "gof\n(goodness-of-fit)",
          x=spine_x - 0.30, y=out_y, h=0.65)
    s.box("stats", "var", "stats\n(per-draw accuracy)",
          x=spine_x + 2.20, y=out_y, h=0.65)
    for nid in ("fit", "gof", "stats"):
        s.connect_line("predict", nid, src_side="bottom", dst_side="top")

    # Plot side-effect
    s.box("plot", "script", "sampling-curve plot\n(if 'do_plot')",
          x=spine_x + 2.60, y=y0 + 2 * gap, h=0.65, w=2.50, font_size=10.0)
    s.connect_line("loop", "plot", src_side="right", dst_side="left")

    s.object_types_panel(
        x=spine_x + 4.85, y=y0 + 4 * gap,
        lines=["fmri_data", "atlas"],
        width=2.10,
    )
    return s
