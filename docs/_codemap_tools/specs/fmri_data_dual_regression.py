"""Code map: fmri_data.dual_regression — FSL-style dual regression on group ICA components."""

from codemap_lib import Slide


def build() -> Slide:
    s = Slide(
        title="CANlab code map: fmri_data.dual_regression",
        description=[
            "Dual regression of group spatial maps onto a 4-D fmri_data object — yields per-image",
            "spatial maps and timecourses (FSL convention; thin wrapper around dual_regression_fsl).",
        ],
    )

    s.use([
        "Call as dual_regression(group_map_obj, data_obj). The first arg holds the group ICA / template maps.",
        "data_obj is automatically resampled into the group_map space before regression.",
        "Returns spatial maps (fmri_data), timecourses (component × timepoint), and z-scored t-maps (fmri_data).",
    ])

    s.notes([
        "Stage 1: regress spatial maps onto each volume to get timecourses. Stage 2: regress those timecourses onto the data to get subject-specific spatial maps.",
        "All work happens inside dual_regression_fsl; this method only does space alignment and dispatch.",
    ], y=3.40, h=1.55)

    s.legend(x=0.30, y=5.55)

    spine_x = 5.00
    y0 = 1.20
    gap = 0.85
    s.spine([
        ("group", "var", "group_map_obj\n(ICA / template maps)"),
        ("dr", "func", "dual_regression( )"),
        ("resample", "func", "resample_space\n(data → group)"),
        ("fsl", "func", "dual_regression_fsl"),
    ], x=spine_x, y0=y0, row_gap=gap, box_h=0.65, bold_ids=["dr"])

    # data_obj on the right
    s.box("data", "var", "data_obj\n(4-D fmri_data)",
          x=spine_x + 2.60, y=y0, h=0.65)
    s.connect_line("data", "dr", src_side="bottomleft", dst_side="topright")

    # Outputs — kept inside slide bounds and clear of Use/Notes column (x ≥ 2.90)
    out_y = y0 + 4 * gap + 0.10
    s.box("smaps", "var", "spatial_maps\n(fmri_data)",
          x=2.90, y=out_y, h=0.65, w=2.80)
    s.box("tc", "var", "timecourses\n[components × T]",
          x=5.95, y=out_y, h=0.65, w=2.80)
    s.box("tmaps", "var", "tmaps\n(z-scored fmri_data)",
          x=9.00, y=out_y, h=0.65, w=2.80)
    for nid in ("smaps", "tc", "tmaps"):
        s.connect_line("fsl", nid, src_side="bottom", dst_side="top")

    s.object_types_panel(
        x=spine_x + 4.85, y=y0,
        lines=["fmri_data", "image_vector"],
        width=2.10,
    )
    return s
