"""Code map: fmri_data.descriptives — summary stats and coverage maps."""

from codemap_lib import Slide


def build() -> Slide:
    s = Slide(
        title="CANlab code map: fmri_data.descriptives",
        description=[
            "Summary statistics for an fmri_data / image_vector object: min/max/percentiles, valid-voxel counts,",
            "and coverage maps showing how many images contribute to each voxel.",
        ],
    )

    s.use([
        "Use as a quick QC check on a freshly loaded dataset before running models.",
        "Returns coverage as fmri_data summary objects you can pass to montage() or surface().",
        "Pass 'noverbose' to suppress the printed summary; 'plotcoverage' to render coverage maps inline.",
    ])

    s.notes([
        "By convention, zero is treated as missing/empty data, not a valid value.",
        "coverage_obj_binned is useful when most images cover most voxels — it highlights the rare gaps.",
    ], y=3.40, h=1.40)

    s.legend(x=0.30, y=5.20)

    spine_x = 5.20
    y0 = 1.30
    gap = 0.95
    s.spine([
        ("indat", "var", "obj.dat\n(voxels × images)"),
        ("desc", "func", "descriptives( )"),
        ("desc_struct", "var", "desc struct"),
    ], x=spine_x, y0=y0, row_gap=gap, box_h=0.65, bold_ids=["desc"])

    # Outputs fanned across the bottom (kept inside slide bounds)
    out_y = y0 + 3 * gap + 0.10
    outs = [
        ("cov_all", "coverage_obj\n(fmri_data)"),
        ("cov_binned", "coverage_obj_binned"),
        ("cov_complete", "coverage_obj_complete"),
        ("nums", "min/max/percentiles\nnonempty / complete"),
    ]
    cx0 = 2.80
    col_w = 2.40
    col_gap = 0.20
    for i, (nid, label) in enumerate(outs):
        s.box(nid, "var", label,
              x=cx0 + i * (col_w + col_gap), y=out_y, h=0.65, w=col_w, font_size=10.0)
        s.connect_line("desc_struct", nid, src_side="bottom", dst_side="top")

    s.object_types_panel(
        x=spine_x + 2.50, y=y0,
        lines=["fmri_data", "image_vector"],
        width=2.10,
    )
    return s
