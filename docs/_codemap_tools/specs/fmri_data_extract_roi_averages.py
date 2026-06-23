"""Code map: fmri_data.extract_roi_averages — ROI / atlas extraction with optional pattern expression."""

from codemap_lib import Slide


def build() -> Slide:
    s = Slide(
        title="CANlab code map: fmri_data.extract_roi_averages",
        description=[
            "Average data within each ROI of a mask or atlas image. Supports contiguous-cluster and atlas-code modes,",
            "and can return mean activity, dot-product / cosine pattern expression, or correlation per ROI.",
        ],
    )

    s.use([
        "Pass an fmri_data plus a mask filename, fmri_mask_image, atlas, or fmri_data with integer codes.",
        "Default 'unique_mask_values' averages each integer code; 'contiguous_regions' uses 0/NaN-bounded blobs.",
        "Add 'pattern_expression' for weighted averages, plus 'cosine_similarity' or 'correlation' to switch metric.",
    ])

    s.notes([
        "Returns a region object (cl) — pass it to montage(), table(), or surface() for downstream use.",
        "When the mask is a filename and the data needs reslicing, the mask is reloaded from disk; manual mask thresholds are not preserved (a feature of map_to_image_space).",
    ], y=3.40, h=1.60)

    s.legend(x=0.30, y=5.40)

    spine_x = 5.30
    y0 = 1.20
    gap = 0.85
    s.spine([
        ("indat", "var", "obj.dat\n(voxels × images)"),
        ("extract", "func", "extract_roi_averages( )"),
        ("resample", "func", "resample_to_image_space\n(if needed)"),
        ("reduce", "func", "average / pattern\nexpression per ROI"),
    ], x=spine_x, y0=y0, row_gap=gap, box_h=0.65, bold_ids=["extract"])

    # Mask input from the right
    s.box("mask", "var", "mask / atlas image\n(file or object)",
          x=spine_x + 2.60, y=y0, h=0.65)
    s.connect_line("mask", "extract", src_side="bottomleft", dst_side="topright")

    # Mode selector (placed on right under mask, out of Use column)
    s.box("mode", "var", "'unique_mask_values' /\n'contiguous_regions'",
          x=spine_x + 2.60, y=y0 + gap - 0.05, h=0.65, w=2.60, font_size=10.0)
    s.connect_line("mode", "extract", src_side="left", dst_side="right")

    # Outputs (kept inside slide, clear of Use/Notes column)
    out_y = y0 + 4 * gap + 0.05
    s.box("cl", "var", "cl\n(region object)",
          x=2.90, y=out_y, h=0.65, w=2.40)
    s.box("roisum", "var", "cl_roisum\n(weighted)",
          x=5.50, y=out_y, h=0.65, w=2.40)
    s.box("dempat", "var", "cl_demeanedpattern\n(mean-centered weights)",
          x=8.10, y=out_y, h=0.65, w=2.80, font_size=10.0)
    s.connect_line("reduce", "cl", src_side="bottomleft", dst_side="top")
    s.connect_line("reduce", "roisum", src_side="bottom", dst_side="top")
    s.connect_line("reduce", "dempat", src_side="bottomright", dst_side="top")

    s.object_types_panel(
        x=spine_x + 5.50, y=y0,
        lines=["fmri_data", "atlas", "region"],
        width=2.10,
    )
    return s
