"""Code map: fmri_data.normalize_gm_by_wm_csf — shift/scale gray matter by WM/CSF references."""

from codemap_lib import Slide


def build() -> Slide:
    s = Slide(
        title="CANlab code map: fmri_data.normalize_gm_by_wm_csf",
        description=[
            "Shift- and scale-normalize gray-matter voxels across subjects, using WM/CSF compartments as references.",
            "Removes a per-subject additive shift from CSF/WM medians and a multiplicative scale from MAD ratios.",
        ],
    )

    s.use([
        "Run on a 2nd-level fmri_data (.dat is [voxels × subjects]).",
        "'do_scale', false applies shift only; 'log_scale', true uses log-regression for the scale model.",
        "Optional 'mask_files' to override canonical GM/WM/CSF masks (must be a 3-element cell, GM/WM/CSF order).",
        "Returns the normalized fmri_data plus a per-subject statistics table appended to obj.metadata_table.",
    ])

    s.notes([
        "Non-GM voxels are copied through unchanged — only GM voxels are normalized.",
        "Masks are resampled to the data space first, so input data does not need to be in canonical mask space.",
    ], y=3.40, h=1.50)

    s.legend(x=0.30, y=5.30)

    spine_x = 5.10
    y0 = 1.20
    gap = 0.80
    s.spine([
        ("indat", "var", "obj.dat\n[voxels × subjects]"),
        ("normalize", "func", "normalize_gm_by_wm_csf( )"),
        ("resample", "func", "resample_space\n(masks → data)"),
        ("estimate", "func", "median (CSF/WM)\n+ MAD ratios"),
        ("apply", "func", "normalize_gm_shift_scale\n(voxel-level)"),
    ], x=spine_x, y0=y0, row_gap=gap, box_h=0.62, bold_ids=["normalize"])

    # Mask files on the right
    s.box("masks", "file", "gray_matter_mask_sparse.img\ncanonical_white_matter.img\ncanonical_ventricles.img",
          x=spine_x + 2.60, y=y0 + 0.50, h=1.20, w=3.50, font_size=10.0)
    s.connect_line("masks", "resample", src_side="left", dst_side="right")

    # Outputs
    out_y = y0 + 5 * gap + 0.10
    s.box("obj_out", "var", "obj_out\n(fmri_data, GM normalized)",
          x=spine_x - 2.40, y=out_y, h=0.65, w=2.80, font_size=10.0)
    s.box("statstab", "var", "statstab\n(per-subject stats table)",
          x=spine_x + 1.00, y=out_y, h=0.65, w=2.80, font_size=10.0)
    s.connect_line("apply", "obj_out", src_side="bottomleft", dst_side="topright")
    s.connect_line("apply", "statstab", src_side="bottomright", dst_side="topleft")

    s.object_types_panel(
        x=spine_x + 4.85, y=y0,
        lines=["fmri_data", "MATLAB table"],
        width=2.10,
    )
    return s
