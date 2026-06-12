"""Code map: fmri_data.extract_gray_white_csf — mean tissue signals + components."""

from codemap_lib import Slide


def build() -> Slide:
    s = Slide(
        title="CANlab code map: fmri_data.extract_gray_white_csf",
        description=[
            "Extract mean values and the top 5 principal-component scores from each of gray, white, and CSF",
            "compartments using canonical CANlab/SPM tissue masks. Useful as nuisance regressors and QC.",
        ],
    )

    s.use([
        "Run on an fmri_data in MNI space (the canonical masks are MNI templates).",
        "Default masks: gray_matter_mask.img, canonical_white_matter.img, canonical_ventricles.img.",
        "Pass 'eval', @fn to swap mean for any function-handle reduction (must accept one arg, NaN-safe).",
        "Useful as nuisance covariates in connectivity / 1st-level / 2nd-level GLMs.",
    ])

    s.notes([
        "WM/CSF masks are conservatively eroded so they are unlikely to bleed signal from gray matter.",
        "l2norms are length-adjusted (divided by sqrt(nvox)) — comparable across compartments of different size.",
    ], y=3.40, h=1.50)

    s.legend(x=0.30, y=5.30)

    spine_x = 5.30
    y0 = 1.30
    gap = 0.95
    s.spine([
        ("indat", "var", "obj.dat\n(MNI space)"),
        ("extract", "func", "extract_gray_white_csf( )"),
        ("loop", "func", "per-tissue reduction\n(mean / pca / 'eval')"),
    ], x=spine_x, y0=y0, row_gap=gap, box_h=0.65, bold_ids=["extract"])

    # Mask files on the right
    s.box("masks", "file", "gray_matter_mask.img\ncanonical_white_matter.img\ncanonical_ventricles.img",
          x=spine_x + 2.60, y=y0 + 0.40, h=1.10, w=3.40, font_size=10.0)
    s.connect_line("masks", "extract", src_side="left", dst_side="right")

    # Outputs (kept inside slide, clear of Use/Notes column)
    out_y = y0 + 3 * gap + 0.15
    s.box("values", "var", "values\n[obs × 3] tissue means",
          x=2.90, y=out_y, h=0.65, w=2.40, font_size=10.0)
    s.box("comp", "var", "components\n[obs × 5] per tissue",
          x=5.40, y=out_y, h=0.65, w=2.40, font_size=10.0)
    s.box("masked", "var", "full_data_objects\n{gray, white, CSF}",
          x=7.90, y=out_y, h=0.65, w=2.40, font_size=10.0)
    s.box("l2", "var", "l2norms\n(length-adjusted)",
          x=10.40, y=out_y, h=0.65, w=2.40, font_size=10.0)
    for nid in ("values", "comp", "masked", "l2"):
        s.connect_line("loop", nid, src_side="bottom", dst_side="top")

    s.object_types_panel(
        x=spine_x + 6.20, y=y0,
        lines=["fmri_data", "image_vector"],
        width=2.10,
    )
    return s
