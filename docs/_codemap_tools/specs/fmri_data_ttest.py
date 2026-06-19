"""Code map: fmri_data.ttest — voxelwise one-sample t-test."""

from codemap_lib import Slide


def build() -> Slide:
    s = Slide(
        title="CANlab code map: fmri_data.ttest",
        description=[
            "Voxelwise one-sample t-test across images of an fmri_data object.",
            "Returns a statistic_image with t, p, and (optionally) thresholded sig vector.",
        ],
    )

    s.use([
        "Run on a 2nd-level fmri_data (one image per subject / contrast).",
        "Pass an optional p-value and threshold type ('unc' / 'fwe' / 'fdr') to threshold immediately.",
        "Output statistic_image plugs straight into threshold(), region(), montage(), surface(), and table().",
    ])

    s.notes([
        "ttest() ignores .X / .Y — it is a true one-sample test on the columns of .dat.",
        "For paired or two-sample tests, build contrast images first and re-run ttest on them.",
    ], y=3.40, h=1.30)

    s.legend(x=0.30, y=5.10)

    spine_x = 4.90
    y0 = 1.20
    gap = 0.85
    s.spine(
        [
            ("nii", "file", "*.nii / *.nii.gz"),
            ("ctor", "func", "fmri_data( )"),
            ("ttest", "func", "ttest( )"),
            ("stats", "var", "statistic_image"),
            ("threshold", "func", "threshold( )"),
            ("region", "func", "region( )"),
        ],
        x=spine_x, y0=y0, row_gap=gap, bold_ids=["ttest"],
    )

    # Branch right off statistic_image → display
    s.box("crf", "func", "canlab_results_\nfmridisplay",
          x=spine_x + 2.20, y=y0 + 3 * gap - 0.05, h=0.55)
    s.connect_line("stats", "crf", src_side="right", dst_side="left")
    s.box("display", "var", "fmridisplay\nobject",
          x=spine_x + 2.20, y=y0 + 4 * gap, h=0.55)
    s.connect_down("crf", "display")

    # Branch right off region → table
    s.box("table", "func", "table( )", x=spine_x + 2.20, y=y0 + 5 * gap)
    s.connect_line("region", "table", src_side="right", dst_side="left")

    # Branch left off region → montage
    s.box("montage", "func", "montage( )", x=spine_x - 2.20, y=y0 + 5 * gap)
    s.connect_line("region", "montage", src_side="left", dst_side="right")

    s.object_types_panel(
        x=spine_x + 4.50, y=1.15,
        lines=["fmri_data", "statistic_image", "region", "fmridisplay"],
        width=2.10,
    )
    return s


if __name__ == "__main__":
    import os
    import sys
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    from codemap_lib import render_pptx_to_png
    here = os.path.dirname(os.path.abspath(__file__))
    pptx_dir = os.path.normpath(os.path.join(here, "..", "..", "code_maps_pptx"))
    png_dir = os.path.normpath(os.path.join(here, "..", "..", "code_maps_png"))
    out = os.path.join(pptx_dir, "fmri_data_ttest_codemap.pptx")
    s = build()
    s.save(out)
    print("Wrote", out)
    print("Wrote", render_pptx_to_png(out, png_dir))
