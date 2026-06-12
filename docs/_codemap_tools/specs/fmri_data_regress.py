"""Code map: fmri_data.regress — voxelwise multiple regression."""

from codemap_lib import Slide


def build() -> Slide:
    s = Slide(
        title="CANlab code map: fmri_data.regress",
        description=[
            "Voxelwise multiple regression on an fmri_data object using a design matrix attached as obj.X.",
            "Returns statistic_image objects with betas, t-statistics, and p-values for each regressor.",
        ],
    )

    s.use([
        "Attach a design matrix to the fmri_data object as obj.X (one column per regressor; intercept added if missing).",
        "Optional name-value pairs control thresholding, robust regression, contrasts, and voxel-varying covariates.",
        "Output statistic_image objects can be passed to threshold(), region(), montage(), or surface().",
    ])

    s.notes([
        "Call as DAT.regress() or regress(DAT, ...). 'robust' switches OLS to robustfit() per voxel — slower.",
        "For group analyses, run per subject then pass beta maps to a 2nd-level call.",
    ], y=3.50, h=1.40)

    s.legend(x=0.30, y=5.30)

    # ---- Main spine (left of center, vertical) ----
    spine_x = 4.80
    row_y = [1.20, 1.95, 2.70, 3.45, 4.20, 4.95, 5.70]

    s.box("nii",          "file", "*.nii / *.nii.gz",       x=spine_x, y=row_y[0])
    s.box("ctor",         "func", "fmri_data( )",            x=spine_x, y=row_y[1])
    s.box("apply_mask",   "func", "apply_mask",              x=spine_x, y=row_y[2])
    s.box("regress",      "func", "regress( )",              x=spine_x, y=row_y[3], bold=True)
    s.box("stats_img",    "var",  "statistic_image",         x=spine_x, y=row_y[4])
    s.box("threshold",    "func", "threshold( )",            x=spine_x, y=row_y[5])
    s.box("region",       "func", "region( )",               x=spine_x, y=row_y[6])

    for a, b in [("nii", "ctor"), ("ctor", "apply_mask"),
                 ("apply_mask", "regress"), ("regress", "stats_img"),
                 ("stats_img", "threshold"), ("threshold", "region")]:
        s.connect_down(a, b)

    # ---- Left-side branch off region: table ----
    s.box("table", "func", "table( )", x=spine_x - 2.20, y=row_y[6])
    s.connect_left("region", "table")

    # ---- Right-side branch: design matrix → obj.X → regress ----
    right_x = spine_x + 2.20
    s.box("design_file", "file", "design matrix\n(behavioral / task)",
          x=right_x, y=row_y[1] - 0.05, h=0.55)
    s.box("obj_X", "var", "obj.X", x=right_x, y=row_y[2])
    s.connect_down("design_file", "obj_X")
    # obj.X feeds regress with a diagonal connector from its bottom-left to regress's top-right
    s.connect_line("obj_X", "regress", src_side="bottomleft", dst_side="topright")
    s.annotate("attach as obj.X", x=right_x + 1.85 + 0.05,
               y=row_y[2] + 0.06, w=1.40, h=0.30)

    # ---- Display branch: stats_img → canlab_results_fmridisplay ----
    s.box("crf", "func", "canlab_results_\nfmridisplay",
          x=right_x, y=row_y[4] - 0.05, h=0.55)
    s.arrow_right(x=spine_x + 1.85 + 0.05, y=row_y[4] + 0.07,
                  w=right_x - (spine_x + 1.85) - 0.10)

    s.box("fmridisplay", "var", "fmridisplay\nobject",
          x=right_x, y=row_y[5], h=0.55)
    s.connect_down("crf", "fmridisplay")
    s.annotate("handles", x=right_x + 1.85 + 0.05,
               y=row_y[5] + 0.10, w=0.95, h=0.30)
    s.box("addblobs", "script", "addblobs / removeblobs",
          x=right_x, y=row_y[6], font_size=10.0)
    s.connect_down("fmridisplay", "addblobs")

    # ---- Object types panel (top-right corner) ----
    s.object_types_panel(
        x=right_x + 1.95,
        y=row_y[0],
        lines=["fmri_data", "image_vector", "statistic_image", "region", "fmridisplay"],
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
    out = os.path.join(pptx_dir, "fmri_data_regress_codemap.pptx")
    s = build()
    s.save(out)
    png = render_pptx_to_png(out, png_dir)
    print(f"Wrote {out}\nWrote {png}")
