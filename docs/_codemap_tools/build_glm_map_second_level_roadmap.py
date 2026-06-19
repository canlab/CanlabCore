"""Graphical roadmap for the second-level (group) glm_map workflow.

Follows the CANlab codemap template via codemap_lib. Builds one widescreen slide
and renders it to PNG with LibreOffice.
    docs/workflows/glm_map_second_level_roadmap.pptx / .png

Run:
    cd docs/_codemap_tools && PYTHONPATH=. python3 build_glm_map_second_level_roadmap.py
"""

from __future__ import annotations

import os

from codemap_lib import Slide, render_pptx_to_png
from pptx.enum.text import PP_ALIGN


def build() -> Slide:
    s = Slide(
        title="CANlab workflow: second-level (group) GLM with glm_map",
        description=[
            "One image per subject (a contrast map) regressed on a static design matrix X at every voxel.",
            "Build the object with fmri_data.regress or the glm_map estimator API; screen, fit (OLS/robust), threshold, visualize.",
        ],
    )

    # ---- Column headers -----------------------------------------------------
    # (left "prepare the data" header omitted: the legend sits top-left, and
    #  the green/blue boxes are self-labeling)
    s._text("BUILD THE OBJECT", x=4.55, y=1.30, w=3.1, h=0.25, size=12.0, bold=True)
    s._text("RESULTS", x=8.45, y=1.30, w=2.9, h=0.25, size=12.0, bold=True)

    # ---- Left: data + preprocessing ----------------------------------------
    s.box("imgs", "var", "Contrast images\n(fmri_data, [voxels x subjects])",
          x=0.35, y=1.60, w=2.7, h=0.62, font_size=10.0)
    s.box("outliers", "func", "outliers('notimeseries')\n-> exclude with get_wh_image",
          x=0.35, y=2.45, w=2.7, h=0.58, font_size=9.5)
    s.box("gwcsf", "func", "extract_gray_white_csf\n-> WM / CSF nuisance cols",
          x=0.35, y=3.20, w=2.7, h=0.52, font_size=9.5)
    s.box("norm", "func", "normalize_gm_by_wm_csf\n(rescale gray matter)",
          x=0.35, y=3.86, w=2.7, h=0.52, font_size=9.5)
    s.box("setX", "var", "dat.X = predictors +\nnuisance + intercept",
          x=0.35, y=4.70, w=2.7, h=0.58, font_size=10.0)
    for b in ("outliers", "gwcsf", "norm"):
        s.connect_line(b, "setX", src_side="bottom", dst_side="top", weight=0.8)
    s.connect_line("imgs", "outliers", src_side="bottom", dst_side="top", weight=0.8)

    s.arrow_right(x=3.20, y=2.55, w=1.10, h=0.40)

    # ---- Center: two build paths -> the object -----------------------------
    cx, cw = 4.55, 3.10
    s.box("regress", "func", "regress(dat, p, type[, 'robust'])",
          x=cx, y=1.60, w=cw, h=0.55, font_size=10.5, bold=True)
    s._text("the quick path", x=cx, y=2.16, w=cw, h=0.2, size=9.0)
    s._text("— or, the estimator API —", x=cx, y=2.40, w=cw, h=0.22, size=9.5, bold=True)
    s.box("api", "func", "glm_map('X', X, 'level', 2)",
          x=cx, y=2.66, w=cw, h=0.48, font_size=10.0)
    s.box("addc", "func", "add_contrasts / create_orthogonal_contrast_set",
          x=cx, y=3.22, w=cw, h=0.48, font_size=9.5)
    s.box("diag", "func", "run_diagnostics  (VIF/cVIF, efficiency)",
          x=cx, y=3.78, w=cw, h=0.48, font_size=9.5)
    s.box("fit", "func", "fit(g, dat[, 'robust'])",
          x=cx, y=4.34, w=cw, h=0.48, font_size=10.0)
    s.connect_line("api", "addc", src_side="bottom", dst_side="top", weight=0.8)
    s.connect_line("addc", "diag", src_side="bottom", dst_side="top", weight=0.8)
    s.connect_line("diag", "fit", src_side="bottom", dst_side="top", weight=0.8)

    s.box("g", "var", "glm_map\nbetas · t · contrasts · diagnostics",
          x=cx, y=5.20, w=cw, h=0.62, font_size=10.5, bold=True)
    s.connect_line("regress", "g", src_side="right", dst_side="top", weight=1.0, kind="bent")
    s.connect_line("fit", "g", src_side="bottom", dst_side="top", weight=1.0)
    s.connect_line("setX", "api", src_side="right", dst_side="left", weight=0.8)

    s.arrow_right(x=7.75, y=5.35, w=0.60, h=0.34)

    # ---- Right: results -----------------------------------------------------
    rx, rw = 8.45, 2.95
    s.box("thr", "func", "threshold(g, p, type, 'k')", x=rx, y=1.60, w=rw, h=0.48, font_size=10.0)
    s.box("tbl", "func", "table(g)  (atlas-labeled)", x=rx, y=2.16, w=rw, h=0.46, font_size=10.0)
    s.box("mtg", "func", "montage(g) / canlab_results_fmridisplay", x=rx, y=2.70, w=rw, h=0.56, font_size=9.5)
    s.box("summ", "func", "summary(g)", x=rx, y=3.34, w=rw, h=0.44, font_size=10.0)
    s.box("tmap", "var", "thresholded t / contrast maps\n(statistic_image) -> region(t)",
          x=rx, y=4.00, w=rw, h=0.58, font_size=9.5)
    s.connect_line("g", "thr", src_side="right", dst_side="left", weight=1.0, kind="bent")
    for b in ("thr", "tbl", "mtg", "summ"):
        s.connect_line(b, "tmap", src_side="bottom", dst_side="top", weight=0.6)

    # ---- object types + legend ---------------------------------------------
    s.object_types_panel(x=11.55, y=1.58,
                         lines=["fmri_data", "glm_map", "statistic_image", "region"],
                         width=1.55)
    s.legend(x=0.30, y=0.16)

    # ---- footer -------------------------------------------------------------
    s._text(
        "Second-level designs are STATIC: X is prespecified (no onsets to convolve). Center continuous predictors "
        "so the intercept maps the group average. Mark covariates of no interest with g.nuisance_columns so diagnostics "
        "report VIFs with and without them. Robust regression down-weights influential subjects voxelwise (vs. hard exclusion).",
        x=0.35, y=6.40, w=12.6, h=0.7, size=9.5, align=PP_ALIGN.LEFT,
    )
    return s


def main() -> int:
    here = os.path.dirname(os.path.abspath(__file__))
    out_dir = os.path.normpath(os.path.join(here, "..", "workflows"))
    os.makedirs(out_dir, exist_ok=True)
    pptx_path = os.path.join(out_dir, "glm_map_second_level_roadmap.pptx")
    slide = build()
    slide.save(pptx_path)
    png_path = render_pptx_to_png(pptx_path, out_dir)
    print(f"[OK] {os.path.basename(pptx_path)} -> {os.path.basename(png_path)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
