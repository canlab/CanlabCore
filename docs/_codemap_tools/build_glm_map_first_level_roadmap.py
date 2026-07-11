"""Graphical roadmap for the first-level (time-series) glm_map workflow.

    docs/workflows/glm_map_first_level_roadmap.pptx / .png

Run:
    cd docs/_codemap_tools && PYTHONPATH=. python3 build_glm_map_first_level_roadmap.py
"""

from __future__ import annotations

import os

from codemap_lib import Slide, render_pptx_to_png
from pptx.enum.text import PP_ALIGN


def build() -> Slide:
    s = Slide(
        title="CANlab workflow: first-level (single-subject) fMRI time-series GLM with glm_map",
        description=[
            "Events -> HRF convolution -> design matrix X, stacked across runs with nuisance regressors.",
            "Import onsets, choose a basis set, screen, fit a BOLD time series (AR errors), threshold, visualize.",
        ],
    )

    s._text("BUILD & SCREEN THE MODEL", x=4.55, y=1.30, w=3.3, h=0.25, size=12.0, bold=True)
    s._text("RESULTS", x=8.55, y=1.30, w=2.9, h=0.25, size=12.0, bold=True)

    # ---- Left: import routes -----------------------------------------------
    s._text("Get event info in:", x=0.35, y=1.46, w=2.7, h=0.22, size=11.0, bold=True)
    s.box("fsl", "func", "import_onsets\n(FSL event table .csv/.xlsx)",
          x=0.35, y=1.72, w=2.75, h=0.56, font_size=9.5)
    s.box("cells", "func", "import_onsets\n(SPM-style cell arrays)",
          x=0.35, y=2.40, w=2.75, h=0.52, font_size=9.5)
    s.box("spmmat", "func", "import_SPM\n(SPM.mat first-level)",
          x=0.35, y=3.04, w=2.75, h=0.52, font_size=9.5)
    s.box("wrap", "func", "glm_map(fmri_glm_design_matrix)",
          x=0.35, y=3.68, w=2.75, h=0.46, font_size=9.5)
    s.box("sim", "var", "simulate data (X·betas+noise)\nor load run time series (fmri_data)",
          x=0.35, y=4.95, w=2.75, h=0.60, font_size=9.0)
    s._text("(feeds fit, below)", x=0.35, y=5.58, w=2.7, h=0.2, size=8.5)

    s.arrow_right(x=3.25, y=2.80, w=1.05, h=0.40)

    # ---- Center: build/screen vertical flow --------------------------------
    cx, cw = 4.55, 3.30
    s.box("build", "func", "build_design  (onsets ⊛ basis set → X)",
          x=cx, y=1.62, w=cw, h=0.52, font_size=10.5, bold=True)
    s.box("basis", "func", "replace_basis_set  (HRF → spline / + derivatives)",
          x=cx, y=2.20, w=cw, h=0.46, font_size=9.5)
    s.box("contr", "func", "add_contrasts / create_orthogonal_contrast_set",
          x=cx, y=2.74, w=cw, h=0.46, font_size=9.5)
    s.box("diag", "func", "run_diagnostics\nVIF/cVIF (± nuisance) · efficiency · HP-filter cutoff",
          x=cx, y=3.28, w=cw, h=0.60, font_size=9.5)
    s.box("fit", "func", "fit(g, timeseries, 'AR', 4 / 'robust')",
          x=cx, y=4.00, w=cw, h=0.48, font_size=10.0)
    for a, b in (("build", "basis"), ("basis", "contr"), ("contr", "diag"), ("diag", "fit")):
        s.connect_line(a, b, src_side="bottom", dst_side="top", weight=0.8)
    s.connect_line("sim", "fit", src_side="right", dst_side="left", weight=0.8, kind="bent")

    s.box("g", "var", "glm_map\nbetas · t · contrasts · diagnostics",
          x=cx, y=4.86, w=cw, h=0.60, font_size=10.5, bold=True)
    s.connect_line("fit", "g", src_side="bottom", dst_side="top", weight=1.0)
    for b in ("fsl", "cells", "spmmat", "wrap"):
        s.connect_line(b, "build", src_side="right", dst_side="left", weight=0.5)

    s.arrow_right(x=7.95, y=5.00, w=0.55, h=0.34)

    # ---- Right: results -----------------------------------------------------
    rx, rw = 8.55, 2.85
    s.box("thr", "func", "threshold(g, p, type, 'k')", x=rx, y=1.62, w=rw, h=0.48, font_size=10.0)
    s.box("mtg", "func", "montage(g) / canlab_results_fmridisplay", x=rx, y=2.18, w=rw, h=0.56, font_size=9.5)
    s.box("tbl", "func", "table(g)  (atlas-labeled)", x=rx, y=2.82, w=rw, h=0.44, font_size=10.0)
    s.box("summ", "func", "summary(g)", x=rx, y=3.34, w=rw, h=0.42, font_size=10.0)
    s.box("tmap", "var", "thresholded t maps per event\n(statistic_image) → region(t)",
          x=rx, y=3.96, w=rw, h=0.56, font_size=9.5)
    s.connect_line("g", "thr", src_side="right", dst_side="left", weight=1.0, kind="bent")
    for b in ("thr", "mtg", "tbl", "summ"):
        s.connect_line(b, "tmap", src_side="bottom", dst_side="top", weight=0.6)

    # ---- object types + legend ---------------------------------------------
    s.object_types_panel(x=11.55, y=1.58,
                         lines=["glm_map", "fmri_glm_design_matrix", "fmri_data", "statistic_image"],
                         width=1.7)
    s.legend(x=0.30, y=0.16)

    # ---- footer -------------------------------------------------------------
    s._text(
        "First-level designs are BUILT from event timing: onsets/durations are convolved with a basis set (canonical HRF by "
        "default) and stacked across runs with per-run intercepts; entered events are flagged of interest, motion / "
        "multiple_regressors as nuisance. Set g.is_timeseries = true to enable AR(p) error models (AR(4) recommended; needs "
        "Econometrics + Signal Processing toolboxes) and HP-filter diagnostics; run_diagnostics recommends a high-pass-filter "
        "cutoff (< 5% variance lost) for the regressors and contrasts.",
        x=0.35, y=6.45, w=12.6, h=0.7, size=9.5, align=PP_ALIGN.LEFT,
    )
    return s


def main() -> int:
    here = os.path.dirname(os.path.abspath(__file__))
    out_dir = os.path.normpath(os.path.join(here, "..", "workflows"))
    os.makedirs(out_dir, exist_ok=True)
    pptx_path = os.path.join(out_dir, "glm_map_first_level_roadmap.pptx")
    slide = build()
    slide.save(pptx_path)
    png_path = render_pptx_to_png(pptx_path, out_dir)
    print(f"[OK] {os.path.basename(pptx_path)} -> {os.path.basename(png_path)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
