"""Graphical roadmap for the ROI / atlas data-extraction workflow.

Follows the CANlab codemap template (docs/canlab_template_codemap.pptx) via
codemap_lib. Builds a single widescreen slide and renders it to PNG with
LibreOffice. Output:
    docs/workflows/ROI_data_extraction_roadmap.pptx
    docs/workflows/ROI_data_extraction_roadmap.png

Run:
    cd docs/_codemap_tools && PYTHONPATH=. python3 build_roi_roadmap.py
"""

from __future__ import annotations

import os

from codemap_lib import Slide, render_pptx_to_png
from pptx.enum.text import PP_ALIGN


def build() -> Slide:
    s = Slide(
        title="CANlab workflow: extracting data from brain regions",
        description=[
            "Choose a method by what you want to summarize: ROI means, multivariate pattern / signature responses,",
            "whole-atlas parcel means, gray / white / CSF compartments, or spheres around coordinates.",
        ],
    )

    # ---- Column headers -----------------------------------------------------
    # (left column header omitted: the legend sits top-left and the green
    #  Data image / Region definition boxes are self-labeling)
    s._text("EXTRACTION METHODS", x=4.60, y=1.30, w=3.0, h=0.25, size=12.0, bold=True)
    s._text("WHAT YOU GET", x=8.30, y=1.30, w=2.9, h=0.25, size=12.0, bold=True)

    # ---- Inputs (left) ------------------------------------------------------
    s.box("data", "var", "Data image\n(fmri_data / image_vector)",
          x=0.35, y=1.62, w=2.55, h=0.66, font_size=10.5)
    s.box("regiondef", "var", "Region definition\nmask · region · atlas · weight map",
          x=0.35, y=2.55, w=2.55, h=0.80, font_size=10.5)
    s._text("Every method takes a data image + a region definition.",
            x=0.35, y=3.45, w=2.6, h=0.5, size=9.5)

    # ways to MAKE a region definition (feed regiondef from below)
    s._text("Make a region definition:", x=0.35, y=3.95, w=2.6, h=0.22, size=10.0, bold=True)
    s.box("sphere", "func", "sphere_mask /\nsphere_roi_tool_2008",
          x=0.35, y=4.22, w=2.55, h=0.52, font_size=9.5)
    s.box("draw", "func", "draw_anatomical_roi",
          x=0.35, y=4.86, w=2.55, h=0.46, font_size=9.5)
    s.box("loadroi", "func", "canlab_load_ROI",
          x=0.35, y=5.44, w=2.55, h=0.46, font_size=9.5)
    for b in ("sphere", "draw", "loadroi"):
        s.connect_line(b, "regiondef", src_side="top", dst_side="bottom", weight=1.0)

    # ---- decorative arrow inputs -> methods ---------------------------------
    s.arrow_right(x=3.05, y=2.12, w=1.30, h=0.42)

    # ---- Methods (center) ---------------------------------------------------
    mx, mw, mh = 4.60, 3.00, 0.66
    ox, ow = 8.30, 2.90
    ys = [1.62, 2.54, 3.46, 4.38, 5.30]
    methods = [
        ("m_roi",   "extract_roi_averages",                 "ROI means\n(region object)"),
        ("m_parc",  "apply_parcellation\n(atlas.extract_data)", "[images × parcels]\nmatrix"),
        ("m_pat",   "apply_mask\n('pattern_expression')",    "Pattern / signature\nresponse"),
        ("m_reg",   "region.extract_data",                  "region.dat\n(no resampling)"),
        ("m_gwc",   "extract_gray_white_csf",               "GM / WM / CSF means\n+ components"),
    ]
    bold_first = True
    for (mid, mlabel, olabel), y in zip(methods, ys):
        s.box(mid, "func", mlabel, x=mx, y=y, w=mw, h=mh, font_size=10.5, bold=bold_first)
        bold_first = False
        oid = mid + "_out"
        s.box(oid, "var", olabel, x=ox, y=y, w=ow, h=mh, font_size=10.0)
        s.connect_line(mid, oid, src_side="right", dst_side="left")

    # imply inputs feed the whole method stack (a couple of representative links)
    s.connect_line("data", "m_roi", src_side="right", dst_side="left", weight=1.0)
    s.connect_line("regiondef", "m_pat", src_side="right", dst_side="left", weight=1.0)

    # ---- object types + legend ---------------------------------------------
    s.object_types_panel(x=11.45, y=1.58,
                         lines=["fmri_data", "image_vector", "atlas", "region"],
                         width=1.65)
    s.legend(x=0.30, y=0.16)   # top-left corner, clear of the centered title

    # ---- footer: spaces & resampling ---------------------------------------
    s._text(
        "Spaces & resampling:  methods assume every image is registered to the same template space (e.g. MNI). "
        "Differing voxel grids are resampled automatically onto the data’s grid (nearest-neighbor for integer "
        "atlases); call resample_space(a, b) to reslice manually. region.extract_data matches mm coordinates with no resampling.",
        x=0.35, y=6.85, w=12.6, h=0.55, size=9.5, align=PP_ALIGN.LEFT,
    )
    return s


def main() -> int:
    here = os.path.dirname(os.path.abspath(__file__))
    out_dir = os.path.normpath(os.path.join(here, "..", "workflows"))
    os.makedirs(out_dir, exist_ok=True)
    pptx_path = os.path.join(out_dir, "ROI_data_extraction_roadmap.pptx")
    slide = build()
    slide.save(pptx_path)
    png_path = render_pptx_to_png(pptx_path, out_dir)
    print(f"[OK] {os.path.basename(pptx_path)} -> {os.path.basename(png_path)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
