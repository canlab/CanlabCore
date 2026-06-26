"""Code map: fmri_data.image_similarity_plot — wedge / polar plot vs basis maps."""

from codemap_lib import Slide


def build() -> Slide:
    s = Slide(
        title="CANlab code map: fmri_data.image_similarity_plot",
        description=[
            "Compute spatial similarity between input images and a set of a priori basis maps,",
            "then render the similarities as a wedge or polar plot. Default metric is point-biserial r.",
        ],
    )

    s.use([
        "Pass an fmri_data with one or more images. Use 'mapset', name to pick a basis set (e.g., 'NPSplus', 'bucknerlab', 'kragelemotion').",
        "Use 'average' to pool across input images with SE bars and statistical tests on each basis.",
        "Pass a grouping vector to test multivariate group differences across the basis maps.",
        "Switch metric with 'cosine_similarity' or 'binary_overlap'.",
    ])

    s.notes([
        "Underlying similarity is computed by canlab_pattern_similarity, which treats zeros as missing in DATA but not in basis maps.",
        "Use this for one-call follow-ups: 'how does my map compare to NPS, SIIPS, the Buckner networks, or Kragel emotion patterns?'",
    ], y=3.40, h=1.60)

    s.legend(x=0.30, y=5.40)

    spine_x = 5.10
    y0 = 1.20
    gap = 0.85
    s.spine([
        ("indat", "var", "fmri_data\n(images to compare)"),
        ("isim", "func", "image_similarity_plot( )"),
        ("simfn", "func", "canlab_pattern_similarity\n(per basis map)"),
        ("plot", "script", "wedge / polar plot"),
    ], x=spine_x, y0=y0, row_gap=gap, box_h=0.65, bold_ids=["isim"])

    # Basis maps load from the right
    s.box("mapset", "file", "basis maps\n(NPSplus, bucknerlab, ...)",
          x=spine_x + 2.60, y=y0 + gap - 0.05, h=0.65, w=2.80, font_size=10.0)
    s.connect_line("mapset", "isim", src_side="left", dst_side="right")

    # Group input (placed on right under mapset)
    s.box("group", "var", "group vector\n(optional)",
          x=spine_x + 2.60, y=y0, h=0.65, w=2.30, font_size=10.0)
    s.connect_line("group", "isim", src_side="left", dst_side="right")

    # Outputs — connected from the BOTTOM of the spine ("plot") to avoid
    # crossing through the wedge/polar-plot box.
    out_y = y0 + 4 * gap + 0.05
    s.box("stats", "var", "stats\n(r, p per basis)",
          x=2.90, y=out_y, h=0.65, w=2.30)
    s.box("hh", "var", "hh, hhfill\n(plot handles)",
          x=5.50, y=out_y, h=0.65, w=2.30)
    s.box("groups", "var", "table_group,\nmultcomp_group",
          x=8.10, y=out_y, h=0.65, w=2.50, font_size=10.0)
    for nid in ("stats", "hh", "groups"):
        s.connect_line("plot", nid, src_side="bottom", dst_side="top")

    s.object_types_panel(
        x=spine_x + 4.30, y=y0,
        lines=["fmri_data", "image_vector", "MATLAB table"],
        width=2.10,
    )
    return s
