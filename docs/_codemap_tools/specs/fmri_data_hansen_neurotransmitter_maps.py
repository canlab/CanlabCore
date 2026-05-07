"""Code map: fmri_data.hansen_neurotransmitter_maps — apply Hansen 2022 PET tracer maps."""

from codemap_lib import Slide


def build() -> Slide:
    s = Slide(
        title="CANlab code map: fmri_data.hansen_neurotransmitter_maps",
        description=[
            "Compare an fmri_data map to the Hansen et al. 2022 (Nat Neurosci) PET neurotransmitter atlas",
            "using cosine similarity on gray-matter-masked images. Returns stats, the maps themselves, and plot handles.",
        ],
    )

    s.use([
        "Pass an fmri_data with one or more images (a single signature, contrast set, or per-subject maps).",
        "Use 'compareGroups' + 'group' to test for between-group differences in tracer associations.",
        "'doAverage' aggregates correlations across input images; 'correlation' switches metric to Pearson r.",
    ])

    s.notes([
        "Internally calls image_similarity_plot with the Hansen mapset and a polar-plot rendering.",
        "Tracer maps are gray-matter-masked first to avoid spurious correlations driven by background voxels.",
    ], y=3.40, h=1.50)

    s.legend(x=0.30, y=5.30)

    spine_x = 5.00
    y0 = 1.20
    gap = 0.85
    s.spine([
        ("indat", "var", "fmri_data\n(map(s) or signature)"),
        ("hansen", "func", "hansen_neurotransmitter_maps( )"),
        ("isim", "func", "image_similarity_plot\n(mapset = 'hansen')"),
    ], x=spine_x, y0=y0, row_gap=gap, box_h=0.65, bold_ids=["hansen"])

    # Hansen tracer maps - prepped .mat distributed with CanlabCore,
    # loaded by name via load_image_set('hansen22'). NOT individual NIfTIs.
    s.box("hansen_files", "file", "hansen22.mat\n(load_image_set)",
          x=spine_x + 2.60, y=y0, h=0.65, w=2.50, font_size=10.0)
    s.connect_line("hansen_files", "isim", src_side="bottomleft", dst_side="topright")

    # Optional group input (placed on right under Hansen input — out of Use column)
    s.box("group", "var", "group\n(optional)",
          x=spine_x + 2.60, y=y0 + gap - 0.05, h=0.65, w=2.30, font_size=10.0)
    s.connect_line("group", "hansen", src_side="left", dst_side="right")

    # Outputs (kept inside slide, clear of Notes column)
    out_y = y0 + 3 * gap + 0.20
    s.box("stats", "var", "stats\n(per-tracer assoc)",
          x=2.90, y=out_y, h=0.65, w=2.30)
    s.box("ntmaps", "var", "ntmaps\n(fmri_data, sorted)",
          x=5.40, y=out_y, h=0.65, w=2.30)
    s.box("plots", "var", "hh, hhfill\n(plot handles)",
          x=7.90, y=out_y, h=0.65, w=2.30)
    s.box("groups", "var", "table_group,\nmultcomp_group",
          x=10.40, y=out_y, h=0.65, w=2.50, font_size=10.0)
    for nid in ("stats", "ntmaps", "plots", "groups"):
        s.connect_line("isim", nid, src_side="bottom", dst_side="top")

    s.object_types_panel(
        x=spine_x + 5.40, y=y0 + 2 * gap + 0.10,
        lines=["fmri_data", "MATLAB table"],
        width=2.30,
    )
    return s
