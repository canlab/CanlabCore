"""Code map: fmri_data.wedge_plot_by_atlas — wedge plots colored by atlas regions."""

from codemap_lib import Slide


def build() -> Slide:
    s = Slide(
        title="CANlab code map: fmri_data.wedge_plot_by_atlas",
        description=[
            "Render one wedge plot per atlas, with one wedge per region. In data mode, wedge size = |parcel mean|;",
            "in 'signature' mode, wedge size is normalized local pattern expression and color encodes valence.",
        ],
    )

    s.use([
        "Pass an image_vector / fmri_data / statistic_image plus 'atlases', {names} to draw one wedge plot per named atlas.",
        "Use 'signature' mode for multivariate patterns (sizes scale with pattern energy, color is uniform-pos / uniform-neg / mixed).",
        "Add 'montage' to render a parcel-colored montage and 'surfaces' for cortical-surface rendering with a legend.",
    ])

    s.notes([
        "Atlases are loaded by name through load_atlas; pass a cell of names like {'canlab2024', 'buckner_networks'}.",
        "output_values_by_region is a struct of per-atlas, per-region values — useful for downstream tabular reports.",
    ], y=3.40, h=1.50)

    s.legend(x=0.30, y=5.30)

    spine_x = 4.85
    y0 = 1.30
    gap = 0.90
    s.spine([
        ("indat", "var", "image_vector / fmri_data\n(map or signature)"),
        ("wpa", "func", "wedge_plot_by_atlas( )"),
        ("ap", "func", "apply_parcellation\n(per atlas)"),
        ("plot", "script", "wedge plot per atlas\n(+ montage / surfaces)"),
    ], x=spine_x, y0=y0, row_gap=gap, box_h=0.65, bold_ids=["wpa"])

    # Atlas list on the right (placed below spine row 2, between input rows)
    s.box("atlases", "file", "atlases\n(load_atlas names)",
          x=spine_x + 3.10, y=y0 + 2 * gap - 0.05, h=0.65, w=2.80, font_size=10.0)
    s.connect_line("atlases", "ap", src_side="left", dst_side="right")

    # Outputs across the bottom — kept inside slide bounds and clear of Use/Notes column
    out_y = y0 + 4 * gap + 0.20
    outs = [
        ("hh", "hh\n(figure handles)"),
        ("vals", "output_values_\nby_region"),
        ("labels", "labels"),
        ("atlas_obj", "atlas_obj\n(per-atlas)"),
        ("table", "atlastab\n(MATLAB table)"),
    ]
    cx0 = 2.90
    col_w = 1.95
    col_gap = 0.10
    for i, (nid, label) in enumerate(outs):
        s.box(nid, "var", label,
              x=cx0 + i * (col_w + col_gap), y=out_y, h=0.85, w=col_w, font_size=9.5)
        s.connect_line("plot", nid, src_side="bottom", dst_side="top")

    s.object_types_panel(
        x=spine_x + 3.10, y=y0,
        lines=["fmri_data", "image_vector", "atlas", "MATLAB table"],
        width=2.60,
    )
    return s
