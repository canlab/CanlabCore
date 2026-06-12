"""Code map: atlas.select_atlas_subset — pick regions out of an atlas by name or index."""

from codemap_lib import Slide


def build() -> Slide:
    s = Slide(
        title="CANlab code map: atlas.select_atlas_subset",
        description=[
            "Carve a subset of regions out of an atlas object — select by label string(s),",
            "by integer index, or both. Optionally flatten the subset into a single binary mask.",
        ],
    )

    s.use([
        "Pass an atlas (from load_atlas) plus a cell of label strings and/or a vector of region IDs.",
        "Use 'flatten' to collapse the subset into one combined mask region.",
        "Use 'deterministic' for probabilistic atlases when winner-take-all boundaries matter more than precision.",
    ])

    s.notes([
        "Probabilistic atlases include voxels with non-zero probability — extracted boundaries can EXCEED the visualized winner-take-all map.",
        "to_extract is a logical vector flagging which original regions ended up in the subset.",
    ], y=3.40, h=1.50)

    s.legend(x=0.30, y=5.40)

    spine_x = 5.40
    y0 = 1.30
    gap = 1.05
    s.spine([
        ("in_atlas", "var", "atlas (full)"),
        ("select", "func", "select_atlas_subset( )"),
        ("out_atlas", "var", "atlas (subset)\n+ to_extract logical"),
    ], x=spine_x, y0=y0, row_gap=gap, box_h=0.65, bold_ids=["select"])

    # Single grouped criteria box on the right
    s.box(
        "criteria",
        "var",
        "selection criteria:\ncell of label strings\ninteger region indices\n'flatten' / 'deterministic'\nalternate field name",
        x=spine_x + 2.50, y=y0 + gap - 0.20,
        w=3.30, h=1.10, font_size=10.0,
    )
    s.connect_line("criteria", "select", src_side="left", dst_side="right")

    # Downstream consumers off out_atlas
    down_y = y0 + 3 * gap + 0.10
    s.box("apply_mask", "func", "apply_mask\n(use as mask)",
          x=spine_x - 2.80, y=down_y, h=0.55, w=1.95)
    s.box("apply_parcel", "func", "apply_parcellation\n(per-parcel means)",
          x=spine_x - 0.40, y=down_y, h=0.55, w=2.50)
    s.box("region_obj", "func", "atlas2region\n(→ region object)",
          x=spine_x + 2.45, y=down_y, h=0.55, w=2.10)
    s.connect_line("out_atlas", "apply_mask", src_side="bottomleft", dst_side="top")
    s.connect_line("out_atlas", "apply_parcel", src_side="bottom", dst_side="top")
    s.connect_line("out_atlas", "region_obj", src_side="bottomright", dst_side="top")

    s.object_types_panel(
        x=spine_x + 2.60, y=y0,
        lines=["atlas", "fmri_data (mask)", "region"],
        width=3.30,
    )
    return s
