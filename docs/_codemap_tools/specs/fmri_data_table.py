"""Code map: fmri_data.table — tabular results for a thresholded image."""

from codemap_lib import Slide


def build() -> Slide:
    s = Slide(
        title="CANlab code map: fmri_data.table",
        description=[
            "Print and return a results table for a thresholded fmri_data / statistic_image — one row per",
            "contiguous cluster, or one row per atlas region when an atlas is supplied.",
        ],
    )

    s.use([
        "Pass a thresholded image (apply threshold() first if needed); without an atlas, table() splits on contiguous clusters.",
        "Add an atlas-class object (e.g., from load_atlas()) as a second argument to subdivide by anatomy.",
        "Use 'k', N to drop clusters smaller than N voxels; 'nosep' to keep positive/negative peaks together; 'nosort' to keep input order.",
        "Returns the labeled region object plus a MATLAB table — pass the region object on to montage() or surface().",
    ])

    s.notes([
        "Without 'nosep', clusters are split into +/- subregions, so the row count may exceed the input region count.",
        "table() on the @region class is the underlying engine; the fmri_data / image_vector wrapper handles atlas dispatch.",
    ], y=3.40, h=1.50)

    s.legend(x=0.30, y=5.30)

    spine_x = 5.10
    y0 = 1.30
    gap = 0.85
    s.spine([
        ("indat", "var", "thresholded\nfmri_data / statistic_image"),
        ("table", "func", "table( )"),
        ("region", "func", "region( ) →\nsubdivide_by_atlas"),
        ("regobj", "var", "region_obj\n(labeled)"),
        ("results", "var", "results_table\n(MATLAB table)"),
    ], x=spine_x, y0=y0, row_gap=gap, box_h=0.65, bold_ids=["table"])

    # Optional atlas input on the right
    s.box("atlas", "var", "atlas object\n(optional)",
          x=spine_x + 2.60, y=y0 + gap - 0.05, h=0.65, w=2.50)
    s.connect_line("atlas", "table", src_side="left", dst_side="right")

    # Options block — placed on right (under atlas) so it doesn't overlap Use column.
    s.box("opts", "var", "options:\n'k', N | 'nosep'\n'nosort' | 'names'",
          x=spine_x + 2.60, y=y0 + 2 * gap - 0.10, h=0.85, w=2.50, font_size=10.0)
    s.connect_line("opts", "table", src_side="left", dst_side="right")

    s.object_types_panel(
        x=spine_x + 2.60, y=y0 + 4 * gap - 0.05,
        lines=["region", "atlas", "MATLAB table"],
        width=2.40,
    )
    return s
