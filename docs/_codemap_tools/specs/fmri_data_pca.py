"""Code map: fmri_data.pca — principal components on the images × voxels matrix."""

from codemap_lib import Slide


def build() -> Slide:
    s = Slide(
        title="CANlab code map: fmri_data.pca",
        description=[
            "Principal-component analysis over an image_vector / fmri_data object.",
            "PCA runs on the images × voxels matrix (covariance, voxel-mean-centered, not z-scored).",
        ],
    )

    s.use([
        "Pass an fmri_data with multiple images (e.g., a time series or a stack of contrasts).",
        "Returns component scores (one per image), eigenmaps (as an fmri_data), and percent variance explained.",
        "Use 'k', N to control the number of components; 'noplot' to suppress the diagnostics figure.",
    ])

    s.notes([
        "Voxels are mean-centered but not standardized — high-variance voxels dominate. For a correlation-style PCA, z-score voxels first.",
        "Mean-centered data can be reconstructed as X = scores * eigenmap_obj.dat'.",
    ], y=3.40, h=1.50)

    s.legend(x=0.30, y=5.20)

    spine_x = 5.10
    y0 = 1.30
    gap = 0.90
    s.spine([
        ("indat", "var", "obj.dat\n(voxels × images)"),
        ("pca", "func", "pca( )"),
        ("svd", "func", "MATLAB pca / svd\n(covariance)"),
    ], x=spine_x, y0=y0, row_gap=gap, box_h=0.65, bold_ids=["pca"])

    # Three outputs branch off
    out_y = y0 + 3 * gap + 0.10
    s.box("scores", "var", "component_scores\n[images × k]",
          x=spine_x - 3.00, y=out_y, h=0.65)
    s.box("eigen", "var", "eigenmap_obj\n(fmri_data, [voxels × k])",
          x=spine_x - 0.40, y=out_y, h=0.65, w=2.40)
    s.box("explained", "var", "explained\n(% variance)",
          x=spine_x + 2.50, y=out_y, h=0.65)
    s.connect_line("svd", "scores", src_side="bottomleft", dst_side="top")
    s.connect_line("svd", "eigen", src_side="bottom", dst_side="top")
    s.connect_line("svd", "explained", src_side="bottomright", dst_side="top")

    # Diagnostic plot off pca itself
    s.box("plot", "script", "diagnostics figure\n(scree + first eigenmaps)",
          x=spine_x + 2.50, y=y0 + gap - 0.05, h=0.65, w=2.50, font_size=10.0)
    s.connect_line("pca", "plot", src_side="right", dst_side="left")

    s.object_types_panel(
        x=spine_x + 5.10, y=y0,
        lines=["fmri_data", "image_vector"],
        width=2.00,
    )
    return s
