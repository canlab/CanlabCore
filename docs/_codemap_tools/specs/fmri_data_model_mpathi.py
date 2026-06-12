"""Code map: fmri_data.model_mpathi — single-pathway multivariate connectivity (MPathI)."""

from codemap_lib import Slide


def build() -> Slide:
    s = Slide(
        title="CANlab code map: fmri_data.model_mpathi",
        description=[
            "Cross-validated PLS regression between one source ROI and one target ROI — the streamlined",
            "single-pathway sibling of model_brain_pathway, focused on the MPathI framework (Kragel et al. 2021 Neuron).",
        ],
    )

    s.use([
        "Pass an fmri_data plus binary source and target masks (each as fmri_data).",
        "'Indices', vec defines CV folds and bootstrap blocks (default: 10-fold).",
        "'Align' fits a hyperalignment model on training data and warps test data into shared space.",
        "'nboot', N runs block bootstrap on V and Z spatial weights.",
    ])

    s.notes([
        "Estimates one directed pathway only (X → Y); for matched/mismatched pathway comparisons use model_brain_pathway.",
        "Conventions: X, Y are [images × voxels]; T, U latent timeseries are [images × 1]; Z, V spatial weights are [voxels × 1].",
    ], y=3.40, h=1.60)

    s.legend(x=0.30, y=5.40)

    spine_x = 5.00
    y0 = 1.20
    gap = 0.85
    s.spine([
        ("indat", "var", "fmri_data\n(images × voxels)"),
        ("mp", "func", "model_mpathi( )"),
        ("cv", "func", "10-fold CV\n(grouped by Indices)"),
        ("pls", "func", "PLS per fold\n(latent TS + weights)"),
        ("boot", "func", "block bootstrap\n(V, Z, p-values)"),
    ], x=spine_x, y0=y0, row_gap=gap, box_h=0.62, bold_ids=["mp"])

    # Mask inputs
    rx = spine_x + 2.60
    s.box("src", "var", "source_mask\n(fmri_data binary)",
          x=rx, y=y0, h=0.65, w=2.50, font_size=10.0)
    s.box("tgt", "var", "target_mask\n(fmri_data binary)",
          x=rx, y=y0 + 0.85, h=0.65, w=2.50, font_size=10.0)
    s.connect_line("src", "mp", src_side="left", dst_side="right")
    s.connect_line("tgt", "mp", src_side="left", dst_side="right")

    # Output
    out_y = y0 + 5 * gap + 0.10
    s.box("stats", "var",
          "stats struct\n(T, U latent TS; Z, V spatial weights;\npathway correlation + bootstrap p-values)",
          x=spine_x - 1.30, y=out_y, h=0.95, w=4.60, font_size=10.0)
    s.connect_line("boot", "stats", src_side="bottom", dst_side="top")

    s.object_types_panel(
        x=rx, y=y0 + 1.80,
        lines=["fmri_data"],
        width=2.30,
    )
    return s
