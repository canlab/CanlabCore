"""Code map: fmri_data.model_brain_pathway — cross-validated PLS connectivity (4 pathways)."""

from codemap_lib import Slide


def build() -> Slide:
    s = Slide(
        title="CANlab code map: fmri_data.model_brain_pathway",
        description=[
            "Model multivariate connectivity between two source ROIs and two target ROIs using cross-validated",
            "PLS. Compares the 4 source × target pathways (matched + 2 mismatched) — Kragel et al. (2021) Neuron.",
        ],
    )

    s.use([
        "Pass an fmri_data plus 4 ROI masks: source_one, source_two, target_one, target_two.",
        "Use 'Indices' to define CV folds (defaults to 10-fold). 'Align' adds hyperalignment to training data.",
        "'nboot', N runs block bootstrap on V/Z weights for voxel-level inference.",
        "Compares matched pathways (X1→Y1, X2→Y2) against mismatched (X1→Y2, X2→Y1) for specificity.",
    ])

    s.notes([
        "PLS yields T (X scores) / U (Y scores) latent timeseries plus V (X-voxel) / Z (Y-voxel) weights.",
        "Higher correlations expected for matched pathways than mismatched if connectivity is region-specific.",
    ], y=3.40, h=1.70)

    s.legend(x=0.30, y=5.50)

    spine_x = 4.80
    y0 = 1.10
    gap = 0.75
    s.spine([
        ("indat", "var", "fmri_data\n(images × voxels)"),
        ("mbp", "func", "model_brain_pathway( )"),
        ("cv", "func", "cross-validation folds"),
        ("simple", "func", "ROI averages →\nsimple correlations"),
        ("pls", "func", "PLS per pathway\n(4 src × tgt combinations)"),
        ("boot", "func", "block bootstrap\n(V, Z weights)"),
    ], x=spine_x, y0=y0, row_gap=gap, box_h=0.62, bold_ids=["mbp"])

    # Source / target mask inputs from the right
    rx = spine_x + 2.50
    s.box("s1", "var", "source_one mask",
          x=rx, y=y0, h=0.55, w=2.40, font_size=10.0)
    s.box("s2", "var", "source_two mask",
          x=rx, y=y0 + 0.62, h=0.55, w=2.40, font_size=10.0)
    s.box("t1", "var", "target_one mask",
          x=rx, y=y0 + 1.24, h=0.55, w=2.40, font_size=10.0)
    s.box("t2", "var", "target_two mask",
          x=rx, y=y0 + 1.86, h=0.55, w=2.40, font_size=10.0)
    for nid in ("s1", "s2", "t1", "t2"):
        s.connect_line(nid, "mbp", src_side="left", dst_side="right")

    # Output
    out_y = y0 + 6 * gap + 0.10
    s.box("stats", "var",
          "stats struct\n(T, U latent TS;  V, Z spatial weights;\npathway correlations + boot p-values)",
          x=spine_x - 1.30, y=out_y, h=0.95, w=4.60, font_size=10.0)
    s.connect_line("boot", "stats", src_side="bottom", dst_side="top")

    s.object_types_panel(
        x=rx, y=y0 + 2.55,
        lines=["fmri_data", "image_vector\n(masks)"],
        width=2.40,
    )
    return s
