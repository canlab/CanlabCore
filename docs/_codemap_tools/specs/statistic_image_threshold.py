"""Code map: statistic_image.threshold — apply (or update) a statistical threshold."""

from codemap_lib import Slide


def build() -> Slide:
    s = Slide(
        title="CANlab code map: statistic_image.threshold",
        description=[
            "Apply (or re-apply) a statistical threshold to a statistic_image, updating the .sig vector",
            "without losing the underlying t / p / effect-size values. Thresholding is reversible.",
        ],
    )

    s.use([
        "Pass a statistic_image returned by ttest(), regress(), or fitlme_voxelwise().",
        "Choose threshold type: 'unc', 'fdr', 'fwe' (needs residuals + df), 'cluster_size_fwe', 'bfr', or a raw value range.",
        "Optional 'k', extent for cluster-extent filtering. Re-call to switch thresholds without recomputing the model.",
    ])

    s.notes([
        "threshold() does NOT change .dat — it only edits .sig and (re)derives p-values for FDR/FWE.",
        "Reverse a threshold by calling threshold(obj, Inf, 'unc') or by reloading from .original_image_data.",
    ], y=3.40, h=1.40)

    s.legend(x=0.30, y=5.30)

    spine_x = 5.20
    y0 = 1.50
    gap = 1.05
    s.spine([
        ("in_stats", "var", "statistic_image\n(t, p, effect)"),
        ("threshold", "func", "threshold( )"),
        ("out_stats", "var", "statistic_image\n(.sig updated)"),
    ], x=spine_x, y0=y0, row_gap=gap, box_h=0.65, bold_ids=["threshold"])

    # Single grouped box of threshold types feeding threshold( )
    s.box(
        "types",
        "var",
        "thresh_type:\n'unc'   |   'fdr'   |   'fwe' (+ df, residuals)\n'cluster_size_fwe'   |   'bfr'   |   raw range",
        x=spine_x + 2.50, y=y0 + gap - 0.20, w=3.80, h=0.95, font_size=10.0,
    )
    s.connect_line("types", "threshold", src_side="left", dst_side="right")

    # 'k' (cluster extent) on the left
    s.box("k_opt", "var", "'k', extent\n(cluster filter)",
          x=spine_x - 2.80, y=y0 + gap - 0.05, h=0.65, w=2.20, font_size=10.0)
    s.connect_line("k_opt", "threshold", src_side="right", dst_side="left")

    # Downstream consumers off out_stats
    cons_y = y0 + 2 * gap + 0.30
    s.box("region", "func", "region( )", x=spine_x - 2.80, y=cons_y)
    s.box("table", "func", "table( )", x=spine_x + 2.50, y=cons_y)
    s.box("display", "func", "montage / surface / orthviews",
          x=spine_x - 0.50, y=cons_y + 0.85, w=2.85, h=0.42, font_size=10.0)
    s.connect_line("out_stats", "region", src_side="bottomleft", dst_side="topright")
    s.connect_line("out_stats", "table", src_side="bottomright", dst_side="topleft")
    s.connect_line("out_stats", "display", src_side="bottom", dst_side="top")

    s.object_types_panel(
        x=spine_x + 2.50, y=cons_y + 1.00,
        lines=["statistic_image", "region"],
        width=2.30,
    )
    return s
