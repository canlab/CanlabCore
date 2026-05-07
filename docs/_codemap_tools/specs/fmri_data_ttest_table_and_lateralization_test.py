"""Code map: fmri_data.ttest_table_and_lateralization_test — Yeo-network ttest + lateralization."""

from codemap_lib import Slide


def build() -> Slide:
    s = Slide(
        title="CANlab code map: fmri_data.ttest_table_and_lateralization_test",
        description=[
            "Average per-image data within left and right hemisphere parcels of the Schaefer/Yeo 17-network atlas,",
            "run a one-sample t-test on each network, and test (R-L) lateralization. Returns a table and wedge plot.",
        ],
    )

    s.use([
        "Pass a 2nd-level fmri_data (one image per subject / contrast).",
        "Renders a network-wedge plot and a montage by default; pass 'nomontage' to skip the brain rendering.",
        "Currently hard-coded for the Schaefer/Yeo 17-network parcellation — extending to other parcellations is straightforward.",
    ])

    s.notes([
        "Calls apply_parcellation under the hood with the Yeo atlas split into left/right hemispheres.",
        "Lateralization test is simply ttest(right_mean - left_mean) per network, alongside the absolute network test.",
    ], y=3.40, h=1.50)

    s.legend(x=0.30, y=5.30)

    spine_x = 4.70
    y0 = 1.30
    gap = 0.90
    s.spine([
        ("indat", "var", "fmri_data\n(2nd-level)"),
        ("tlat", "func", "ttest_table_and_\nlateralization_test( )"),
        ("ap", "func", "apply_parcellation\n(Yeo 17 × L/R)"),
        ("tt", "func", "ttest per network\n+ ttest(R - L)"),
    ], x=spine_x, y0=y0, row_gap=gap, box_w=2.85, box_h=0.70, bold_ids=["tlat"])

    # Atlas reference (on right)
    s.box("yeo", "file", "Schaefer / Yeo 17-net\nL/R atlas",
          x=spine_x + 3.40, y=y0 + 2 * gap - 0.05, h=0.70, w=2.50, font_size=10.0)
    s.connect_line("yeo", "ap", src_side="left", dst_side="right")

    # Outputs across the bottom — clear of Use/Notes column (x ≥ 2.90)
    out_y = y0 + 4 * gap + 0.20
    s.box("roi", "var", "roi_table\n(MATLAB table: means + ttest stats)",
          x=2.90, y=out_y, h=0.70, w=3.40, font_size=10.0)
    s.box("subj", "var", "subj_dat\n[images × networks]",
          x=6.50, y=out_y, h=0.70, w=2.80, font_size=10.0)
    s.box("plots", "script", "wedge plot + montage",
          x=9.50, y=out_y, h=0.70, w=2.80, font_size=10.0)
    for nid in ("roi", "subj", "plots"):
        s.connect_line("tt", nid, src_side="bottom", dst_side="top")

    s.object_types_panel(
        x=spine_x + 3.40, y=y0,
        lines=["fmri_data", "atlas", "MATLAB table"],
        width=2.50,
    )
    return s
