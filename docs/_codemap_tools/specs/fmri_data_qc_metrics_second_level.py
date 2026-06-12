"""Code map: fmri_data.qc_metrics_second_level — QC metrics for a 2nd-level dataset."""

from codemap_lib import Slide


def build() -> Slide:
    s = Slide(
        title="CANlab code map: fmri_data.qc_metrics_second_level",
        description=[
            "Compute scale-invariant quality metrics for a 2nd-level fmri_data (e.g., subject-level beta / contrast images).",
            "Returns group-level summaries, per-image metrics, and tissue-compartment signals.",
        ],
    )

    s.use([
        "Pass an fmri_data with one image per subject (or contrast).",
        "Designed to be comparable across studies — metrics are scale-invariant.",
        "Pass 'noverbose' to suppress the printed summary.",
    ])

    s.notes([
        "Internally extracts gray/white/CSF signals via extract_gray_white_csf and computes Mahalanobis-style outlier scores.",
        "Use group_metrics for a quick at-a-glance check; individual_metrics for per-subject troubleshooting.",
    ], y=3.40, h=1.50)

    s.legend(x=0.30, y=5.30)

    spine_x = 5.00
    y0 = 1.20
    gap = 0.85
    s.spine([
        ("indat", "var", "fmri_data\n(2nd-level: subj × img)"),
        ("qc", "func", "qc_metrics_second_level( )"),
        ("inner", "func", "extract_gray_white_csf\n+ similarity / Mahalanobis"),
    ], x=spine_x, y0=y0, row_gap=gap, box_h=0.65, bold_ids=["qc"])

    # Outputs across the bottom — kept inside slide bounds and clear of the
    # Notes column (x ≥ 2.90, x_max ≤ 13.00).
    out_y = y0 + 3 * gap + 0.20
    outs = [
        ("group", "group_metrics\n(study-level)"),
        ("indiv", "individual_metrics\n(per image)"),
        ("vals", "values\n(raw measures)"),
        ("gwcsf", "gwcsf, gwcsfmean\n(tissue signals)"),
        ("l2", "gwcsf_l2norm\n(length-adjusted)"),
    ]
    cx0 = 2.90
    col_w = 1.95
    col_gap = 0.10
    for i, (nid, label) in enumerate(outs):
        s.box(nid, "var", label,
              x=cx0 + i * (col_w + col_gap), y=out_y, h=0.65, w=col_w, font_size=10.0)
        s.connect_line("inner", nid, src_side="bottom", dst_side="top")

    s.object_types_panel(
        x=spine_x + 2.40, y=y0,
        lines=["fmri_data", "MATLAB struct"],
        width=2.10,
    )
    return s
