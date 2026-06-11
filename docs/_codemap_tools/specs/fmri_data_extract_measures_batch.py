"""Code map: fmri_data.extract_measures_batch — aggregate canonical measures from a dataset."""

from codemap_lib import Slide


def build() -> Slide:
    s = Slide(
        title="CANlab code map: fmri_data.extract_measures_batch",
        description=[
            "Run a battery of pattern- and network-based extractors on an fmri_data and bundle the results",
            "into a single DAT struct — Mahalanobis QC, RMSSD, GM/WM/CSF tables, signature responses, parcellations.",
        ],
    )

    s.use([
        "Pass an fmri_data — typically a single subject's preprocessed time series, or a stack of contrast images.",
        "Output DAT has one field per measure; pull whatever subset is needed for downstream analysis.",
        "Use as a standardized 'first-pass extractor' so multiple studies have comparable measure inventories.",
    ])

    s.notes([
        "DAT.mahalanobis.wh_outlier_uncorr is a useful logical vector of outlier images.",
        "Signature fields (npsplus, kragelemotion, kragel18, pain_pdm) require their masks to be on the path (Neuroimaging_Pattern_Masks).",
    ], y=3.40, h=1.50)

    s.legend(x=0.30, y=5.30)

    spine_x = 4.90
    y0 = 1.20
    gap = 0.85
    s.spine([
        ("indat", "var", "fmri_data\n(subject TS or contrasts)"),
        ("batch", "func", "extract_measures_batch( )"),
    ], x=spine_x, y0=y0, row_gap=gap, box_h=0.65, bold_ids=["batch"])

    # Sub-extractors arranged as a single grouped panel — avoids the visual
    # mess of fan-in arrows from a 2-row grid into a single sink.
    sub_y0 = y0 + 2 * gap - 0.10
    subs_label = (
        "internal extractors:\n"
        "• outliers( )                              → DAT.mahalanobis\n"
        "• rmssd_movie                       → DAT.rmssd\n"
        "• extract_gray_white_csf      → DAT.gray_white_csf_table\n"
        "• apply_all_signatures (NPS+) → DAT.npsplus / kragel*\n"
        "• apply_parcellation              → DAT.PARCELS\n"
        "• image_similarity_plot          → DAT.pain_pdm"
    )
    s.box("subs", "func", subs_label,
          x=4.20, y=sub_y0, w=5.30, h=2.00, font_size=10.0)
    s.connect_line("batch", "subs", src_side="bottom", dst_side="top")

    # DAT output (single arrow from the extractor block down)
    out_y = sub_y0 + 2.20
    s.box("DAT", "var", "DAT struct\n(one field per extractor + provenance)",
          x=spine_x - 0.40, y=out_y, h=0.75, w=3.80, font_size=10.0)
    s.connect_line("subs", "DAT", src_side="bottom", dst_side="top")

    s.object_types_panel(
        x=spine_x + 4.20, y=y0,
        lines=["fmri_data", "MATLAB struct", "MATLAB table"],
        width=2.40,
    )
    return s
