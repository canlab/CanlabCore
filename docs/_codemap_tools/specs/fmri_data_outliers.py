"""Code map: fmri_data.outliers — multi-criterion outlier detection on a 4-D image set."""

from codemap_lib import Slide


def build() -> Slide:
    s = Slide(
        title="CANlab code map: fmri_data.outliers",
        description=[
            "Detect outlier images using DVARS (RMSSD), spatial MAD, Mahalanobis distances on covariance / correlation,",
            "high missing-data fraction, and (optionally) framewise displacement. Returns indicator vectors and tables.",
        ],
    )

    s.use([
        "Run on a 4-D fmri_data — either a time series or a stack of subject contrasts.",
        "For 2nd-level (between-subject) data, pass 'notimeseries' to skip RMSSD/DVARS.",
        "Pass 'fd' with a movement matrix (rotations 1-3, translations 4-6) to add framewise displacement.",
    ])

    s.notes([
        "outlier_indicator_table can be added to a design matrix as nuisance regressors.",
        "Mahalanobis on covariance picks up scale outliers; on correlation it picks up pattern outliers.",
    ], y=3.40, h=1.40)

    s.legend(x=0.30, y=5.20)

    spine_x = 5.30
    y0 = 1.20
    gap = 0.85
    s.spine([
        ("indat", "var", "obj.dat\n(voxels × images)"),
        ("outliers", "func", "outliers( )"),
    ], x=spine_x, y0=y0, row_gap=gap, box_h=0.65, bold_ids=["outliers"])

    # Outputs row directly below outliers (kept inside slide bounds, x ≥ 2.90)
    out_y = y0 + 2 * gap + 0.10
    s.box("uncorr", "var", "est_outliers_uncorr",
          x=2.90, y=out_y, h=0.55, w=2.80, font_size=10.0)
    s.box("corr", "var", "est_outliers_corr",
          x=5.95, y=out_y, h=0.55, w=2.80, font_size=10.0)
    s.box("tables", "var", "outlier_tables\n(+ indicator regressors)",
          x=9.00, y=out_y, h=0.55, w=3.20, font_size=10.0)

    s.connect_line("outliers", "uncorr", src_side="bottom", dst_side="top")
    s.connect_line("outliers", "corr", src_side="bottom", dst_side="top")
    s.connect_line("outliers", "tables", src_side="bottom", dst_side="top")

    # Criterion subfunctions panel — placed BELOW the outputs as a footer so
    # no arrow needs to cross it.
    crits = [
        "RMSSD / DVARS",
        "spatial MAD",
        "Mahalanobis (cov)",
        "Mahalanobis (corr)",
        "missing > 25%",
        "framewise_displacement (optional)",
    ]
    crit_label = (
        "criterion subfunctions called inside outliers( ):\n"
        + "      ".join(f"• {c}" for c in crits[:3]) + "\n"
        + "      ".join(f"• {c}" for c in crits[3:])
    )
    s.box("crit_block", "func", crit_label,
          x=2.90, y=out_y + 0.85, h=1.05, w=9.30, font_size=10.0)

    s.object_types_panel(
        x=2.90, y=out_y + 2.10,
        lines=["fmri_data", "image_vector", "MATLAB table"],
        width=3.40,
    )
    return s
