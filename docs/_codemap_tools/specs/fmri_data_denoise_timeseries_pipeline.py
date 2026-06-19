"""Code map: fmri_data.denoise_timeseries_pipeline — full nuisance-regression pipeline."""

from codemap_lib import Slide


def build() -> Slide:
    s = Slide(
        title="CANlab code map: fmri_data.denoise_timeseries_pipeline",
        description=[
            "Build a comprehensive nuisance regression model and residualize a 4-D fmri_data time series.",
            "Combines 24 motion regressors, run intercepts, spike/outlier indicators, CSF, and HP-filter components.",
        ],
    )

    s.use([
        "Pass an fmri_data, the TR, the high-pass cutoff (sec), and the path to the realignment movement file (rp_*.txt).",
        "Use 'images_per_run', vec to add per-run intercepts (or pre-fill obj.images_per_session).",
        "'pca_denoise', true adds an extra PCA-based denoising step on top of the regression.",
        "'save', true writes the denoised object next to obj.fullpath as <name>_denoised.mat.",
    ])

    s.notes([
        "Nuisance regressors are stored on the output's .covariates and itemized in .metadata_table.",
        "The function delegates the actual regression to canlab_connectivity_preproc; this wrapper assembles the design matrix.",
    ], y=3.40, h=1.50)

    s.legend(x=0.30, y=5.30)

    spine_x = 4.80
    y0 = 1.10
    gap = 0.62
    s.spine([
        ("indat", "var", "fmri_data\n(4-D time series)"),
        ("s1", "func", "1. load mvmt → 24 regressors"),
        ("s2", "func", "2. add run-intercept columns"),
        ("s3", "func", "3. spike / outlier indicators"),
        ("s4", "func", "4. extract CSF signal"),
        ("s5", "func", "5. HP-filter discrete cosine basis"),
        ("s6", "func", "6. (optional) PCA denoise"),
        ("s7", "func", "7. residualize via regression"),
    ], x=spine_x, y0=y0, row_gap=gap, box_w=3.05, box_h=0.50, bold_ids=[])

    # Inputs from the right (kept clear of the Use/Notes column)
    in_x = spine_x + 3.30
    s.box("mvmt", "file", "rp_*.txt\n(realignment params)",
          x=in_x, y=y0 + gap * 0.5, h=0.55, w=2.50, font_size=10.0)
    s.box("tr", "var", "TR (sec)",
          x=in_x, y=y0 + gap * 1.5 + 0.20, h=0.45, w=2.50, font_size=10.0)
    s.box("hp", "var", "HP_cutoff_sec",
          x=in_x, y=y0 + gap * 2.5 + 0.40, h=0.45, w=2.50, font_size=10.0)
    s.connect_line("mvmt", "s1", src_side="left", dst_side="right")
    s.connect_line("tr", "s5", src_side="left", dst_side="right")
    s.connect_line("hp", "s5", src_side="left", dst_side="right")

    # Output
    out_y = y0 + 8 * gap + 0.10
    s.box("out", "var",
          "obj_denoised\n(.dat = residuals, .covariates = nuisance,\n.metadata_table = regressor index)",
          x=spine_x - 0.50, y=out_y, h=0.95, w=4.00, font_size=10.0)
    s.connect_line("s7", "out", src_side="bottom", dst_side="top")

    s.object_types_panel(
        x=in_x, y=y0 + gap * 5,
        lines=["fmri_data"],
        width=2.10,
    )
    return s
