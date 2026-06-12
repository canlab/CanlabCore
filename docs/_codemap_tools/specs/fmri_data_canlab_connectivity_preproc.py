"""Code map: fmri_data.canlab_connectivity_preproc — connectivity-prep pipeline."""

from codemap_lib import Slide


def build() -> Slide:
    s = Slide(
        title="CANlab code map: fmri_data.canlab_connectivity_preproc",
        description=[
            "Prepare a 4-D fmri_data for connectivity analysis: regress nuisance covariates,",
            "remove WM/CSF signal, bandpass-filter, residualize, windsorize, and extract ROI averages.",
        ],
    )

    s.use([
        "Run on a denoised 4-D time series fmri_data; outputs the residualized data and per-ROI signals.",
        "Pass nuisance regressors (motion, spikes) and optional ROI masks to extract from.",
        "Filter step uses conn_filter (FFT-based bandpass); pass cutoffs in Hz.",
        "Optional 'additional GLM' step runs a separate model with extra regressors and returns betas.",
    ])

    s.notes([
        "Filtering is applied to BOTH data and covariates before regression (Lindquist 2019 orthogonalization).",
        "Extracting ROI averages first then filtering gives the same answer as filtering then averaging.",
    ], y=3.40, h=1.50)

    s.legend(x=0.30, y=5.40)

    spine_x = 4.90
    y0 = 1.10
    gap = 0.65
    s.spine([
        ("indat", "var", "fmri_data\n(4-D time series)"),
        ("step1", "func", "1. nuisance + linear trend"),
        ("step2", "func", "2. remove WM / CSF"),
        ("step3", "func", "3. bandpass (conn_filter)"),
        ("step4", "func", "4. residualize via regression"),
        ("step5", "func", "5. (optional) additional GLM"),
        ("step6", "func", "6. windsorize"),
        ("step7", "func", "7. extract ROI averages"),
    ], x=spine_x, y0=y0, row_gap=gap, box_w=2.85, box_h=0.50, bold_ids=[])

    # Outputs to the right
    out_x = spine_x + 3.20
    s.box("dat_out", "var", "preprocessed_dat\n(residualized fmri_data)",
          x=out_x, y=y0 + 5 * gap - 0.10, h=0.65, w=2.80, font_size=10.0)
    s.box("roi_val", "var", "roi_val\n(ROI signals or pattern expr)",
          x=out_x, y=y0 + 6 * gap, h=0.65, w=2.80, font_size=10.0)
    s.box("beta_dat", "var", "beta_dat / beta_roi_val\n(optional GLM)",
          x=out_x, y=y0 + 7 * gap + 0.05, h=0.65, w=2.80, font_size=10.0)
    s.connect_line("step7", "dat_out", src_side="topright", dst_side="left")
    s.connect_line("step7", "roi_val", src_side="right", dst_side="left")
    s.connect_line("step5", "beta_dat", src_side="right", dst_side="left", kind="straight")

    # Mask inputs from upper right
    s.box("masks", "file", "canonical WM / CSF masks\nROI masks (optional)",
          x=out_x, y=y0, h=0.70, w=2.80, font_size=10.0)
    s.connect_line("masks", "step2", src_side="left", dst_side="right")

    s.object_types_panel(
        x=out_x, y=y0 + 1.10,
        lines=["fmri_data", "image_vector"],
        width=2.40,
    )
    return s
