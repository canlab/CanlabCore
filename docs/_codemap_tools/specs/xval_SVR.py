"""Code map: xval_SVR — within or between-person SVR with repeated CV and bootstrap inference."""

from codemap_lib import Slide


def build() -> Slide:
    s = Slide(
        title="CANlab code map: xval_SVR",
        description=[
            "Cross-validated support vector regression using fitrsvm. The SVR analogue of xval_SVM —",
            "subject-grouped folds, optional nested CV hyperparameter optimization, bootstrap inference on weights.",
        ],
    )

    s.use([
        "Pass S = xval_SVR(X, Y, id, ...). X is observations × features; Y is the continuous outcome; id keeps subjects together.",
        "Toggle 'doOptimize', true for Bayesian hyperparameter optimization (5-fold inner CV).",
        "'doRepeats', N repeats outer CV with new random splits to stabilize accuracy estimates.",
        "'highdimensional', true switches to fitrlinear for voxel-scale brain data (no hyperopt).",
    ])

    s.notes([
        "Linear regressors only; nonlinear models are commented-out hooks in the source.",
        "Optimizing hyperparameters AND repeating CV gives nested CV — slow but unbiased.",
        "Default loss is mean squared error; correlation between Y and yfit is the headline accuracy metric.",
    ], y=3.40, h=1.80)

    s.legend(x=0.30, y=5.60)

    spine_x = 4.80
    y0 = 1.10
    gap = 0.65
    s.spine([
        ("X", "var", "X (obs × features)"),
        ("xs", "func", "xval_SVR( )"),
        ("hold", "func", "stratified holdout\n(grouped by id)"),
        ("sanity", "func", "fit overall model\n(sanity-check fit)"),
        ("opt", "func", "(optional) nested CV\nhyperparam optimization"),
        ("cv", "func", "fitrsvm per fold\n→ yfit, dist_from_hyperplane"),
        ("repeats", "func", "(optional) repeat CV\nwith new splits"),
        ("boot", "func", "bootstrap weights\n→ FDR-corrected P-values"),
        ("S", "var", "S struct\n(yfit, MSE, r, weights, p-values)"),
    ], x=spine_x, y0=y0, row_gap=gap, box_h=0.50, bold_ids=["xs"])

    rx = spine_x + 2.80
    s.box("Y", "var", "Y (continuous outcome)",
          x=rx, y=y0, h=0.55, w=2.60, font_size=10.0)
    s.box("id", "var", "id (subject grouping)",
          x=rx, y=y0 + 0.65, h=0.55, w=2.60, font_size=10.0)
    s.connect_line("Y", "xs", src_side="left", dst_side="right")
    s.connect_line("id", "xs", src_side="left", dst_side="right")

    s.box("plot", "script", "scatter + residual plots",
          x=rx, y=y0 + 1.40, h=0.55, w=2.60, font_size=10.0)
    s.connect_line("cv", "plot", src_side="right", dst_side="left")

    s.object_types_panel(
        x=rx, y=y0 + 2.10,
        lines=["MATLAB struct", "fitrsvm model"],
        width=2.50,
    )
    return s
