"""Code map: xval_SVM — within-person SVM classification with repeated CV and bootstrap inference."""

from codemap_lib import Slide


def build() -> Slide:
    s = Slide(
        title="CANlab code map: xval_SVM",
        description=[
            "Repeated-measures SVM classification using fitcsvm. Selects subject-grouped stratified folds,",
            "optionally runs nested CV for hyperparameter optimization, and bootstraps voxel weights with FDR-corrected P-values.",
        ],
    )

    s.use([
        "S = xval_SVM(X, Y, id, ...).",
        "X = obs × features, Y = binary, id groups subjects.",
        "'doOptimize' true → Bayesian hyperopt (5-fold inner CV).",
        "'doRepeats', N → repeat outer CV with new splits.",
        "'highdimensional' true → fitclinear (voxel-scale).",
    ])

    s.notes([
        "ROC and single-interval CV accuracy can differ; ROC_plot picks an optimal threshold.",
        "Linear classifiers only; nonlinear hooks are commented out.",
        "Hyperopt + repeated CV = nested CV. Slow but unbiased.",
    ], y=3.65, h=1.55)

    s.legend(x=0.30, y=5.60)

    spine_x = 4.80
    y0 = 1.10
    gap = 0.65
    s.spine([
        ("X", "var", "X (obs × features)"),
        ("xs", "func", "xval_SVM( )"),
        ("hold", "func", "stratified holdout\n(grouped by id)"),
        ("sanity", "func", "fit overall model\n(sanity-check accuracy)"),
        ("opt", "func", "(optional) nested CV\nhyperparam optimization"),
        ("cv", "func", "fitcsvm per fold\n→ yfit, scores, probability"),
        ("repeats", "func", "(optional) repeat CV\nwith new splits"),
        ("boot", "func", "bootstrap weights\n→ FDR-corrected P-values"),
        ("S", "var", "S struct\n(yfit, accuracy, weights, p-values)"),
    ], x=spine_x, y0=y0, row_gap=gap, box_h=0.50, bold_ids=["xs"])

    # Other inputs
    rx = spine_x + 2.80
    s.box("Y", "var", "Y (binary outcome)",
          x=rx, y=y0, h=0.55, w=2.50, font_size=10.0)
    s.box("id", "var", "id (subject grouping)",
          x=rx, y=y0 + 0.65, h=0.55, w=2.50, font_size=10.0)
    s.connect_line("Y", "xs", src_side="left", dst_side="right")
    s.connect_line("id", "xs", src_side="left", dst_side="right")

    # Plot
    s.box("plot", "script", "score plot + ROC_plot",
          x=rx, y=y0 + 1.40, h=0.55, w=2.50, font_size=10.0)
    s.connect_line("cv", "plot", src_side="right", dst_side="left")

    s.object_types_panel(
        x=rx, y=y0 + 2.10,
        lines=["MATLAB struct", "fitcsvm model"],
        width=2.50,
    )
    return s
