"""Code map: xval_classify — k-fold cross-validated discriminant classification."""

from codemap_lib import Slide


def build() -> Slide:
    s = Slide(
        title="CANlab code map: xval_classify",
        description=[
            "Stratified k-fold cross-validated linear discriminant classification using MATLAB's fitcdiscr.",
            "Optionally optimizes Delta / Gamma / DiscrimType hyperparameters via nested CV.",
        ],
    )

    s.use([
        "Pass observation matrix X [N × M] and integer class label vector labels [N × 1].",
        "Use 'id', vec to keep all observations from the same subject in the same fold (leave-whole-subject-out).",
        "Set 'optimizeHyperparameters', true for nested-CV Bayesian hyperparameter search.",
        "Returns S struct with per-fold models, predictions, true labels, and accuracies.",
    ])

    s.notes([
        "Fold assignment: with 'id', uses xval_stratified_holdout_leave_whole_subject_out; otherwise stratified_holdout_set.",
        "Default optimization uses 'auto' hyperparameter set with expected-improvement-plus acquisition.",
    ], y=3.40, h=1.50)

    s.legend(x=0.30, y=5.30)

    spine_x = 4.90
    y0 = 1.20
    gap = 0.78
    s.spine([
        ("X", "var", "X [N × M]"),
        ("xc", "func", "xval_classify( )"),
        ("folds", "func", "stratified holdout\n(grouped by id if given)"),
        ("opt", "func", "(optional) nested CV\nhyperparam optimization"),
        ("fit", "func", "fitcdiscr per fold"),
        ("S", "var", "S struct\n(models, predictions, accuracy)"),
    ], x=spine_x, y0=y0, row_gap=gap, box_h=0.60, bold_ids=["xc"])

    # Other inputs
    rx = spine_x + 2.60
    s.box("labels", "var", "labels [N × 1]",
          x=rx, y=y0, h=0.55, w=2.30, font_size=10.0)
    s.box("id", "var", "id (optional)\nleave-subj-out grouping",
          x=rx, y=y0 + 0.65, h=0.65, w=2.50, font_size=10.0)
    s.box("nfolds", "var", "nFolds (default 5)",
          x=rx, y=y0 + 1.45, h=0.55, w=2.30, font_size=10.0)
    s.connect_line("labels", "xc", src_side="left", dst_side="right")
    s.connect_line("id", "xc", src_side="left", dst_side="right")
    s.connect_line("nfolds", "xc", src_side="left", dst_side="right")

    # Plot
    s.box("plot", "script", "confusion matrix\nper fold (if doplot)",
          x=rx, y=y0 + 4 * gap, h=0.65, w=2.50, font_size=10.0)
    s.connect_line("fit", "plot", src_side="right", dst_side="left")

    s.object_types_panel(
        x=rx, y=y0 + 2.20,
        lines=["MATLAB struct", "fitcdiscr model"],
        width=2.50,
    )
    return s
