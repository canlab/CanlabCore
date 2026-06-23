"""Code map: fmri_data.predict — cross-validated multivariate prediction."""

from codemap_lib import Slide


def build() -> Slide:
    s = Slide(
        title="CANlab code map: fmri_data.predict",
        description=[
            "Train and cross-validate a multivariate predictive model from brain images.",
            "Returns CV error, a stats struct (weight map, fold predictions, error metrics), and per-fold raw outputs.",
        ],
    )

    s.use([
        "Attach predictors as obj.dat (voxels × images) and outcomes as obj.Y.",
        "Choose an algorithm with 'algorithm_name' (cv_pcr, cv_lassopcr, cv_svr, cv_svm, cv_pls, ...).",
        "k-fold CV is stratified on outcome by default; pass a vector of fold IDs to use custom holdouts.",
        "Set 'bootweights' for bootstrap inference on voxel weights; 'numcomponents' for PCA-based methods.",
    ])

    s.notes([
        "Each fold trains on (k-1) folds and tests on the held-out one. Algorithm dispatch happens by function name — drop a new train/test .m on the path to add a model.",
        "stats.weight_obj is a statistic_image holding the model's voxel weights; pass it through threshold(), region(), montage() like any other map.",
    ], y=3.40, h=2.10)

    s.legend(x=0.30, y=5.85)

    spine_x = 4.90
    y0 = 1.20
    gap = 0.85
    s.spine([
        ("indat", "var", "obj.dat\n(voxels × images)"),
        ("predict", "func", "predict( )"),
        ("folds", "func", "stratified holdout\n(k folds)"),
        ("algo", "func", "<algorithm_name>\ne.g. cv_lassopcr"),
        ("stats", "var", "stats struct\n+ cverr, optout"),
    ], x=spine_x, y0=y0, row_gap=gap, box_h=0.62, bold_ids=["predict"])

    # obj.Y joins predict from the right
    s.box("Y", "var", "obj.Y\n(outcome)",
          x=spine_x + 2.50, y=y0, h=0.62)
    s.connect_line("Y", "predict", src_side="bottomleft", dst_side="topright")

    # Optional bootweights branch off algo
    s.box("boot", "func", "bootweights\n(resample → P, FDR)",
          x=spine_x + 2.50, y=y0 + 3 * gap - 0.05, h=0.62)
    s.connect_line("algo", "boot", src_side="right", dst_side="left")

    # Outputs from stats
    out_y = y0 + 5 * gap
    s.box("weight_obj", "var", "stats.weight_obj\n(statistic_image)",
          x=spine_x - 2.40, y=out_y, h=0.62)
    s.box("yfit", "var", "stats.Y_xval, yfit_xval\n(per-fold predictions)",
          x=spine_x + 2.40, y=out_y, h=0.62)
    s.connect_line("stats", "weight_obj", src_side="bottomleft", dst_side="topright")
    s.connect_line("stats", "yfit", src_side="bottomright", dst_side="topleft")

    s.object_types_panel(
        x=spine_x + 4.85, y=y0,
        lines=["fmri_data", "statistic_image", "predictive_model"],
        width=2.10,
    )
    return s
