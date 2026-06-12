"""Code map: fmri_data.rsa_regression — representational similarity regression."""

from codemap_lib import Slide


def build() -> Slide:
    s = Slide(
        title="CANlab code map: fmri_data.rsa_regression",
        description=[
            "Construct a representational dissimilarity matrix (RDM) from spatial covariance across images,",
            "then regress the RDM on a binary task / construct design and bootstrap statistics on each regressor.",
        ],
    )

    s.use([
        "Pass an fmri_data, a binary design matrix (one column per condition / construct), and a 'study' grouping vector.",
        "Distance metric: 'cosine' (default), 'correlation', 'euclidean', 'average_Euclidean'.",
        "Bootstrapping resamples WITHIN study blocks to test generalization across construct/study combinations.",
    ])

    s.notes([
        "Designed for assessing across-construct or across-study generalizability (Kragel et al. 2018, Nat Neurosci).",
        "Currently expects one image per subject (no repeated-measures support yet).",
    ], y=3.40, h=1.50)

    s.legend(x=0.30, y=5.30)

    spine_x = 5.00
    y0 = 1.20
    gap = 0.85
    s.spine([
        ("indat", "var", "fmri_data\n(.dat = images × voxels)"),
        ("rsa", "func", "rsa_regression( )"),
        ("rdm", "func", "build RDM\n(spatial similarity)"),
        ("regress", "func", "regress RDM ~ design"),
        ("boot", "func", "bootstrap within studies"),
    ], x=spine_x, y0=y0, row_gap=gap, box_h=0.62, bold_ids=["rsa"])

    # Design and study inputs from the right
    s.box("design", "var", "design\n(binary [obs × k])",
          x=spine_x + 2.60, y=y0, h=0.65)
    s.box("study", "var", "study\n(integer block IDs)",
          x=spine_x + 2.60, y=y0 + 0.85, h=0.65)
    s.connect_line("design", "rsa", src_side="left", dst_side="right")
    s.connect_line("study", "rsa", src_side="left", dst_side="right")

    # Output
    out_y = y0 + 5 * gap + 0.10
    s.box("stats", "var",
          "stats struct\n(per-regressor betas, bootstrap CIs and p-values)",
          x=spine_x - 1.30, y=out_y, h=0.95, w=4.60, font_size=10.0)
    s.connect_line("boot", "stats", src_side="bottom", dst_side="top")

    s.object_types_panel(
        x=spine_x + 2.60, y=y0 + 1.80,
        lines=["fmri_data", "image_vector"],
        width=2.30,
    )
    return s
