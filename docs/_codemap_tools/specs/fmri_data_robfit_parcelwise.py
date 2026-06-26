"""Code map: fmri_data.robfit_parcelwise — robust regression at the parcel level."""

from codemap_lib import Slide


def build() -> Slide:
    s = Slide(
        title="CANlab code map: fmri_data.robfit_parcelwise",
        description=[
            "Robust multiple regression on per-parcel mean activations from an atlas.",
            "Produces t-maps, atlas-labeled tables, montages/surfaces, and per-subject diagnostics in one OUT struct.",
        ],
    )

    s.use([
        "Pass a 2nd-level fmri_data with .X already filled (one row per subject, one column per regressor).",
        "Use 'names', {...} to label regressors; 'noplot' / 'plotdiagnostics' to control figure output.",
        "Underlying parcellation defaults to a CANlab atlas — pass 'atlas', atlas_obj to override.",
    ])

    s.notes([
        "Running per-parcel cuts compute time and reduces multiple-comparison burden vs voxelwise robust regression.",
        "OUT contains the t-map (statistic_image), tables, diagnostics struct, and per-image / per-parcel intermediates.",
    ], y=3.40, h=1.50)

    s.legend(x=0.30, y=5.30)

    spine_x = 4.90
    y0 = 1.20
    gap = 0.80
    s.spine([
        ("indat", "var", "fmri_data\n(images, .X = design)"),
        ("rfp", "func", "robfit_parcelwise( )"),
        ("ap", "func", "apply_parcellation"),
        ("rob", "func", "robustfit per parcel"),
        ("diag", "func", "coverage / scaling /\nglobal-signal diagnostics"),
    ], x=spine_x, y0=y0, row_gap=gap, box_h=0.62, bold_ids=["rfp"])

    # Atlas input
    s.box("atlas", "var", "atlas object\n(load_atlas)",
          x=spine_x + 2.50, y=y0, h=0.65)
    s.connect_line("atlas", "rfp", src_side="bottomleft", dst_side="topright")

    # OUT struct — single arrow from the last spine node (diag) to avoid the
    # crossing-through-diag arrow that the previous rob→OUT line caused.
    out_y = y0 + 5 * gap + 0.10
    s.box("OUT", "var",
          "OUT struct\n(t-map statistic_image, results tables,\nmontages/surfaces, diagnostics)",
          x=spine_x - 1.40, y=out_y, h=0.95, w=4.80, font_size=10.0)
    s.connect_line("diag", "OUT", src_side="bottom", dst_side="top")

    s.object_types_panel(
        x=spine_x + 2.50, y=y0 + gap,
        lines=["fmri_data", "atlas", "statistic_image", "MATLAB table"],
        width=2.40,
    )
    return s
