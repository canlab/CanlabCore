"""Code map: fmri_data.annotate_binary_results_map — annotate a binary map against
gradients, neurochemistry, networks, and topics."""

from codemap_lib import Slide


def build() -> Slide:
    s = Slide(
        title="CANlab code map: fmri_data.annotate_binary_results_map",
        description=[
            "Compare a binary results map to a battery of reference maps — connectivity gradients, transcriptomic gradients,",
            "PET neurotransmitter maps, Neurosynth topics, and resting-state networks — and report the associations.",
        ],
    )

    s.use([
        "Pass an fmri_data containing a thresholded / binary map (e.g., from threshold(t, .005)).",
        "Renders kernel-density histograms and matrix plots; correlates the map against each reference set.",
        "Useful as a one-call follow-up after a primary analysis to characterize *what* the active regions tend to overlap.",
    ])

    s.notes([
        "Reference maps are loaded from companion files via load_image_set / load_atlas — Neuroimaging_Pattern_Masks must be on the path.",
        "Outputs are returned as a single RESULTS struct with one field per reference set.",
    ], y=3.40, h=1.50)

    s.legend(x=0.30, y=5.30)

    spine_x = 4.70
    y0 = 1.20
    gap = 0.95
    s.spine([
        ("indat", "var", "binary fmri_data\n(thresholded map)"),
        ("annotate", "func", "annotate_binary_results_map( )"),
        ("compare", "func", "kernelDensitySmoothedHistogram\n+ correlation per reference set"),
    ], x=spine_x, y0=y0, row_gap=gap, box_w=2.85, box_h=0.65, bold_ids=["annotate"])

    # Reference sets fan in from the right
    refs = [
        ("ref1", "Principal gradient\n(transmodal vs unimodal)"),
        ("ref2", "Allen transcriptomic\ngradients"),
        ("ref3", "Hansen Neuromaps\nPET tracers"),
        ("ref4", "Neurosynth\ntopics & terms"),
        ("ref5", "Yeo/Buckner\nresting-state nets"),
    ]
    rx = spine_x + 3.10
    for i, (nid, label) in enumerate(refs):
        s.box(nid, "var", label, x=rx, y=y0 + i * 0.62, h=0.55, w=2.80, font_size=10.0)
        s.connect_line(nid, "annotate", src_side="left", dst_side="right")

    s.annotate("reference map sets", x=rx, y=y0 + 5 * 0.62 + 0.05, w=2.80, h=0.25, size=10.0)

    # Output
    s.box("results", "var", "RESULTS struct\n(one field per ref set)",
          x=spine_x, y=y0 + 3 * gap + 0.10, h=0.65, w=2.85, font_size=10.0)
    s.connect_line("compare", "results", src_side="bottom", dst_side="top")

    s.object_types_panel(
        x=rx, y=y0 + len(refs) * 0.62 + 0.50,
        lines=["fmri_data", "MATLAB struct"],
        width=2.50,
    )
    return s
