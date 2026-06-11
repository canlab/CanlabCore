"""Code map: fmri_data.jackknife_similarity — leave-one-out similarity vs median reference."""

from codemap_lib import Slide


def build() -> Slide:
    s = Slide(
        title="CANlab code map: fmri_data.jackknife_similarity",
        description=[
            "For each image j, compute similarity between image j and the voxelwise median of the remaining N-1.",
            "Yields a robust per-image agreement score with the dataset's central tendency.",
        ],
    )

    s.use([
        "Pass an image_vector / fmri_data with multiple images, or a numeric [voxels × images] matrix.",
        "Choose 'similarity_metric': 'correlation', 'cosine_similarity', 'dot_product', 'euclidean'.",
        "low_agreement returns indices of images falling below an automatic threshold — useful for QC.",
    ])

    s.notes([
        "Missing-data handling follows canlab_compute_similarity_matrix conventions (zeros treated as missing in data).",
        "Median reference (rather than mean) keeps the metric robust to outlier images.",
    ], y=3.40, h=1.50)

    s.legend(x=0.30, y=5.30)

    spine_x = 5.10
    y0 = 1.20
    gap = 0.85
    s.spine([
        ("indat", "var", "obj.dat\n(voxels × images)"),
        ("jk", "func", "jackknife_similarity( )"),
        ("loop", "func", "for j = 1..N:\nref_j = median(dat[:, ~j])"),
        ("simfn", "func", "canlab_compute_similarity_matrix\n(per metric)"),
    ], x=spine_x, y0=y0, row_gap=gap, box_h=0.65, bold_ids=["jk"])

    # Outputs
    out_y = y0 + 4 * gap + 0.10
    s.box("sim", "var", "sim_values\n[N × 1]",
          x=spine_x - 3.30, y=out_y, h=0.65, w=2.20)
    s.box("d", "var", "d (effect size\nvs. permutation)",
          x=spine_x - 0.85, y=out_y, h=0.65, w=2.40)
    s.box("low", "var", "low_agreement\n(indices)",
          x=spine_x + 1.80, y=out_y, h=0.65, w=2.20)
    s.box("nvox", "var", "Nvox\n(per pair)",
          x=spine_x + 4.20, y=out_y, h=0.65, w=2.00)
    for nid in ("sim", "d", "low", "nvox"):
        s.connect_line("simfn", nid, src_side="bottom", dst_side="top")

    s.object_types_panel(
        x=spine_x + 2.70, y=y0,
        lines=["fmri_data", "image_vector"],
        width=2.10,
    )
    return s
