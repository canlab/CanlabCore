"""Code map: fmri_data.fitlme_voxelwise — voxelwise mixed-effects models via fitlme."""

from codemap_lib import Slide


def build() -> Slide:
    s = Slide(
        title="CANlab code map: fmri_data.fitlme_voxelwise",
        description=[
            "Run the same MATLAB fitlme model at every voxel of an fmri_data, returning beta, t,",
            "and contrast t-maps as statistic_image objects. Spec is either a formula string or a role table.",
        ],
    )

    s.use([
        "obj.dat = [voxels × observations]; tbl is a MATLAB table with one row per observation (no Y column needed).",
        "Pass either a fitlme formula ('Y ~ 1 + Run + (1+Run|Subject)') OR a role table with Name/Role columns.",
        "Optional p-threshold + 'unc'/'fdr' triggers immediate thresholding; 'C', M evaluates fixed-effect contrasts.",
    ])

    s.notes([
        "Per-voxel fitlme is slow — expect minutes per 1000s of voxels. Use 'mask' or pre-mask the object to constrain.",
        "Role-table mode auto-builds the formula: fixed effects + random slopes per 'mixed' role, grouped by 'group'.",
    ], y=3.40, h=1.60)

    s.legend(x=0.30, y=5.40)

    spine_x = 5.10
    y0 = 1.20
    gap = 0.80
    s.spine([
        ("indat", "var", "obj.dat\n[voxels × obs]"),
        ("fitlme", "func", "fitlme_voxelwise( )"),
        ("formula", "func", "build / parse formula\n(string or role table)"),
        ("loop", "func", "per-voxel fitlme"),
    ], x=spine_x, y0=y0, row_gap=gap, box_h=0.65, bold_ids=["fitlme"])

    # Tbl input from the right
    s.box("tbl", "var", "tbl (MATLAB table)\n[obs × covariates]",
          x=spine_x + 2.60, y=y0, h=0.65, w=2.80, font_size=10.0)
    s.connect_line("tbl", "fitlme", src_side="bottomleft", dst_side="topright")

    # Spec input
    s.box("spec", "var", "spec:\nformula string OR role table",
          x=spine_x + 2.60, y=y0 + gap - 0.05, h=0.65, w=2.80, font_size=10.0)
    s.connect_line("spec", "fitlme", src_side="left", dst_side="right")

    # Outputs
    out_y = y0 + 4 * gap + 0.10
    s.box("betas", "var", "out.beta\nstatistic_image",
          x=spine_x - 2.80, y=out_y, h=0.65)
    s.box("t", "var", "out.t\nstatistic_image",
          x=spine_x - 0.30, y=out_y, h=0.65)
    s.box("con", "var", "out.con_t\n(if 'C' supplied)",
          x=spine_x + 2.20, y=out_y, h=0.65)
    for nid in ("betas", "t", "con"):
        s.connect_line("loop", nid, src_side="bottom", dst_side="top")

    s.object_types_panel(
        x=spine_x + 4.85, y=y0 + 2 * gap,
        lines=["fmri_data", "statistic_image", "MATLAB table"],
        width=2.10,
    )
    return s
