"""Phase 3: hyperlink documented method names in the per-class methods pages.

For each method that now has an individual_functions/<name>.md page, replace
the bare-backtick mention in the table row of each *_methods.md file with a
hyperlink. Uses precise regex matching so that only the FIRST table cell
column (the method name) is hyperlinked, leaving the @class column and the
description column unchanged.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

DOCS = Path("/Users/f003vz1/Documents/GitHub/CanlabCore/docs")

# Mapping: method name -> individual_functions/<file>.md basename
METHOD_TO_PAGE: dict[str, str] = {
    # @fmri_data
    "regress": "fmri_data_regress",
    "ttest": "fmri_data_ttest",
    "predict": "fmri_data_predict",
    "fitlme_voxelwise": "fmri_data_fitlme_voxelwise",
    "robfit_parcelwise": "fmri_data_robfit_parcelwise",
    "rsa_regression": "fmri_data_rsa_regression",
    "dual_regression": "fmri_data_dual_regression",
    "evaluate_spatial_scale": "fmri_data_evaluate_spatial_scale",
    "model_brain_pathway": "fmri_data_model_brain_pathway",
    "model_mpathi": "fmri_data_model_mpathi",
    "annotate_binary_results_map": "fmri_data_annotate_binary_results_map",
    "canlab_connectivity_preproc": "fmri_data_canlab_connectivity_preproc",
    "denoise_timeseries_pipeline": "fmri_data_denoise_timeseries_pipeline",
    "extract_measures_batch": "fmri_data_extract_measures_batch",
    "extract_roi_averages": "fmri_data_extract_roi_averages",
    "normalize_gm_by_wm_csf": "fmri_data_normalize_gm_by_wm_csf",
    "ttest_table_and_lateralization_test": "fmri_data_ttest_table_and_lateralization_test",
    # @image_vector (inherited by fmri_data, statistic_image, atlas)
    "pca": "fmri_data_pca",
    "outliers": "fmri_data_outliers",
    "descriptives": "fmri_data_descriptives",
    "apply_parcellation": "fmri_data_apply_parcellation",
    "extract_gray_white_csf": "fmri_data_extract_gray_white_csf",
    "qc_metrics_second_level": "fmri_data_qc_metrics_second_level",
    "jackknife_similarity": "fmri_data_jackknife_similarity",
    "hansen_neurotransmitter_maps": "fmri_data_hansen_neurotransmitter_maps",
    "image_similarity_plot": "fmri_data_image_similarity_plot",
    "table_of_atlas_regions_covered": "fmri_data_table_of_atlas_regions_covered",
    "wedge_plot_by_atlas": "fmri_data_wedge_plot_by_atlas",
    # @atlas
    "select_atlas_subset": "atlas_select_atlas_subset",
}

# Methods that have a class-specific override page rather than the shared one.
# `table` lives in @image_vector / @statistic_image / @region — the codemap
# documents the @image_vector base, so all subclasses link to fmri_data_table.md.
TABLE_PAGE = "fmri_data_table"

# `threshold` lives in @atlas / @image_vector / @statistic_image. Codemap is
# the @statistic_image variant — only that class's row should link.
THRESHOLD_PAGE = "statistic_image_threshold"

# Files to update. (path, special_case_dict where key=method, value=which page or None to skip)
FILES = [
    "fmri_data_methods.md",
    "image_vector_methods.md",
    "statistic_image_methods.md",
    "atlas_methods.md",
    "region_methods.md",
]

# A table row looks like:    | `methodname` | `@class` | description |
# We want to replace just the first cell.
ROW_PATTERN = re.compile(
    r"^(\| )(`(?P<name>[a-zA-Z_][a-zA-Z0-9_]*)`)( \| )",
    re.MULTILINE,
)


def link_for(method: str, file_basename: str) -> str | None:
    """Return the individual_functions/<page>.md basename, or None to skip."""
    if method == "table":
        return TABLE_PAGE
    if method == "threshold":
        # Only link in statistic_image_methods.md; @atlas.threshold is a different method.
        if file_basename == "statistic_image_methods.md":
            return THRESHOLD_PAGE
        if file_basename == "image_vector_methods.md":
            return THRESHOLD_PAGE
        if file_basename == "fmri_data_methods.md":
            return THRESHOLD_PAGE
        return None
    return METHOD_TO_PAGE.get(method)


def relink_file(path: Path) -> int:
    text = path.read_text()
    n = 0

    def repl(m: re.Match) -> str:
        nonlocal n
        name = m.group("name")
        page = link_for(name, path.name)
        if page is None:
            return m.group(0)
        n += 1
        return f"{m.group(1)}[`{name}`](individual_functions/{page}.md){m.group(4)}"

    new = ROW_PATTERN.sub(repl, text)
    if new != text:
        path.write_text(new)
    return n


def main() -> int:
    total = 0
    for fname in FILES:
        path = DOCS / fname
        if not path.exists():
            print(f"  (skip, missing) {fname}")
            continue
        n = relink_file(path)
        print(f"  {fname}: {n} link(s) added")
        total += n
    print(f"\nTotal links added: {total}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
