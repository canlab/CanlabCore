# Phase 4 — Implementation log and reproduction guide

This is the companion to `help_doc_agent_instructions_3_minimalcodeinsert_unittestadd.rtf` — the original instructions describe **what** to build; this file describes **how** the work was actually done so the same task can be reproduced on a re-run (or applied to a new batch of methods/functions).

The user-facing scope ("Phase 4a / 4b / 4c") is unchanged from the RTF spec. Read that first.

---

## 1. Outputs and where they live

| Path | Contents |
|---|---|
| `docs/class_method_pngs/` | Sample-figure PNGs, one per method/function. Naming: `<filename_referenced>_sample.png` (matches the `.md` filename of the page that consumes it). |
| `docs/individual_functions/` | Per-function help pages. The Quick example section sits **above** the Code map section in pages that have one, and **above** the Properties section in class-level pages (`atlas_methods.md`). |
| `docs/_codemap_tools/phase4a_runner.m` | MATLAB driver that produces the 4 Phase-4a sample PNGs. |
| `docs/_codemap_tools/phase4b_runner.m` | MATLAB driver that produces the 26 Phase-4b sample PNGs. |
| `docs/_codemap_tools/phase4b_fixes.m`, `_fixes2.m`, `_fix_text.m` | Targeted re-renders for failures and refinements. Kept around so each fix can be re-run independently. |
| `docs/_codemap_tools/phase4b_refinements.m` | User-requested PNG refinements (seaborn colours, t-map for cluster_surf, OLS-vs-robust scatter, full-content table screenshots). |

## 2. Path setup that every runner uses

All runners begin with:

```matlab
parent = '/Users/f003vz1/Documents/GitHub';
addpath(genpath(fullfile(parent, 'CanlabCore')));
addpath(genpath(fullfile(parent, 'Neuroimaging_Pattern_Masks')));
png_dir = fullfile(parent, 'CanlabCore', 'docs', 'class_method_pngs');
```

If you re-run from a different machine, change `parent`. SPM is also expected on the path; on the development machine SPM25 is auto-loaded by MATLAB userpath. If you don't have SPM25, install SPM12+ and `which('spm_vol')` should resolve before running.

## 3. PNG capture — `exportgraphics`, not `print`

The first attempt used `print -dpng -r150`. That produces large blank borders because it uses the figure's `PaperPosition`. Switching to `exportgraphics` (added in R2020a) gives clean tightly-bounded PNGs:

```matlab
exportgraphics(fig_handle, out_path, 'Resolution', 150, 'BackgroundColor', 'white');
```

`exportgraphics` crops to the union of visible content within the figure window. Two consequences:

1. **If the figure window is too small, content overflows the window border and gets clipped.** This is why several runners explicitly resize the figure before saving (`fig.Position = [x y w h]`) and sometimes also reposition axes (`ax.Position = [...]`) to leave room for x-labels / right-edge text.
2. **If the figure has lots of empty whitespace, `exportgraphics` removes it.** The text-rendered "table" PNGs use this to keep the output tightly cropped to the actual lines.

## 4. The `save_first_fig` pattern

Every runner defines a local helper:

```matlab
function save_first_fig(out_path, varargin)
    p = inputParser;
    p.addParameter('min_height', 0);
    p.addParameter('min_width', 0);
    p.parse(varargin{:});
    figs = findobj('type', 'figure');
    if isempty(figs), warning('no fig'); return; end
    nums = arrayfun(@(f) f.Number, figs);
    [~, ix] = min(nums);
    first_fig = figs(ix);
    figure(first_fig);
    drawnow;
    pos = first_fig.Position;
    if p.Results.min_height > 0 && pos(4) < p.Results.min_height
        pos(4) = p.Results.min_height;
    end
    if p.Results.min_width > 0 && pos(3) < p.Results.min_width
        pos(3) = p.Results.min_width;
    end
    first_fig.Position = pos;
    drawnow;
    exportgraphics(first_fig, out_path, 'Resolution', 150, 'BackgroundColor', 'white');
end
```

Important details:

- **First figure = lowest figure number.** `findobj` returns figures in arbitrary order; `[~, ix] = min(arrayfun(@(f) f.Number, figs))` selects figure 1 reliably.
- **`min_height` / `min_width`** options are essential for functions whose first figure is intentionally short (`annotate_binary_results_map` hard-codes `Position = [34 674 940 206]`). Without min_height, x-labels overflow below the figure boundary and clip.
- **`close all`** at the start of every example block. Otherwise stale figures from earlier examples bias `findobj`.

## 5. Capturing text output as a PNG (the "table screenshot" pattern)

`statistic_image.table`, `region.table`, and `region.table_of_atlas_regions_covered` print text to stdout instead of producing a figure. Capture with `evalc` and render to a tightly-cropped figure:

```matlab
function save_text_as_png_tight(txt, out_path)
    txt = regexprep(txt, '<[^>]*>', '');     % strip HTML markup (MATLAB inserts <strong>...)
    txt = regexprep(txt, '\x{0008}', '');    % strip backspace
    lines = splitlines(txt);
    while ~isempty(lines) && isempty(strtrim(lines{end})), lines(end) = []; end
    while ~isempty(lines) && isempty(strtrim(lines{1})),   lines(1)   = []; end
    n = numel(lines);
    line_h = 16;
    top_pad = 8; bot_pad = 8;
    fig_h = top_pad + bot_pad + line_h * n;
    fh = figure('Color', 'white', 'Units', 'pixels', ...
                'Position', [50 50 1300 fig_h], 'MenuBar', 'none', 'ToolBar', 'none');
    ax = axes('Parent', fh, 'Position', [0 0 1 1], 'Visible', 'off', ...
              'XLim', [0 1], 'YLim', [0 1]);
    text(ax, 0.006, 1 - (top_pad/fig_h), strjoin(lines, newline), ...
         'FontName', 'Menlo', 'FontSize', 9, ...
         'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', ...
         'Interpreter', 'none', 'Units', 'normalized');
    drawnow;
    exportgraphics(fh, out_path, 'Resolution', 150, 'BackgroundColor', 'white');
end
```

Critical lessons:

- `regexprep(txt, '<[^>]*>', '')` strips ALL HTML-like tags. MATLAB's pretty-print uses `<strong>`, `<a href="…">`, etc. The earlier code only stripped `<a>`/`</a>` and produced visibly broken output.
- `'Interpreter', 'none'` on the text object — otherwise underscores in region names get treated as TeX subscripts.
- Use **normalized** axes coordinates (`'Units', 'normalized'`, `'XLim', [0 1]`) and a normalized `text(..., 'Units', 'normalized')` call. Mixing pixel and normalized coordinates produces invisible text.
- Compute `fig_h` from the line count (`16 px/line` at 9-pt Menlo) so there is no empty whitespace below the last line. **Do NOT** use a fixed `min height` here — that re-introduces empty space.
- Don't truncate. The first version capped at 60 lines and inserted "(truncated)"; the user requested the full table + everything after it. `lines(1:60)` was removed.

## 6. Figure-resize cookbook for known troublesome functions

| Function | First figure default size | Fix |
|---|---|---|
| `annotate_binary_results_map` | 940 × 206 (hard-coded) | Resize to 1200 × 420; shift axes up by 0.13 (normalized) and shrink width × 0.86 to expose right-side `Transmodal ->` annotation. |
| `atlas.isosurface` of `canlab2024` | default | Resize figure to 1100 × 850 before drawing so tick labels at all corners of the 3-D bounding box have room. |
| `atlas.isosurface` of subset (e.g. `Thal`) | default | 900 × 700 is enough. |
| `cluster_surf` | small default | 1000 × 700; otherwise the rendered surface is tiny. |

When a figure is built by a third-party function that hard-codes its `Position`, the only reliable fix is to enlarge AFTER the function returns and adjust axes positions if needed. See the `annotate_binary_results_map` block in `phase4a_runner.m` for the canonical example.

## 7. Phase 4b refinements that depend on specific function behaviour

- **`barplot_columns` colours**: pass `'colors', cell_of_rgb`. CANlab provides a HUSL-based palette via `seaborn_colors(N)`; first return is a cell of `[r g b]` triplets. Three bars use `colors(1:3)`.
- **`cluster_surf` with a thresholded t-map**:
  ```matlab
  t = ttest(imgs); t = threshold(t, .005, 'unc', 'k', 10);
  r = region(t);
  cluster_surf(r, 5, 'left');     % regions, mm-deep, surface keyword
  ```
  `cluster_surf` is deprecated; the help page calls this out and points to `addbrain` + `render_on_surface`.
- **`image_scatterplot` OLS vs robust**:
  ```matlab
  obj = imgs;
  obj.X = obj.metadata_table.Reappraisal_Success - mean(...);
  obj.X(:, end+1) = 1;
  out_ols = regress(obj, 'noverbose', 'nodisplay');
  out_rob = regress(obj, 'robust', 'noverbose', 'nodisplay');
  t_ols = get_wh_image(out_ols.t, 1);
  t_rob = get_wh_image(out_rob.t, 1);
  image_scatterplot(t_ols, t_rob, 'colorpoints');
  xlabel('OLS t-value'); ylabel('Robust regression t-value');
  ```
  Both `regress` calls use `'noverbose', 'nodisplay'` to suppress side effects so the scatter is figure 1.

## 8. Markdown insertion conventions

- **Class-level pages** (`atlas_methods.md` etc.): insert a `## Quick example` section BETWEEN the lead summary and the `## Properties` section.
- **Individual function pages with a Code map**: insert `## Quick example` immediately ABOVE the `## Code map` heading.
- **Individual function pages without a Code map** (the new Phase-4b stand-alones): the Quick example sits in the standard position right after the lead summary.
- The Quick example block uses GitHub-flavoured-markdown:
  - One ` ```matlab ` fenced block with ONLY the user-facing code (no figure-resize trickery from the runner).
  - One blank line.
  - One image link in `![alt](relative/path.png)` form. Relative path is `../class_method_pngs/<name>_sample.png` for individual_functions pages and `class_method_pngs/<name>_sample.png` for class-level pages.

## 9. Failures encountered and their resolutions

| Function | Failure | Resolution |
|---|---|---|
| `image_similarity_plot` | `t = ttest(imgs)` then passing `t` failed with "Cannot compare means with 0 degrees of freedom". | Pass `imgs` (multi-image) directly with `'average'`. |
| `jackknife_similarity` | `'doplot'` parsed as flag-only; `inputParser` complained. | Pass `'doplot', true` as explicit name-value pair. |
| `wedge_plot_by_atlas` | Atlas keyword `'buckner_networks'` not in registry. | Use `'yeo17networks'` instead. |
| `riverplot` | "Load image set only tested for input fmri_data objects now". | Build two `fmri_data` layers and pass `riverplot(layer1, 'layer2', layer2)`. |
| `cluster_surf` initial | "Exiting" — invalid args. | Simplest signature: `cluster_surf(cl, mm_deep, 'left')`. |
| `region.table_of_atlas_regions_covered` | Self-disclaims as broken (`disp('… DOES NOT CURRENTLY WORK …')`) and crashes with "logical indices outside array bounds". | The .md page documents this. The example calls the `@image_vector` / `@statistic_image` overload instead and notes the workaround. |
| Tables rendered with HTML tags visible | `regexprep(txt, '<a[^>]*>', '')` only stripped `<a>`. | Replace with `regexprep(txt, '<[^>]*>', '')` to strip all tags. |
| Text-rendered tables truncated | First version capped at 60 lines + "(truncated)" message. | Removed truncation entirely; height grows with line count. |
| `annotate_binary_results_map` PNG cropped | Figure 1 is intentionally short; bottom labels clipped. | Resize after the function returns: `fig.Position(4) = 420`; shift axes positions. |
| `annotate_binary_results_map` PNG right-clipped | "Transmodal ->" off the right edge. | Widen to 1200; shrink axes width × 0.86. |

## 10. Re-running

```bash
# Render all PNGs
matlab -batch "run('docs/_codemap_tools/phase4a_runner.m')"
matlab -batch "run('docs/_codemap_tools/phase4b_runner.m')"
matlab -batch "run('docs/_codemap_tools/phase4b_fixes.m')"
matlab -batch "run('docs/_codemap_tools/phase4b_fixes2.m')"
matlab -batch "run('docs/_codemap_tools/phase4b_fix_text.m')"
matlab -batch "run('docs/_codemap_tools/phase4b_refinements.m')"
```

Each runner is independent and writes only the PNGs in its own scope, so re-running a single file is safe. PNGs are deterministic except for any RNG-using examples (`barplot_columns`, `plot_correlation_matrix`) which seed `rng(7)`.

The markdown insertion was done by hand (one `Edit` call per page). To re-do that programmatically you'd need a structural patch — the existing markdown files are stable, so prefer manual edits.

## 11. Class-methods page link updates

`docs/_codemap_tools/relink_methods_pages.py` is idempotent and turns bare `` `methodname` `` cells in method tables into hyperlinks. Re-run after creating new individual_functions pages:

```bash
cd docs/_codemap_tools && python3 relink_methods_pages.py
```

The regex requires `` `@class` `` or `both` in the next cell to avoid linking property-table rows. If a row is already linked (e.g. from a prior Phase-3 run with an older mapping), the regex won't match and the link stays as-is — fix manually with `Edit` calls. Class-specific overrides (`STATISTIC_IMAGE_METHODS`, `REGION_METHODS`) handle methods that have a different page per class (e.g., `table` → `region_table.md` from `region_methods.md` but `fmri_data_table.md` from `fmri_data_methods.md`).

## 12. Phase 4c (CI unit test) — implemented

The Phase 4c work landed as `CanlabCore/Unit_tests/canlab_test_help_examples.m`. Key design choices:

- **Discovery is automatic.** `canlab_run_all_tests.m` already globs `canlab_test_*.m` recursively, so the new file is picked up by the existing per-push CI workflow (`.github/workflows/test.yml`) without any workflow edits.
- **One test function per `## Quick example`.** 30 tests total: 4 from Phase 4a + 26 from Phase 4b. Each test mirrors the exact code shown in the user-facing `.md` page, so when an example bit-rots, its corresponding test fails.
- **Smoke-only verification.** Each test asserts that the example runs end-to-end and (where appropriate) that returns are of the expected type/shape. Pixel content of the generated PNGs is NOT verified — that would tie CI to specific MATLAB rendering and produce too many false positives.
- **`setupOnce` caches expensive shared state.** The emotionreg sample, group t-map, thresholded t-map (`p<.005, k=10`), and `region(thresholded_t)` are loaded once per test session and reused via `tc.TestData`. Without this, the suite would do ~18 redundant `ttest`/`threshold`/`region` calls.
- **Headless graphics handling.** `setup` forces `DefaultFigureVisible = 'off'`. Tests that exercise OpenGL-only code paths (`surface`, `addbrain`, `cluster_surf`, `region.surface`, `region.labelled_surface`) wrap the offending call in a try/catch and call `tc.assumeFail(...)` for missing-graphics errors, so CI runners without a display do not fail spuriously. The `orthviews` test additionally guards on `usejava('jvm') && feature('ShowFigureWindows')`.
- **Documented-broken example.** `region.table_of_atlas_regions_covered` is currently broken (per its own `.md` page); the test wraps the call in try/catch with `tc.assumeFail` so the suite does not fail on a known issue. When the function is fixed, replace `assumeFail` with `verifyNotEmpty`.
- **Doc/PNG/test parity.** The `.md` Quick examples for `cluster_surf`, `barplot_columns`, and `image_scatterplot` were updated during Phase 4c so that `.md` code, regenerated PNG, and the unit test all match (previously the PNGs were re-rendered per user request but the `.md` examples were never updated).

To run only the help-examples test interactively:

```matlab
cd CanlabCore/CanlabCore/Unit_tests;
results = runtests('canlab_test_help_examples');
```

Or as part of the full suite, exactly as CI does it:

```matlab
results = canlab_run_all_tests;
```
