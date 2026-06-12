# `fmridisplay` methods, organized by area

The `fmridisplay` class is a CANlab display container for anatomical underlays
plus one or more activation maps rendered as montages, surfaces, and / or SPM
orthviews. Unlike `fmri_data` / `atlas`, it is **not** a subclass of
`image_vector` — it is a standalone class whose job is to hold the graphics
state of a multi-panel brain figure and to expose methods for adding,
removing, and styling blobs and points on it. The default underlay is
`fmriprep20_template.nii.gz` (MNI152NLin2009cAsym); other underlays can be
supplied via the `'overlay'` constructor argument (e.g., the legacy SPM8
single-subject T1, or the Keuken 2014 enhanced underlay used prior to 2024).

A typical workflow is `o2 = fmridisplay; o2 = montage(o2, ...); o2 =
addblobs(o2, region_obj, ...);` then `legend(o2)`. All methods listed below
are defined on the `fmridisplay` class itself — there is no inheritance.
Use `methods(fmridisplay)` in MATLAB for the live list, and
`help fmridisplay.montage` / `help fmridisplay.render_blobs` for the full
list of slice-display and rendering options.

## Properties

| Property | Description |
|---|---|
| `overlay` | Filename of the anatomical underlay image used for slice montages |
| `SPACE` | Sampling-space struct returned by `define_sampling_space(spm_vol(overlay))` |
| `activation_maps` | Cell array of activation-map structs added by `addblobs` (graphics handles, color info, etc.) |
| `montage` | Cell array of montage structs, one per registered montage, each with `axis_handles`, slice info, etc. |
| `surface` | Cell array of registered surface plot structs (handles + metadata) |
| `orthviews` | Cell array of registered SPM-orthviews plot structs |
| `history` | Cell array: names of methods applied to this object, in order |
| `history_descrip` | Description string for the history field (legacy / introspection) |
| `additional_info` | Free-form struct for attaching arbitrary metadata (legacy / extension slot) |

## Display and visualization

Adding, removing, and styling brain panels and overlays. Most of these
require a graphics environment (won't fully work in headless `matlab -batch`).

| Method | From | One-liner |
|---|---|---|
| `montage` | `fmridisplay` | Create a slice montage (axial / saggital / coronal) and register it on the object |
| `surface` | `fmridisplay` | Add cortical surface(s) to the figure and register them on the object |
| `addblobs` | `fmridisplay` | Render a region/cluster object as blobs on every registered montage and surface |
| `addthreshblobs` | `fmridisplay` | Add blobs from a `statistic_image` at multiple thresholds (single color) |
| `addpoints` | `fmridisplay` | Plot point coordinates (with optional text/markers) on registered montages |
| `removeblobs` | `fmridisplay` | Remove all rendered blobs (and their legends) from the object |
| `removepoints` | `fmridisplay` | Delete plotted point handles from every montage |
| `legend` | `fmridisplay` | Build a colorbar legend for the object's activation maps |
| `title_montage` | `fmridisplay` | Add a title to a registered montage (by montage index) |
| `transparency_change` | `fmridisplay` | Scale the transparency (AlphaData) of all rendered blobs |
| `zoom_in_on_regions` | `fmridisplay` | Zoom each per-region axis in a montage onto its cluster |
| `activate_figures` | `fmridisplay` | Bring the object's figures (or a montage subset) to the front |
