# onavg (OpenNeuro Average) registration spheres

Vendored spherical meshes for the **onavg** cortical surface template, used by
`fmri_surface_data/resample_surface` to resample data to/from onavg.

- `onavg_sphere_fsLR_{lh,rh}_41k.mat` — onavg den-41k (40,962 vertices/hemi)
- `onavg_sphere_fsLR_{lh,rh}_10k.mat` — onavg den-10k (10,242 vertices/hemi)

Each `.mat` holds `vertices` (`single`, N×3) and `faces` (`uint32`, M×3).

These are the **`space-fsLR`** onavg registration spheres: onavg vertices expressed
in the fs_LR-aligned spherical frame. That is the same common frame the bundled
fsaverage↔fs_LR "deformed" sphere uses, so onavg resamples to fsaverage / fs_LR
with the same native spherical-interpolation engine (no FreeSurfer / Connectome
Workbench). onavg is an equal-area tessellation (not nested with the FreeSurfer
icosahedra), so onavg resampling always uses interpolation.

## Source and license

Downloaded from TemplateFlow `tpl-onavg`
(https://github.com/templateflow/tpl-onavg), files
`tpl-onavg_space-fsLR_hemi-{L,R}_den-{41k,10k}_sphere.surf.gii`, converted to
`.mat`. License: **CC0** (public domain).

onavg was built from 1,031 participants across 30 OpenNeuro datasets, with vertex
locations optimized so cortical regions are sampled at uniform resolution.

Reference: Feilong M, Guo J (Jiahui), Gobbini MI, Haxby JV (2024). A cortical
surface template with vertex-level, equal-area sampling. *Nature Methods* 21,
2069–2078. https://doi.org/10.1038/s41592-024-02346-y

Additional group statistics (sulcal depth, curvature, vertex area) are available
via GIN (https://gin.g-node.org/neuroboros/core) and are not vendored here.
