# `fmri_mask_image` methods, organized by area

`fmri_mask_image` is a legacy subclass of `image_vector` that holds a
single binary or weighted mask image plus the bookkeeping needed to
track which space it lives in and whether it has been resampled to match
another image. It inherits all of `image_vector`'s storage (`dat`,
`volInfo`, `removed_voxels`, `fullpath`, `history`) and methods.

Most modern CANlab code accepts an `fmri_data` or `image_vector` object
wherever a mask is needed (legacy; many newer methods accept `fmri_data`
or `image_vector` as a mask), so explicit `fmri_mask_image` objects are
mainly used inside class constructors (e.g. `fmri_data`,
`fmri_timeseries`) and a few resampling pipelines. Type
`methods(my_obj)` in MATLAB for the live list on any instance.

## Properties

Inherits all `image_vector` properties (`dat`, `volInfo`,
`removed_voxels`, `fullpath`, `image_names`, `files_exist`, `history`,
etc.). Adds:

| Property | Description |
|---|---|
| `volInfo_descrip` | Text description of the volume / in-mask info |
| `space_defining_image_name` | File name of the image whose space the mask is currently in |

## Methods specific to `@fmri_mask_image`

| Method | From | One-liner |
|---|---|---|
| `fmri_mask_image` | `@fmri_mask_image` | Constructor: load a mask file (or copy from an `image_vector`); `'implicit'` derives an implicit mask |
| `resample_to_image_space` | `@fmri_mask_image` | Reslice the mask onto another image's space using `scn_map_image` (re-reads from disk) |

For everything else (display, masking other objects, write, etc.), see
the methods inherited from `image_vector` documented in
[`fmri_data` methods](fmri_data_methods.md).
