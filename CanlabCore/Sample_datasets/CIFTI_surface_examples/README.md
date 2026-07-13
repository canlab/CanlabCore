# CIFTI surface / grayordinate sample data

Example CIFTI-2 surface data used by the `fmri_surface_data` object documentation
and walkthrough (see `docs/workflows/fmri_surface_data_howto.md` and
`docs/fmri_surface_data_methods.md`). Loads with:

```matlab
s = fmri_surface_data(which('S1200.MyelinMap_MSMAll.32k_fs_LR.dscalar.nii'));
surface(s);
```

| File | Type | Contents |
|---|---|---|
| `S1200.MyelinMap_MSMAll.32k_fs_LR.dscalar.nii` | `.dscalar` | HCP S1200 group-average cortical myelin map (T1w/T2w), 59,412 cortical grayordinates (fs_LR-32k, medial wall excluded). A continuous map — good for surface visualization. |

## Additional surface data (Neuroimaging_Pattern_Masks)

CanlabCore ships only this one small continuous map. **Surface atlases,
grayordinate parcellations, and other CIFTI maps used by parts of the walkthrough
live in the companion `Neuroimaging_Pattern_Masks` repository** (a standard CANlab
dependency; add it with `canlab_toolbox_setup`). Once it is on the path you can
load e.g.:

```matlab
atl = fmri_surface_data(which('Gordon333.32k_fs_LR_Tian_Subcortex_S2.dlabel.nii'));
grd = fmri_surface_data(which('transcriptomic_gradients.dscalar.nii'));
```

See the inventory table in `docs/fmri_surface_data_methods.md` (section
"Surface data in Neuroimaging_Pattern_Masks").

## Provenance and license

- **HCP S1200 group myelin map** — from the Human Connectome Project S1200 release
  (`Q1-Q6_RelatedParcellation210.MyelinMap_BC_MSMAll_2_d41_WRN_DeDrift.32k_fs_LR.dscalar.nii`),
  distributed under the [HCP Open Access Data Use Terms](https://www.humanconnectome.org/study/hcp-young-adult/document/wu-minn-hcp-consortium-open-access-data-use-terms)
  (redistributable for research/education with acknowledgment; WU-Minn HCP
  Consortium, funded by the NIH Blueprint for Neuroscience Research).

This is a group-average template file (no individual-subject data). If you
redistribute CanlabCore, retain this README with its attribution.
