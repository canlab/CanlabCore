%% Working with surface and grayordinate data: the fmri_surface_data object
%
% This walkthrough demonstrates the common operations on |fmri_surface_data|,
% the CANlab object for cortical-surface and grayordinate (HCP CIFTI) data. It
% is the surface analogue of |fmri_data|: data are stored flat as a
% [grayordinates x maps] matrix in |.dat|, with a |brain_model| describing the
% surface vertices and subcortical voxels, so the familiar CANlab method names
% (mean, threshold, apply_mask, surface, write, ...) work on surface data.
%
% Everything here runs natively in MATLAB -- no gifti toolbox, FieldTrip,
% Connectome Workbench, or FreeSurfer is required at runtime.
%
% Reference: docs/fmri_surface_data_methods.md
%
% To run: make sure CanlabCore (with Surface_tools/ and @fmri_surface_data/) and
% Neuroimaging_Pattern_Masks are on your path (canlab_toolbox_setup), then step
% through the cells below.

%% Executive summary
%
% A complete tour in a few commands:
%
%     s   = fmri_surface_data(which('transcriptomic_gradients.dscalar.nii')); % load CIFTI
%     surface(s, 'which_image', 1);          % render on the native fs_LR surface
%     r   = reconstruct_image(s);            % dense vertex arrays + subcortex volume
%     vol = to_fmri_data(s);                 % subcortex -> fmri_data (writeable .nii)
%     s2  = vol2surf(fmri_data('some.nii')); % project a volume onto the surface
%     back= surf2vol(s2);                    % ...and back to a volume
%     write(s, 'out.dscalar.nii');           % write CIFTI natively
%
% Now step through it in detail.

%% 1. Load a CIFTI grayordinate file
%
% The constructor auto-detects CIFTI (.dscalar/.dtseries/.dlabel.nii) and GIFTI
% (.surf/.func/.shape/.label.gii) by extension and reads them natively.

ciftifile = which('transcriptomic_gradients.dscalar.nii');   % ships with Neuroimaging_Pattern_Masks
s = fmri_surface_data(ciftifile);

disp(s.intent)                 % 'dscalar'
disp(s.surface_space)          % 'fsLR_32k'
fprintf('%d grayordinates x %d maps\n', size(s.dat,1), size(s.dat,2));

%% 2. Inspect the geometry (brain_model)
%
% brain_model mirrors the CIFTI BrainModels: one entry per cortical hemisphere
% (surface vertices) and per subcortical structure (voxels). It plays the role
% volInfo plays for fmri_data.

cellfun(@(m) sprintf('%-18s %4s  %d rows', m.struct, m.type, m.count), ...
    s.brain_model.models, 'UniformOutput', false)'

% The inherited volInfo describes ONLY the subcortical voxel sub-block:
disp(s.volInfo)

%% 3. Render on the native cortical surface
%
% surface() loads the bundled mesh matching the object's surface_space (here
% fs_LR-32k) and colors vertices directly -- no resampling. The default is a
% 4-panel figure (left/right x lateral/medial). The medial wall renders gray.

surface(s, 'which_image', 1);
drawnow, snapnow;

% Pick a different colormap range or map:
% surface(s, 'which_image', 2, 'clim', [-4 4]);

%% 4. Reconstruct to dense surfaces and a subcortical volume
%
% reconstruct_image returns per-hemisphere dense vertex arrays (medial wall =
% NaN) and a 3-D/4-D subcortical volume. .dat is always full (no replace_empty
% needed for fmri_surface_data).

r = reconstruct_image(s);
fprintf('cortex_left  : %d vertices x %d maps\n', size(r.cortex_left,1), size(r.cortex_left,2));
fprintf('subcortex vol: %s\n', mat2str(size(r.volume)));

%% 5. Export the subcortical part to an fmri_data object (and a .nii)
%
% to_fmri_data returns the volumetric grayordinates as a standard fmri_data in
% MNI space, which you can montage, threshold, or write to NIfTI.

vol = to_fmri_data(s);
fprintf('subcortex fmri_data: %d voxels x %d maps\n', size(vol.dat,1), size(vol.dat,2));
% vol.fullpath = fullfile(pwd, 'subctx.nii'); write(vol);     % uncomment to write

%% 6. Project a volume onto the surface (vol2surf) and back (surf2vol)
%
% vol2surf samples any MNI152 volume (fmri_data/statistic_image) onto the
% fsaverage cortical surface using the vendored CBIG registration-fusion warps.
% surf2vol is its self-consistent inverse.

img = load_image_set('emotionreg');     % 30 contrast images (fmri_data)
m   = mean(img);                        % a single mean volume
ssurf = vol2surf(m);                    % -> fmri_surface_data, fsaverage_164k
surface(ssurf);                         % render the projected map
drawnow, snapnow;

backvol = surf2vol(ssurf);              % fsaverage surface -> MNI fmri_data
fprintf('surf2vol -> fmri_data with %d cortical voxels\n', size(backvol.dat,1));

%% 7. Render on ANY existing surface (e.g. an addbrain MNI surface)
%
% Pass an addbrain keyword (or your own patch handles). If the surface is not
% the object's native mesh, the data is projected through a volume automatically
% (resampling), reusing image_vector.render_on_surface.

create_figure('mni surface');
surface(ssurf, 'mni_surface', 'left');  % render onto an addbrain MNI cortical surface
drawnow, snapnow;

%% 8. Threshold and find contiguous clusters on the mesh
%
% threshold zeros sub-threshold grayordinates (raw-value). reparse_contiguous
% labels contiguous clusters using the cortical mesh graph (cortex) and 26-voxel
% connectivity (subcortex), storing labels in brain_model.cluster.

st = threshold(s, 1.0, 'positive');
[st, ncl] = reparse_contiguous(st, 'which_image', 1);
fprintf('Found %d contiguous clusters above threshold\n', ncl);

%% 8b. Parcellate with a surface atlas
%
% apply_parcellation averages each map within the parcels of a surface atlas (a
% .dlabel object on the same space, or an integer key vector). Background / medial
% wall (key 0) is excluded. surface_region summarizes contiguous clusters.

atlasfile = which('Gordon333.32k_fs_LR_Tian_Subcortex_S2.dlabel.nii');
if ~isempty(atlasfile)
    atl = fmri_surface_data(atlasfile);                  % 91k surface atlas
    if compare_space(s, atl) == 0
        [parcel_means, parcel_labels] = apply_parcellation(s, atl);
        fprintf('Parcellated into %d parcels (%d maps)\n', size(parcel_means,2), size(parcel_means,1));
    end
end

% Clusters of the thresholded map as region-like structs:
reg = surface_region(st, 'which_image', 1);
fprintf('surface_region: %d clusters; sizes: %s\n', numel(reg), mat2str([reg.numVox]));

%% 9. Other analysis-style operations
%
% mean() averages across maps; apply_mask() keeps a subset of grayordinates
% (zeroing the rest -- no shrinking, since grayordinate data is already compact).

mn = mean(s);                                   % single-map mean
keep = s.dat(:,1) > 0;                          % a simple grayordinate mask
smasked = apply_mask(s, keep);                  % zero out non-positive grayordinates

%% 9b. Group analysis: cat, ttest, regress, predict
%
% Combine per-subject surface maps with cat (or [a, b, ...]) into a single
% object, then run the same analyses as fmri_data. Here we simulate a few
% "subjects" by projecting individual emotionreg contrast images to the surface.

subj = cell(1, 5);
for i = 1:5
    subj{i} = vol2surf(get_wh_image(img, i));   % each subject's contrast on the surface
end
group = cat(subj{:});                            % one object, 5 maps
fprintf('group object: %d grayordinates x %d maps\n', size(group.dat,1), size(group.dat,2));

tmap = ttest(group);                             % grayordinate-wise t-test
surface(tmap, 'clim', [-4 4]);
drawnow, snapnow;

% Cross-validated prediction of an outcome (here a synthetic one):
group.Y = (1:5)';
[cverr, stats] = predict(group, 'algorithm_name', 'cv_lassopcr', 'nfolds', 5, 'verbose', 0);
% surface(stats.weight_obj);   % render the predictive weight map

% OLS regression onto a design matrix:
group.X = [ones(5,1), (1:5)'];
b = regress(group);             % betas in b.dat; t/p in b.additional_info.statistic

%% 10. Write CIFTI / GIFTI natively
%
% write() dispatches on the output extension. CIFTI map metadata (scalar names,
% label tables, series info) and the subcortical affine are preserved.

outdir = tempdir;
write(s,     fullfile(outdir, 'gradients_copy.dscalar.nii'));   % CIFTI-2
write(ssurf, fullfile(outdir, 'emotionreg_surf.func.gii'));     % GIFTI (cortex)
fprintf('Wrote CIFTI and GIFTI to %s\n', outdir);

% Round-trip check:
s_reloaded = fmri_surface_data(fullfile(outdir, 'gradients_copy.dscalar.nii'));
fprintf('Round-trip max abs diff: %g\n', max(abs(double(s.dat(:)) - double(s_reloaded.dat(:)))));

%% 11. Quick QC plot
%
% plot() shows a value histogram, per-map mean/sd, coverage, and a mean-map
% surface render.

plot(s);
drawnow, snapnow;

%% See also
%
% docs/fmri_surface_data_methods.md      - full method/option reference
% docs/fmri_surface_data_design_plan.md  - design rationale and roadmap
% methods(fmri_surface_data)             - the full method list
