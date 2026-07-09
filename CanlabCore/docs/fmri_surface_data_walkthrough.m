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
% To run: you only need CanlabCore on your path (canlab_toolbox_setup). The core
% sections use a small CIFTI map that ships with CanlabCore and the built-in
% 'emotionreg' dataset. A few sections that show subcortical grayordinates and
% atlas parcellation use files from the companion Neuroimaging_Pattern_Masks
% repository; those sections are guarded with `if ~isempty(which(...))` so the
% script still runs top to bottom without it.

%% Executive summary
%
% A complete tour in a few commands:
%
%     s     = fmri_surface_data(which('S1200.MyelinMap_MSMAll.32k_fs_LR.dscalar.nii'));
%     surface(s);                            % render on the native fs_LR surface
%     r     = reconstruct_image(s);          % dense per-hemisphere vertex arrays
%     ssurf = vol2surf(ttest(load_image_set('emotionreg')));  % volume result -> surface
%     back  = surf2vol(ssurf);               % ...and back to an MNI fmri_data volume
%     write(s, 'out.dscalar.nii');           % write CIFTI natively
%
% Now step through it in detail.

%% 1. Load a CIFTI surface file
%
% The constructor auto-detects CIFTI (.dscalar/.dtseries/.dlabel.nii) and GIFTI
% (.surf/.func/.shape/.label.gii) by extension and reads them natively. Here we
% load the HCP group cortical myelin map bundled with CanlabCore.

ciftifile = which('S1200.MyelinMap_MSMAll.32k_fs_LR.dscalar.nii');   % ships with CanlabCore
s = fmri_surface_data(ciftifile);

disp(s.imagetype)                 % 'dscalar'
disp(s.surface_space)          % 'fsLR_32k'
fprintf('%d grayordinates x %d maps\n', size(s.dat,1), size(s.dat,2));

%% 2. Render on the native cortical surface
%
% surface() loads the bundled mesh matching the object's surface_space (here
% fs_LR-32k) and colors vertices directly -- no resampling. The default is a
% 4-panel figure (left/right x lateral/medial). The medial wall renders gray.

surface(s, 'pos_colormap', hot(256));
drawnow, snapnow;

% Options: pick a map ('which_image'), color range ('clim'), or mesh ('surftype',
% e.g. 'midthickness' for fs_LR). For example:
% surface(s, 'clim', [1 1.8], 'surftype', 'midthickness');

%% 3. The grayordinate model (brain_model)
%
% brain_model mirrors the CIFTI BrainModels: one entry per cortical hemisphere
% (surface vertices) and per subcortical structure (voxels). It plays the role
% volInfo plays for fmri_data. The myelin map above is cortex-only; a full "91k"
% grayordinate file (from Neuroimaging_Pattern_Masks) also has subcortical
% volume structures.

fprintf('myelin map models:\n');
cellfun(@(m) sprintf('  %-14s %-4s %6d rows', m.struct, m.type, m.count), ...
    s.brain_model.models, 'UniformOutput', false)'

sg = [];   % a full grayordinate object (cortex + subcortex), if available
gfile = which('Gordon333.32k_fs_LR_Tian_Subcortex_S2.dlabel.nii');   % from NPM
if ~isempty(gfile)
    sg = fmri_surface_data(gfile);
    types = cellfun(@(m) m.type, sg.brain_model.models, 'UniformOutput', false);
    fprintf('\n91k grayordinate atlas: %d grayordinates, %d models (%d surface + %d volume)\n', ...
        size(sg.dat,1), numel(sg.brain_model.models), nnz(strcmp(types,'surf')), nnz(strcmp(types,'vox')));
    % The inherited volInfo describes ONLY the subcortical voxel sub-block:
    disp(sg.volInfo)
else
    fprintf('\n(Install Neuroimaging_Pattern_Masks to see subcortical grayordinates.)\n');
end

%% 4. Reconstruct to dense surfaces (and a subcortical volume)
%
% reconstruct_image returns per-hemisphere dense vertex arrays (medial wall =
% NaN); for a grayordinate object with subcortex it also returns a 3-D/4-D
% volume. .dat is always full (no replace_empty needed for fmri_surface_data).

r = reconstruct_image(s);
fprintf('cortex_left : %d vertices x %d maps (medial wall = NaN)\n', ...
    size(r.cortex_left,1), size(r.cortex_left,2));

if ~isempty(sg)
    rg = reconstruct_image(sg);
    fprintf('subcortex volume: %s\n', mat2str(size(rg.volume)));
end

%% 5. Export the subcortical part to an fmri_data object (and a .nii)
%
% to_fmri_data returns the volumetric (subcortical) grayordinates as a standard
% fmri_data in MNI space, which you can montage, threshold, or write to NIfTI.
% (This needs a grayordinate object with subcortex.)

if ~isempty(sg)
    vol = to_fmri_data(sg);
    fprintf('subcortex fmri_data: %d voxels x %d maps\n', size(vol.dat,1), size(vol.dat,2));
    % vol.fullpath = fullfile(pwd, 'subctx.nii'); write(vol);   % uncomment to write NIfTI
end

%% 6. Project a volume onto the surface (vol2surf) and back (surf2vol)
%
% vol2surf samples any MNI152 volume (fmri_data/statistic_image) onto the
% fsaverage cortical surface using the vendored CBIG registration-fusion warps.
% surf2vol is its self-consistent inverse.

img   = load_image_set('emotionreg');   % 30 contrast images (fmri_data)
t     = ttest(img);                      % volumetric statistic_image (group t-map)
ssurf = vol2surf(t);                     % -> fmri_surface_data, fsaverage_164k
surface(ssurf, 'clim', [-6 6]);          % render the projected map
drawnow, snapnow;

backvol = surf2vol(ssurf);               % fsaverage surface -> MNI fmri_data
fprintf('surf2vol -> fmri_data with %d cortical voxels\n', size(backvol.dat,1));

%% 7. Render on ANY existing surface (e.g. an addbrain MNI surface)
%
% Pass an addbrain keyword (or your own patch handles). If the surface is not
% the object's native mesh, the data is projected through a volume automatically
% (resampling), reusing image_vector.render_on_surface.

surface(ssurf, 'mni_surface', 'left');   % render onto an addbrain MNI cortical surface
drawnow, snapnow;

%% 8. Threshold and find contiguous clusters on the mesh
%
% threshold zeros sub-threshold grayordinates (raw value; add 'k', N for a
% cluster-extent threshold). reparse_contiguous labels contiguous clusters using
% the cortical mesh graph (cortex) and 26-voxel connectivity (subcortex).

st = threshold(ssurf, 3, 'positive', 'k', 20);   % t > 3, clusters >= 20 grayordinates
[st, ncl] = reparse_contiguous(st, 'which_image', 1);
fprintf('Found %d contiguous clusters above threshold\n', ncl);

% Summarize the clusters as region-like structs (centroid, size, mean value):
reg = surface_region(st, 'which_image', 1);
if ~isempty(reg)
    fprintf('surface_region: %d clusters; sizes: %s\n', numel(reg), mat2str([reg.numVox]));
end

%% 9. Parcellate with a surface atlas
%
% apply_parcellation averages each map within the parcels of a surface atlas (a
% .dlabel object on the same grayordinate space, or an integer key vector).
% Background / medial wall (key 0) is excluded. The data and atlas must be on the
% same grayordinate space; here we make a demo map on the atlas's own space.

if ~isempty(sg)
    demo = sg;                                  % same grayordinate space as the atlas
    demo.dat = single(sqrt(double(sg.dat)));    % a continuous map on those grayordinates
    demo.imagetype = 'dscalar';
    [parcel_means, parcel_labels] = apply_parcellation(demo, sg);
    fprintf('Parcellated into %d parcels; first parcel "%s" mean = %.3f\n', ...
        size(parcel_means,2), parcel_labels{1}, parcel_means(1,1));
end

%% 10. Other data operations: mean, apply_mask
%
% mean() averages across maps; apply_mask() keeps a subset of grayordinates
% (zeroing the rest -- no shrinking, since grayordinate data is already compact).

mn      = mean(ssurf);                    % single-map mean across maps
keep    = ssurf.dat(:,1) > 0;             % a simple grayordinate mask (positive t)
smasked = apply_mask(ssurf, keep);        % zero the non-positive grayordinates

%% 11. Group analysis: cat, ttest, regress, predict
%
% Combine per-map/per-subject surface objects with cat (or [a, b, ...]) into one
% object, then run the same analyses as fmri_data. Here we build a small "group"
% by projecting individual emotionreg contrast images to the surface.

subj = cell(1, 5);
for i = 1:5
    subj{i} = vol2surf(get_wh_image(img, i));   % each subject's contrast on the surface
end
group = cat(subj{:});                            % one object, 5 maps
fprintf('group object: %d grayordinates x %d maps\n', size(group.dat,1), size(group.dat,2));

tmap = ttest(group);                             % grayordinate-wise one-sample t-test
% surface(tmap, 'clim', [-4 4]);                 % (render the group t-map)

% OLS regression onto a design matrix (set X BEFORE calling regress):
group.X = [ones(5,1), (1:5)'];                   % intercept + a linear predictor
b = regress(group);                              % betas in b.dat; t/p in additional_info.statistic

% Cross-validated multivariate prediction (set Y first):
group.Y = (1:5)';
[cverr, stats] = predict(group, 'algorithm_name', 'cv_lassopcr', 'nfolds', 5, 'verbose', 0);
fprintf('predict cverr = %.3f; weight map is an fmri_surface_data (%dx%d)\n', ...
    cverr, size(stats.weight_obj.dat,1), size(stats.weight_obj.dat,2));

%% 12. Write CIFTI / GIFTI natively
%
% write() dispatches on the output extension. CIFTI map metadata (scalar names,
% label tables, series info) and the subcortical affine are preserved.

outdir = tempdir;
write(s,     fullfile(outdir, 'myelin_copy.dscalar.nii'));   % CIFTI-2
write(ssurf, fullfile(outdir, 'emotionreg_surf.func.gii'));  % GIFTI (cortex)
fprintf('Wrote CIFTI and GIFTI to %s\n', outdir);

% Round-trip check:
s_reloaded = fmri_surface_data(fullfile(outdir, 'myelin_copy.dscalar.nii'));
fprintf('Round-trip max abs diff: %g\n', max(abs(double(s.dat(:)) - double(s_reloaded.dat(:)))));

%% 13. Quick QC plot
%
% plot() shows a value histogram, per-map mean/sd, coverage, and a mean-map
% surface render.

plot(s);
drawnow, snapnow;

%% See also
%
% docs/fmri_surface_data_methods.md      - full method/option reference
% docs/workflows/fmri_surface_data_howto.md - the workflow how-to (this content as markdown)
% docs/fmri_surface_data_design_plan.md  - design rationale and roadmap
% methods(fmri_surface_data)             - the full method list
