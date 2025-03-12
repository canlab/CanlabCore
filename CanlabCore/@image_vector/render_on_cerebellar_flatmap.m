function giftiname = render_on_cerebellar_flatmap(target_obj, varargin)
% Render an object with a single image onto a standard cerebellar flat map.
%
% - requires SPM, the SUIT toolbox, and some standard templates on your Matlab path.
% - saves a Gifti .gii image of your map, called <image_name>.cblm_surf.gii
%
% Diedrichsen, J. & Zotow, E. (2015). Surface-based display of volume-averaged cerebellar data. PLoS One, 7, e0133402.
% https://www.diedrichsenlab.org/imaging/suit_flatmap.htm
%
% :Inputs:
%
% :Optional Inputs:
%   **'target_obj':**
%        A CANlab fmri_data object you want to visualize
%
%   **'newfigure':** [numeric scalar]
%        Create a new figure; default = uses existing axes
%
%   **'color_map':** [matrix]
%        An n-by-3 colormap matrix used for rendering.
%        Default = colormap_tor([0 0 1], [1 1 0], [0.5 0.5 0.5], [0 0.5 1], [1 0.5 0]).
%
% Some examples:
% 
% Pain-related activation from 2021 Spisak Placebo meta N = 603 maps
% fname = which('full_pain_g_pperm_FWE05.nii.gz');
% pain603 = fmri_data(fname);
% render_on_cerebellar_flatmap(pain603)
%
% target_obj = load_image_set('nps');
% render_on_cerebellar_flatmap(target_obj)
%
% transgrad = cifti_read('transcriptomic_gradients.dscalar.nii');
% r = cifti_struct_2_region_obj(cifti_struct, 'which_image', 3); % 3rd gradient
% subctx_fmri_data_obj = region2fmri_data(r);
% render_on_cerebellar_flatmap(subctx_fmri_data_obj, 'color_map', colormap('summer'))
% render_on_cerebellar_flatmap(subctx_fmri_data_obj)


% Some other examples:
% fname = which('full_pla_g_pperm_tfce_FWE05.nii.gz');
% placebo603 = fmri_data(fname);
% 
% fname = which('full_pla_rrating_pperm_tfce_FWE05.nii.gz');
% placebocorr603 = fmri_data(fname);
% 
% gunzip(pain603.fullpath);
% pain603fname = strrep(pain603.fullpath, '.gz', '');
% 
% gunzip(placebo603.fullpath);
% placebo603fname = strrep(placebo603.fullpath, '.gz', '');
% 
% gunzip(placebocorr603.fullpath);
% placebocorr603fname = strrep(placebocorr603.fullpath, '.gz', '');

%%
% -------------------------------------------------------------------------
% DEFAULTS AND INPUTS
% -------------------------------------------------------------------------

% Parse variable inputs using inputParser
ARGS = parse_inputs(varargin{:});

fn = fieldnames(ARGS);
for i = 1:length(fn)
    eval([fn{i}, ' = ARGS.(fn{i});']);
end

% -------------------------------------------------------------------------
% Turn target object into a temporary .nii file
% -------------------------------------------------------------------------

if size(target_obj.dat, 2) > 1
    error('Target object should contain only a single image')
end

target_obj.fullpath = fullfile(pwd, 'tmp_target_nii.nii');
write(target_obj, 'overwrite');

targetniifile = fullfile(pwd, 'tmp_target_nii.nii');
fprintf('created %s\n', targetniifile)

% Check for SUIT toolbox
% ------------------------------
suitname = which('suit_reslice_dartel.m');
if isempty(suitname)
    disp('Can''t find suit_reslice_dartel.m')
    error('Add SUIT folder to path')
end

% Check for SPM toolbox
% ------------------------------
if isempty(which('spm_dartel_norm.m')), error('Start SPM to find spm_dartel_norm.m'), end

% Find SUIT in correct place and add to path
% Note: SUIT must be in spm12/toolbox/suit/flatmap' for suit_map2surf to work
% ------------------------------
spmtoolboxdir = fileparts(fileparts(which('spm_dartel_norm.m')));
suitdir = fullfile(spmtoolboxdir, 'suit');
if ~isdir(suitdir), error('Add suit toolbox to spm12/toolbox/suit'), end
g = genpath(suitdir); addpath(g);

% Check for target .nii file
% ------------------------------
targetniifile = which(targetniifile);       % get full path
if isempty(targetniifile)
    disp('Can''t find target .nii file on path')
    error('Check path and file names.')
end

[dd, ff, ee] = fileparts(targetniifile);

% Set output names
% ------------------------------

wdtargetniifile = fullfile(dd, ['wd' ff ee]);

giftiname = fullfile(dd, [ff '.cblm_surf.gii']);

% Find templates for flat map transformation
% -------------------------------------------------------------------------

affinename = which('Affine_MNI152NLin6Asym_T1_1mm_seg1.mat');
if isempty(affinename)
    disp('Can''t find Affine_MNI152NLin6Asym_T1_1mm_seg1.mat')
    error('Add Canlab Neuroimaging_Pattern_Masks github subfolders to path')
end

flowname = which('u_a_MNI152NLin6Asym_T1_1mm_seg1.nii');
if isempty(flowname)
    disp('Can''t find u_a_MNI152NLin6Asym_T1_1mm_seg1.nii')
    error('Add Canlab Neuroimaging_Pattern_Masks github subfolders to path')
end

cerebmaskfile = which('c_MNI152NLin6Asym_T1_1mm_pcereb.nii');
if isempty(cerebmaskfile)
    disp('Can''t find c_MNI152NLin6Asym_T1_1mm_pcereb.nii')
    error('Add Canlab Neuroimaging_Pattern_Masks github subfolders to path')
end

% Some code to create templates above
% ------------------------------
% Normalize template to SUIT cerebellar template - define mapping to SUIT space
% suit_normalize_dartel(struct('subjND', struct('gray', {{'MNI152NLin6Asym_T1_1mm_seg1.nii'}}, 'white', {{'MNI152NLin6Asym_T1_1mm_seg2.nii'}}, 'isolation', {{'c_MNI152NLin6Asym_T1_1mm_pcereb.nii'}})))

% Apply normalization : Reslice the target image, already in MNI space, in SUIT space

% See also:
% apply_spm_warp(mvg_img0, fxd_img0, pre_affine_mat, warp_img, post_affine_mat, out_img, interp)

% Reslice the target image
% ------------------------------

suit_input_struct = struct('subj', struct('affineTr', {{affinename}}, ...
    'flowfield', {{flowname}}, ...
    'resample', {{targetniifile}}, ...
    'mask', {{cerebmaskfile}}));

% creates wdtargetniifile
suit_reslice_dartel(suit_input_struct)
fprintf('created %s\n', wdtargetniifile)

% Map to cerebellar surface
% ------------------------------

C = [];
C.cdata = suit_map2surf(wdtargetniifile, 'space','SUIT', 'stats', @nanmean);

% Save GIFTI gii file, if requested
% ------------------------------
if ~isempty(giftiname)

    save(gifti(C), giftiname)
    disp('Writing GIFTI .gii image:')
    fprintf('created %s\n', giftiname)
    
end


% Plot the map
% ------------------------------

if newfigure
create_figure('Cerebellar flatmap')
end

% cmap = colormap_tor([.5 0 1], [1 1 0]);

mymax = prctile(abs(C.cdata), 95);
suit_plotflatmap(C.cdata, 'cmap', color_map, 'cscale',[-mymax mymax]);

% f = gifti(giiname)
% mymax = prctile(abs(f.cdata), 95);
% suit_plotflatmap(f.cdata, 'cmap', cmap, 'cscale',[-mymax mymax]);

% Clean up
% ------------------------------

delete(targetniifile)
delete(wdtargetniifile)


end % main function



% ------------------------------------------------------------------------
% Subfunction: parse_inputs
% ------------------------------------------------------------------------
function ARGS = parse_inputs(varargin)
% parse_inputs parses optional input arguments.
%
% :Usage:
% ::
%     ARGS = parse_inputs(optional_name_value_pairs)
%
% :Optional Inputs:
%
%   **'newfigure':** [numeric scalar]
%        Create a new figure; default = uses existing axes
%
%   **'color_map':** [matrix]
%        An n-by-3 colormap matrix used for rendering.
%        Default = colormap_tor([0 0 1], [1 1 0], [0.5 0.5 0.5], [0 0.5 1], [1 0.5 0]).
%
% ------------------------------------------------------------------------
p = inputParser;

addParameter(p, 'newfigure', false, @(x) islogical(x) && isscalar(x));

defaultColorMap = colormap_tor([0 0 1], [1 1 0], [0.5 0.5 0.5], [0 0.5 1], [1 0.5 0]);
addParameter(p, 'color_map', defaultColorMap, @(x) isnumeric(x) && size(x,2)==3);

addParameter(p, 'verbose', false, @(x) islogical(x) && isscalar(x));

parse(p, varargin{:});
ARGS = p.Results;

end