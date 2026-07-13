function p = canlab_cbig_warp_path(name)
% canlab_cbig_warp_path Resolve a vendored CBIG registration-fusion warp file path.
%
% :Usage:
% ::
%     p = canlab_cbig_warp_path('lh_ras')
%
% Returns the absolute path to a CBIG Registration-Fusion (RF-ANTs) MNI152<->
% fsaverage mapping file vendored under CanlabCore/Cifti_plotting/. These warps
% (MIT-licensed) back the native vol2surf / surf2vol mappers -- no FreeSurfer or
% Connectome Workbench binary is required.
%
% :Inputs:
%   **name:** one of
%       'lh_ras' / 'rh_ras' : per-fsaverage-vertex MNI152 RAS coords (3x163842),
%                             variable 'ras' (MNI152 orig -> fsaverage).
%
% :Outputs:
%   **p:** absolute path to the .mat warp file (errors if not found).
%
% :See also: vol2surf, surf2vol, fmri_surface_data

% CanlabCore root = parent of Surface_tools (this file's folder)
root = fileparts(fileparts(mfilename('fullpath')));
base = fullfile(root, 'Cifti_plotting', 'CBIG_registration_fusion_surf2vol_vol2surf', ...
    'standalone_scripts_for_MNI_fsaverage_projection', 'final_warps_FS5.3');

switch lower(name)
    case 'lh_ras'
        p = fullfile(base, 'lh.avgMapping_allSub_RF_ANTs_MNI152_orig_to_fsaverage.mat');
    case 'rh_ras'
        p = fullfile(base, 'rh.avgMapping_allSub_RF_ANTs_MNI152_orig_to_fsaverage.mat');
    otherwise
        error('canlab_cbig_warp_path:badname', 'Unknown warp name: %s', name);
end

if exist(p, 'file') ~= 2
    error('canlab_cbig_warp_path:notfound', ...
        ['CBIG warp not found: %s\nThe registration-fusion warps should ship under ' ...
         'CanlabCore/Cifti_plotting/CBIG_registration_fusion_surf2vol_vol2surf/.'], p);
end
end
