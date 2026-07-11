function [hot, warm, ok] = canlab_get_dpsp_hot_warm()
%CANLAB_GET_DPSP_HOT_WARM Load the DPSP single-subject Hot/Warm sample maps.
%
%   [hot, warm, ok] = canlab_get_dpsp_hot_warm()
%
%   Returns the single-subject Hot and Warm condition maps from
%   Sample_datasets/DPSP_pain_rejection_participant_maps as fmri_data objects.
%   Used by the predictive_model / xval_* tests so the load semantics live in
%   one place. ok is false (and hot/warm are []) if the sample files are not
%   on the path, so callers can assume/skip gracefully.

sample_dir = fullfile(fileparts(fileparts(which('fmri_data'))), ...
    'Sample_datasets', 'DPSP_pain_rejection_participant_maps');
hot_file  = fullfile(sample_dir, 'DPSP_single_subject_images_hot.mat');
warm_file = fullfile(sample_dir, 'DPSP_single_subject_images_warm.mat');

ok = exist(hot_file, 'file') == 2 && exist(warm_file, 'file') == 2;
hot = [];
warm = [];

if ok
    H = load(hot_file);
    W = load(warm_file);
    hot  = H.single_subject_images_hot;
    warm = W.single_subject_images_warm;
end
end
