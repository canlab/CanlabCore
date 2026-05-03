% canlab_glm_publish_subject_levels Publish subject-level SPM design reviews via scn_spm_design_check.
%
% :Usage:
% ::
%
%     % Run from a directory containing subject-level SPM analyses
%     canlab_glm_publish_subject_levels
%
% Child script of canlab_glm_publish. Runs scn_spm_design_check in a
% publishable manner (intended for use with MATLAB's publish() so the
% review can be rendered as HTML/PDF).
%
% Behavior:
%
%   - If the current working directory itself contains an SPM.mat file,
%     it is treated as the single subject-level analysis to review.
%   - Otherwise the script searches one directory level deep for any
%     subdirectories that contain an SPM.mat file and reviews each.
%   - For every located analysis, scn_spm_design_check(..., 'events_only')
%     is called and snapnow is issued so that figures appear in the
%     published output.
%
% A diary log of the review is written to
% model_review_work_log.txt in the current working directory.
%
% :Inputs:
%
%   None. This file is a script and operates on the current working
%   directory.
%
% :Outputs:
%
%   No MATLAB outputs. Side effects:
%
%   - model_review_work_log.txt log file in pwd.
%   - Figures from scn_spm_design_check captured by snapnow when the
%     script is run via MATLAB's publish().
%
% :See also:
%   - canlab_glm_publish_group_levels
%   - canlab_glm_subject_levels
%   - scn_spm_design_check
%   - publish

% child script of canlab_glm_publish
% runs scn_spm_design_check in publishable manner
% runs in directory containing subject level SPM analyses


z = '_________________________________________________';
diaryname = fullfile('model_review_work_log.txt');


if exist(fullfile(pwd,'SPM.mat'),'file')
    sublevs{1} = pwd;
else
    % figure out which directories are lower level analyses
    % criterion: they contain an SPM.mat file
    spmfiles = filenames(fullfile('*','SPM.mat'),'absolute');
    sublevs = cellfun(@fileparts,spmfiles,'UniformOutput',false);
end

for i = 1:numel(sublevs)
    diary(diaryname), fprintf('%s\n%s\n%s\n', z, sublevs{i}, z), diary off
    
    try
        scn_spm_design_check(sublevs{i}, 'events_only');
        snapnow
    catch exc
        diary(diaryname)
        disp(getReport(exc,'extended'))
%         fprintf('Design review FAILED for: %s.\n', sublevs{i})
        diary off
    end    
end

close all
