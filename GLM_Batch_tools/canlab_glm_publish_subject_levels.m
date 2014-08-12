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
