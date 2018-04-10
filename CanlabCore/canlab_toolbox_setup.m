function canlab_toolbox_setup
% Help:
% In matlab, change to the directory with the CANlab tools repository
% Then run canlab_toolbox_setup
% This assumes you have all CANlab repositories in the same master folder
%
% CanlabCore
% CANlab_help_examples
% FMRI_simulations
% MasksPrivate
% MediationToolbox
% Neuroimaging_Pattern_Masks
% RobustToolbox

main_repository_dir = pwd;

% Assume main dir is 2 levels above this, in master dir
main_repository_dir = fileparts(fileparts(main_repository_dir));

% Find these files, which will tell us where toolbox folders are:
key_files_to_find = {'fmri_data.m' 'a2_second_level_toolbox_check_dependencies.m' ...
    'i_density.m' 'power_figure3_num_comparisons.m' 'apply_nps.m' ...
    'mediation_brain.m' 'apply_all_signatures.m' 'robust_results_batch.m' 'power_calc.m'};

% Find each file
% Get enclosing folder
% Add with subfolders to path

for i = 1:length(key_files_to_find)
    
    
    [status,result] = system(['find ' main_repository_dir ' -name "' key_files_to_find{i} '"']);
    
    if ~status && ~isempty(result)
        
        % Get enclosing folder
        mydir = fileparts(fileparts(result(1, :)));
        
        % Add with subfolders
        g = genpath(mydir);
        addpath(g)
        
        fprintf('Added to path with subfolders: %s\n', mydir);
        
    else
        
        fprintf('Did not find: %s\n', key_files_to_find{i});
        
    end
    
end % for

savepath

%%


