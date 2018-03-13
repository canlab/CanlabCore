% Help:
% In matlab, change to the directory with the CANlab tools repository
% Then run this script.
%
% CanlabCore
% CANlab_help_examples
% FMRI_simulations
% MasksPrivate
% MediationToolbox
% Neuroimaging_Pattern_Masks
% RobustToolbox

main_repository_dir = pwd;

% Find these files, which will tell us where toolbox folders are:
key_files_to_find = {'fmri_data.m' 'a2_second_level_toolbox_check_dependencies.m' ...
                     'i_density.m' 'power_figure3_num_comparisons.m' 'apply_nps.m' ...
                     'mediation_brain.m' 'apply_all_signatures.m' 'robust_results_batch.m'};

% Find each file
% Get enclosing folder
% Add with subfolders to path

for i = 1:length(key_files_to_find)
    
    if ~status
        
        [status,result] = system(['find ' main_repository_dir ' -name "' key_files_to_find{i} '"']);
        
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


