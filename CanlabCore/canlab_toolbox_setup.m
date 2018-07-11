function canlab_toolbox_setup
% canlab_toolbox_setup
% ------------------------------------------------------------------------
% This function is run with no input arguments. It checks for a number of
% CANlab repositories, and gives you the option to download them from
% Github, add them to the Matlab path with subfolders, and save the Matlab
% path. Run it from the directory ABOVE your individual repositories. e.g.,
% if you have a local directory for your Github repositories, run it from
% there. This function assumes you have all CANlab repositories in the same master
% folder.
%
% The following repositories are included. This is a subset of the
% repositories available, but contains must of the publicly useful ones:
% ------------------------------------------------------------------------
% CANlab core tools - Object-oriented interactive analysis and tools used in other functions
% CANlab help examples - Example scripts and html output, second-level analysis batch scripts
% MKDA meta-analysis toolbox - Coordinate-based meta-analysis of neuroimaging data
% FMRI_simulations - Misc tools for simulating and testing power and image data
% Private masks and signatures - Selected models shared on request
% M3 Mediation toolbox - Single- and Multi-level mediation tools
% Neuroimaging_Pattern_Masks - meta-analysis masks, signatures, parcellations
% Robust regression toolbox - Voxel-wise robust regression for neuroimaging data
% CanlabPrivate - Private repository for lab members and collaborators
%
% July 2018, Tor Wager

main_repository_dir = pwd;

% If you drag and drop the script in, Matlab changes to the folder with the script automatically,
% then changes back.  If so, assume main dir is 2 levels above this.
% If you have it on your path and run it by typing the script name only, it
% will not change directories.  
[dd, ff] = fileparts(main_repository_dir);
if strcmp(ff, 'CanlabCore') && exist(fullfile(pwd, 'canlab_toolbox_setup.m'), 'file')
    % We are in CanlabCore subdirectory
    
    main_repository_dir = fileparts(fileparts(main_repository_dir));
end

dashes = '-----------------------------------------';
fprintf('\n%s\nChecking and/or installing CANlab tools in:\n%s\n%s\n', dashes, main_repository_dir, dashes);


repo_descrip = {'CANlab core tools - Object-oriented interactive analysis and tools used in other functions' ...
    'CANlab help examples - Example scripts and html output, second-level analysis batch scripts' ...
    'MKDA meta-analysis toolbox - Coordinate-based meta-analysis of neuroimaging data' ...
    'FMRI_simulations - Misc tools for simulating and testing power and image data' ...
    'Private masks and signatures - Selected models shared on request' ...
    'M3 Mediation toolbox - Single and Multi-level mediation tools' ...
    'Neuroimaging_Pattern_Masks - meta-analysis masks, signatures, parcellations' ...
    'Robust regression toolbox - Voxel-wise robust regression for neuroimaging data' ...
    'CanlabPrivate - Private repository for lab members and collaborators'};

% Find these files on disk, which will tell us where toolbox folders are:

key_files_to_find = {'fmri_data.m' ...
    'a2_second_level_toolbox_check_dependencies.m' ...
    'i_density.m' ...
    'power_figure3_num_comparisons.m' ...
    'apply_nps.m' ...
    'mediation_brain.m' ...
    'apply_all_signatures.m' ...
    'robust_results_batch.m' ...
    'power_calc.m'};

group_url = 'https://github.com/canlab/';

repo_names = {'CanlabCore' ...
    'CANlab_help_examples' ...
    'Canlab_MKDA_MetaAnalysis' ...
    'FMRI_simulations' ...
    'MasksPrivate' ...
    'MediationToolbox' ...
    'Neuroimaging_Pattern_Masks' ...
    'RobustToolbox' ...
    'CanlabPrivate'};

n_repos = length(key_files_to_find);

% ------------------------------------------------------------------------
% List repos
% ------------------------------------------------------------------------

disp('Looking for CANlab repositories and adding to Matlab path:');
disp(' ');

for i = 1:n_repos
    
    fprintf('%s\t%s\n', repo_descrip{i}, repo_names{i});
    
end

fprintf('\n%s\n', dashes);

% ------------------------------------------------------------------------
% Find each file
% Get enclosing folder
% Add with subfolders to path
% ------------------------------------------------------------------------

was_missing = true(1, n_repos);

for i = 1:n_repos
    
    fprintf('\n%s\n', dashes);

    fprintf('\nRepository: %s\nSite:%s\n', repo_descrip{i}, repo_names{i});
    
    myrepodir = fullfile(main_repository_dir, repo_names{i});
    
    % Future update: 
    % Check if already on path. If so, add with subfolders to be sure all
    % are added. <could use regexp to find base dir and use that>
    % ------------------------------------------------------------------------
    %
    %     checkfile = which([repo_names{i} filesep '@fmri_data' filesep 'fmri_data.m']);
    %     toolboxdir = fileparts(fileparts(checkfile));
    %
    %     if ~exist(toolboxdir, 'dir')
    %
    %         disp('Cannot find CANlabCore toolbox');
    %         disp('You need the CANlab CanlabCore Github repository on your path, with subfolders, to run these scripts.')
    %
    %     else
    %         disp('Found CANlabCore toolbox');
    %         g = genpath(toolboxdir);
    %         addpath(g);
    %     end

    disp('Finding repo: Running command:')
    disp(['find ' myrepodir ' -name "' key_files_to_find{i} '"']);
    
    [status,result] = system(['find ' myrepodir ' -name "' key_files_to_find{i} '"']);
    
    if ~status && ~isempty(result)
        
        was_missing(i) = false;
        
        % Get enclosing folder
        % mydir = fileparts(fileparts(result(1, :)));
        
        % Add with subfolders
        g = genpath(myrepodir);
        addpath(g)
        
        fprintf('Found and added to path with subfolders: %s\n', myrepodir);
        
    else
        
        fprintf('Did not find: %s\n', key_files_to_find{i});
        
    end
    
end % for

fprintf('\n%s\n', dashes);

% ------------------------------------------------------------------------
% Possibly get repositories that are missing
% ------------------------------------------------------------------------

if ~any(was_missing)
    fprintf('\nAll CANlab repositories previously installed. Saving path\n');
    savepath
else
    fprintf('\nMissing repositories:\n');
    disp(char(repo_descrip{was_missing}));
    
    doinstall = input(sprintf('\nTry to download from Github and install in %s?\nPress 1 for yes, 0 for no: ', main_repository_dir));
    
    if doinstall
        
        for i = 1:n_repos
            
            if was_missing(i)
                
                fprintf('\n%s\n', dashes);

                canlab_clone_github_repository('repo', repo_names{i}, 'noverbose');
                
            end
            
        end
        
    end % doinstall
    
    disp('Finished install. Saving path.');
    savepath

end % any missing

% remove .git dirs
strip_git_dirs;



end % function
%%


