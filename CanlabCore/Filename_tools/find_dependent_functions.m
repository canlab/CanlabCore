function foundFiles = find_dependent_functions(targetFunc, varargin)
% find_dependent_functions - Recursively search folders for files calling a target function.
% - Prints the names of files and returns them in foundFiles
% - Prints line numbers and lines of text in each file
% - Designed by default for CANlab toolbox folders, but can be used for any
% text files and target strings (targetFunc)
%
% :Usage:
% ::
%     foundFiles = find_dependent_functions(targetFunc);
%     foundFiles = find_dependent_functions(targetFunc, folders);
%
% :Inputs:
%
%   **targetFunc:** function name to search for, e.g., apply_parcellation
%        Note: do not include .m in the name in order to fund function
%        calls, which do not include .m.
%        More generally, any text string you want to search for.
%
%   **folders:** (Optional) Cell array of folder paths to search. If not provided or empty,
%        the function uses find_canlab_toolbox_paths to automatically locate CANlab toolbox
%        directories.
%
% :Outputs:
%
%   **foundFiles:** Cell array of full file paths (strings) to MATLAB files (.m) that contain
%        calls to the target function.
%
% :Examples:
% ::
%     % Automatically search CANlab toolbox paths for files that call "apply_atlas":
%     foundFiles = find_dependent_functions('apply_atlas');
%
%     % Specify folders manually:
%     folders = {'/Users/username/Documents/GitHub/CanlabCore'};
%     foundFiles = find_dependent_functions('apply_atlas', folders);
%
% This function recursively searches the specified folders (or the CANlab toolbox paths
% if folders is empty) for MATLAB files that call the target function "apply_atlas".
% For each file that contains the target function call, the function prints the file
% name (without extension), line numbers, and the text of each line where the function
% is called (similar to the Unix grep command).
%
% :References:
%   CANlab Core Tools GitHub Repository: https://github.com/canlab
%
% :See also:
%   dir, fileread, strsplit

% Author: Tor Wager
% Date: March 2025
% License: GNU General Public License v3 or later

    % Optional input: folders
    if nargin < 2 || isempty(varargin{1})
        folders = find_canlab_toolbox_paths;
    else
        folders = varargin{1};
    end

    foundFiles = {};

    % Loop through each folder and search recursively for .m files
    for i = 1:length(folders)

        fileList = dir(fullfile(folders{i}, '**', '*.m'));

        for j = 1:length(fileList)

            fullFileName = fullfile(fileList(j).folder, fileList(j).name);
            if ~isfile(fullFileName)
                continue;
            end

            fileContents = fileread(fullFileName);

            if contains(fileContents, targetFunc)
                foundFiles{end+1, 1} = fullFileName;
            end

        end
    end

    % Display the files that call the target function
    fprintf('Files that call "%s":\n', targetFunc);
    disp(foundFiles);

    % For each found file, print line numbers and text of lines containing targetFunc
    for i = 1:length(foundFiles)
        fname = foundFiles{i};
        [~, fname_no_ext] = fileparts(fname);
        fid = fopen(fname, 'r');
        if fid == -1
            warning('Could not open file: %s', fname);
            continue;
        end
        lineNum = 0;
        while ~feof(fid)
            lineText = fgetl(fid);
            lineNum = lineNum + 1;
            if contains(lineText, targetFunc)
                fprintf('%s:%d: %s\n', fname_no_ext, lineNum, lineText);
            end
        end
        fclose(fid);
    end
end

%% Subfunction: find_canlab_toolbox_paths
function folders = find_canlab_toolbox_paths
% find_canlab_toolbox_paths - Automatically locate CANlab toolbox directories.
%
% :Usage:
% ::
%     folders = find_canlab_toolbox_paths();
%
% :Outputs:
%
%   **folders:** Cell array of folder paths where key CANlab toolbox files are located.
%
% This helper function searches for several key CANlab toolbox files using which()
% and extracts their parent folder (two levels above the file) as candidate paths.
%
% :Examples:
% ::
%     folders = find_canlab_toolbox_paths();
%     disp(folders);
%
% :See also:
%   which, strsplit, fullfile


    key_files_to_find = { 'fmri_data.m', ...
                          'a2_second_level_toolbox_check_dependencies.m', ...
                          'i_density.m', ...
                          'power_figure3_num_comparisons.m', ...
                          'apply_nps.m', ...
                          'mediation_brain.m', ...
                          'apply_all_signatures.m', ...
                          'robust_results_batch.m', ...
                          'power_calc.m' };
    folders = {};
    for i = 1:length(key_files_to_find)
        filePath = which(key_files_to_find{i});
        if ~isempty(filePath)
            parts = strsplit(filePath, filesep);

            % Remove the last two parts to get the parent folder
            targetFolder = fullfile(filesep, parts{1:end-2});
            fprintf('Searching CANlab toolbox folder: %s\n', targetFolder);
            folders{end+1} = targetFolder;
        end
    end

    disp(' ')

end