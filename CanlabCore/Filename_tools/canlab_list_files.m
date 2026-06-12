function [names, isfile, isdir, hasduplicates] = canlab_list_files(varargin)
% canlab_list_files Build a list of file paths by concatenating subfolder names with optional cell expansion.
%
% :Usage:
% ::
%
%     [names, isfile, isdir, hasduplicates] = canlab_list_files(dir1, dir2, dir3, ...)
%
% Lists files in a series of subfolders, in a particular order. Enter
% as many dirs as you want; they are treated as nested subfolders and
% concatenated in order with fullfile.
%
% If any dir argument is a cell array of names, this function will
% expand it: it will produce one output entry for every combination of
% that cell with the other dir strings (in order). This is convenient
% for building per-subject file lists.
%
% :Inputs:
%
%   **dir1, dir2, ...:**
%        Each argument is either a string (treated as a single
%        subfolder/filename appended to all current paths) or a cell
%        array of strings (treated as alternatives, each one expanding
%        the current list of paths). Strings are concatenated in order
%        with fullfile.
%
% :Outputs:
%
%   **names:**
%        Cell array of constructed paths.
%
%   **isfile:**
%        Numeric/logical column vector; whether each entry in names is
%        an existing file (per exist(name, 'file')).
%
%   **isdir:**
%        Numeric/logical column vector; whether each entry in names is
%        an existing directory (per exist(name, 'dir')).
%
%   **hasduplicates:**
%        Logical scalar; true if any entries in names are
%        duplicates of one another.
%
% :Examples:
% ::
%
%     subjects = {'nsf1' 'nsf2' 'NSF909'};
%     names = canlab_list_files(subjects, 'Structural', 'SPGR', 'wspgr.img');
%     names = canlab_list_files(pwd, subjects, 'Structural', 'SPGR', 'wspgr.img');
%
%     % Return list of whether each is an existing file or dir
%     [names, isfile, isdir, hasduplicates] = ...
%         canlab_list_files(pwd, subjects, 'Structural', 'SPGR', 'wspgr.img');
%
%     % Convert the resulting cell array to a char string matrix
%     names = char(names{:});
%
% :See also:
%   - filenames
%   - dir
%   - fullfile

d = varargin;

names = {};

for i = 1:length(d)
    
    n = length(names);
    
    if iscell(d{i})
        
        if isempty(names)
            names = d{i};
            if diff(size(names)) % is row
                names = names';
            end
            continue
        end
        
        for j = 1:n  % for each cell in original names...
            
            namebase = names{j};
            
            for k = 1:length(d{i})
                % Expand: For every cell in original names, append a list of
                % names for each cell in d{i}
                
                newnames{k, 1} = fullfile(namebase, d{i}{k});
            end
            
            if j == 1
                names = newnames;
            else
                names = [names; newnames];
            end
            
        end
        
    else
        % d{i} is a string, not a list; simply append it to each dir in
        % names
        
        if isempty(names)
            names = d(i);
            continue
        end
        
        for k = 1:length(names)
            names{k, 1} = fullfile(names{k}, d{i});
        end
        
    end % end what to do with d{i}
    
    
end % for

for i = 1:length(names)
    
    names{i} = strtrim(names{i});
    isfile(i, 1) = exist(names{i}, 'file');
    isdir(i, 1) = exist(names{i}, 'dir');
    
end

nunique = size(unique(char(names{:}), 'rows'), 1);

hasduplicates = nunique < length(names);

end % function

