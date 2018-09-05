function [names, isfile, isdir, hasduplicates] = canlab_list_files(varargin)
% [names, isfile, isdir, hasduplicates]  = canlab_list_files(dir1, dir2, dir3, etc.)
%
% Lists files in a series of subfolders in a particular order.
% Enter as many dirs as you want, which are subfolders.
% If a dir is a cell array of names, this function will look in all combos
% of those cells with other dir strings, in order.
% dir1...n strings are concatenated in order.
%
% Returns:
% Vector of whether each name in names is an existing file, and whether
% each is an existing directory
%
% e.g.,
% subjects = {'nsf1' 'nsf2' 'NSF909'};
% names = canlab_list_files(subjects, 'Structural', 'SPGR', 'wspgr.img')
% names = canlab_list_files(pwd, subjects, 'Structural', 'SPGR', 'wspgr.img')
% 
% Return list of whether each is an existing file or dir
% [names, isfile, isdir, hasduplicates] = canlab_list_files(pwd, subjects, 'Structural', 'SPGR', 'wspgr.img')
%
% Make string matrix:
% names = char(names{:});

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

