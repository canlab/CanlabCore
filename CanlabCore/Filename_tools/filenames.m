% A utility function that uses ls or dir on a pattern and returns a list of
% filenames
%
% :Usage:
% ::
%
%     files = FILENAMES(pattern, ['absolute'], ['cellstr'|'char'], ['sortField', fieldNums], ['sortFieldSeparator', separator_string], ['verbose'])
%
% :Inputs:
%
%   **'absolute':**
%        forces filenames to have an absolute path from root - default
%        is to return based upon the pattern (i.e., a relative pattern will
%        produce relative paths, an absolute pattern will return absolute
%        paths
%
%   **'cellstr':**
%        returns filenames as a cell array of strings, a cellstr -
%        this is the default behavior
%
%   **'char':**
%        return filenames as a character array
%
%   **'verbose':**
%        display commands
%
% :PARAMS:
%
%   **'sortField':**
%        if set, sort the filenames numerically based on the fields specified -
%        type "man sort" in a terminal window for more info
%
%   **'sortFieldSeparator':**
%        string to determine what is considered a field
%        in the filename - defaults to whitespace
%
% :Examples:
% ::
%
%    % returns a cellstr of dcm files, sorted by the 3rd field, as
%    % separated by the '-' character
%    dcm_files = filenames('*.dcm', 'sortField', 3, 'sortFieldSeparator', '-');
%
% ..
%    tor wager
% ..

function files = filenames(pattern, varargin)

    returnChar = 0;
    mustBeAbsolute = 0;
    sortField = [];
    sortFieldSeparator = [];
    verbose = 0;

    for i=1:length(varargin)
        if(ischar(varargin{i}))
            switch(varargin{i})
                case 'char'
                    returnChar = 1;
                case 'cellstr'
                    returnChar = 0;
                case 'absolute'
                    mustBeAbsolute = 1;
                case 'sortField'
                    sortField = varargin{i+1};
                case 'sortFieldSeparator'
                    sortFieldSeparator = varargin{i+1};
                case 'verbose'
                    verbose = 1;
            end
        end
    end

    if(mustBeAbsolute && ~isAbsolutePath(pattern))
        pattern = fullfile(pwd(), pattern);
    end

    isPCandIsRelative = pcIsRelativePath(pattern);
    command = shellCommand(pattern, isPCandIsRelative, sortField, sortFieldSeparator);
    if(verbose), disp(command); end
    [status, output] = system(command);
    outputString = java.lang.String(deblank(output));

    if(outputString.indexOf('Parameter format not correct') ~= -1)
        error('Incorrect format for "%s". Are you sure you''re using the right kind of slashes?', pattern);
    elseif(outputString.indexOf('No such file or directory') ~= -1 || outputString.indexOf('File Not Found') ~= -1)
        if(returnChar), files = []; else files = {}; end
    else
        files = char(outputString.split('\n'));
        if(isPCandIsRelative)
            basepath = fileparts(pattern);
            if(~isempty(basepath))
                files = [repmat([basepath filesep()], size(files,1), 1) files];
            end
        end
        if(~returnChar)
            files = cellstr(files);
        end
    end
    files = removeBlankLines(files);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Subfunctions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function command = shellCommand(pattern, isPCandIsRelative, sortField, sortFieldSeparator)
    if(isunix())
        command = sprintf('ls -1d %s', escapeForShell(pattern));
        if(~isempty(sortField))
            command = [command ' | ' sortCommand(sortField, sortFieldSeparator)];
        end
    elseif(ispc())
        if(~isempty(sortField))
            error('Cannot currently sort on PCs. No sort command.');
        end
        whFileSeps = strfind(pattern, filesep());
        whWildcards = strfind(pattern, '*');
        % if there's a wildcard before the last file separator, we have a problem
        if(length(whFileSeps) > 0 && length(whWildcards) > 0 && whWildcards(1) < whFileSeps(end))
            error('There appears to be a wildcard before the last file separator in "%s", which Windows cannot handle, unfortunately.', pattern);
        else
            if(isPCandIsRelative)
                command = sprintf('dir /B %s', escapeForShell(pattern));
            else
                command = sprintf('dir /B /S %s', escapeForShell(pattern));
            end
        end
    else
        error('What kinda system you got?');
    end
end

function command = sortCommand(sortField, sortFieldSeparator)
    commandParams = {};
    if(~isempty(sortFieldSeparator))
        commandParams{end+1} = sprintf('-t ''%s''', sortFieldSeparator);
    end

    for i=1:length(sortField)
        commandParams{end+1} = [' -k ' num2str(sortField(i))];
    end
    command = sprintf('sort -n %s', implode(commandParams, ' '));
end

function isPCandIsRelative = pcIsRelativePath(pattern)
    isPCandIsRelative = ispc() && ~isAbsolutePath(pattern);
end

function isAbsolute = isAbsolutePath(pattern)
    isAbsolute = 0;
    if(ispc() && (~isempty(strfind(pattern, '\\')) || ~isempty(strfind(pattern, ':\'))))
        isAbsolute = 1;
    elseif(isunix())
        location = strfind(pattern, '/');
        if(~isempty(location) && any(location == 1))
            isAbsolute = 1;
        end
    end
end

function files = removeBlankLines(files)
    if(iscellstr(files))
        wh_empty = cellfun(@isempty, files);
        files(wh_empty) = [];
    elseif(ischar(files))
        files = char(removeBlankLines(cellstr(files)));
    end
end
