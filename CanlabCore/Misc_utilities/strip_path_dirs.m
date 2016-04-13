function strip_path_dirs(patterns)
% Removes directories from the Matlab path, based on the rgex patterns passed in.
%
% :Usage:
% ::
%
%     function strip_path_dirs(regexes)

    if(ischar(patterns))
        patterns = cellstr(patterns);
    end

    %matches = regexp(path(), sprintf('[^%s]+%s[^%s]+', pathsep, patterns{i}, pathsep), 'match');
    %matches = regexp(path(), sprintf('([^%s]+)', pathsep()), 'match');
    for i=1:length(patterns)
        matches = regexp(path(), sprintf('[^%s]*%s[^%s]*', pathsep, patterns{i}, pathsep), 'match');
        rmpath(implode(matches, pathsep));
    end

%     for i=1:length(bad_dirs)
%         rmpath(implode(bad_dirs{i}, pathsep));
%     end
% 
    savepath();
end
