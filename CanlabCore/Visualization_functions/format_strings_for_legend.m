function legendnames = format_strings_for_legend(legendnames)
% Takes a string matrix or cell array of strings (legendnames)
% Strips out characters, _, ^, .nii, .img, and .hdr
% Returns cell array with one string per cell
%
% See also: makelegend

if isempty(legendnames), return, end

if ~iscell(legendnames)
    legendnames = cellstr(legendnames);
end

legendnames = cellfun(@(x) strrep(x, '_', ' '), legendnames, 'UniformOutput', false);
legendnames = cellfun(@(x) strrep(x, '.nii', ''), legendnames, 'UniformOutput', false);
legendnames = cellfun(@(x) strrep(x, '.img', ''), legendnames, 'UniformOutput', false);
legendnames = cellfun(@(x) strrep(x, '.hdr', ''), legendnames, 'UniformOutput', false);
legendnames = cellfun(@(x) strrep(x, '^', ' '), legendnames, 'UniformOutput', false);


end % function