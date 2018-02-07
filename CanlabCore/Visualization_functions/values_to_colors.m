function mycolors = values_to_colors(myvalues, value_limits, startcolor, endcolor)
% Turn a series of data values into a cell array of mapped colors for each value
%
% mycolors = values_to_colors(myvalues, startrange, endrange, startcolor, endcolor)

% get colormap values
cm = colormap_tor(startcolor, endcolor);

k = size(cm, 1);

valence_vals = linspace(value_limits(1), value_limits(2), k); % to match with valence to get colors

if ~iscolumn(myvalues), myvalues = myvalues'; end

% get closest color-mapped value for each pattern_valence value
diffs = bsxfun(@minus, myvalues, valence_vals);
[~, wh] = min(abs(diffs), [], 2);

% get colors for these values
mycolors = cm(wh, :);

% turn into cell
mycolors = mat2cell(mycolors, ones(size(mycolors, 1), 1), 3)';

end % function