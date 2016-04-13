function [vec, outnames] = contrast_code(vec)
% :Usage:
% ::
%
%     [vec,outnames] = contrast_code(vec)
%
% Changes values to 1, -1, or 0 for contrast coding

outnames = [];

if iscell(vec) && ~isempty(vec)
    [indic, names]  = string2indicator(vec);
    
    for i = 2:length(names)
        vec = indic(:, i - 1) - indic(:, i);
        outnames{i - 1} = [names{i - 1} ' - ' names{i}];
    end
    
elseif ~isempty(vec)
    wh = find(vec > 0);

    vec(wh) = 1;

    wh = find(vec < 0);

    vec(wh) = -1;

end

return

