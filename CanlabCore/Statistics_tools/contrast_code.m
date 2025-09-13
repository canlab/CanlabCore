function [vec, outnames] = contrast_code(vec)
% :Usage:
% ::
%
%     [vec,outnames] = contrast_code(vec)
%
% Attempts to change values to 1, -1, or 0 for contrast coding
% Success depends on inputs. Check output

outnames = [];

if iscell(vec) && ~isempty(vec)
    [indic, names]  = string2indicator(vec);
    
    for i = 2:length(names)
        vec = indic(:, i - 1) - indic(:, i);
        outnames{i - 1} = [names{i - 1} ' - ' names{i}];
    end
    
elseif ~isempty(vec)

    if all(vec == 0 | vec == 1) % dummy codes

        vec(vec == 0) = -1;

    else

    wh = find(vec > 0);

    vec(wh) = 1;

    wh = find(vec < 0);

    vec(wh) = -1;

    end

end

return

