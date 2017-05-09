function riverplot_toggle_lines(ribbons)
%
% riverplot_set_ribbon_property(ribbons, propertyname, propertyval)
%
% e.g., riverplot_set_ribbon_property(ribbons, 'LineStyle', 'none')

if ~iscell(ribbons) || isempty(ribbons) || all(all(cellfun(@isempty, ribbons)))
    return % nothing to do
end

allribbons = cat(1, ribbons{:});

% all ribbon patches
allpatches = cat(1, allribbons(:).patchh);


for i = 1:length(allribbons)
    line1h{i} = allribbons(i).line1.h;
    line2h{i} = allribbons(i).line2.h;
    
end

line1h = cat(1, line1h{:});
line2h = cat(1, line2h{:});

visibility = get(line1h(1), 'Visible');

switch visibility
    case 'off'
        set(allpatches, 'LineStyle', '-');
        set(line1h, 'Visible', 'on');
        set(line2h, 'Visible', 'on');
        
    case 'on'
        set(allpatches, 'LineStyle', 'none');
        set(line1h, 'Visible', 'off');
        set(line2h, 'Visible', 'off');
        
    otherwise
        error('unknown visibility...check code');
end

end % function