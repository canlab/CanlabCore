function [o2, title_handle] = title_montage(o2, wh_montage, title_str)
% Add a title to a montage figure registered in an fmridisplay object o2
%
% [o2, title_handle] = title_montage(o2, wh_montage, title_str)
%
% wh_montage is an integer that specifies the image axis to add the title to
% it is usually 5 for standard canlab 2020 double-row montages.

title_handle = [];

wh_axis = floor(length(o2.montage{wh_montage}.axis_handles) ./ 2);
wh_axis = max([wh_axis 1]);

axis_han = o2.montage{wh_montage}.axis_handles(wh_axis);

  
if ~isempty(title_str)
    title_str = strrep(title_str, '_', ' ');
end

% Save title
o2.montage{wh_montage}.title = title_str;

% Return if its not a real handle (maybe figure was closed?)
if ~ishandle(axis_han)
  return
end
    
% Display it
title_handle = title(axis_han, title_str);

set(title_handle, 'FontSize', 18);

end % function 