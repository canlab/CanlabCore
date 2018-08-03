function [o2, title_handle] = title_montage(o2, wh_montage, title_str)
% o2 = title_montage(o2, wh_montage, title_str)

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

end % function 