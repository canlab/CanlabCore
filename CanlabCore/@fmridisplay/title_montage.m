function o2 = title_montage(o2, wh_montage, title_str)
% o2 = title_montage(o2, wh_montage, title_str)

wh_axis = floor(length(o2.montage{wh_montage}.axis_handles) ./ 2);
wh_axis = max([wh_axis 1]);

axis_han = o2.montage{wh_montage}.axis_handles(wh_axis);

title_str = strrep(title_str, '_', ' ');

% Save title
o2.montage{wh_montage}.title = title_str;

% Display it
title(axis_han, title_str);

end % function 