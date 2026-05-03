function [o2, title_handle] = title_montage(o2, wh_montage, title_str)
% title_montage Add a title to a montage registered in an fmridisplay object.
%
% Adds a title to a montage figure registered in an fmridisplay object o2.
% Underscores in title_str are replaced with spaces before display.
%
% :Usage:
% ::
%
%     [o2, title_handle] = title_montage(o2, wh_montage, title_str)
%
% :Inputs:
%
%   **o2:**
%        An fmridisplay object with one or more montages attached.
%
%   **wh_montage:**
%        An integer that specifies the image axis (montage index) to
%        add the title to. It is usually 5 for standard canlab 2020
%        double-row montages.
%
%   **title_str:**
%        Title string to display. Underscores are replaced with spaces.
%
% :Outputs:
%
%   **o2:**
%        The fmridisplay object, with .montage{wh_montage}.title set
%        to title_str.
%
%   **title_handle:**
%        Handle to the title graphics object, or [] if the target axis
%        handle was invalid (e.g., figure was closed).
%
% :See also:
%   - fmridisplay
%   - montage

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