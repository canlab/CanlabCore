function out = get_wh_image(dat, wh)
% For an image_vector with multiple images (cases, contrasts, etc.), select one.
%
% :Usage:
% ::
%
%    function obj_out = get_wh_image(obj1, wh)
%
% :Inputs:
%
%   **obj1:**
%        An image_vector object
%
%   **wh:**
%        An integer indicating which image
%
% :Examples:
% ::
%
%    my_image_vector = get_wh_image(dat, 3) %to get 3rd image
%
% ..
%    Programmer Notes
%    Created 3/3/16 by Yoni Ashar
% ..

% check that wh is in range

if wh <= 0 || wh > size(dat.dat, 2)
    error('No image at index %i', wh);
end

out = dat;

fnames = fieldnames(dat);

for f=fnames'
    
    field = char(f);
    if isequal( size(out.(field)) , size(dat.dat) ) % this is a field to be selected down
        
        out.(field) = out.(field)(:,wh);
        
    end
end
