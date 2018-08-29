function is_empty = isempty(obj)
% Tests if an image_vector (or fmri_data, statistic_image, atlas) object is empty.
%
% :Usage:
% ::
%
%     is_empty = isempty(obj)
%
% :Inputs:
%
%   **obj:**
%        An image_vector object or one of its subclasses
%
% :Outputs:
%
%   **is_empty:**
%        Logical value (1/0)
%
% :See also:
%   - methods(image_vector)
%   - methods(fmri_data)
%


is_empty = isempty(obj.dat) || isempty(obj.volInfo) || all(obj.dat(:) == 0);

end