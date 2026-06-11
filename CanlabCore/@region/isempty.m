function is_empty = isempty(obj)
% isempty Test whether a region object is empty.
%
% Returns true when the region object array has length 0, or when its
% first element has an empty XYZ field. Useful as a guard before
% iterating over regions.
%
% :Usage:
% ::
%
%     is_empty = isempty(obj)
%
% :Inputs:
%
%   **obj:**
%        A region-class object array.
%
% :Outputs:
%
%   **is_empty:**
%        Logical scalar; true if obj is empty in the sense above.
%
% :See also:
%   - region

is_empty = length(obj) == 0 || isempty(obj(1).XYZ) || length(obj(1).XYZ) == 0;


end

