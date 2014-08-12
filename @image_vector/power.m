function c = power(obj, b)
% function s = power(obj1, obj2)
%
% Implements the power (^) operator on image_vector objects across voxels.
%
% Examples:
% c = dat1^2;
%
% Programmer Notes
% Created 3/14/14 by Luke Chang


%Check if image_vector object
if ~isa(obj,'image_vector')
    error('Input Data is not an image_vector object')
end

c = obj;
c.dat = obj.dat.^b;


end % function