function c = power(obj, b)
% Implements the power (^) operator on image_vector objects across voxels.
%
% :Examples:
% ::
%
%    c = dat1^2;
%
% ..
% Programmer Notes:
% Created 3/14/14 by Luke Chang
% ..

error('This method is deprecated.  Use image_math instead.');


%Check if image_vector object
if ~isa(obj,'image_vector')
    error('Input Data is not an image_vector object')
end

c = obj;
c.dat = obj.dat.^b;


end % function
