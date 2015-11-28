function c = minus(obj1,obj2)
% Implements the minus (-) operator on image_vector objects across voxels.
% Requires that each object has an equal number of columns and voxels
%
% ..
%    Programmer Notes:
%    Created 3/14/14 by Luke Chang
% ..

error('This method is deprecated.  Use image_math instead.');

%Check if image_vector object
if ~isa(obj1,'image_vector') || ~isa(obj2,'image_vector')
    error('Input Data is not an image_vector object')
end

%Check number of rows
if size(obj1.dat,1)~=size(obj2.dat,1)
    error('number of voxels is different between objects.')
end

%Check number of columns
if size(obj1.dat,2)~=size(obj2.dat,2)
    error('number of observations is different between objects.')
end

c = obj1;
c.dat = obj1.dat - obj2.dat;


end % function
