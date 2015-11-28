function c = horzcat(varargin)
% Implements the horzcat ([a b]) operator on image_vector objects across voxels.
% Requires that each object has an equal number of columns and voxels
%
% :Usage:
% ::
%
%    function s = horzcat(varargin)
%
% :Examples:
% ::
%
%    c = [dat1 dat2];
%
% ..
%    Programmer Notes
%    Created 3/14/14 by Luke Chang
% ..

error('This method is deprecated.  Use image_math instead.');

dat = [];
for i = 1:nargin
    %Check if image_vector object
    if ~isa(varargin{i}, 'image_vector')
        error('Input Data is not an image_vector object')
    end
    
    dat = [dat, varargin{i}.dat];
end

c = varargin{1};
c.dat = dat;
