function c = horzcat(varargin)
% Implements the horzcat ([a b]) operator on image_vector objects across voxels.
% Requires that each object has an equal number of columns and voxels
%
% :Usage:
% ::
%
%    function s = horzcat(varargin)
%
% :Example:
% ::
%
%    c = [dat1 dat2];
%
% ..
%    Programmer Notes
%    Created 3/14/14 by Luke Chang for image_vector; updated for fmri_data 8/2015 Yoni Ashar
% ..

hasX = 1; X = [];
hasY = 1; Y = [];
% check whether ALL inputs have X, Y values

for i = 1:nargin
    if isempty(varargin{i}.X), hasX=0;end % if any is missing, set to 0
    if isempty(varargin{i}.Y), hasY=0;end % if any is missing, set to 0
end

dat = [];
for i = 1:nargin
    %Check if fmri_data object
    if ~isa(varargin{i}, 'fmri_data')
        error('Input Data is not an image_vector object')
    end
    
    dat = [dat, varargin{i}.dat];
    
    if hasX, X = [X; varargin{i}.X]; end
    
    if hasY, Y = [Y; varargin{i}.Y]; end
end

c = varargin{1};
c.dat = dat;
if hasX, c.X = X; end
if hasY, c.Y = Y; end
