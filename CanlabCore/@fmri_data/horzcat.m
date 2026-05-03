function c = horzcat(varargin)
% horzcat Implements the horzcat ([a b]) operator on fmri_data objects.
%
% :Usage:
% ::
%
%     c = horzcat(dat1, dat2, ...)
%     c = [dat1 dat2 ...]            % equivalent
%
% Concatenates fmri_data objects horizontally across the image (column)
% dimension. Each input object's .dat is expected to have the same number
% of voxels (rows); the resulting object has the union of images in
% column order. If every input has a non-empty .X (predictors), the .X
% matrices are vertically concatenated; same for .Y (outcomes). If any
% input is missing X or Y, those fields are dropped from the result.
%
% :Inputs:
%
%   **varargin:**
%        Two or more fmri_data objects with matching voxel counts. All
%        inputs must be of class fmri_data; passing any non-fmri_data
%        argument raises an error.
%
% :Outputs:
%
%   **c:**
%        An fmri_data object with .dat = [dat1.dat, dat2.dat, ...]. If
%        every input had a populated .X / .Y, those are stacked
%        vertically and attached.
%
% :Examples:
% ::
%
%     c = [dat1 dat2];
%     c = horzcat(dat1, dat2, dat3);
%
% :See also:
%   - cat (concatenation method on fmri_data)
%   - get_wh_image (subset of images)
%
% ..
%    Programmer Notes:
%    Created 3/14/14 by Luke Chang for image_vector; updated for fmri_data
%    8/2015 Yoni Ashar.
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
