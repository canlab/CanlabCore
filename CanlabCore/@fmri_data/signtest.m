function [out, statimg] = signtest(dat, varargin)
% Sign test for each voxel of an fmri_data object
% returns voxel-wise statistic images.
%
% :Usage:
% ::
%
%     [out, statimg] = signtest(dat, [p-val threshold], [thresh_type])
%
% :Inputs:
%
%   **dat:**
%        Should be an fmri_data object with .dat field containing voxels x observations matrix
%
% :Optional inputs: in [  ] above are
%   **p-value threshold:**
%        string indicating threshold type (see help statistic_image.threshold for options)
%
% :Outputs:
%   **out:**
%        is a structure of information about the sign test
%
%   **statimg:**
%        is a statistic_image object that can be thresholded and
%        plotted/imaged.  statimg.dat contains signed direction values,
%       .p contains p-values 
%
% ..
% c Tor Wager, 2011
% ..
%
% See also: fmri_data.regress
%

inputargs = {.001, 'uncorrected'};
% default options for thresholding

if ~isempty(varargin)
    for i = 1:length(varargin)
        inputargs{i} = varargin{i};
    end
end
    
n = size(dat.dat, 2);

% display info about regression
fprintf('sign test > Y: %3.0f voxels. X: %3.0f obs, %3.0f tests\n', size(dat.dat, 1), n);


% model fit
tic

out = signtest_matrix(dat.dat');

toc

% save results
statimg = statistic_image;
statimg.type = 'signtest';
statimg.p = out.p;
statimg.ste = NaN .* out.p;
statimg.N = out.n;
statimg.dat = out.direction;
statimg.dat_descrip = sprintf('Sign test direction: Pos for > 0, Neg for < 0');
statimg.volInfo = dat.volInfo;
statimg.removed_voxels = dat.removed_voxels;
statimg.removed_images = false;  % this image does not have the same dims as the original dataset

statimg = threshold(statimg, inputargs{:});

% view results
orthviews(statimg);

end % function
