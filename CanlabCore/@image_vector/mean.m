function [m, varargout] = mean(obj, varargin)
% Mean across a set of images. Returns a new image_vector object.
% Creates an image_vector object with mean values for each voxel (cols)
% across images (rows) of an image_vector (e.g., fmri_data) object.
%
% :Usage:
% ::
%
%    function [m, imagemeans, voxelmeans]  = mean(obj, [optional args])
%
% m is an image_vector object whose data contains the mean values.
%
% - Average available valid data in each voxel. Some images may have
% missing data for some voxels.
% 
% :Optional Inputs:
%   - 'write', followed by file name
%   - 'path', followed by location for file (default = current directory)
%   - 'orthviews' -> show orthviews for this image, same as orthviews(m)
%   - 'histogram' -> show histogram for this image, same as histogram(m)
%   - 'plot' -> do both
%
% :Examples:
% ::
%
%    % If sdat is an fmri_data object with multiple images,
%    m = mean(sdat, 'plot', 'write', anatmeanname, 'path', maskdir);
%

% Programmers' notes:
% - Tor: Average available valid data in each voxel (June 2017)

fname = [];
fpath = pwd;
dohist = 0;
doorth = 0;
doplot = 0;

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            % functional commands
            case 'write', fname = varargin{i+1}; varargin{i+1} = [];
            case 'path', fpath = varargin{i+1}; varargin{i+1} = [];
            case 'plot', dohist = 1; doorth = 1;
                
            case {'hist', 'histogram'}, dohist = 1;
            case {'orth', 'orthviews'}, doorth = 1;
                
            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end


% get missing values and replace with NaNs, so we average only valid data in each
% voxel
wh = obj.dat == 0; % NaNs are already OK | isnan(obj.dat);
obj.dat(wh) = NaN;

% return output in the same format as input object

if isa(obj, 'fmri_data')
    m = image_vector('dat', nanmean(obj.dat', 1)', 'volInfo', obj.mask.volInfo);
    m = fmri_data(m);
    m.mask = obj.mask;
    
else
    m = image_vector('dat', mean(obj.dat', 1, 'omitnan')', 'volInfo', obj.volInfo);
end

% Not completed for statistic_image
% if isa(obj, 'statistic_image')
%     m = statistic_image(m);
% end

m.removed_voxels = obj.removed_voxels;

if doplot || doorth
    orthviews(m);
end

if doplot || dohist
    create_figure('image_histogram');
    histogram(m);
end

if ~isempty(fname)
    fullp = fullfile(fpath, fname);
    fprintf('Writing mean image to disk')
    
    m.filename = fname;
    m.fullpath = fullp;
    
    write(m);
end

% Optional outputs:
% Row means, Column means, double-centered object

if nargout > 1
    
    imagemeans = mean(obj.dat);
    varargout{1} = imagemeans;
    
end

if nargout > 2
    voxelmeans = nanmean(obj.dat, 2);
    varargout{2} = voxelmeans;
end

end % function
