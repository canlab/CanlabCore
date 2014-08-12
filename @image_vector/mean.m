function m = mean(obj, varargin)
% function m = mean(obj, [optional args])
%
% Create an image_vector object with mean values for each voxel (cols)
% across images (rows) of an fmri_data object.
%
% m is an image_vector object whose data contains the mean values.
%
% Options are:
% 'write', followed by file name
% 'path', followed by location for file (default = current directory)
% 'orthviews' -> show orthviews for this image, same as orthviews(m)
% 'histogram' -> show histogram for this image, same as histogram(m)
% 'plot' -> do both
%
% Examples:
% If sdat is an fmri_data object with multiple images,
% m = mean(sdat, 'plot', 'write', anatmeanname, 'path', maskdir);
%

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

if isa(obj, 'fmri_data')
    m = image_vector('dat', mean(obj.dat', 1)', 'volInfo', obj.mask.volInfo);
else
    m = image_vector('dat', mean(obj.dat', 1)', 'volInfo', obj.volInfo);
end

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


end % function