function m = mean(obj, varargin)
% mean Average across maps/images of an fmri_surface_data, returning one map.
%
% :Usage:
% ::
%     m = mean(obj)
%     m = mean(obj, 'omitnan')
%
% Averages obj.dat across the image (column) dimension and returns a new
% single-map fmri_surface_data carrying the same geometry (via rebuild_like).
% Surface analogue of fmri_data.mean; no volume rebuild.
%
% :Optional Inputs:
%   **'omitnan':** ignore NaNs in the average (default includes them).
%
% :Outputs:
%   **m:** fmri_surface_data with one map (the mean across images).
%
% :See also: fmri_surface_data, rebuild_like

flag = 'includenan';
if any(strcmpi(varargin, 'omitnan')), flag = 'omitnan'; end

mdat = mean(double(obj.dat), 2, flag);
m = rebuild_like(obj, mdat);
m.image_names = {'mean'};
m.history{end+1} = sprintf('mean across %d maps', size(obj.dat, 2));
end
