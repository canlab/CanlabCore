function icadat = ica(obj, varargin)
% ica Spatial ICA decomposition of grayordinate data.
%
% :Usage:
% ::
%     icadat = ica(obj, [number of ICs])
%
% Runs independent component analysis on obj.dat (mirroring image_vector.ica) by
% delegating to a proxy and remapping the spatial component maps back to an
% fmri_surface_data, so the components can be surface()-rendered.
%
% DEPENDENCY: like image_vector.ica, this requires icatb_fastICA (the GIFT /
% GroupICA "icatb" toolbox) on the path. This is the one fmri_surface_data method
% that is not fully self-contained; all other methods need no external toolbox.
%
% :Inputs:
%   **obj:** fmri_surface_data.
%   **nic:** (optional) number of components to save (passed to image_vector.ica).
%
% :Outputs:
%   **icadat:** fmri_surface_data whose maps are the spatial IC maps.
%
% :See also: predict, image_vector.ica, rebuild_like, fmri_surface_data

proxy = as_fmri_data_proxy(obj);
comp_obj = ica(proxy, varargin{:});

if isa(comp_obj, 'image_vector') && size(comp_obj.dat, 1) == size(obj.dat, 1)
    icadat = rebuild_like(obj, double(comp_obj.dat));
    icadat.image_names = arrayfun(@(k) sprintf('IC%d', k), 1:size(icadat.dat,2), ...
        'UniformOutput', false)';
    icadat.history{end+1} = sprintf('ica: %d spatial components', size(icadat.dat,2));
else
    icadat = comp_obj;   % unexpected shape: pass through
end
end
