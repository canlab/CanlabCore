function [cverr, stats, optout] = predict(obj, varargin)
% predict Cross-validated multivariate prediction on grayordinate data.
%
% :Usage:
% ::
%     [cverr, stats, optout] = predict(obj, 'algorithm_name', 'cv_lassopcr', 'nfolds', 5)
%
% Runs CANlab's cross-validated prediction (fmri_data.predict) on surface/
% grayordinate data. The algorithm operates on obj.dat and obj.Y exactly as for
% volumetric data, so all algorithms and options are supported unchanged; this
% method delegates to fmri_data.predict via a proxy (each grayordinate treated as
% a feature) and then remaps the geometry-bearing weight map back to an
% fmri_surface_data (so you can surface()/write() it).
%
% Set the outcome to predict in obj.Y before calling (as with fmri_data).
%
% :Inputs:
%   **obj:** fmri_surface_data with obj.Y set (one value per map/observation).
%
% :Optional Inputs:
%   All fmri_data.predict options ('algorithm_name', 'nfolds', 'error_type', ...).
%
% :Outputs:
%   **cverr:**  cross-validated error/accuracy (as fmri_data.predict).
%   **stats:**  stats struct; stats.weight_obj is remapped to an
%               fmri_surface_data carrying the object's geometry.
%   **optout:** optional algorithm outputs (passed through).
%
% :Examples:
% ::
%     obj.Y = behavioral_scores(:);
%     [err, stats] = predict(obj, 'algorithm_name', 'cv_lassopcr', 'nfolds', 5);
%     surface(stats.weight_obj);     % render the predictive weight map
%
% :See also: regress, fmri_data.predict, rebuild_like, fmri_surface_data

if isempty(obj.Y)
    error('fmri_surface_data:predict:noY', ...
        'Set obj.Y (one outcome value per map/observation) before calling predict.');
end

proxy = as_fmri_data_proxy(obj);

[cverr, stats, optout] = predict(proxy, varargin{:});

% Remap the predictive weight map back to a surface object
if isfield(stats, 'weight_obj') && ~isempty(stats.weight_obj) ...
        && size(stats.weight_obj.dat, 1) == size(obj.dat, 1)
    stats.weight_obj = rebuild_like(obj, double(stats.weight_obj.dat));
    stats.weight_obj.image_names = {'predictive weights'};
end
end
