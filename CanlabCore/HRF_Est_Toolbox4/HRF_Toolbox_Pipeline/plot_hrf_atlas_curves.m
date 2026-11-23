function fig = plot_hrf_atlas_curves(varargin)
%PLOT_HRF_ATLAS_CURVES Deprecated alias for PLOT_HRF_CURVES.
%
% This function was renamed to plot_hrf_curves once it grew to plot
% signatures and image-set networks in addition to atlas regions. This
% thin shim forwards all arguments unchanged so existing scripts keep
% working. Prefer plot_hrf_curves in new code.
%
% See also: plot_hrf_curves.

fig = plot_hrf_curves(varargin{:});
end
