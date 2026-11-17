function ax = plot_hrf_results(results, varargin)
%PLOT_HRF_RESULTS Backward-compatible wrapper for plot_hrf_by_condition.
%
% ax = plot_hrf_results(results, 'Model', 'sfir', 'Conditions', {'pain'})
%
% New code should call plot_hrf_by_condition directly. This wrapper keeps
% existing examples working while using the clearer labeling and SE behavior.

ax = plot_hrf_by_condition(results, varargin{:});
end
