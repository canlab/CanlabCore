function obj = threshold(obj, varargin)
% Re-threshold the statistic maps of a glm_map object (no refitting).
%
% Applies a statistical threshold to the fitted t and/or contrast_t maps by
% delegating to statistic_image.threshold. The underlying statistic values
% are preserved; only the significance mask changes.
%
% :Usage:
% ::
%
%     obj = threshold(obj, pval, thresh_type, varargin)
%
% :Inputs:
%
%   **obj:**
%        A fitted glm_map object (obj.is_fitted == true).
%
%   Remaining inputs are passed through unchanged to
%   statistic_image.threshold, e.g. (.001, 'unc', 'k', 10) or (.05, 'fdr').
%
% :Optional Inputs:
%
%   **'which_map', 'contrast' | 't' | 'both':**
%        Which map(s) to threshold. Default 'both' (t and contrast_t when
%        present).
%
% :Outputs:
%
%   **obj:**
%        glm_map with the selected statistic_image maps re-thresholded.
%
% :Examples:
% ::
%
%     g = threshold(g, .001, 'unc', 'k', 10);
%     g = threshold(g, .05, 'fdr', 'which_map', 'contrast');
%
% :See also:
%   - statistic_image.threshold
%
% ..
%    Programmers' notes:
%    2026 - Initial implementation (thin delegation).
% ..

if ~obj.is_fitted
    error('glm_map:NotFitted', 'Object is not fitted; run fit(obj, data) before thresholding.');
end

% Pull out the which_map keyword if present
which_map = 'both';
wh = find(strcmpi(varargin, 'which_map'));
if ~isempty(wh)
    which_map = varargin{wh(1) + 1};
    varargin(wh(1):wh(1) + 1) = [];
end

do_t   = ismember(lower(which_map), {'t', 'both'});
do_con = ismember(lower(which_map), {'contrast', 'both'});

if do_t && ~isempty(obj.t)
    obj.t = threshold(obj.t, varargin{:});
end

if do_con && ~isempty(obj.contrast_t)
    obj.contrast_t = threshold(obj.contrast_t, varargin{:});
end

obj.history{end + 1} = 'threshold: re-thresholded statistic maps';

end % threshold
