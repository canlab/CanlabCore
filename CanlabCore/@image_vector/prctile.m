function threshold_values = prctile(obj, percentile_values)
% prctile Calculate thresholds in image data values for specified percentile values.
%
% Cumulative distribution: percentile_values(i) % of data values are at
% or below threshold_values(i). Thresholds for percentiles are estimated
% across the entire image matrix (4-D). Excludes zero-valued elements of
% the matrix.
%
% :Usage:
% ::
%
%     threshold_values = prctile(obj, percentile_values)
%
% :Inputs:
%
%   **obj:**
%        An image_vector / fmri_data / statistic_image object.
%
%   **percentile_values:**
%        Vector of percentile values (0-100) at which to compute
%        thresholds.
%
% :Outputs:
%
%   **threshold_values:**
%        Vector of data thresholds, the same length as
%        percentile_values, such that percentile_values(i) % of data
%        values are at or below threshold_values(i).
%
% :Examples:
% ::
%
%     perc_vals = [.1 1 5 95 99 99.9];
%     obj = load_image_set('emotionreg');
%     thr = prctile(obj, perc_vals);
%
% :See also:
%   - prctile
%   - threshold
%
% ..
%    Tor Wager, May 2018
% ..

mydat = double(obj.dat(:)); 

mydat(mydat == 0) = [];
threshold_values = prctile(mydat, percentile_values);

end % function


