function threshold_values = prctile(obj, percentile_values)
% Cumulative distribution: Calculate thresholds in data values for a series of percentile values you specify
% - percentile_values(i) % of data values are at or below threshold_values(i) 
% - Thresholds for percentiles are estimated across the entire image matrix (4-D)
% - Excludes zero-valued elements of matrix
%
% Example:
%
% perc_vals = [.1 1 5 95 99 99.9]; 
% obj = load_image_set('emotionreg');
% thr = prctile(obj, perc_vals);
%
% Tor Wager, May 2018

mydat = double(obj.dat(:)); 

mydat(mydat == 0) = [];
threshold_values = prctile(mydat, percentile_values);

end % function


