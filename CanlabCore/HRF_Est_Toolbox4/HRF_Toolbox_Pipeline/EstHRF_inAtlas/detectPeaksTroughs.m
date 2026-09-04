function [peaks, troughs] = detectPeaksTroughs(waveform, shouldPlot)
    % Detect Peaks and Troughs from a waveform using Matlab Signal
    % Processing Toolbox
    % Michael Sun, Ph.D.
    % - Takes a vector timeseries e.g., generated from EstimateHRF_inAtlas
    % or hrf_fit_one_voxel().
    % - shouldPlot: 1: yes, plot; 0: no, don't plot
    %
    % *Usage:
    % ::
    %    [peaks, troughs] = detectPeaksTroughs(waveform, 1)

    % [tc, roi_val, maskdat]=canlab_connectivity_preproc(fmri_d, PREPROC_PARAMS.R, 'hpf', .018, PREPROC_PARAMS.TR, 'average_over', 'no_plots');

    if nargin < 2
        shouldPlot = false;
    end

    % Initialize
    peaks = struct([]);
    troughs = struct([]);
    auc_values = [];


    % Detecting Peaks:
    % [d_pos, l_pos, pos_width, pos_prominences] = findpeaks(waveform, 'MinPeakHeight', waveform(1), 'MinPeakDistance', 1);
    % threshold_prominence = 0.5 * max(pos_prominences);
    % 
    % % Detecting Troughs:
    % [d_neg, l_neg, neg_width, neg_prominences] = findpeaks(-waveform, 'MinPeakHeight', waveform(1), 'MinPeakDistance', 1);
    % threshold_prominence = 0.5 * max(neg_prominences);

    % Detecting Peaks:
    [d_pos, l_pos, ~, pos_prominences] = findpeaks(waveform, 'MinPeakHeight', waveform(1), 'MinPeakDistance', 1);
    if d_pos > 0
        peak_intervals=diff(l_pos);
        threshold_prominence = .5*max(pos_prominences);
        [d_pos, l_pos, pos_width, pos_prominences] = findpeaks(waveform, 'MinPeakHeight', waveform(1), 'MinPeakDistance', 1, 'MinPeakProminence', threshold_prominence, 'Annotate', 'extents');
    end

    % Detecting Troughs
    [d_neg, l_neg, ~, neg_prominences] = findpeaks(-waveform, 'MinPeakHeight', waveform(1), 'MinPeakDistance', 1);
    if d_neg > 0
        trough_intervals=diff(l_neg);
        threshold_prominence = .5*max(neg_prominences);
        [d_neg, l_neg, ~, neg_prominences] = findpeaks(-waveform, 'MinPeakHeight', waveform(1), 'MinPeakDistance', 1, 'MinPeakProminence', threshold_prominence, 'Annotate', 'extents');
    end    

    % fields = {'height', 'time_to_peak', 'width', 'start_time', 'end_time', 'half_height', 'AUC'};
    % Initialize with default values
    default_struct = struct('height', NaN, 'time_to_peak', NaN, 'width', NaN, 'start_time', NaN, 'end_time', NaN, 'half_height', NaN, 'AUC', NaN);
    peaks = repmat(default_struct, 1, length(l_pos));
    troughs = repmat(default_struct, 1, length(l_neg));

    % Optional Plotting
    if shouldPlot
        % figure;
        plot(waveform);
        hold on;
    end

    % Process Peaks
    for i = 1:length(l_pos)
        [start_point, end_point, auc] = processRegion(l_pos(i), waveform, d_pos(i), shouldPlot, 'ro', [1, 0.6, 0.6], false);
        peaks(i) = storeFeatures(d_pos(i), l_pos(i), start_point, end_point, auc);
    end

    % Process Troughs
    for j = 1:length(l_neg)
        [start_point, end_point, auc] = processRegion(l_neg(j), -waveform, d_neg(j), shouldPlot, 'bo', [0.6, 0.6, 1], true);
        troughs(j) = storeFeatures(-d_neg(j), l_neg(j), start_point, end_point, auc);
    end

    % hold off;

end

%% HELPER FUNCTIONS

function [start_point, end_point, auc] = processRegion(loc, waveform, height, shouldPlot, marker, fillColor, isTrough)
    start_point = find(waveform(1:loc) <= 0, 1, 'last') + 1;
    if isempty(start_point)
        start_point = 1;
    end
    
    end_point = find(waveform(loc:end) <= 0, 1, 'first') + loc - 2;
    if isempty(end_point)
        end_point = length(waveform);
    end

    % Ensure start_point is less than end_point
    if start_point > end_point
        temp = start_point;
        start_point = end_point;
        end_point = temp;
    end

    region_data = waveform(start_point:end_point);
    auc = abs(trapz(region_data));

    if shouldPlot
        x = start_point:end_point;
        y = waveform(x);
        x_fill = [x, fliplr(x)];
        y_fill = [y', zeros(1, numel(y))];
        
        if isTrough
            y_fill = -y_fill;  % Flip the y-values for troughs
            height = -height;  % Make the height negative for troughs
        end
        
        fill(x_fill, y_fill, fillColor, 'EdgeColor', 'none');
        plot(loc, height, marker);
        text(loc, height, sprintf('AUC: %.2f', auc));
    end
end


function features = storeFeatures(height, time_to_peak, start_time, end_time, auc)

    features.height = height;
    features.time_to_peak = time_to_peak;
    features.width = end_time - start_time;
    features.start_time = start_time;
    features.end_time = end_time;
    features.half_height = height / 2;
    features.AUC = auc;

end

