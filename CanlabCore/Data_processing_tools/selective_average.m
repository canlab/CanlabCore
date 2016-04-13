function [averages, stderrs, data, indices] = selective_average(y, onsets, varargin)
% Purpose: Get a selective average of values of data vector y, given
% onsets specified in onsets.  Onsets can be fractional; in this case,
% linear interpolation is used.
%
% :Usage:
% ::
%
%    [averages, stderrs, data, indices] = selective_average(y, onsets, varargin)
%
% :Inputs:
%
%   **y:**
%        is a data vector to get selective averages from.
%        It should be a column vector.
%
%   **onsets:**
%        should be a cell array, with one cell per condition
%        each cell should contain a column vector of onset times in SAMPLES (same
%        resolution as y; e.g., in TRs, if y is an fMRI time series.
%
% :Optional Inputs:
%
%   **t:**
%        followed by number of time points following onset to use;
%        default is 20
%
%   **plot:**
%        plot results. 
%
%   **baseline:**
%        followed by vector of which time points are baseline
%        values; will subtract from each
%
% :Outputs:
%
%   **data:**
%        indices, averages, stderrs:  Cell vectors, one cell per condition
%
%   **data, indices:**
%        time points (observations) x trials (onsets)
%
%   **indices:**
%        Cell vector, one cell per condition; time points (observations) x trials (onsets)
%
%
% :Examples:
% ::
%
%    onsets = {[1 10 30 80]'  [20 60 90]'}; y = (1:120)';
%    [averages, stderrs, data, indices] = selective_average(y, onsets, 't', 20)
%
%    V = spm_vol(EXPT.FILES.im_files{1});
%    y = spm_get_data(V, [10 10 10 1]');
%    [averages, stderrs, data, indices] = selective_average(y, onsets2(1), 't', 20, 'plot');
%
% ..
%    Tor Wager, Dec 2007
%    Minor update: June 2009
% ..

    % ..
    %    Optional inputs
    % ..
    t = 20;
    doplot = 0;
    basepts = [];
    
    % plot colors
    colors = {'ro-' 'go-' 'bo-' 'yo-' 'co-' 'mo-' 'r^:' 'g^:' 'b^:' 'y^:' 'c^:' 'm^:'};

        
    for i = 1:length(varargin)
        if ischar(varargin{i})
            switch varargin{i}

                % functional commands
                case {'t', 'timepoints'}, t = varargin{i+1};
                case 'plot', doplot = 1;
                    
                case {'baseline' 'basepts' 'basepoints'}, basepts = varargin{i+1};
                  
                case {'color' 'colors'}, colors = varargin{i+1};
                    
                otherwise, warning(['Unknown input string option:' varargin{i}]);
            end
        end
    end


    % ---------------------------------------------------------------------
    % Initialize vars
    % ---------------------------------------------------------------------
    % onsets should be a cell array, with one cell per condition
    % each cell should contain a column vector of onset times in SAMPLES
    n_conditions = length(onsets);

    indices = cell(1, n_conditions);
    data = cell(1, n_conditions);
    averages = cell(1, n_conditions);
    stderrs = cell(1, n_conditions);


    % ---------------------------------------------------------------------
    % Run
    % ---------------------------------------------------------------------
    for i = 1:n_conditions

        ons = onsets{i};

        [indices{i}, n_onsets{i}] = get_indices(ons, t);

        data{i} = get_data(y, indices{i}, t, n_onsets{i});

        [averages{i} stderrs{i}] = get_averages(data{i}, basepts, t);

    end


    if doplot
        
        plot_results(averages, stderrs, t, colors);
        
    end

end  % END MAIN FUNCTION




% Sub-functions


function [indices, n_onsets] = get_indices(ons, t)

    n_onsets = size(ons, 1);

    to_add = (1:t)' - 1;    % elements to add to form indices

    % matrix of values to add to each onset
    add_mtx = to_add(:, ones(1, n_onsets));

    % matrix form of onsets
    ons_mtx = ons'; ons_mtx = ons_mtx(ones(t, 1), :);

    indices = ons_mtx + add_mtx;

end




function data = get_data(y, indices, varargin)

    if all(indices(:) == round(indices(:)))
        % easy, all indices are integers; just get data
        data = y(indices);

    else
        t = varargin{1};
        n_onsets = varargin{2};
        
        % we have fractions of indices; linear interpolation
        % should return same as above for integer indices
        yest = interp1((1:length(y))', y, indices(:), 'linear');

        data = reshape(yest, t, n_onsets);
    end

end




function [avg, stderr] = get_averages(data, basepts, t)

    % Adjust data, if asked for: subtract baseline points
    % -----------------------------------------------------
    if ~isempty(basepts)

% %         fprintf('Subtracting baseline timepoints:')
% %         fprintf(' %3.0f', basepts);
% %         fprintf('\n')

        basemean = nanmean(data(basepts, :), 1);
        basemean = basemean(ones(t, 1), :);
        data = data - basemean;

    end
    
    % Now data should be time points (observations) x trials (onsets)
    avg = nanmean(data, 2);
    stderr = ste(data')';

end


function plot_results(averages, stderrs, t, varargin)
    
    to_add = (1:t)' - 1;    % elements to add to form indices
        
    if ~isempty(varargin), colors = varargin{1}; end
    
    %create_figure('Selective Average Plot');
    
    for i = 1:length(averages)
    
        if ischar(colors{i})
            plot(to_add, averages{i}, colors{i}, 'LineWidth', 2);
        else
            plot(to_add, averages{i}, 'Color', colors{i}, 'LineWidth', 2);
        end

    end
    
    xlabel('Time from onset (samples)');
    ylabel('Average response level');
    
end
