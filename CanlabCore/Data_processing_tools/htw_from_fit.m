function [h, t, w, auc, w_times, halfh] = htw_from_fit(hrf, b, dt, varargin)
% Estimates height, time to peak,  width, and area under the curve of a fitted response
%
% :Usage:
% ::
%
%     h, t, w, auc, w_times, halfh] = htw_from_fit(hrf, b, dt, [optional arguments])
%
% ..
%    tor wager, Sept 08
% ..
%
% :Inputs:
%
%   **hrf:**
%        a matrix of columns that form a linear basis set for an
%           event type in an fMRI design, t time points x k basis
%           functions
%   **b:**
%        betas associated with the columns of hrf, k x 1
%   **dt:**
%        the sampling resolution of hrf (in seconds)
% 
% :Optional Inputs:
%
%   **plot:**
%        make plot
%
%   **verbose:**
%        verbose output
%
%   **startval:**
%        followed by the starting value in sec within which to calculate peak
%
%   **endval:**
%        followed by the ending value in sec within which to calculate peak
%
%   **colors:**
%        followed by a cell array, for example, {'r'} or {[1 0 0]}
%
% This function is essentially the same as fir2htw2.m, but the main
% differences are:
%   1) It imposes a time constraint on the peak amplitude automatically,
%      which is constrained to be between 4 seconds and 12 seconds
%      (endval, which was hconstraint) by default. YOU MAY WANT TO CHANGE hconstraint
%      depending on whether you're expecting delayed hemodynamic responses.
%      This requires input of the sampling resolution (e.g., dt)
%
%   2) This function will automatically create fitted responses, given a
%      basis set and betas (parameters).  This is different from fir2htw2.m,
%      which takes the fitted response as input.
%
%   3) The method for getting width (w) has been changed to work better
%      for multi-modal (multi-peak) responses.
%
% Otherwise, the algorithm is the same.
% See Lindquist and Wager, 2007, for simulations that use a version of
% this method to estimate HRFs using different kinds of models.
%
% :Notes on scaling:
% The scaling of the amplitude depends on the scaling of the hrf basis
% set, which (in SPM) depends on the time resolution.  You should at
% least use an hrf basis set with the same scaling for all subjects in
% a group analysis.  The amplitude of the fitted response is
% interpreted as the amplitude of the unit "impulse response," assuming
% that the hrf you enter here is the same as the impulse response
% function you used to create the design matrix.  In SPM5, higher-res
% impulse response functions are normalized by their positive sum, and
% the higher the time resolution, the lower the amplitude of the unit
% HRF.  (The scaling of regressors in the SPM design matrix isn't
% affected, because the hrf basis functions are convolved with a boxcar
% that also depends on the time resolution.  The bottom line is that if
% all your subjects have the same scaling, you should be fine.  And,
% secondly, the amplitudes that come out of this function reflect the
% scaling of the HRF you put in and are for an impulse response, NOT
% for an "event," and so the scaling here would not be expected to
% match up with the amplitudes on a plot that you'd get from a
% selective average time-course plot, unless you adjust by multiplying
% by the number of elements in SPMs "hi-res" onset boxcar to adjust.
%
% :For example:
% With "zero-duration" events, an hrf input scaled to reflect "event response" amplitudes
% might look something like this: (***may not be exactly right because
% i think dt is in sec)
% figure; plot(conv(SPM.xBF.bf(:, 1), my_ons))
% my_ons = ones(1, TR ./ SPM.xBF.dt .* SPM.Sess.U(1).dur(1));
%
% If you have epochs and want "epoch response" amplitude, you have to consider that as well.
% If your durations are specified in TRs, and all durations are the
% same:
% TR = SPM.xY.RT;
% my_ons = ones(1, TR ./ SPM.xBF.dt .* SPM.Sess.U(1).dur(1));
%
% minh = min height
%
% This version uses turning points (zero gradient) to find the largest
% "hump" in the data and the time it occurs.
%
% :Examples:
% ::
%
%    % Load and SPM mat file and use the basis set stored in that, and use
%    % that as the hrf.  Generate some arbitrary combos to test different shapes:
%    % cd('my spm directory')
%    % load SPM
%    [h, t, w, auc] = htw_from_fit(SPM.xBF.bf, [1 .4 .4]', SPM.xBF.dt, 'plot', 'verbose');
%
%    for i = 1:20
%        [h, t, w, auc] = htw_from_fit(SPM.xBF.bf, randn(3, 1), SPM.xBF.dt, 'plot'); pause(1.5);
%    end
%
%    % Generate an SPM basis set at a lower resolution, and try that:
%    bf = spm_get_bf(struct('name', 'hrf (with time and dispersion derivatives)', 'length', 30, 'dt', 1));
%    for i = 1:20
%       [h, t, w, auc] = htw_from_fit(bf.bf, randn(3, 1), bf.dt, 'plot'); h, t, w, auc, pause(1.5)
%    end
%

% ..
%    Variable input arguments
% ..

    doplot = 0;
    startval_sec = 4;
    endval_sec = 12;
    verbose = 0;
    colors = {'ro-'};

    for i = 1:length(varargin)
        if ischar(varargin{i})
            switch varargin{i}

                % reserved keywords
                case 'plot', doplot = 1;
                case 'verbose', verbose = 1;

                    % optional inputs
                case 'startval', startval_sec = varargin{i+1}; % the ending value in sec within which to calculate peak
                case 'endval', endval_sec = varargin{i+1}; % the ending value in sec within which to calculate peak

                case 'colors', colors = varargin{i+1};

                otherwise, warning(['Unknown input string option:' varargin{i}]);
            end
        end
    end

    % Set starting and ending values in samples (depending on dt, time res of hrf)
    % --------------------------------------------------------

    startval_samples = round(startval_sec ./ dt);
    endval_samples = round(endval_sec ./ dt);

    [rows, cols] = size(hrf);

    if startval_samples < 1
        error('Illegal starting value in seconds.');
    end

    if endval_samples > rows
        error('There do not seem to be enough samples in the HRF for the desired amplitude estimation window.')
    end

    if verbose
        fprintf('Number of basis functions: %3.0f\n', cols);
        fprintf('Length of HRF in seconds : %3.1f seconds\n', rows .* dt);

        fprintf('Time resolution is : %3.4f seconds per sample in hrf\n', dt);
        fprintf('Computing amplitude from max/min between : %3.2f to %3.2f seconds, %3.0f to %3.0f samples\n', startval_sec, endval_sec, startval_samples, endval_samples);
    end

    % Make sure data format is correct, betas are column vector

    if size(b, 1) ~= size(hrf, 2)
        warning('htw:BadInput', 'Betas should be column vectors...trying transpose...');
        b = b'; 
    end

    % get fitted response
    % --------------------------------------------------------

    fit = hrf * b;

    % for matrix
    % ***still working
    ntime = size(fit, 1);
    ntests = size(fit, 2);
    
    isbad = all(fit - fit(ones(ntime, 1), :) == 0);
    isbad = isbad | any(fit - fit(ones(size(fit, 1), 1), :) == NaN);

    % exit if empty
    if all(isbad)  %isempty(fit) || all(fit==0) || all(isnan(fit))
        h = NaN; t = h; w = h; w_times = [NaN NaN]; halfh = NaN;
        warning('Empty, nan, or all zero values!')
        return
    end

    % exit if all vals are same
    if all(fit == mean(fit)), h = NaN; t = h; w = h; w_times = [NaN NaN]; halfh = NaN;
        warning('All fitted values are the same!')
        return
    end

    dat_for_amp_est = fit(startval_samples:endval_samples);

    % --------------------------------------------------------
    % find turning points in data (where approx. derivative is zero)
    % --------------------------------------------------------
    tmp = diff(dat_for_amp_est);
    turnpts = find([0; diff(sign(tmp))]);

    if isempty(turnpts)
        % monotonic increasing or decreasing function, can't define
        turnpts = [find(dat_for_amp_est == min(dat_for_amp_est)) find(dat_for_amp_est == max(dat_for_amp_est))];
    end


    tpdat = abs(dat_for_amp_est(turnpts));

    t = turnpts(tpdat == max(tpdat));       % time of max absolute turning point
    t = t(1);
    h = dat_for_amp_est(t);                             % the value at max

    % turnpts is in samples from startval_samples on; add to adjust to
    % samples from time zero
    t = t + startval_samples - 1;

    % --------------------------------------------------------
    % calculations for width
    % --------------------------------------------------------
    halfh = .5 * h;                         % 1/2 the max to min distance


    % First time point before t that crosses half-height (increasing for
    % activations, decreasing for deactivations)
    % ------------------------------------------
    crossings = diff(fit > halfh);

    eligible_times = find(crossings(1:t) == sign(h));
    if isempty(eligible_times)
        w_times(1) = NaN;
    else
        w_times(1) = eligible_times(end);

        % interpolate to get more accurate estimate (relevant for low sampling res)
        y1 = fit(w_times(1));
        y2 = fit(w_times(1) + 1);
        m = y2 - y1;          % slope, x reduces to 1 (x2 - x1 = 1)

        w_times(1) = w_times(1) + (halfh - y1) ./ m; % solve y = mx + b for x* given m; y1 = b

    end


    % First time point after t that crosses half-height (decreasing for activations, increasing for deactivations)
    % ------------------------------------------
    eligible_times = find(crossings(t:end) == -sign(h));
    if isempty(eligible_times)
        w_times(1) = NaN;
    else
        w_times(2) = t + eligible_times(1);


        % interpolate to get more accurate estimate (relevant for low sampling
        % res)
        y1 = fit(w_times(2));
        y2 = fit(w_times(2) + 1);
        m = y2 - y1;                      % slope, x reduces to 1 (x2 - x1 = 1)
        w_times(2) = w_times(2) + (halfh - y1) ./ m;     % solve y = mx + b for x* given m; y1 = b
    end

    w = w_times(2) - w_times(1);        % width in elements

    
    
    % Finish up and plot
    % --------------------------------------------------------

    % Add area under curve measure
    auc = sum(dat_for_amp_est);


    % re-express times in sampling units (seconds or TRs, depending on what dt
    % represents)

    t = t .* dt;
    w_times = w_times * dt;
    w = w * dt;
    
    
    if doplot
        minh = 0; %mean(fit(1:2));

        create_figure('HTW Estimates', 1, 2);
        plot(hrf); title('Basis set');

        ylims = get(gca, 'YLim');
        h1 = drawbox(startval_samples, endval_samples - startval_samples, ylims(1), ylims(2) - ylims(1), [.9 .9 .9]);
        set(h1, 'EdgeColor', 'none');
        text(startval_samples, ylims(2) - .05 * (ylims(2) - ylims(1)), 'Eligible area for peak', 'FontSize', 16);

        plot(hrf, 'LineWidth', 2);
        xlabel('Samples');

        
        secs = (1:rows) .* dt;
        
        subplot(1, 2, 2);
        plot(secs, fit, 'k', 'LineWidth', 3);
        
        ylims = get(gca, 'YLim');
        h1 = drawbox(startval_samples.* dt, endval_samples.* dt - startval_samples.* dt, ylims(1), ylims(2) - ylims(1), [.9 .9 .9]);
        set(h1, 'EdgeColor', 'none');
        
        plot(secs, fit, 'k', 'LineWidth', 3);
        plot([startval_samples : endval_samples] * dt, dat_for_amp_est, 'g', 'LineWidth', 1);

        title('Fitted response and summary statistics')
        
        try

            hold on
            h1 = arrow([t minh],[t h+minh],'Length',12);        % height arrow
            h2 = arrow([0 h+minh],[t h+minh],'Length',12);      % delay arrow
            h3 = arrow([w_times(1) halfh+minh],[w_times(2) halfh+minh],'Length',12);
            h4 = arrow([w_times(2) halfh+minh],[w_times(1) halfh+minh],'Length',12);

            h5 = text(t + .1 * w, halfh+minh + .4*(halfh+minh),'h','FontSize',16,'FontWeight','bold');
            h6 = text(t - .5 * w, h+minh - .2*(halfh+minh),'t','FontSize',16,'FontWeight','bold');
            h7 = text(t - .5 * w, halfh+minh - .2*(halfh+minh),'w','FontSize',16,'FontWeight','bold');

            set(h1,'Color',colors{1}(1))
            set(h2,'Color',colors{1}(1))
            set(h3,'Color',colors{1}(1))
            set(h4,'Color',colors{1}(1))
            set(h5,'Color',colors{1}(1))
            set(h6,'Color',colors{1}(1))
            set(h7,'Color',colors{1}(1))

        catch
            warning('Error drawing arrows!')
        end
    end

    xlabel('Sec (Samples * dt)');

    drawnow



end


% OLD WAY OF GETTING W
% %     % first half
% %     x2 = max(wh(1),1);    % first above halfh
% %     x1 = max(wh(1)-1,1);  % first above-half - 1
% %     y2 = dat_for_amp_est(x2);
% %     y1 = dat_for_amp_est(x1);
% %     m = y2 - y1;          % slope, x reduces to 1 (x2 - x1 = 1)
% %     if m == 0
% %         w_times(1) = x1;    % exact match
% %     else
% %         w_times(1) = x1 + (halfh - y1) ./ m; % solve y = mx + b for x* given m; y1 = b
% %     end
% %
% %     % second half
% %     x1 = min(wh(end),hconstraint);      % last one above halfh
% %     x2 = min(wh(end)+1,hconstraint);    % + 1
% %     y2 = dat_for_amp_est(x2);  y1 = dat_for_amp_est(x1);
% %     m = y2 - y1;                      % slope, x reduces to 1 (x2 - x1 = 1)
% %     if m == 0
% %         w_times(2) = x1;    % exact match
% %     else
% %         w_times(2) = x1 + (halfh - y1) ./ m;     % solve y = mx + b for x* given m; y1 = b
% %     end
