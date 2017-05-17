function [h,t,w,w_times,halfh, auc] = fir2htw2(b,varargin)
% Estimates height, time to peak, and width of FIR response
%
% :Usage:
% ::
%
%     [h,t,w,w_times,halfh, auc] = fir2htw2(b,[hconstraint],[doplot],[colors cell])
%
% :Inputs:
%
%   **b:**
%        beta/estimate series for hemodynamic response curve
%
%   hconstraint:**
%        max time in samples that can be considered the peak
%        (default = last sample)
%
%   doplot:**
%        flag for plot, 1/0
%
%   colors:**
%        cell vector of colors for plot
%
% ..
%    tor wager, 2/14/05
% ..
%
% minh = min height
%
% This version uses turning points (zero gradient) to find the largest
% "hump" in the data and the time it occurs.
%
% :Example:
% ::
%
%    hrf = spm_hrf(.5); hrf = hrf ./ max(hrf); hrf = hrf + .1 * randn(length(hrf), 1);
%    create_figure('hrf'); plot(hrf);
%    [h,t,w,w_times,halfh, auc] = fir2htw2(hrf, [], 1);
%

    if size(b,2) < length(b), b = b'; end

    if ~isempty(varargin) && ~isempty(varargin{1})
        hconstraint = varargin{1};
    else
        hconstraint = length(b);
    end

    if hconstraint > length(b)
        warning('fir2htw is trying to use more betas than there are!  limiting.')
        hconstraint = length(b);
    end

    if length(varargin) > 1,
        doplot = varargin{2};
    else
        doplot = 1;
    end


    % exit if empty
    if isempty(b) || all(b==0) || all(isnan(b))
        h = NaN; t = h; w = h; w_times = [NaN NaN]; halfh = NaN;
        warning('fir2htw:  Empty, nan, or all zero values!')
        return
    end

    % exit if all vals are same
    if all(b == mean(b)), h = NaN; t = h; w = h; w_times = [NaN NaN]; halfh = NaN;
        warning('fir2htw:  all values are the same!')
        return
    end



    colors = {'ro-'};
    if length(varargin) > 2, colors = varargin{3}; end


    % find turning points in data
    tmp = diff(b);
    turnpts = find([0 diff(sign(tmp))]);

    if isempty(turnpts)
        % monotonic increasing or decreasing function, can't define
        turnpts = [find(b==min(b)) find(b==max(b))];
    end



    dat = b(1:hconstraint) - mean(b(1:2));   %- b(1); % ; - mean(b(1:2));          % deviations from first point (baseline)

    turnpts(find(turnpts > length(dat))) = [];  % eliminate extra turn points after hconst

    if isempty(turnpts)
        % monotonic increasing or decreasing function, can't define
        turnpts = [find(dat==min(dat)) find(dat==max(dat))];
    end

    tpdat = abs(dat(turnpts));

    t = turnpts(tpdat == max(tpdat));       % time of max absolute turning point
    t = t(1);
    h = dat(t);                             % the value at max


    % calculations for width
    halfh = .5 * h;                         % 1/2 the max to min distance

    %d = distance(halfh,dat);

    if h > 0
        wh = find(dat > halfh);
    else
        wh = find(dat < halfh);
    end

    if isempty(wh), wh = NaN; end

    % first half
    x2 = max(wh(1),1);    % first above halfh
    x1 = max(wh(1)-1,1);  % first above-half - 1
    y2 = dat(x2);
    y1 = dat(x1);
    m = y2 - y1;          % slope, x reduces to 1 (x2 - x1 = 1)
    if m == 0
        w_times(1) = x1;    % exact match
    else
        w_times(1) = x1 + (halfh - y1) ./ m; % solve y = mx + b for x* given m; y1 = b
    end

    % second half
    x1 = min(wh(end),hconstraint);      % last one above halfh
    x2 = min(wh(end)+1,hconstraint);    % + 1
    y2 = dat(x2);  y1 = dat(x1);
    m = y2 - y1;                      % slope, x reduces to 1 (x2 - x1 = 1)
    if m == 0
        w_times(2) = x1;    % exact match
    else
        w_times(2) = x1 + (halfh - y1) ./ m;     % solve y = mx + b for x* given m; y1 = b
    end


    w = w_times(2) - w_times(1);        % width in elements

    % Add area under curve measure
    if nargout > 5, auc = sum(dat); end


    if doplot
        minh = mean(b(1:2));
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



end


% OLD EXTRA STUFF

% % range = max(b(1:hconstraint)) - min(b(1:hconstraint));
% %
% % % find all values within some tolerance of the half-max height
% % % if not enough points found to determine width, then interpolate until we
% % % get values.
% % % it would be better to interpolate!
% % tolval = .02; wh = []; resampleval = 1;
% %
% % tol = tolval * range;   % tolerance of tolval % of range
% % wh = find(d <= tol);    % find all points within tolerance
% %
% % while length(find(wh>t)) < 1 | length(find(wh<t)) < 1
% %     resampleval = resampleval + 1;
% %
% %     b2 = resample(b(1:hconstraint),resampleval,1); % interpolate
% %
% %     d = distance(halfh,b2);
% %     wh = find(d <= tol);    % find all points within tolerance
% %     wh = wh ./ resampleval;  % convert back to original time units
% %
% %     if resampleval > 10, warning('Could not find width!'), break, end
% % end
% %
% % % now find nearest elements to peak that are at max height
% % d_from_t = wh - t;
% % w_times = [max(d_from_t(d_from_t < 0)) min(d_from_t(d_from_t > 0))];
% %
% %
% %
% % if length(w_times) > 1
% %     w = w_times(2) - w_times(1);        % width in elements
% %     w_times = t + w_times;              % time indices of elements before and after peak at 1/2 max height
% % else
% %     w_times = [NaN NaN];
% %     w = NaN;                            % width is undefined; 1st or last tp was max
% % end


