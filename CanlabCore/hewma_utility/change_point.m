function [cp,tm,b,t,wh,oocvec] = change_point(t,tt,Mode,tthresh)
% :Usage:
% ::
%
%     [cp,tm,b,t,wh,ooc_vector] = change_point(t,tt,Mode,tthresh)
%
%     [cp,baseline mean,max t value,time of max t] = change_point(Zm,sterr,tt,tthresh)
%
% :Inputs:
%
%   **tt:**
%        is number of baseline timepoints
%
%   **t:**
%        is group t-value timeseries
% 
%
%   **Mode: 'thresh' or 'time':**
%        if 'thresh': tthresh = threshold for significant t-values.
%
%        change_point finds the first significant supra-threshold t-value and
%        looks back in time to estimate the change point.
% 
%        if Mode = 'time'
%
%        change_point finds the estimated change-point for process determined
%        to be out of control at time tthresh
%
% :Outputs:
%
%   **wh:**
%        indices of out of control points
%
%   **ooc_vector:**
%        indicator vector for ooc points
%

    switch Mode
        case 'thresh'
            % find t-values above threshold
            % see timeseries_mc_pvalue to get threshold
        
            oocvec = abs(t) > abs(tthresh);
            wh = find(oocvec);
            wh(wh <= tt) = [];      % eliminate values in baseline period
            
            if isempty(wh),
                b = [];
            else
                b = wh(1);
            end
        
        case 'time'
            b = tthresh;
            wh = b;
            oocvec = [];
        otherwise
            error('Mode must be thresh or time');
    end
    
    % tm = tmax, b = time of max t
    %[tm, b] = max(abs(t(tt+1:end)));              % Calculate maximum absolute t-value
    % b = b+tt;                                     % max t time
    
    tm = t(b);                                      % put sign back in
    
    % zero-crossing

    if ~isempty(b)
        [a,cp] = max((sign(t(1:b)) ~= sign(tm)).*(1:b));  % change point, last zero crossing
    else
        cp = [];
    end
      
return
