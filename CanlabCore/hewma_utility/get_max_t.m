function [t,tm,b] = get_max_t(Zpop,sterr,tt)
% :Usage:
% ::
%
%     [t,tm,b] = get_max_t(Zpop,sterr,tt)
%
% :Outputs:
%
%   **t:**
%        t-value timeseries
%
%   **tm:**
%        max t-value (abs)
%
%   **b:**
%        time (index) of max t-value
%

mu = mean(Zpop(1:tt));                              % population mean
t = (Zpop - mu) ./ sterr;
[tm, b] = max(abs(t(tt+1:end)));              % Calculate maximum absolute t-value
b = b+tt;                                     % max t time
tm = t(b);                                      % put sign back in        
    
return
