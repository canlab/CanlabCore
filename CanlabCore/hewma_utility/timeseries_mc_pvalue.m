function [Zcor,p,tthresh,q,C] = timeseries_mc_pvalue(varZm,lambda,Wm,sesq,df,tm)
% Uses Monte Carlo simulation on a timeeseries to get the null hypothesis
% max t-value over time and establish an appropriate statistical threshold
%
% :Usage:
% ::
%
%     [Zcor,p,tthresh,q,C] = timeseries_mc_pvalue(varZm,lambda,Wm,sesq,df,tm)
%
% :Inputs:
%
%   **varZm:**
%        within subjects variances, subj x time
%
%
%   **lambda:**
%        smoothing parameter in EWMA
%
%   **Wm:**
%        individual case (subject) weights
%
%   **sesq:**
%        squared standard error between subjects (variance of group error estimate)
%
%   **df:**
%        degrees of freedom between subjects (est.)
%
%   **tm:**
%        max t stat for group over time
%
% ..
%    NOTE: TOR CHANGED INPUT TO ASSUME THAT WE SHOULD ENTER VARWI + VARBETWEEN
% ..

C = get_correlationmat(varZm,lambda,Wm,sesq); % Calculate correlation matrix


%%%%% Calculate p-value corresponding to max t-statistic


rep = 5000;                                          % Number of repetitions in Monte-Carlo integration 



q = mvtrnd(C,df,rep);                                 % Simulate random variables according to multivariate T-distribution

q = max(abs(q)')';                                   % max t-value across t-timeseries, + or - (2-tailed)

p = sum(q >= abs(tm)) ./ rep;                        % p of getting tm or greater max t-value by chance

%Zcor = (tm - mean(q)) ./ std(q);                     % convert to z-score corrected for mult. comps.
Zcor = sign(tm) .* norminv(1-p);

%ind = sum((q > repmat(tm,rep,T)),2)/T;                % Count number of dimensions for which (q_i > thres) for i=1,...T
%p = sum((ind > 0))/rep;                               % Probability that (q_i > thres) atleast one time

try
    tthresh = prctile(q,95);                         % 2-tailed because q is 2-sided maximum
catch
    disp('Error getting t-threshold: No stats toolbox?');
    tthresh = [];
end
