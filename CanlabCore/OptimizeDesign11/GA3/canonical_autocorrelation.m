function [xc,Vi] = canonical_autocorrelation(TR,n,varargin)
% [xc,Vi] = canonical_autocorrelation(TR,n,[a])
% Tor Wager, 2 / 23 / 04
%
% exponential autocorrelation function of the form xc = e ^ (-ax)
% Optional input "a" is an exponent; higher = more white noise, low =
% colored
%
% TR = sampling rate (time per data point)
% n = number of samples in timeseries
%
% xc = vector containing autocorrelation coefficients
% Vi = autocorrelation matrix

f = inline('A * exp(-a * x)','A','a','x');  % exponential function
xa = 0:TR:60;

if length(varargin) > 0
    xc = f(1,varargin{1},xa); 
else
    xc = f(.9475,.5323,xa);    % params at 3T from vnl experiment, n = 10
    xc = xc ./ max(xc);
end

% autocorrelation
Vi = getV('make',xc,n);

return
