function [yi, xi] = downsample_scnlab(y, orig_samprate, new_samprate, varargin)
% Uses linear interpolation to resample a vector from one sampling
% rate to another
%
% :Usage:
% ::
%
%     [yi, xi] = downsample_scnlab(y, orig_samprate, new_samprate, [doplot])
%
% :Inputs:
%
%   **orig_samprate:**
%        sampling rate in Hz
%
%   **new_samprate:**
%        desired sampling rate in Hz
%
% ..
%    I prefer this to matlab's downsample.m because linear interpolation is
%    robust and does fewer weird things to the data.
%
%     Tor Wager, June 2009
% ..
%
% :Example:
% ::
%
%    % Downsample a 100 Hz signal to a scanning TR of 2 sec
%    % signal at 100 Hz, sample to low-freq TR of 0.5 hz (2 sec TR)
%    % every 100 / TR = 100/.5 = 200 samples
%
%    [yi, xi] = downsample_scnlab(y, 100, .5)
%

doplot = 0;
if length(varargin) > 0, doplot = varargin{1}; end

downsamplerate = orig_samprate ./ new_samprate;

nobs = length(y);


% observations in resampled output
nobs_out = round(nobs ./ downsamplerate);


x = 1:nobs;

xi = linspace(0, nobs, nobs_out);


yi = interp1(x, y, xi, 'linear', 'extrap');

if doplot

    create_figure('Downsample plot');
    plot(y)
    hold on; plot(xi, yi, 'r')

end

end
