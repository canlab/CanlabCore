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
%   **doplot:**
%        visualize the plot. Default: off
% 
%   **method:**
%        Interpolation method. Default: linear
%        'linear' (default) | 'nearest' | 'next' | 'previous' | 'pchip' | 
%        'cubic' | 'v5cubic' | 'makima' | 'spline'
% 
%   **extrap:**
%        Extrapolation strategy, specified as 'extrap' or a real scalar value.
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
%    [yi, xi] = downsample_scnlab(y, 100, .5, 'doplot', 'method', 'spline')
%
% Programmers notes:
% 8/25/2023 - Byeol
% Add more method options for interp1 function.


doplot = 0;
method = 'linear';
extrapolation = 'extrap';

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            % functional commands
            case {'doplot'}
                doplot = 1;
            case {'method'}
                method = varargin{i+1};
            case {'extrap'}
                extrapolation = varargin{i+1};
        end
    end
end

downsamplerate = orig_samprate ./ new_samprate;

nobs = length(y);


% observations in resampled output
nobs_out = round(nobs ./ downsamplerate);


x = 1:nobs;

xi = linspace(0, nobs, nobs_out);


yi = interp1(x, y, xi, method, extrapolation);

if doplot

    create_figure('Downsample plot');
    plot(y)
    hold on; plot(xi, yi, 'r')

end

end
