function [mag, freq, h, h2, peak_freq] = fft_calc(dat, TR, varargin)
% Simple function to calculate the FFT power of a 
% data vector (dat) as a function of frequency,
% given a sample-to-sample repetition time (TR)
%
% :Usage:
% ::
%
%     [mag, freq, line_handle, nyquist_line_handle, peak_freq] = fft_calc(dat, TR)
%
% Example:
% Create and plot a sin wave:
% ------------------------------------------------------------------
% Fs = 1000;            % Sampling frequency, e.g., Hz (samples/sec)                    
% T = 1/Fs;             % Sampling period (time between samples, e.g., sec)
% L = 1000;             % Length of signal (in samples)
% t = (0:L-1)*T;        % Time vector
% 
% % Sine wave parameters
% theta = 0;            % Phase, in radians
% F = 10;               % Frequency of sin wave (cycles/sec)
% 
% y = sin(2 * pi * F * t + theta);
% 
% create_figure('sin wave');
% plot(t, y)
% xlabel('Time (sec)')
%
% Downsample the sin wave at new frequency Fs and plot signal and FFT:
% ------------------------------------------------------------------
% Fs = 20;  % New sampling frequency
% downsample_by = round(1000/Fs);
% 
% figure;
% plot(t, y, 'k-'); hold on
% set(gcf, 'Position', [0 644        1358         231])
% 
% plot(t(1:downsample_by:end), y(1:downsample_by:end), '.-', 'Color', [1 .3 .6], 'LineWidth', 2);
% Nyquist = Fs/2;
% title(sprintf('Sampling rate = %3.1f, Nyquist limit = %3.2f', Fs, Nyquist))
% 
% figure; hold on;
% [mag, freq, line_han] = fft_calc(y, 1/1000);
% set(line_han, 'Color', 'k', 'LineWidth', 2, 'LineStyle', '-');
% 
% [mag, freq, line_han] = fft_calc(y(1:downsample_by:end), 1/Fs);
% set(line_han, 'Color', [1 .3 .6], 'LineWidth', 2, 'LineStyle', '-');
% 
% set(gca, 'XLim', [0 (max(Fs, 1.2 * F))]);



n = length(dat);

if ~isempty(varargin)
    if ischar(varargin{1}) && strcmp(varargin{1}, 'plot')
    else
        doplot = varargin{1};
    end
end

% from matlab example
% power = abs(Y(1:floor(n/2))).^2;
% nyquist = 1/2;
% freq = (1:n/2)/(n/2)*nyquist

nyq = 1 ./ (2 * TR);  % Sampling rate frequency / 2. TR = sample-to-sample time, 1/Fs

timepts = floor(n ./ 2);

freq = (0:timepts-1)/timepts * nyq;
%freq = linspace(0, nyq, timepts)';

mag = fft(dat); %real(abs(fft(dat)));
mag = abs(mag(1:timepts)) .^ 2;  % power

mag = mag ./ sum(mag);

h = plot(freq, mag, 'LineWidth', 2); 
xlabel('Frequency (Hz)')
title('Frequency domain')

% nyquist = Fs ./ 2;  % Any signal faster than Nyquist freq will be aliased
hold on
h2 = plot_vertical_line(nyq);
set(h2, 'LineStyle', ':')

peak_freq = freq(mag == max(mag));

end

