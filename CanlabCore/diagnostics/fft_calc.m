function [myfft, freq] = fft_calc(dat, TR)
% Simple function to calculate the FFT power of a 
% data vector (dat) as a function of frequency,
% given a sample-to-sample repetition time (TR)
%
% :Usage:
% ::
%
%     [myfft, freq] = fft_calc(dat, TR)
%

n = length(dat);

% from matlab example
% power = abs(Y(1:floor(n/2))).^2;
% nyquist = 1/2;
% freq = (1:n/2)/(n/2)*nyquist

nyq = 1 ./ (2 * TR);

timepts = floor(n ./ 2);

freq = (0:timepts-1)/timepts * nyq;
%freq = linspace(0, nyq, timepts)';

myfft = fft(dat); %real(abs(fft(dat)));
myfft = abs(myfft(1:timepts)) .^ 2;  % power

myfft = myfft ./ sum(myfft);
