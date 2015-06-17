function noise = noisevector(vlength,xcfunction,noisevariance)
% function noise = noisevector(vlength,xcfunction,noisevariance)
%
% returns a noise vector of 'vlength' samples 
% with characteristics of the autocorrelation function that you specify
%   (must be same sampling rate as timeseries!)
%
% the vector is normalized such that the variance is what you specify
% and the mean is zero
%
% vlength: vector length
% xcfunction: any cross-correlation (i.e., noise autocorrelation) function
%   - also called a periodogram
%   - a 1/f function is generally good for fMRI data
% noisevariance: desired variance of noise
%
% 2/11/01 Tor Wager

if isempty(xcfunction), xcfunction = [1 0];,end

% define series of random "shocks".  Last one is time 1, paramount one is time 2, etc.
%       pad with zeros, because 1st shock has no history to influence it.
xclength = size(xcfunction,2);
a = randn(1,vlength); a(end + 1:end + xclength) = 0;
% the noise is the series of shocks multiplied by the autocorrelation function
% so the current value of the system = weighted sum of past shocks, up to size of autocorr function
for i = 1:vlength
    noise(i) = dot(a(end-i-xclength+1:end-i),xcfunction');
end


% make sure that the variance is one and the mean is zero

noise = noise - mean(noise);

noise = noise / sqrt(var(noise));

% now make the variance match your estimate of the variance

noise = noise * sqrt(noisevariance);

%disp(['Noise variance is ' num2str(var(noise))])
%figure;plot(noise)

return



