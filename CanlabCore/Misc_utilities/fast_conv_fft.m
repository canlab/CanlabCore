function y = fast_conv_fft(hrf,x,varargin)
% :Usage:
% ::
%
%     y = fast_conv_fft(hrf,x)
%
% Much faster than conv or using matrix multiplication.
%
% :Inputs:
%
%   **hrf:**
%        should be length of x
%
% :Examples:
% ::
%
%    y = fast_conv_fft(hrf,x);
%    y = fast_conv_fft(hrf,x3,'deconv');
%
% ..
%    tor wager, jan 07
% ..


dodeconv = 0;
if length(varargin)
    for i = 1:length(varargin)
        switch varargin{i}
            case {'deconv','deconvolve'};
                dodeconv = 1;
            otherwise
                error('Unrecognized input string.');

        end
    end
end

% make sure HRF is long enough
if length(hrf) < length(x)
    nel = length(x) - length(hrf);
    hrf = [hrf; zeros(nel,1)];
end


len = length(x);
nz = round(len ./ 2);
z = zeros(nz,1);


if dodeconv
    % deconvolve

    % x is a truncated, convolved timeseries
    % we want to deconvolve to find y given the kernel hrf
    % because x is truncated, we need to pad x with the values
    % that we would get if we convolved [x; zeros(kl,1)] with
    % the hrf.

    kl = max(find(hrf ~= 0));  % kernel length
        
    % pad with mean to avoid truncation artifact at end
    % and zeros to avoid aliasing (?) at beginning
    % this is an approx. solution, which we use in padding in next step
    m = x(end-1:-1:1); %repmat(mean(x), 10*sum(hrf ~= 0),1);
    zm = zeros(size(m));

    fftx = fft([x;m]);
    ffthrf = fft([hrf;zm]);

    y = real(ifft(fftx ./ ffthrf));
    y = y(1:len);
   
    % now we use our estimate of stim. series y to get ending values that
    % have been truncated from x
    % final values are inaccurate, so let's taper down...
    
    w = linspace(1,0,4)';
    
    % taper towards mean
    endy = y(end-kl+1:end);
    lastvals = w .* endy(end-3:end);
    endy(end-3:end) = lastvals;
    endy = [endy; zeros(kl,1)];
    
    % get what truncated vals should look like
    endvals = fast_conv_fft(hrf(1:kl),endy);
    extra = endvals(end-kl+1:end);
    
    % re-do deconv
    zm = zeros(size(extra));
    
    fftx = fft([x;extra]);
    ffthrf = fft([hrf;zm]);

    y = real(ifft(fftx ./ ffthrf));

    % not needed...the above is pretty good
    % taper because we still have the problem
% %     y = y(1:len);
% %     
% %     
% %     endvals = w .* y(end-3:end) + (1-w) .* mean(y(1:end-4));
% %     y(end-3:end) = endvals;
    
else
    % convolve
% pad with zeros to avoid edge effects/aliasing
fftx = fft([x;z]);
ffthrf = fft([hrf;z]);

    y = real(ifft(fftx .* ffthrf));
end

y = y(1:len);

end
