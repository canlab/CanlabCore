function noise = make1overf(length, sampling_rate)

% function noise = make1overf(length, sampling_rate)

%

% returns a noise vector of 'length' smaples 

% with 1/f characteristics, sampled at 'sampling_rate' (in Hz)

% the vector is normalized such that the variance is one and the 

% mean is zero



%length = 1200;

%sampling_rate = 0.5;



% Generate 100 vectors of flat noise.

noise = 0.5 - rand(length,100);

fnoise = fft(noise);



% Make a sqrt(1/f) filter. The power spectrum is 1/f

% (Must account for sampling rate, and make it in radians)

w = [1:length/2  length/2:-1 : 1];

w = ( w * sampling_rate/length);

w = 2*pi*w;

filter = sqrt( 1./ w );



for i=1:100

   fnoise(:,i) = fnoise(:,i) .* filter';

end



fnoise  = sum(fnoise,2)/100;

noise = real(ifft(fnoise));



% make sure that the variance is one and the mean is zero

noise = noise - mean(noise);

noise = noise / sqrt(var(noise));



%subplot 311, plot(filter)

%subplot 312, plot(noise)

%subplot 313, plot(abs(fnoise))





return



