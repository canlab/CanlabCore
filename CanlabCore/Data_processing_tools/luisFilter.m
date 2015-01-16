function outdata = luisFilter(data,  TR, cutoff, verbose)
% function outdata = luisFilter(data,  TR, cutoff [, verbose])
%
% Luis Hernandez. University of Michigan.  Last Edit 9/13/01
%
% This is a low pass FIR filter using the Parks Mclellan design algorithm
% The phase introduced by the filter is linear and is un-done by 
% filtering the data again, backwards
% if you want to see what it does to the data, use verbose=1
% if just want to filter the data, don't use the argument at all.
%
   % close all
 
    %%%%%%  These are the important lines of the code   %%%%%%%%%%%%%
    cutoff = cutoff / (2*TR)
    b = remez(10, [0 cutoff-0.05  cutoff+0.05  1], [1 1 0 0]);
    outdata = filtfilt(b,1, data);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
if (nargin==4)
    freqz(b)
    figure
 
    % Read the data from file
%    data = load(data_file);
%    data = data(:,2);
    len = size(data, 1);
    time = [1:len] * TR;
    subplot(2,1,1) , plot(time, abs(data));
    xlabel('sec.')
    
    % Fourier Transform the data
    fdata = fftshift(fft(data));   
    scalefactor = max(abs(fdata)); 
    
    %frequency range:
    f = [1:len]';
    f = (f - len/2 ) * 1/(TR * len);
    whos
    % Plotting the data before the filter:
     
    subplot(2,1,1) , hold on,plot(time, abs(data))
    subplot(2,1,2), plot(f, abs(fdata));
    subplot(2,1,2), axis([0  1/(2*TR)  0  scalefactor]);
    xlabel('Hz.');
      
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Make  Fourier Domain Filter by hand
    %ffilt = ones(len, 1);
    %co  = len/2  +  cutoff * (TR * len);
    %ffilt(co:end) = exp(-((f(co:end) - cutoff).^2) / 0.001 );
    %co  = len/2  - cutoff * (TR * len);
    %ffilt(1:co)   = exp( -((f(1:co) + cutoff).^2) / 0.001 );
        
    %subplot(2,1,2) , hold on , plot(f, ffilt*scalefactor/2 ,'g');
    %whos
 
    % Apply the filter to the data and IFT it
    %fdata  = fdata .* ffilt;
    %data = ifft(fdata);
   
    
    
    
     % Fourier Transform the data
    outfdata = fftshift(fft(outdata));
    
    
    % Let's look at the impulse response of the filter:
    impulse= zeros(size(data));
    
    impulse(size(data,1)/2) = 0.1 *scalefactor;
    impulse_response = filtfilt(b,1,impulse);
    
    freq_response = fftshift(fft(impulse_response));
    freq_impulse = fftshift(fft(impulse));
    
    freq_gain = scalefactor/2 + 20*log10( abs(freq_response).^2 ./ abs(freq_impulse) ) ;
    
    %Plotting the reponse of the filter
    subplot(2,1,1) , hold on,plot(time, abs(outdata), 'r'), title ('Time Domain'), legend('Before', 'After')
    subplot(2,1,2), hold on, plot(f, abs(outfdata), 'r'), title ('Frequency Domain');
    
    subplot(2,1,2), hold on, plot(f, freq_gain, 'g'),legend('Before','After','Filter Response');
 
    %figure;
    %plot(f, freq_gain, 'g');
end
 
return
