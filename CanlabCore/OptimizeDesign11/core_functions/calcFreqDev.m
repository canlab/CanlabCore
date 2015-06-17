function [maxFreqDev,omega,freqs] = calcFreqDev(stimList,conditions,freqConditions)
% function [maxFreqDev,omega,freqs] = calcFreqDev(stimList,conditions,freqConditions)
%
%
% conditions is integer vector of conditions of interest
% freqConditions is vector of propotions of all conditions

			nc = size(conditions,2);

            for i = 1:nc
         		freqs(i) = sum(stimList == conditions(i));
      	    end
			
			% not necessary if forcing sum to 1
			% freqs =  freqs ./ size(stimList,1);
			
			% force frequencies to sum to 1 b/c otherwise rest intervals will influence frequencies
            freqs = freqs ./ sum(freqs); 
			
			maxFreqDev = max(abs(freqConditions(1:nc) - freqs));
			
			omega = 1 - maxFreqDev;
			
return
			