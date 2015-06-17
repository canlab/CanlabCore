function [model,delta] = getPredictors(stimList, HRF,nsess,varargin)
% function [model,delta] = getPredictorsOddEven(stimList, HRF, nsess,[downsample factor])
% 
% Build predictors and delta functions, given a condition function or delta
% function and either a convolution matrix or vector.
%
% IMPORTANT: The intercept in this function is added.
%
% stimList: condition function or delta function
% HRF:      hemodynamic response function, or convolution matrix (columns
%           are HRF)
% [downsamp]    takes every nth element of the design matrix
% 
%
%inputs: 
% 1     a col. vector of stimulus conditions OR a delta function matrix
% 2	    an HRF vector sampled at the frequency of the stimulus vector, OR
%	    a convolution matrix H
% 3     Number of sessions (must be a divisor of the stimList length)
%
%outputs: 
% 1     a n x 2 matrix of regressors (cols) for each condition
% 2     a n x k delta matrix with onsets
%
% Tor Wager, last modified 2/22/04 to center predictors
%            modified to optionally take H convolution matrix as input
%            in place of HRF vector.  See convmtx.m (Buracas,mseq toolbox)
% 
% stimList can be condition function e.g., [1 3 2 4 3 2 1]' or 
% delta matrix (n x k), n samples and k conditions, e.g., [1 0 0 0 1 0 1]'
% 
% Resampling: the default N in matlab resample has built-in antialiasing,
% but may not be good for fmri designs!  The appropriate downsampling 
% is expected to be res*TR (res is units of samples/s), but we use 0
% because the model will depend on the analysis method used, and this is
% the most veridical approach.  With N = 0, every ith sample is used, where
% i is the downsampling factor you input.  Popular choices are 16*TR (for
% onsets2delta.m), using the SPM default res of 16.
% Delta is NOT resampled.
%
% example: TR = 2, 16 samples per second in hi-res delta dhr
% [tmp,d] = downsample_delta(dhr,16*2); X=getPredictors(d,hrf);

model = [];

% build odd/even session list
div = size(stimList,1) ./ nsess; 
if div ~= round(div), error('The stimList length must be divisible by nsess!'),end
tmp = repmat([ones(div,1); 2*ones(div,1)],ceil(nsess./2),1);
odde = tmp(1:size(stimList,1));


if min(size(stimList)) > 1 % delta matrix
    
    delta = stimList;
    odde = repmat(odde,1,size(delta,2));
    delta = [delta delta];
    odde = [odde==1 odde==2];
    delta = delta .* odde;
    
    
    if min(size(HRF)) == 1
        for i = 1:size(delta,2)
            model(:,i) = conv(delta(:,i), HRF);
        end
    end
    
else
    
    for i = 1:max(stimList(:,1)) % condition function

        delta(:,i) = (stimList == i) & odde==1;
        d2(:,i) = (stimList == i) & odde==2;
        
        if min(size(HRF)) == 1
            model(:,i) = conv(delta(:,i), HRF);
            model2(:,i) = conv(d2(:,i), HRF);
        end

    end
    
    delta = [delta d2];
    model = [model model2];
end

if min(size(HRF)) > 1 % convolution matrix
    model = HRF * delta;
end
    
model = model(1:size(stimList,1),:);      	% eliminate extra values

% downsample, if necessary
if length(varargin) > 1
    model = model(1:varargin{1}:end,:); % equivalent to resample(X,1,varargin{1},0)
end


model = model - repmat(mean(model),size(model,1),1);    % center predictors

return

