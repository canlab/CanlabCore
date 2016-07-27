function [model,hires_model,delta] = designvector2model(stimList,ISI,HRF,TR,numsamps,nonlinthreshold,S,varargin)   
% function [model,hires_model,delta] = designvector2model(stimList,ISI,HRF,TR,numsamps,dosaturation,S,[opt] hox)
%
% transforms a design vector (final "expressed" form) to model matrix
% hi-res model at .1 s sampling rate
%
% Tor Wager, 11.17.01
%
% Functions called:
% sampleInSeconds.m
% getPredictors.m
% resample.m    (toolbox; free version available from Sourceforge.net)
% modelSaturation.m
% 

% if using variable ISI or TR using hox genes,
% this segment gets the new parameters.
% ------------------------------------------------------
if ~isempty(varargin)    % do hox
    % compute new number of samples, and extend stimList if necessary
    numSecs = numsamps * TR;
    
    hox = varargin{1};
    if hox(1) > 0, ISI = hox(1); end
    if length(hox) > 1, if hox(2) > 0, TR = hox(2); end, end
    
    numStim = numsamps * TR / ISI;
    while numSecs > size(stimList,1)*ISI, stimList = [stimList; stimList];  end
    stimList = stimList(1:ceil(numStim));
    
    % does not work yet with variable TR.
end



model = sampleInSeconds(stimList,ISI);
[model,delta] = getPredictors(model,HRF);
hires_model = model;

model = resample(model,10,TR*100);

if ~isempty(nonlinthreshold),model = modelSaturation(model,nonlinthreshold); end		% saturation (nonlinear responses)
    
if size(model,1) < numsamps, model = [model; model(end, :)]; end
if size(model,1) > numsamps, model = model(1:numsamps,:); end
if size(delta,1) > numsamps*TR*10, delta = delta(1:numsamps*TR*10,:); end  

if ~isempty(S)
		model = S * model;                                                  			% temporal smoothing and HP/LP filter
end
	
model(:,end+1) = 1;																		% add intercept

return