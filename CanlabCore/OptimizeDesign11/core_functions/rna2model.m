function [model,hires_model,delta] = rna2model(stimList,ISI,HRF,TR,numsamps,nonlinthreshold,S)   
% function [model,hires_model] = rna2model(stimList,ISI,HRF,TR,numsamps,dosaturation,S)
%
% transforms a design vector (final "expressed" form) to model matrix
% hi-res model at .1 s sampling rate
%
% Tor Wager, 11.17.01

model = sampleInSeconds(stimList,ISI);
[model,delta] = getPredictors(model,HRF);
hires_model = model;

model = resample(model,10,TR*100);

if ~isempty(nonlinthreshold),model = modelSaturation(model,nonlinthreshold);,end		% saturation (nonlinear responses)
        
if size(model,1) > numsamps, model = model(1:numsamps,:);,end
        
if ~isempty(S)
		model = S * model;                                                  			% temporal smoothing and HP/LP filter
end
	
model(:,end+1) = 1;																		% add intercept
hires_model(:,end+1) = 1;																% add intercept

return