function [model,delta] = getRandfxPredictors(times, HRF)

% function [model,delta] = getPredictors(stimList, HRF)

% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ % 

% GETPREDICTORS: BUILD REGRESSORS / PREDMATRIX FOR A CONDITION LIST

% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ %

%input: a col. vector of stimulus conditions

%	    an HRF vector sampled at the frequency of the stimulus vector.

%output: a n x 2 matrix of regressors (cols) for each condition

model = [];

for i = 1:size(times,1)
   delta(:,i) = (stimList i);
   model(:,i) = conv(delta(:,i), HRF);

end

model = model(1:size(stimList),:);      	% eliminate extra values

return

