function data = perturb_classes(data,labels,SDRatio)
%PERTURB_CLASSES Given a set of data, labels for the data, and a ratio,
% perturb each class of data independently by a vector, each component of
% which is sampled from a normal random variable with mean 0 and standard
% deviation given by the ratio times the full range of that component
% across the whole dataset.
%
% data = perturb_classes(data, labels, SDRatio)
%
% Parameters
% ----------
% data: double array of size (R, C)
%     A set of R datapoints, each of which has C measurements.
%
% labels: integer array of size (R, 1)
%     A column of class labels for the data.
%
% SDRatio: double (optional, default 0.05)
%     The proportion of the range the standard deviation of each random
%     variable will be.
%
% Returns
% -------
% data: double array of size (R, C)
%     A new dataset with each class independently perturbed.
%
%   AUTHORSHIP
%   Math Lead & Primary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Provided by the Herzenberg Lab at Stanford University 
%   License: BSD 3 clause
if nargin < 3
    SDRatio = 0.05;
end

range = max(data) - min(data);
SD = SDRatio*range;
subsetIds=unique(labels);
nIds = length(subsetIds);

for i = 1:nIds
    id = subsetIds(i);
    
    perturbation = normrnd(0,SD);
    data(labels==id,:) = data(labels==id,:) + perturbation;
end

end

