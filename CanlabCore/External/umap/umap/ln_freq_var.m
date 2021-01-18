function elimIdxs = ln_freq_var(labels,mean)
%LN_FREQ_VAR Given a set of class labels and a desired mean, multiply the
% incidence of each class label by a frequency sampled from a random
% variable with a log-normal distribution and the given mean (frequency
% capped at 1).
%
% elimIdxs = ln_freq_var(labels,mean)
%
% Parameters
% ----------
% labels: integer vector
%     A set of class labels
%
% mean: double (optional, default 0.5)
%     The mean of the variable with log-normal distribution that determines
%     the frequency of each class
%
% Returns
% -------
% elimIdxs: logical vector, same size as "labels"
%     A logical array indicating the indices of the labels to be
%     eliminated.
%
%   AUTHORSHIP
%   Math Lead & Primary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Provided by the Herzenberg Lab at Stanford University 
%   License: BSD 3 clause

if nargin < 2
    normMean = -log(2)-1/2; %lognrnd(normMean, 1) has mean 1/2.
else
    normMean = log(mean)-1/2; %lognrnd(normMean, 1) has mean "mean".
end

subsetIds=unique(labels);
nIds = length(subsetIds);
elimIdxs=false(size(labels));

for i = 1:nIds
    id = subsetIds(i);
    
    freq = min(lognrnd(normMean, 1), 1);
    idxs= labels==id;
    nIdxs=nnz(idxs);
    nEliminating=floor((1-freq)*nIdxs);
    if nEliminating > 0
        elimIdxs(find(idxs,nEliminating))=true;
    end
end
end

