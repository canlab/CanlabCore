%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Math Lead & Secondary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Provided by the Herzenberg Lab at Stanford University 
%   License: BSD 3 clause

function [cnt, unsupervisedIdxs]=detectUnsupervised(umap, inData, ...
    sduLimit, parameterLimit)
supervisors=umap.supervisors;
subsetIds=unique(supervisors.labels);
subsetIds=subsetIds(subsetIds>0);
nSubsets=length(subsetIds);

a=zeros(nSubsets, size(inData,1));
for i=1:nSubsets
    r=umap.raw_data(supervisors.labels==subsetIds(i),:);
    means_=mean(r);
    stds_=std(r);
    B=(abs(inData-means_(1,:)))./stds_(1,:);
    a(i,:)=(sum(B>sduLimit, 2)>parameterLimit)';
end
S=sum(a);
unsupervisedIdxs=S==nSubsets;
cnt=sum(unsupervisedIdxs);
end
        