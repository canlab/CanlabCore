function [synData,synLabels] = generate_synthetic_set(data,labels,varargin)
%
%   AUTHORSHIP
%   Math Lead & Primary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Provided by the Herzenberg Lab at Stanford University 
%   License: BSD 3 clause
WARNING_THRESHOLD = 100;

p=parseArguments();
parse(p,varargin{:});
args=p.Results;
randomize = args.randomize;
exact_proportion = args.exact_proportion;
sz = args.size;
targetProps = args.targetProps;

[R,C] = size(data);
if size(labels, 1) ~= R
    error(['Number of labels must equal number of data points! There were ' num2str(R) ' data points entered and ' num2str(size(labels, 1)) ' labels entered!']);
end

subsetIds=unique(labels);
nIds = length(subsetIds);
idxs=false(R,nIds);
nIdxs = zeros(1, nIds);
means = zeros(nIds,C);
covariances=zeros(C,C,nIds);

for i = 1:nIds
    id = subsetIds(i);
    idxs(:,i)= labels==id;
    nIdxs(i)=nnz(idxs(:,i));
    means(i,:) = mean(data(idxs(:,i),:));
    covariances(:,:,i) = cov(data(idxs(:,i),:));
end

if ~randomize
    rng default;
end

if isempty(targetProps)
    targetProps = nIdxs/R;
else
    targetProps = targetProps/sum(targetProps);
end

if exact_proportion
    synSize = round(targetProps*sz);
    sz = sum(synSize);
else
    synSize = mnrnd(sz, targetProps);
end

if any(synSize < WARNING_THRESHOLD)
    warning(['At least one of your classes in the synthetic dataset will have fewer than ' num2str(WARNING_THRESHOLD) ' data points. Proceed with caution']);
end

synData = zeros(sz, C);
synLabels = zeros(sz, 1);
partialSizes=cumsum(synSize);


synData(1:partialSizes(1),:) = mvnrnd(means(1,:), covariances(:,:,1),synSize(1));
synLabels(1:partialSizes(1)) = subsetIds(1);

for i = 1:nIds-1
    synData((partialSizes(i)+1):partialSizes(i+1), :) = mvnrnd(means(i+1,:), covariances(:,:,i+1),synSize(i+1));
    synLabels((partialSizes(i)+1):partialSizes(i+1)) = subsetIds(i+1);
end
    

function p=parseArguments()
        p = inputParser;
        addParameter(p,'randomize',false,@(x) islogical(x));
        addParameter(p,'exact_proportion',true,@(x) islogical(x));
        addParameter(p, 'size', 1e5, @(x) isnumeric(x) && x > 0 );
        addParameter(p, 'targetProps', [], @(x) isnumeric(x) && all(x >= 0 ) && all (x <= 1));
end

end

