function [CE, NSO] = plot_objectives(csv_file_or_data, varargin)
%PDIST_PLOT Given a high-dimensional dataset, plot_objectives calculates
% the objective functions (cross entropy and the negative sampling
% objective) used by UMAP to evaluate the quality of a low-dimensional
% embedding for each embedding dimension between 2 and the original data
% dimension.
%
% [CE, NSO] = plot_objectives(csv_file_or_data)
%
% Required input
% ----------
% csv_file_or_data:
%   This is either 
%   A) a char array identifying a CSV text file containing the data 
%      to be reduced. 
%   B) the actual data to be reduced; a numeric matrix.
%
% Optional input
% ----------
% 'verbose': char array
%     controls verbosity of output in the same manner as the argument for
%     run_umap.m.
%
% Returns
% -------
% CE: double array of size (1, maxD-1)
%     The cross entropy of the final UMAP embedding for each embedding
%     dimension from 2 to maxD.
%
% NSO: double array of size (1, maxD-1)
%     The value of the negative sampling objective function of the final
%     UMAP embedding for each embedding dimension from 2 to maxD.
%
%   AUTHORSHIP
%   Math Lead & Primary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Secondary Developer: Stephen Meehan <swmeehan@stanford.edu>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Provided by the Herzenberg Lab at Stanford University 
%   License: BSD 3 clause
%

MAX_DIMENSIONS = 50;
MAX_SAMPLE_SIZE = 1e4;

p=parseArguments();
parse(p,varargin{:});
args=p.Results;
verbose = args.verbose;

if ischar(csv_file_or_data)
    if ~exist(csv_file_or_data, 'file')
        csv_file_or_data = UmapUtil.RelocateExamples(csv_file_or_data);
    end
    t=readtable(csv_file_or_data, 'ReadVariableNames', true);
    inData=table2array(t);
else
    inData = csv_file_or_data;
end

if size(inData, 1) > MAX_SAMPLE_SIZE
    warning(['Looks like your dataset is too large. Trimming it down to the first ' num2str(MAX_SAMPLE_SIZE) 'rows...']);
    inData = inData(1:MAX_SAMPLE_SIZE,:);
end

maxD = min(MAX_DIMENSIONS, size(inData, 2));

CE = zeros(1, maxD-1);
NSO = zeros(1, maxD-1);

for d = 2:maxD
    [reduction,umap] = run_umap(inData,'n_components', d, 'verbose', verbose);
    
    CE(d-1) = cross_entropy(reduction, reduction, umap.head, umap.tail, umap.weights, umap.a, umap.b, true);
    NSO(d-1) = neg_sampling_objective(reduction, reduction, umap.head, umap.tail, umap.weights, umap.a, umap.b, umap.negative_sample_rate, true);
end

Cfig = figure('Name', 'Cross Entropy vs. Embedding Dimension');
Cax=axes('Parent', Cfig);
plot(Cax, 2:maxD, CE);
xlabel(Cax, 'Dimension');
ylabel(Cax, 'Cross Entropy');

Nfig = figure('Name', 'Negative Sampling Objective vs. Embedding Dimension');
Nax=axes('Parent', Nfig);
plot(Nax, 2:maxD, NSO);
xlabel(Nax, 'Dimension');
ylabel(Nax, 'Negative Sampling Objective');

    function p=parseArguments()
        p = inputParser;
        DEFAULT_VERBOSE = 'none';
        EXPECTED_VERBOSE = {'graphic','text','none'};
        addParameter(p,'template_file',[], @(x) ischar(x) || isa(x, 'UMAP'));
        addParameter(p,'label_column',0,@(x) isnumeric(x) && x>0);
        addParameter(p,'verbose',DEFAULT_VERBOSE,...
            @(x) any(validatestring(x,EXPECTED_VERBOSE)));
    end
end