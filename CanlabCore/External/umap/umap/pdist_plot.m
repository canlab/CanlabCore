function [corr, embedding] = pdist_plot(csv_file_or_data, varargin)
%PDIST_PLOT Given a high-dimensional dataset, pdist_plot produces a density plot
% comparing pairwise distance between sample points in high dimension to
% pairwise distance between the points in a low-dimensional representation
% produced by UMAP.
%
% Ideally, one will see a high correlation between these two variables.
%
% corr = pdist_plot(csv_file_or_data)
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
% 'template_file':
%     identifies a .mat file with a saved instance of the UMAP class that
%     run_umap previously produced.
%
% 'label_column': integer
%     identifies the column in the input data matrix which contains numeric
%     identifiers to label the data for UMAP supervision mode.
%
% 'sample_size': integer
%     number of sample points to use to create the pairwise distance plot.
%     Because the computation time is O(n^2) in this input, the maximum
%     accepted value is 10,000.
% 
% 'verbose': char array
%     controls verbosity of output in the same manner as the argument for
%     run_umap.m.
%
% 'vary_dims': boolean
%     if true, pdist_plot will repeat itself by runing UMAP to produce a
%     d-dimensional embedding for all d between 2 and the dimensionality of
%     the original data.
%
% Returns
% -------
% corr: double
%     The correlation coefficient of the two pairwise distance variables:
%     high dimension and low dimension.
%
%   AUTHORSHIP
%   Math Lead & Primary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Secondary Developer: Stephen Meehan <swmeehan@stanford.edu>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Provided by the Herzenberg Lab at Stanford University 
%   License: BSD 3 clause
%

MAX_DIMENSIONS = 50;

p=parseArguments();
parse(p,varargin{:});
args=p.Results;
sample_size = args.sample_size;
label_column = args.label_column;
template_file = args.template_file;
verbose = args.verbose;
vary_dims = args.vary_dims;
supervised = label_column > 0;
template_applied = ~isempty(template_file);
if ischar(csv_file_or_data)
    csv_file_or_data=UmapUtil.RelocateExamples(csv_file_or_data);
    t=readtable(csv_file_or_data, 'ReadVariableNames', true);
    inData=table2array(t);
else
    inData = csv_file_or_data;
end

[n_rows, n_dims] = size(inData);
sample_size = min(sample_size, n_rows);

sample_inds = sort(randperm(n_rows, sample_size));

sample_hi = inData(sample_inds, :);
if template_applied
    embedding = run_umap(csv_file_or_data, 'template_file', template_file, 'verbose', verbose);
elseif supervised
    sample_hi(:, label_column) = [];
    embedding = run_umap(csv_file_or_data, 'label_column', label_column, 'verbose', verbose);
else
    embedding = run_umap(csv_file_or_data, 'verbose', verbose);
end

sample_lo = embedding(sample_inds, :);

hi_dists = pdist(sample_hi);
lo_dists = pdist(sample_lo);

data = [hi_dists' lo_dists'];

ff=figure('Name', ['2-dimensional pairwise distance vs. ' num2str(n_dims) '-dimensional pairwise distance']);
ax=axes('Parent', ff);
ProbabilityDensity2.Draw(ax,data);
xlabel(ax, [num2str(n_dims) '-dimensional pairwise distance (original data)']);
ylabel(ax, '2-dimensional pairwise distance (UMAP embedding)');
% %scatter(hi_dists, lo_dists, 1, 'filled');
% L = lsline;
% L.Color = 'r';

R_MAT = corrcoef(hi_dists, lo_dists);
corr = R_MAT(1, 2);

if ~strcmpi(verbose, 'none')
    disp(['The correlation coefficient here is ' num2str(corr) '!']);
end

if vary_dims && n_dims > 2
    maxD = min(n_dims, MAX_DIMENSIONS);
    corr_vec = zeros(1,maxD-1);
    corr_vec(1) = corr;
    embedding_array = cell(1,maxD-1);
    embedding_array{1} = embedding;
    
    for d = 3:maxD
        if template_applied
            embedding = run_umap(csv_file_or_data, 'template_file', template_file, 'verbose', verbose, 'n_components', d);
        elseif supervised
            embedding = run_umap(csv_file_or_data, 'label_column', label_column, 'verbose', verbose, 'n_components', d);
        else
            embedding = run_umap(csv_file_or_data, 'verbose', verbose, 'n_components', d);
        end
        
        embedding_array{d-1} = embedding;
    
        sample_lo = embedding(sample_inds, :);

        lo_dists = pdist(sample_lo);

        data = [hi_dists' lo_dists'];

        ff=figure('Name', [num2str(d) '-dimensional pairwise distance vs. ' num2str(n_dims) '-dimensional pairwise distance']);
        ax=axes('Parent', ff);
        ProbabilityDensity2.Draw(ax,data);
        xlabel(ax, [num2str(n_dims) '-dimensional pairwise distance (original data)']);
        ylabel(ax, [num2str(d) '-dimensional pairwise distance (UMAP embedding)']);

        R_MAT = corrcoef(hi_dists, lo_dists);
        corr_vec(d-1) = R_MAT(1, 2);

        if ~strcmpi(verbose, 'none')
            disp(['The correlation coefficient here is ' num2str(corr) '!']);
        end
    end
    Rfig = figure('Name', 'Correlation Coefficient vs. Embedding Dimension');
    ax=axes('Parent', Rfig);
    plot(ax, 2:maxD, corr_vec);
    xlabel(ax, 'Dimension');
    ylabel(ax, 'Correlation coefficient of pairwise distance plot');
    
    corr = corr_vec;
    embedding = embedding_array;
end

function p=parseArguments()
        p = inputParser;
        DEFAULT_SAMPLE_SIZE = 5e3;
        DEFAULT_VERBOSE = 'none';
        EXPECTED_VERBOSE = {'graphic','text','none'};
        addParameter(p,'template_file',[], @(x) ischar(x) || isa(x, 'UMAP'));
        addParameter(p,'label_column',0,@(x) isnumeric(x) && x>0);
        addParameter(p,'sample_size',DEFAULT_SAMPLE_SIZE,@(x) isnumeric(x) && x>1);
        addParameter(p,'verbose',DEFAULT_VERBOSE,...
            @(x) any(validatestring(x,EXPECTED_VERBOSE)));
        addParameter(p,'vary_dims',false,@(x) islogical(x));
end
end