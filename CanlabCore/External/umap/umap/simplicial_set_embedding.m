function [embedding, method, head, tail, graph_data] = ...
    simplicial_set_embedding(data, graph, n_components, initial_alpha, ...
    a, b, gamma, varargin)
%SIMPLICIAL_SET_EMBEDDING Perform a fuzzy simplicial set embedding, using a
% specified initialisation method and then minimizing the fuzzy set cross
% entropy between the 1-skeletons of the high and low dimensional fuzzy
% simplicial sets.
% 
% [embedding, method] = SIMPLICIAL_SET_EMBEDDING(data, graph, n_components,
% initial_alpha, a, b, gamma)
% 
% Parameters
% ----------
% data: array of size (n_samples, n_features)
%     The source data to be embedded by UMAP.
% 
% graph: sparse matrix
%     The 1-skeleton of the high dimensional fuzzy simplicial set as
%     represented by a graph for which we require a sparse matrix for the
%     (weighted) adjacency matrix.
% 
% n_components: double
%     The dimensionality of the euclidean space into which to embed the data.
% 
% initial_alpha: double
%     Initial learning rate for the SGD.
% 
% a: double
%     Parameter of differentiable approximation of right adjoint functor
% 
% b: double
%     Parameter of differentiable approximation of right adjoint functor
% 
% gamma: double
%     Weight to apply to negative samples.
% 
% negative_sample_rate: double (optional, default 5)
%     The number of negative samples to select per positive sample
%     in the optimization process. Increasing this value will result
%     in greater repulsive force being applied, greater optimization
%     cost, but slightly more accuracy.
% 
% n_epochs: double (optional, default 0)
%     The number of training epochs to be used in optimizing the
%     low dimensional embedding. Larger values result in more accurate
%     embeddings. If 0 is specified a value will be selected based on
%     the size of the input dataset (200 for large datasets, 500 for small).
% 
% init: char array (optional, default 'spectral')
%     How to initialize the low dimensional embedding. Options are:
%         * 'spectral': use a spectral embedding of the fuzzy 1-skeleton
%         * 'random': assign initial embedding positions at random.
%         * An array of initial embedding positions.
% 
% random_state: boolean (optional, default true)
%     If true, MATLAB's RNG will be set to default for reproducibility.
% 
% metric: char array (optional, default 'euclidean')
%     The metric used to measure distance in high dimensional space; used if
%     multiple connected components need to be laid out. Currently not
%     implemented in MATLAB.
% 
% dist_args: cell array (optional, default {})
%     Keyword arguments to be passed to the metric function; used if
%     multiple connected components need to be laid out.
% 
% verbose: boolean (optional, default false)
%     Whether to report information on the current progress of the algorithm.
% 
% Returns
% -------
% embedding: array of size (n_samples, n_components)
%     The optimized of "graph" into an "n_components"-dimensional
%     Euclidean space.
%
%   AUTHORSHIP
%   Math Lead & Primary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Secondary Developer: Stephen Meehan <swmeehan@stanford.edu>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Provided by the Herzenberg Lab at Stanford University 
%   License: BSD 3 clause
        
    p=parseArguments();
    parse(p,varargin{:});
    args=p.Results;
    probability_bin_limit = args.probability_bin_limit;
    epoch_reports = args.epoch_reports;
    progress_callback = args.progress_callback;
    method = args.method;
    verbose = args.verbose;
    min_dist=args.min_dist;
    init = args.init;
    n_epochs = args.n_epochs;
    negative_sample_rate = args.negative_sample_rate;
    random_state = args.random_state;
    
    if random_state
        rng default;
    end

    [n_rows, n_vertices] = size(graph);
    if ~issparse(graph)
        graph = sparse(graph);
    end
    if isempty(n_epochs)
        if n_rows <= 10000
            n_epochs = 500;
        else
            n_epochs = 200;
        end
    end
    C=size(data, 2);
    [~, ~, graph_data] = find(graph);
    if ischar(init) && strcmpi(init, UMAP.INIT_RANDOM)
        embedding = -10 + 20*rand(n_rows, n_components);
    elseif ischar(init) && (strcmpi(init, UMAP.INIT_SPECTRAL) || ...
            strcmpi(init, UMAP.INIT_EIGEN))
        embedding=spectral_layout_binned(data, init, probability_bin_limit, n_components);
        if isempty(embedding)
            graph=remove_sparse(graph, @(a)lt(a, max(graph_data)/n_epochs));
            %tic 
            embedding = spectral_layout(graph, n_components);
            %toc
        end

        expansion = 10 / max(max(embedding));
        embedding = (embedding * expansion) + 0.0001*randn(n_rows, n_components);

        if strcmpi(init, UMAP.INIT_EIGEN)
            head = [];
            tail = [];
            graph_data = [];
            return;
        end
    else %Otherwise, init is the desired initial embedding.
        embedding = init;
    end

    [head, tail, graph_data] = find(graph);
    epochs_per_sample = make_epochs_per_sample(graph_data);

    debugTiming('Cost of embedding with eigen variables -->' )
    [embedding, method] = choose_optimize_layout(embedding, embedding, head, tail,...
        n_epochs, n_vertices, epochs_per_sample, a, b, gamma, initial_alpha,...
        negative_sample_rate, verbose, method, progress_callback, ...
        epoch_reports, random_state, min_dist, [], C);
    debugTiming('Cost of stochastic gradient descent--> ' );
    embedding=double(embedding);
    
    function p=parseArguments(varargin)
        p = inputParser;
        addParameter(p,'probability_bin_limit', 20*4096);
        addParameter(p,'epoch_reports', 0);
        addParameter(p,'progress_callback', []);
        addParameter(p,'method','Java');
        addParameter(p,'verbose', false, @islogical);
        addParameter(p,'dist_args', [], @isnumeric);
        addParameter(p,'metric', 'euclidean');
        addParameter(p,'init', 'spectral');
        addParameter(p,'n_epochs', []);
        addParameter(p,'negative_sample_rate', 5);
        addParameter(p,'random_state', true);
        addParameter(p,'min_dist', .3);
    end
end