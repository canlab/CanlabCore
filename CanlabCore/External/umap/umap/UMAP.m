classdef UMAP < handle
%UMAP Uniform Manifold Approximation and Projection. Finds a low
% dimensional embedding of the data that approximates an underlying
% manifold.
% 
% Parameters
% ----------
% n_neighbors: double (optional, default 15)
%     The size of local neighborhood (in terms of number of neighboring
%     sample points) used for manifold approximation. Larger values result
%     in more global views of the manifold, while smaller values result in
%     more local data being preserved. In general values should be in the
%     range 2 to 100.
% 
% n_components: integer (optional, default 2)
%     The dimension of the space to embed into. This defaults to 2 to
%     provide easy visualization, but can reasonably be set to any integer
%     value in the range 2 to 100.
% 
% metric: string or function (optional, default 'euclidean')
%     The metric to use to compute distances in high dimensional space. If
%     a string is passed, it must match a valid predefined metric. For now,
%     valid string metrics include:
%         * euclidean (or l2)
%         * manhattan (or l1)
%         * chebyshev (or linf)
%         * correlation
%         * cosine
%         * hamming
%         * jaccard
%         * mahalanobis
%         * minkowski
%         * seuclidean
% 
% n_epochs: integer (optional)
%     The number of training epochs to be used in optimizing the low
%     dimensional embedding. Larger values result in more accurate
%     embeddings. If 0, a value will be selected based on the size of the
%     input dataset (200 for large datasets, 500 for small).
% 
% learning_rate: double (optional, default 1)
%     The initial learning rate for the embedding optimization.
% 
% init: string (optional, default 'spectral')
%     How to initialize the low dimensional embedding. Options are:
%         * 'spectral': use a spectral embedding of the fuzzy 1-skeleton
%         * 'random': assign initial embedding positions at random.
%         * An array of initial embedding positions.
% 
% min_dist: double (optional, default 0.1)
%     The effective minimum distance between embedded points. Smaller
%     values will result in a more clustered/clumped embedding where nearby
%     points on the manifold are drawn closer together, while larger values
%     will result on a more even dispersal of points. The value should be
%     set relative to the "spread" value, which determines the scale at
%     which embedded points will be spread out.
% 
% spread: double (optional, default 1)
%     The effective scale of embedded points. In combination with
%     "min_dist" this determines how clustered/clumped the embedded points
%     are.
% 
% set_op_mix_ratio: double (optional, default 1)
%     Interpolate between (fuzzy) union and intersection as the set
%     operation used to combine local fuzzy simplicial sets to obtain a
%     global fuzzy simplicial sets. Both fuzzy set operations use the
%     product t-norm. The value of this parameter should be between 0 and
%     1; a value of 1 will use a pure fuzzy union, while 0 will use a pure
%     fuzzy intersection.
% 
% local_connectivity: integer (optional, default 1)
%     The local connectivity required -- i.e. the number of nearest
%     neighbors that should be assumed to be connected at a local level.
%     The higher this value the more connected the manifold becomes
%     locally. In practice this should be not more than the local intrinsic
%     dimension of the manifold.
% 
% repulsion_strength: double (optional, default 1)
%     Weighting applied to negative samples in low dimensional embedding
%     optimization. Values higher than one will result in greater weight
%     being given to negative samples.
% 
% negative_sample_rate: integer (optional, default 5)
%     The number of negative samples to select per positive sample in the
%     optimization process. Increasing this value will result in greater
%     repulsive force being applied, greater optimization cost, but
%     slightly more accuracy.
% 
% transform_queue_size: double (optional, default 4)
%     For transform operations (embedding new points using a trained model)
%     this will control how aggressively to search for nearest neighbors.
%     Larger values will result in slower performance but more accurate
%     nearest neighbor evaluation.
% 
% a: double (optional)
%     More specific parameters controlling the embedding. If empty these
%     values are set automatically as determined by "min_dist" and
%     "spread".
%
% b: double (optional)
%     More specific parameters controlling the embedding. If empty these
%     values are set automatically as determined by "min_dist" and
%     "spread".
% 
% randomize: boolean (optional, default false)
%     If false, MATLAB's RNG will be set to default for reproducibility.
% 
% metric_kwds: cell array (optional)
%     Arguments to pass on to the metric, such as the "p" value for
%     Minkowski distance. If empty then no arguments are passed on.
% 
% dist_args: 1 or more doubles required by metric.  For example, P for
%     minkowski, Q for Mahalanobis or S for SEuclidean
%
% target_n_neighbors: integer (optional, default -1)
%     The number of nearest neighbors to use to construct the target
%     simplicial set. If set to -1 use the "n_neighbors" value.
% 
% target_metric: string or callable (optional, default 'categorical')
%     The metric used to measure distance for a target array is using
%     supervised dimension reduction. By default this is 'categorical'
%     which will measure distance in terms of whether categories match or
%     are different. Furthermore, if semi-supervised is required target
%     values of -1 will be treated as unlabelled under the 'categorical'
%     metric. If the target array takes continuous values (e.g. for a
%     regression problem) then metric of 'l1' or 'l2' is probably more
%     appropriate.
% 
% target_metric_kwds: cell array (optional)
%     Keyword argument to pass to the target metric when performing
%     supervised dimension reduction. If empty then no arguments are passed
%     on.
% 
% target_weight: double (optional, default 0.5)
%     weighting factor between data topology and target topology. A value
%     of 0 weights entirely on data, a value of 1 weights entirely on
%     target. The default of 0.5 balances the weighting equally between
%     data and target.
% 
% verbose: boolean (optional, default false)
%     Controls verbosity of logging.
%
%   ALGORITHMS
%   UMAP is the invention of Leland McInnes, John Healy and 
%   James Melville at Canada's Tutte Institute for Mathematics and 
%   Computing.  See https://umap-learn.readthedocs.io/en/latest/
%
%   AUTHORSHIP
%   Math Lead & Primary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Secondary Developer: Stephen Meehan <swmeehan@stanford.edu>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Provided by the Herzenberg Lab at Stanford University 
%   License: BSD 3 clause
        
    properties(Constant)
        VERSION='2.0';
        INIT_EIGEN='eigen only';
        INIT_SPECTRAL='spectral';
        INIT_RANDOM='random';
        EIGEN_LIMIT=2*4096;
        METRIC_DICT=containers.Map({'precomputed', 'euclidean', 'l2', 'manhattan', 'l1',...
        'taxicab', 'cityblock', 'seuclidean', 'standardised_euclidean',...
        'chebychev', 'linfinity', 'linfty', 'linf', 'minkowski',...
        'mahalanobis', 'cosine', 'correlation', 'hamming', 'jaccard',...
        'spearman'}, {'precomputed', 'euclidean', 'euclidean', 'cityblock', 'cityblock',...
        'cityblock', 'cityblock', 'seuclidean', 'seuclidean', 'chebychev',...
        'chebychev', 'chebychev', 'chebychev', 'minkowski', 'mahalanobis',...
        'cosine', 'correlation', 'hamming', 'jaccard', 'spearman'});
        DEFAULT_A_1 = 1.576943570285991;
        DEFAULT_A_2 = 0.992174408960354;
        DEFAULT_B_1 = 0.895060673920128;
        DEFAULT_B_2 = 1.112255768443176;
        
        REDUCTION_TYPES='(ub=basic, us=supervised, ubt=template, ust=supervised template)';
        REDUCTION_BASIC='ub';
        REDUCTION_SUPERVISED='us';
        REDUCTION_TEMPLATE='ubt';
        REDUCTION_SUPERVISED_TEMPLATE='ust';
    end
    
    properties(SetAccess=private)
        rawMeans;
        rawStds;
        xLimit;
        yLimit;
        zLimit;
        supervisors;
    end
    
    
    properties
        pythonTemplate;
        usingPython;
        dimNames;
        epoch_reports;
        probability_bin_limit=48*4096;
        progress_callback;
        method='Java'; %OR C OR OR 'C vectorized',  OR 'MatLab' OR 'MatLab Vectorized'
        n_neighbors=15
        n_components=2
        metric='euclidean'
        n_epochs
        learning_rate=1
        init=UMAP.INIT_SPECTRAL
        min_dist=0.1
        spread=1
        set_op_mix_ratio=1
        local_connectivity=1
        repulsion_strength=1
        negative_sample_rate=5,
        transform_queue_size=4
        a
        b
        head
        tail
        weights
        randomize=false
        target_n_neighbors=-1
        target_metric='categorical'
        target_metric_kwds
        target_weight=0.5
        verbose=false
        
        initial_alpha=1
        embedding
        raw_data
        sparse_data=false
        small_data=true
        graph
        search_graph
        distance_func
        dist_args
        IncludeTies
        BucketSize
        NSMethod
        knn_indices
        knn_dists
        parameter_names
        
        
        %WHEN should UMAP invoke nearest neighbor descent 
        %   instead of nearest_neighbors which invokes knnsearch?
        nn_descent_min_rows= 50000 % invoke if data rows>=70000, 0 means NEVER
        nn_descent_min_cols=11 %invoke if data columns >=11 
        nn_descent_max_neighbors=45 % invoke if n_neighbors<-30 
        nn_descent_transform_queue_size=1.35;
    end
    methods
        function U = UMAP(varargin)
    
            if nargin > 0
                p=parseArguments();
                parse(p,varargin{:});
                args=p.Results;
                
                U.probability_bin_limit = args.probability_bin_limit;          
                U.method = args.method;
                U.n_neighbors = args.n_neighbors;
                U.n_components = args.n_components;
                U.metric = args.metric;
                U.n_epochs = args.n_epochs;
                U.learning_rate = args.learning_rate;
                U.init = args.init;
                U.min_dist = args.min_dist;
                U.spread = args.spread;
                U.set_op_mix_ratio = args.set_op_mix_ratio;
                U.local_connectivity = args.local_connectivity;
                U.repulsion_strength = args.repulsion_strength;
                U.negative_sample_rate = args.negative_sample_rate;
                U.transform_queue_size = args.transform_queue_size;
                U.randomize = args.randomize;
                U.target_n_neighbors = args.target_n_neighbors;
                U.target_metric = args.target_metric;
                U.target_weight = args.target_weight;
                U.verbose = args.verbose;

                U.initial_alpha = args.initial_alpha;
                U.sparse_data = args.sparse_data;
                U.small_data = args.small_data;
            end
        end
        
        function method=setMethod(this, method)
            if ~strcmpi(method, 'C') && ~strcmpi(method, 'MatLab')...
                    && ~strcmpi(method, 'C vectorized')  ...
                    && ~strcmpi(method, 'MatLab vectorized')...
                    && ~strcmpi(method, 'MATLAB experimental') ...
                    && ~strcmpi(method, 'MATLAB experimental 2') ...
                    && ~strcmpi(method, 'C++') ...
                    && ~strcmpi(method, 'Mex')
                method='Java';
            end
            this.method=method;
        end
        
        function validate_parameters(U)
            if U.set_op_mix_ratio < 0 || U.set_op_mix_ratio > 1
                error('set_op_mix_ratio must be between 0 and 1');
            end
            if U.repulsion_strength < 0
                error('repulsion_strength cannot be negative');
            end
            if U.min_dist > U.spread
                error('min_dist must be less than or equal to spread');
            end
            if U.min_dist < 0
                error('min_dist must be greater than 0');
            end
            if ~ischar(U.init) && ~isfloat(U.init)
                error('init must be a string or array');
            end
            if ischar(U.init) && ~strcmpi(U.init, UMAP.INIT_SPECTRAL) ...
                    && ~strcmpi(U.init, UMAP.INIT_RANDOM) ...
                    && ~strcmpi(U.init, UMAP.INIT_EIGEN)
                error('string init values must be "spectral" or "random"');
            end
            if isfloat(U.init) && size(U.init, 2) ~= U.n_components
                error('init array must match n_components value');
            end
            if ~ischar(U.metric) && ~isa(U.metric, 'function_handle')
                error('metric must be string or callable');
            end
            if U.negative_sample_rate < 0
                error('negative sample rate must be positive');
            end
            if U.initial_alpha < 0
                error('learning_rate must be positive');
            end
            if U.n_neighbors < 2
                error('n_neighbors must be greater than 2');
            end
            if U.target_n_neighbors < 2 && U.target_n_neighbors ~= -1
                error('target_n_neighbors must be greater than 2');
            end
            if floor(U.n_components) ~= U.n_components
                error('n_components must be an integer');
            end
            if U.n_components < 1
                error('n_components must be greater than 0');
            end
            if U.n_components ~= 2 && ~(strcmpi(U.method, 'MatLab') ...
                    || strcmpi(U.method, 'MatLab vectorized')...
                    || strcmpi(U.method, 'Java')...
                    || strcmpi(U.method, 'C++')...
                    || strcmpi(U.method, 'Mex'))
                error('The C method currently only supports reducing to 2 dimensions');
            end
            if ~isempty(U.n_epochs) && (U.n_epochs <= 10 || (floor(U.n_epochs) ~= U.n_epochs))
                error('n_epochs must be a positive integer larger than 10');
            end
        end
        
        function prepareForTemplate(this, ax)
            this.rawMeans=mean(this.raw_data);
            this.rawStds=std(this.raw_data);
            this.progress_callback=[];
            this.graph=[];
            if ~isempty(this.embedding) && (nargin<2 || isempty(ax))
                mx=max(this.embedding);
                mn=min(this.embedding);
                this.xLimit=[mn(1) mx(1)];
                this.yLimit=[mn(2) mx(2)];
                if size(this.embedding, 2)>2
                    this.zLimit=[mn(3) mx(3)];
                end
            elseif ~isempty(ax)
                this.xLimit=xlim(ax);
                this.yLimit=ylim(ax);
                if this.n_components>2
                    this.zLimit=zlim(ax);
                end
            end
        end
        
        function adjustLims(this, ax, data)
            if isprop(this, 'xLimit') && ~isempty(this.xLimit)
                if this.n_components<3
                    xlim(ax, this.xLimit);
                    ylim(ax, this.yLimit);
                elseif ~isempty(this.zLimit)
                    xlim(ax, this.xLimit);
                    ylim(ax, this.yLimit);
                    zlim(ax, this.zLimit);
                else
                    return;
                end
                if nargin>2 && ~isempty(data)
                    Gui.StretchLims(ax, data, .04);
                end
            end 
        end

        function clearLimits(this)
            this.xLimit=[];
            this.yLimit=[];
        end
        
        function ok=report(this, log)
            ok=false;
            if this.verbose
                if isequal('function_handle', class(this.progress_callback))
                    if ~feval(this.progress_callback, log)
                        return;
                    end
                else
                    disp(log);
                end
            end
            ok=true;
        end
        
        function U = fit(U, X, y)
            if nargin < 3
                y = [];
            end
            
            if ~isequal('function_handle', class(U.metric)) 
                U.metric = U.METRIC_DICT(lower(U.metric));
            end
            if isa(U.metric, 'function_handle')
                U.distance_func = U.metric;
            elseif isKey(U.METRIC_DICT, U.metric)
                 U.metric = U.METRIC_DICT(U.metric);
                 U.distance_func = U.metric;
            elseif strcmpi(U.metric, 'precomputed')
                warning('Using precomputed metric; transform will be unavailable for new data');
            else
                error('Metric is neither callable, nor a recognised string');
            end
            
            X_rows = size(X, 1);
            
            U.raw_data = X;

            if isempty(U.a) || isempty(U.b)
                if U.spread == 1 && U.min_dist == 0.1
                    U.a = U.DEFAULT_A_1;
                    U.b = U.DEFAULT_B_1;
                elseif U.spread == 1 && U.min_dist == 0.3
                    U.a = U.DEFAULT_A_2;
                    U.b = U.DEFAULT_B_2;
                else
                    [U.a, U.b] = find_ab_params(U.spread, U.min_dist);
                end
            end

            U.initial_alpha = U.learning_rate;

            validate_parameters(U);

            if U.verbose
                disp(str(U));
            end

            if X_rows <= U.n_neighbors
                if X_rows == 1
                    U.embedding = zeros(1, U.n_components);
                    return
                end

                warning('n_neighbors is larger than the dataset size; truncating to size(X, 1) - 1');
                U.n_neighbors = X_rows - 1;
            end

            if issparse(X)
                U.sparse_data = true;
            end
            
            if islogical(U.randomize) && ~U.randomize
                rng default;
            elseif isnumeric(U.randomize) && ~isempty(U.randomize)
                rng(U.randomize);
            else
                rng;
            end
            
            [R,C]=size(U.raw_data);
            strR=CommaFormat(R);
            strC=CommaFormat(C);
            tic;
            tag=sprintf(' (%s x %s values)', strR, strC);
            if X_rows < 4096
                U.report(['  Computing fuzzy simplicial set' tag]);
                U.small_data = true;
                if ~strcmpi(U.metric, 'precomputed')
                    dmat = squareform(pdist(X, U.metric, U.dist_args));
                else
                    dmat = X;
                end
                U.graph = fuzzy_simplicial_set(dmat, U.n_neighbors, 'precomputed',...
                    'set_op_mix_ratio', U.set_op_mix_ratio,...
                    'local_connectivity', U.local_connectivity);
            else
                U.small_data = false;
                if KnnFind.CanAccelerate(U, size(X))
                    if ~U.report(['  Approximating ' num2str(U.n_neighbors) ...
                            ' neighbors for ' strR ' data points'])
                        U=[];
                        return;
                    end
                    if ~islogical(U.randomize)
                        U.randomize=U.randomize~=0;
                    end
                    [U.knn_indices, U.knn_dists]=KnnFind.Approximate(X, [], ...
                        'metric', U.metric, 'K', U.n_neighbors, ...
                        'dist_args', U.dist_args, 'randomize', U.randomize,...
                        'progress_callback', U.progress_callback);
                    if isempty(U.knn_indices)
                        U=[]; %NnDescent supports progress callbacks!
                        return;
                    end
                else
                    if ~U.report(['  Determining ' num2str(U.n_neighbors) ...
                            ' neighbors for ' strR ' data points'])
                        U=[];
                        return;
                    end
                    
                    [U.knn_indices, U.knn_dists] = nearest_neighbors(X, U);
                end
                if ~U.report(['  Computing fuzzy simplicial set' tag])
                    U=[];
                    return;
                end
                U.graph = fuzzy_simplicial_set(X, U.n_neighbors, U.metric,...
                    'knn_indices', U.knn_indices, 'knn_dists', ...
                    U.knn_dists, 'set_op_mix_ratio', U.set_op_mix_ratio,...
                    'local_connectivity', U.local_connectivity);
            end
            
            debugTiming('Cost of fuzzy simplicial (knnsearch) -->' )
            
            if ~isempty(y)
                if ~U.report(['  Fitting the supervisory labels ' tag])
                    U=[];
                    return;
                end
                if strcmpi(U.target_metric, 'categorical')
                    if U.target_weight < 1
                        far_dist = 2.5 * (1 / (1 - U.target_weight));
                    else
                        far_dist = 1e12;
                    end
                    U.graph=categorical_simplicial_set_intersection(...
                        U.graph, y, 1, far_dist);
                else
                    if U.target_n_neighbors == -1
                        U.target_n_neighbors = U.n_neighbors;
                    end

                    if size(y, 1) < 4096
                        ydmat = squareform(pdist(y, U.target_metric, ...
                            U.target_metric_kwds));
                        target_graph = fuzzy_simplicial_set(ydmat, ...
                            U.target_n_neighbors,...
                            'precomputed', 'set_op_mix_ratio', 1,...
                            'local_connectivity', 1);
                    else
                        target_graph = fuzzy_simplicial_set(y, ...
                            U.target_n_neighbors,...
                            U.target_metric,...
                            'set_op_mix_ratio', 1, 'local_connectivity', 1);
                        
                    end
                    U.graph = general_simplicial_set_intersection(U.graph, target_graph, U.target_weight);
                    U.graph = reset_local_connectivity(U.graph);
                end
                debugTiming('Cost of categorigal simplicial (supervising) -->' )
            end
            if ~U.report(['  Computing embedding' tag ])
                    U=[];
                    return;
            end
            
            U.search_graph = (U.graph > 0) + speye(size(U.graph,1));
            
            [U.embedding, U.method, U.head, U.tail, U.weights] = simplicial_set_embedding(...
                U.raw_data, U.graph, U.n_components,...
                U.initial_alpha, U.a, U.b, U.repulsion_strength, ...
                'negative_sample_rate', U.negative_sample_rate, ...
                'n_epochs', U.n_epochs, 'init', U.init, 'metric', U.metric,...
                'dist_args', U.dist_args, 'verbose', U.verbose, ...
                'method', U.method, 'progress_callback', ...
                U.progress_callback, 'epoch_reports', U.epoch_reports,...
                'probability_bin_limit', U.probability_bin_limit, ...
                'random_state', ~U.randomize, ...
                'min_dist', U.min_dist);
        end
        
        function X_new = fit_transform(U, X, y)
        
            if nargin < 3
                y = [];
            end

            U = fit(U, X, y);
            if nargout>0
                if ~isempty(U)
                    X_new = U.embedding;
                else
                    X_new=[];
                end
            end
        end
        
        function matched = compare_new_knn(U)
            
            if isempty(U.raw_data)
                disp('UMAP object has no data!');
            elseif isempty(U.embedding)
                disp('No dimension reduction has been done!');
            end
            
            orig_indices = knnsearch(U.raw_data,U.raw_data,'K',U.n_neighbors,'Distance',U.metric);
            new_indices = knnsearch(U.embedding,U.embedding,'K',U.n_neighbors,'Distance',U.metric);
            
            matched = 0;
            for i = 1:size(orig_indices,1)
                matched = matched + nnz(ismember(orig_indices(i,:), new_indices(i,:)));
            end
               
            matched = matched/numel(orig_indices);
        end

        function X_new = transform(U, X)
        
            if size(U.embedding, 1) == 1
                error('Transform unavailable when model was fit with only a single data sample.');
            end

            if U.sparse_data
                error('Transform not available for sparse input.');
            elseif strcmpi(U.metric, 'precomputed')
                error('Transform of new data not available for precomputed metric.');
            end

            if islogical(U.randomize) && ~U.randomize
                rng default;
            elseif isnumeric(U.randomize) && ~isempty(U.randomize)
                rng(U.randomize);
            else
                rng;
            end
            
            [n_samples, C] = size(X);
            strR=CommaFormat(n_samples);
            strC=CommaFormat(C);
            
            tag=sprintf(' (%s x %s values)', strR, strC);
            embeddingCount = size(U.raw_data, 1);
            
            if KnnFind.CanAccelerate(U, size(X), U.search_graph)
                if U.verbose
                    if ~U.report(['  Approximating ' num2str(U.n_neighbors) ...
                            ' template neighbors for ' strR ' data points'])
                        X_new=[];
                        return;
                    end
                end
                
                cb=[];
                %cb=U.progress_callback; % prevents parallel run
                if ~islogical(U.randomize)
                    U.randomize=U.randomize~=0;
                end
                [indices, dists]=KnnFind.Approximate(U.raw_data, X, ...
                    'metric', U.metric, 'K', U.n_neighbors, ...
                    'dist_args', U.dist_args, 'randomize', U.randomize,...
                    'X_search_graph', U.search_graph, ...
                    'nn_descent_transform_queue_size', U.nn_descent_transform_queue_size);
                if isempty(indices)
                    X_new=[];%NnDescent supports progress callbacks!
                    return;
                end
            else
                if ~U.report(['  Determining ' num2str(U.n_neighbors) ...
                        ' template neighbors for ' strR ' data points'])
                    X_new=[];
                    return;
                end
                [indices, dists] = KnnFind.Determine(U.raw_data,X, U);
            end
            if U.verbose
                if isequal('function_handle', class(U.progress_callback))
                    if ~feval(U.progress_callback, ...
                            'Incrementally embedding template')
                        X_new=[];
                        return;
                    end
                else
                    disp('Incrementally embedding template...');
                end
            end
            if ~isempty(U.supervisors)
                U.supervisors.knnIndices=indices;
                U.supervisors.supervisingMean=mean(U.raw_data(:));
                U.supervisors.supervisedMean=mean(X(:));
            end
            adjusted_local_connectivity = max(0, U.local_connectivity - 1);
            [sigmas, rhos] = smooth_knn_dist(dists, U.n_neighbors, 'local_connectivity', adjusted_local_connectivity, 'same_set', false);

            [rows, cols, vals] = compute_membership_strengths(indices, dists, sigmas, rhos, false);

            weight_graph = sparse(rows, cols, vals, n_samples, embeddingCount);

            norms = sum(abs(weight_graph),2);
            norms(norms == 0) = 1;
            
            vals = reshape(vals, [U.n_neighbors, n_samples])';
            cols = reshape(cols, [U.n_neighbors, n_samples])';
            init_weights = vals./full(norms);
            
            
            init_embedding = init_transform(cols, init_weights, U.embedding);

            if isempty(U.n_epochs) || isnan(U.n_epochs)
                if n_samples <= 10000
                    n_trans_epochs = 100;
                else
                    n_trans_epochs = 30;
                end
            else
                n_trans_epochs = floor(U.n_epochs / 3);
            end

            [~,~,graph_data]=find(weight_graph);
            weight_graph=remove_sparse(weight_graph, @(a)lt(a, max(graph_data)/n_trans_epochs));
            
            [U.head, U.tail, U.weights] = find(weight_graph);

            epochs_per_sample = make_epochs_per_sample(U.weights);
            [X_new, ~] = choose_optimize_layout(init_embedding, U.embedding, U.head,...
                U.tail, n_trans_epochs, embeddingCount, epochs_per_sample, U.a, ...
                U.b, U.repulsion_strength, U.initial_alpha,...
                U.negative_sample_rate, U.verbose, U.method, ...
                U.progress_callback, U.epoch_reports, ~U.randomize, ...
                U.min_dist, [], C);
        end
        
        function X_new = transform2(U, X)
            [n_samples, C] = size(X);
            if size(U.embedding, 1) == 1
                error('Transform unavailable when model was fit with only a single data sample.');
            end

            if U.sparse_data
                error('Transform not available for sparse input.');
            elseif strcmpi(U.metric, 'precomputed')
                error('Transform of new data not available for precomputed metric.');
            end

            if islogical(U.randomize) && ~U.randomize
                rng default;
            elseif isnumeric(U.randomize) && ~isempty(U.randomize)
                rng(U.randomize);
            else
                rng;
            end
            
            if U.verbose
                if isequal('function_handle', class(U.progress_callback))
                    if ~feval(U.progress_callback, 'Initializing template')
                        X_new=[];
                        return;
                    end
                else
                    disp('Initializing template...');
                end
            end            
            
            embeddingCount = size(U.raw_data, 1);
            
            [indices, dists] = knnsearch(U.raw_data,X,'K',U.n_neighbors,'Distance',U.metric);
            adjusted_local_connectivity = max(0, U.local_connectivity - 1);
            [old_sigmas, old_rhos] = smooth_knn_dist(dists, U.n_neighbors, 'local_connectivity', adjusted_local_connectivity, 'same_set', false);
            [rows, cols, vals] = compute_membership_strengths(indices, dists, old_sigmas, old_rhos, false);
            weight_graph = sparse(rows, cols, vals, n_samples, embeddingCount);
            norms = sum(abs(weight_graph),2);
            norms(norms == 0) = 1;
            old_cols = reshape(cols, [U.n_neighbors, n_samples])';
            vals = reshape(vals, [U.n_neighbors, n_samples])';
            old_weights = vals./full(norms);
            
            XY = vertcat(X, U.raw_data);
            move_point = vertcat(true(n_samples, 1), false(embeddingCount, 1));

            [indices, dists] = knnsearch(XY,X,'K',U.n_neighbors,'Distance',U.metric);
            
            [sigmas, rhos] = smooth_knn_dist(dists, U.n_neighbors, 'local_connectivity', U.local_connectivity);

            [rows, cols, vals] = compute_membership_strengths(indices, dists, sigmas, rhos, true);

            weight_graph = sparse(rows, cols, vals, n_samples, n_samples+embeddingCount);
            
            if U.verbose
                if isequal('function_handle', class(U.progress_callback))
                    if ~feval(U.progress_callback, ...
                            'Incrementally embedding template')
                        X_new=[];
                        return;
                    end
                else
                    disp('Incrementally embedding template...');
                end
            end
            
            init_embedding = init_transform(old_cols, old_weights, U.embedding);
            XY_low = vertcat(init_embedding, U.embedding);

            if isempty(U.n_epochs) || isnan(U.n_epochs)
                if n_samples <= 10000
                    n_trans_epochs = 200;
                else
                    n_trans_epochs = 100;
                end
            else
                n_trans_epochs = floor(U.n_epochs / 2);
            end

            [~,~,graph_data]=find(weight_graph);
            weight_graph=remove_sparse(weight_graph, @(a)lt(a, max(graph_data)/n_trans_epochs));
            
            [U.head, U.tail, U.weights] = find(weight_graph);

            epochs_per_sample = make_epochs_per_sample(U.weights);
            %X_new = trans2_optimize_layout(embedding, XY_low, head,...
            %    tail, n_epochs, embeddingCount, epochs_per_sample, U.a, ...
            %    U.b, U.repulsion_strength, U.initial_alpha,...
            %    U.negative_sample_rate, U.verbose, U.min_dist, move_point);

            X_new = choose_optimize_layout(init_embedding, XY_low, U.head,...
                U.tail, n_trans_epochs, embeddingCount, epochs_per_sample, U.a, ...
                U.b, U.repulsion_strength, U.initial_alpha,...
                U.negative_sample_rate, U.verbose, U.method, ...U.method, ...
                U.progress_callback, U.epoch_reports, ~U.randomize, ...
                U.min_dist, move_point, C);

        end

        function char = str(U)
            if isequal('function_handle', class(U.distance_func)) 
                dstr='custom';
            else
                dstr=mat2str(U.distance_func);
            end
            if isequal('function_handle', class(U.metric)) 
                mstr='custom';
            else
                mstr=mat2str(U.metric);
            end
            if ischar(U.init)
                char = ['UMAP(method=' U.method...
                    ', n_neighbors=' mat2str(U.n_neighbors)...
                    ', n_components=' mat2str(U.n_components)...
                    ', metric=' mstr...
                    ', n_epochs=' mat2str(U.n_epochs)...
                    ', learning_rate=' mat2str(U.learning_rate)...
                    ', init=' U.init...
                    ', min_dist=' mat2str(U.min_dist)...
                    ', spread=' mat2str(U.spread)...
                    ', set_op_mix_ratio=' mat2str(U.set_op_mix_ratio)...
                    ', local_connectivity=' mat2str(U.local_connectivity)...
                    ', repulsion_strength=' mat2str(U.repulsion_strength)...
                    ', negative_sample_rate=' mat2str(U.negative_sample_rate)...
                    ', transform_queue_size=' mat2str(U.transform_queue_size)...
                    ', a=' mat2str(U.a)...
                    ', b=' mat2str(U.b)...
                    ', randomize=' mat2str(U.randomize)...
                    ', target_n_neighbors=' mat2str(U.target_n_neighbors)...
                    ', target_metric=' mat2str(U.target_metric)...
                    ', target_metric_kwds=' mat2str(U.target_metric_kwds)...
                    ', target_weight=' mat2str(U.target_weight)...
                    ', verbose=' mat2str(U.verbose)...
                    ', initial_alpha=' mat2str(U.initial_alpha)...
                    ', sparse_data=' mat2str(U.sparse_data)...
                    ', small_data=' mat2str(U.small_data)...
                    ', distance_func=' dstr...
                    ', dist_args=' mat2str(U.dist_args)];
            else
                char = ['UMAP(method=' U.method...
                    ', n_neighbors=' mat2str(U.n_neighbors)...
                    ', n_components=' mat2str(U.n_components)...
                    ', metric=' mstr...
                    ', n_epochs=' mat2str(U.n_epochs)...
                    ', learning_rate=' mat2str(U.learning_rate)...
                    ', init=custom input'...
                    ', min_dist=' mat2str(U.min_dist)...
                    ', spread=' mat2str(U.spread)...
                    ', set_op_mix_ratio=' mat2str(U.set_op_mix_ratio)...
                    ', local_connectivity=' mat2str(U.local_connectivity)...
                    ', repulsion_strength=' mat2str(U.repulsion_strength)...
                    ', negative_sample_rate=' mat2str(U.negative_sample_rate)...
                    ', transform_queue_size=' mat2str(U.transform_queue_size)...
                    ', a=' mat2str(U.a)...
                    ', b=' mat2str(U.b)...
                    ', randomize=' mat2str(U.randomize)...
                    ', target_n_neighbors=' mat2str(U.target_n_neighbors)...
                    ', target_metric=' mat2str(U.target_metric)...
                    ', target_metric_kwds=' mat2str(U.target_metric_kwds)...
                    ', target_weight=' mat2str(U.target_weight)...
                    ', verbose=' mat2str(U.verbose)...
                    ', initial_alpha=' mat2str(U.initial_alpha)...
                    ', sparse_data=' mat2str(U.sparse_data)...
                    ', small_data=' mat2str(U.small_data)...
                    ', distance_func=' dstr...
                    ', dist_args=' mat2str(U.dist_args)];
            end
        end
        function supervisors=setSupervisors(this, labels, labelMap, ax)
            supervisors=Supervisors(labels, labelMap, this.embedding, ...
                ax, size(this.raw_data,2));
            this.supervisors=supervisors;
        end
    end
end

function p=parseArguments(varargin)
    p = inputParser;
    addParameter(p,'probability_bin_limit',20*4096);          
    addParameter(p,'method','Java');
    addParameter(p,'n_neighbors',15);
    addParameter(p,'n_components',2);
    addParameter(p,'metric','euclidean');
    addParameter(p,'n_epochs',[]);
    addParameter(p,'learning_rate',1);
    addParameter(p,'init',UMAP.INIT_SPECTRAL);
    addParameter(p,'min_dist',0.1);
    addParameter(p,'spread',1);
    addParameter(p,'set_op_mix_ratio',1);
    addParameter(p,'local_connectivity',1);
    addParameter(p,'repulsion_strength',1);
    addParameter(p,'negative_sample_rate',5);
    addParameter(p,'transform_queue_size',4);
    addParameter(p,'randomize',false);
    addParameter(p,'target_n_neighbors',-1);
    addParameter(p,'target_metric','categorical');
    addParameter(p,'target_weight',0.5);
    addParameter(p,'verbose',false, @islogical);
    addParameter(p,'dist_args',[], @isnumeric);
    addParameter(p,'initial_alpha',1);
    addParameter(p,'sparse_data',false);
    addParameter(p,'small_data',true);
end
       