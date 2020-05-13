% function mse = mlpcr_cv_pred(C,X,Y,...)
%
% Performs efficient kfold cross validation and returns hold out set 
%   predictions, and associated statistics. Centers outcome before fitting.
%
% Input ::
%
%   X           - n x p matrix of features
%
%   Y           - n x 1 matrix of outcomes
%
%   C           - n x k matrix specifying k-fold membership. 0 for training
%                   set, 1 for test set
%
%   opts        - Arguments you want passed on to mlpcr. 
%
% Output ::
%
%   pred        - n x 1 vector of out of fold predictions of Y
%
%   STATS       - structure containing information on CV. Currently has MSE
%                   and kfolds used
%
% Notes ::
%
% mlpcr_cv_pred by default automatically precompute X variance 
% decomposition. This only needs to be done once and can be efficiently 
% updated across folds by simply recentering the decomposition for each 
% fold.
%
% WARNING: by default mlpcr_cv_pred won't even merge weights across levels, 
% because the out of fold variance decomposition is always available. 
% However this results in optimistic MSE estimates, because computing
% merged weights is a source of error (often the solution is nonunique, and
% only approximate). Add {'mergedWeights',true} to 'opts' argument to
% override this behavior. Future versions may alter the default, but user
% specified options to 'opts' will always be respected. Using non-merged
% weights may be a valid way of speeding up hyperparameter optimization if
% the bias introduced by foregoing merger is independent of hyperparameter
% choice. This seems like a plausible assumption, but should be tested.
%
% Additional speed ups are possible. Future versions should also perform 
% MLPCA before calling MLPCR to avoid recomputation MLPCA across folds. 
% While the variance decompositions will not differ across crossvalidation 
% folds (except for addative shifts which can be corrected by recentering), 
% PC and CM will differ, but efficient methods exist for recomputing these
% from the full dataset's PCs and CMs without rerunning the full mlpca
% algorithm. These methods are known as "rank 1 downdating" methods (as 
% opposed to 'updating'). These may result in an order of magnitude speed
% boost.
%
% If calling mlpcr_cv_pred as part of an outer cross validation 
% loop even better performance can be achieved by precomputing variance
% components (and in the future MLPCA) before calling mlpcr_cv_pred and 
% passing them as lvlOpts:'varcomps', lvlOpts:'CM' and lvlOpts:'PC' 
% arguments to mlpcr (see example 2 below). If this route is taken make 
% sure to precompute the full eigendecomposition of X rather than a partial 
% decomposition.
%
% Example 1 ::
%
%   Cross validate MLPCR with random subject intercept (default) and slopes, 
%	   using maximal degrees of freedom both within and between 
%	   subjects:
%
%   % prep mlpcr options
%   sid     % a vector with a unique entry for each unique subject. e.g. 
%           % subjid = [1 1 1 2 2 2 3 3 3 4 4 4] suggests the first three
%           % rows of X and Y belong to subject 1, the next three to subject
%           % 2, etc. They do not need to be ordered. Here there are 3
%	        % subject "blocks", one designated by 1s, another by 2s, etc.
%   n_subj = 4; % there are 3 unique subjects identified in subjid
%   fit_lme_options = {'FitMethod','REML','CovariancePattern','isotropic'}
%   mlpcrOpts = {'topLvl', {ones(length(Y),1),n_subj-1}, ...
%                    'withinSubj', {sid,2}, fit_lme_options};
%
%   % prep C to prevent fragmenting of subjects across training/test sets
%   n = length(sid); % num obs
%   uniq_sid = unique(sid); % list of unique subj labels
%   nn = length(uniq_sid); % num subj
%   kfolds = 10;
%   cv = cvpartition(nn,'Kfold',kfolds);
%   C = zeros(n,kfolds);
%   for i = 1:kfolds
%       test_sid = uniq_sid(cv.test(i));
%       test_idx = ismember(sid,test_sid);
%       C(test_idx,i) = 1;
%   end
%
%   mse = mlpcr_cv_pred(C,X,Y,mlpcrOpts{:})
%
%
% Example 2 ::
% 
%   Cross validate MLPCR with random subject intercept (default) and slopes, 
%	   using maximal degrees of freedom both within and between 
%	   subjects, using precomputed variance decomposition and PCA. Reuse
%	   variance decomposition and PCA to compute final predictive weight
%	   map. 11 MLPCR maps will be computed (10 folds + final map) but MLPCA
%	   and variance decomposition are only performed once, providing
%	   potentially dramatic speed and memory savings.
%
%   % metadata
%   sid     % a vector with a unique entry for each unique subject. e.g. 
%           % subjid = [1 1 1 2 2 2 3 3 3 4 4 4] suggests the first three
%           % rows of X and Y belong to subject 1, the next three to subject
%           % 2, etc. They do not need to be ordered. Here there are 3
%	        % subject "blocks", one designated by 1s, another by 2s, etc.
%   n_subj = 4; % there are 3 unique subjects identified in subjid
%
%   % precompute variance decomposition and PCA   
%   pcaLvlOpts = {'topLvl', {ones(length(Y),1),Inf}, 'withinSubj', {sid,Inf}};
%   vcomps = get_nested_var_comps(X,pcaLvlOpts);
%   [PC,CW] = mlpca(X,pcaLvlOpts);
%
%   % construct mlpcr call, note inclusion of lvlOpts:'PC', lvlOpts:'CM' 
%   %    and lvlOpts:'varcomps' args at each level.
%   fit_lme_options = {'FitMethod','REML','CovariancePattern','isotropic'}
%   mlpcrOpts = {'topLvl', ...
%                   {ones(length(Y),1), n_subj-1, 'PC', PC{1}, 'CM', CM{1}, 'varcomps', vcomps{1}}, ...
%               'withinSubj', ...
%                   {sid, 2, 'PC', PC{2}, 'CM', CM{2}, 'varcomps', vcomps{2}}, ...
%               'fitlmeoptions', fit_lme_options};
%
%   % prep C to prevent fragmenting of subjects across training/test sets
%   n = length(sid); % num obs
%   uniq_sid = unique(sid); % list of unique subj labels
%   nn = length(uniq_sid); % num subj
%   kfolds = 10;
%   cv = cvpartition(nn,'Kfold',kfolds);
%   C = zeros(n,kfolds);
%   for i = 1:kfolds
%       test_sid = uniq_sid(cv.test(i));
%       test_idx = ismember(sid,test_sid);
%       C(test_idx,i) = 1;
%   end
%
%   mse = mlpcr_cv_pred(C,X,Y,mlpcrOpts{:})
%
%   [w,Intercept,lme] = mlpcr(X,Y,mlpcrOpts{:});
%
% References ::
%
%  House G, Street G (1995). "The Efficient cross-validation of 
%  principal components applied to principal component regression". 
%  Statistics and Computing (5), 227-235.
%
% Dependencies ::
%
%   mlpcr                   (required by this script)
%   get_nested_var_comps    (required by this script and mlpcr)
%   mlpca                   (required by this script and mlpcr)
%   get_cntrng_mat          (required by mlpca and get_nested_var_comps)

function [pred, STATS] = mlpcr_cv_pred(C,X,Y,varargin)
    kfolds = size(C,2);
    n = size(X,1);
    if n ~= length(Y)
        error('X does not have as many rows as Y has elements');
    end
    
    mlpcrOpts = {'dependentJob','fitlmeoptions'};
    n_d = {};
    dLabels = {};
    id = {};
    vcomps = {};
    CM = {};
    SCW = {};
    vcomps = {};
    for i = 1:length(varargin)
        if ischar(varargin{i}) && ~ismember(varargin{i},mlpcrOpts)
            dLabels{end+1} = varargin{i}; % an arbitrary name to use as a label
            levelOpts = varargin{i+1};
            id{end+1} = levelOpts{1}(:); % block identifying vector
            n_d{end+1} = levelOpts{2}; % number of dims to retain in mlpcaif length(levelOpts) > 2
            aux_levelOpts = levelOpts(3:end);
            for j = 1:length(aux_levelOpts)
                if ischar(aux_levelOpts{j})
                    switch aux_levelOpts{j}
                        case 'PC'
                            CM{end+1} = aux_levelOpts{j+1};
                        case 'CW'
                            SCW{end+1} = aux_levelOpts{j+1};
                        case 'varcomp'
                            vcomps{end+1} = aux_levelOpts{j+1};
                    end
                end
            end
        end
    end
    n_lvls = length(dLabels);
    
    STATS = struct('mse',var(Y,1),'kfolds',kfolds);

    sse = zeros(kfolds,1);
    pred = cell(kfolds,1);
    
    
    % construct MLPCA arguments from available information
    % Note this is used both by PCA and varcomps calls.
    pcaOpts = {};
    for i = 1:length(dLabels)
        pcaOpts = {pcaOpts{:}, dLabels{i}, id{i}, n_d{i}};
    end
    
    newVcomps = 0;
    if isempty(vcomps)
        vcomps = get_nested_var_comps(X,pcaOpts{:});
        newVcomps = 1;
    end
    
    % candidate for parallelization after PCA and varcomp steps are
    % excluded (otherwise those lead to large memory overhead)
    parfor i = 1:kfolds
        train_idx = find(C(:,i) == 0);
        test_idx = find(C(:,i) == 1);
        if length(train_idx) + length(test_idx) ~= n
            error('C matrix must contain only zeros and ones');
        end
        
        X_train = X(train_idx,:);
        Y_train = Y(train_idx);
        arg_train = varargin;
        lvl = 1;
        for j = 2:length(arg_train)
            if iscell(arg_train{j})
                if newVcomps && ~ismember(arg_train{j-1},mlpcrOpts)
                    arg_train{j}{end+1} = 'varcomp';
                    arg_train{j}{end+1} = vcomps{lvl};
                    
                    lvl = lvl + 1;
                end
                for k = 1:length(arg_train{j}) % this will probably need to be updated to not break PCA downdating whenever that gets implemented
                    if size(arg_train{j}{k},1) == n % if it's n x m, drop test rows
                        arg_train{j}{k} = arg_train{j}{k}(train_idx,:);
                    else
                        if size(arg_train{j}{k},2) == n % if it's m x n, drop test columns
                            arg_train{j}{k} = arg_train{j}{k}(:,train_idx);
                        end
                    end
                end
                
                % recenter varcomps w.r.t. this CV fold
                vcIdx = find(strcmp(arg_train{j},'varcomp'));
                if ~isempty(vcIdx)
                    arg_train{j}{vcIdx+1} = ...
                        arg_train{j}{vcIdx+1} - repmat(mean(arg_train{j}{vcIdx+1}),length(train_idx),1);
                end
            end
        end
        
        
        [w, offset] = mlpcr(X_train,Y_train,arg_train{:});

        X_test = X(test_idx,:) - mean(X_train); % MLPCR will have automatically removed X_train mean during training, so recenter accordingly here.
        Y_test = Y(test_idx);
        pred{i} = X_test*w{end} + offset{end};
        
        sse(i) = sum((pred{i} - Y_test).^2);
    end
    
    STATS.mse = sum(sse)/n;
    
    p = zeros(n,1);
    for i = 1:kfolds
        test_idx = find(C(:,i) == 1);
        p(test_idx) = pred{i};
    end
    pred = p;
end
