% [weights, Intercept, lme] = mlpcr(X,Y,options)
%
% Performs multilevel principal components regression using ANOVA principles
% to linearly separate between and within level variance. 
%
% X is separated into between and within level variance, PCA is run on each 
% level of X separately, components and loadings are rotated to ensure 
% all levels are orthogonal, and resultant loadings are used in a linear 
% mixed effects model to predict Y. The resultant model fit is returned
% (lme), and coefficients are projected back into voxel space. Separate
% weight maps will be provided for each level evaluated.
%
% Level specific weight maps can be used individually on novel data, or
% a sum level specific weights is also returned. For best results subtract 
% out mean parameter values in training data from all test data (e.g. 
% subtract mean voxel map across training brains before applying predictive
% "signatures" to new data). The difference should be slight though.
%
% With multiple observations per level in test data additional more
% substantial improvements are possible by noting that maps for one level
% should not predict lower level variance, so for between subject maps
% should not predict within subject variance. Any within subject variance
% predicted by between subject maps is noise and can be subtracted out. See
% example 1.
%
% Multiple levels can be specified, although the first level specified must 
% be the fixed effects level. Subsequent levels will be nested in the order 
% they're specified. So for instance if you have pain ratings, nested within
% stimulus level, nested within subjects, then first provide subject
% information using 'lvlName',{wiBlocks,n_dims} optional parameters. Then
% provide information for within subject blocks, and finally provide 
% information for within stimulus blocks. This is only a description of the 
% order in which information needs to be provided for nesting to be treated 
% properly. See USAGE section for details on what level specific information 
% is required.
%
% USAGE ::
%
%   weights         (lvl+1) x 1 cell array of weight maps. One for each level
%                       evaluated plus weights{end} which is a sum of all the
%                       others
%
%   Intercept       Model intercepts. I{end} is the mean outcome value. Each 
%			            variance-component-specific weight map will have its 
%			            own associated intercept. These should be ~0, but are
%			            computed anyway as a sanity check.
%
%
%   X               n x p predictor variable matrix
%
%   Y               n x 1 outcome vector
%
%  'lvlName'        A string specifying level name. Must be followed by
%                       cell of level options. Level options cell must have a
%			            vector of random effects block labels. Random effect 
%                       labels must in turn be followed by a dimension number
%                       specifying the number of principal components to
%                       retain for the PCA decomposition of that levels 
%			            variance components. See below regarding the latter. 
%                   Note that the name of the variable doesn't matter. It's
%                       not actually used for anything. Mostly just a bookend
%                       telling the script where to look for block vectors 
%                       and to help the user keep track of inputs.
%
%   lvlOpts         A cell array containing options to apply to the grouping
%			            level in question. Must contain at least two elements:
%			            A random effects label vector (see lvlOpts:blocks) 
%                       and a scalar value indicating number of PCA 
%                       dimensions to retain at that level.
%
%   lvlOpts:blocks  Mandatory first element of lvlOpts cell. n x 1 numeric vector
%                      specifying block membership. For fixed effects provide
%                      a unit vector (i.e. for the first/topmost level).
%
%   lvlOpts:n_dims  Must follow lvlOpts:blocks. Specifies number of principle
%                      components to retain for LME fitting.
%
%   lvlOpts:'covs'  Optional string argument. Must be followed by an n x m matrix 
%			           of confounds to control for during lme model fitting. 
%			           Ommit if no relevant covariates are available.
%
%   lvlOpts:'nointercept'
%                   By default level will be fit with intercept. Specify 
%			           'nointercept' to fit without an intercept.
%
%   lvlOpts:'PC'    Principal component coefficients returned by MLPCA for 
%                       the level in question. Must be followed by an 
%                       p x d, d >= n_dims matrix with component loadings. 
%                       d > n_dims will be silently dropped.
%                   An associated 'CW' argument must also be provided. 
%                   Finally these arguments must either be specified 
%                       for all levels or no levels.
%                   Note: In most circumstances the user should not specify
%                       this argument. MLPCR will automatically perform 
%                       MLPCA. The option exists for use in relatively
%                       sophisticated applications such as efficient CV 
%                       algorithms and (potentially) kernel PCA methods.
%
%   lvlOpts:'CW'    Principal component scores returned by MLPCA for the 
%                       level in question. Must be followed by a n x d, 
%                       d >= n_dims matrix with component weights. 
%                       d > n_dims will be silently dropped. 
%                   An associated 'PC' argument must also be provided. 
%                   Finally these arguments must either be specified 
%                       for all levels or no levels.
%                   Note: In most circumstances the user should not specify
%                       this argument. MLPCR will automatically perform 
%                       MLPCA. The option exists for use in relatively
%                       sophisticated applications such as efficient CV 
%                       algorithms and (potentially) kernel methods.
%
%   lvlOpts:'varcomp'
%                   Variance fraction of X that belongs to this level. Must
%                       be followed by n x m matrix. If specified for this
%                       level must be specified for all levels.
%                   Note: In most circumstances the user should not specify
%                       this argument. MLPCR will automatically perform 
%                       variance decomposition. One circumstances where it
%                       may be helpful to precompute this is in CV
%                       algorithms that call MLPCR. Use
%                       get_nested_var_comps() to compute these values.
%
%   'fitlmeoptions' Must be followed by a cell array of options to pass to 
%                       lmefit. See "help fitlme" for details.
%
%   'dependentJob'  MLPCA can be very memory intensive. If running mlpcr in
%                       parallel it can be a good idea to only run one mlpca
%                       at a time. Specify an array of jobs that should
%                       wait for this job to complete before running mlpca.
%                       jobIDs in array should correspond to "labindices".
%                       use labRecieve(<this jobs labindex>) in those jobs
%                       on the other end to get them going upon receipt of 
%                       this jobs message. Note: PCA is automatically 
%                       parallelized by matlab, but if using compute labs
%                       PCA will run serially instead, and hence slow down.
%                       Parallelization of MLPCR will only improve speed if
%                       what's gained by running multiple parallel fitlme()s
%                       outweighs what's lost from serializing PCA. Invoke
%                       multithreadWorkers() after starting a parpool for
%                       best results.
%
%   'verbose'           Print timing and formula information (useful for
%                       realtime sanity checks for model fits that take 
%                       time)
%
% EXAMPLES ::
%
%   Example 1:
%
%   Mixed effects PCR with random subject intercept (default) and slopes, 
%	   using maximal degrees of freedom both within and between 
%	   subjects:
%
%   subjid; % a vector with a unique entry for each unique subject. e.g. 
%           % subjid = [1 1 1 2 2 2 3 3 3 4 4 4] suggests the first three
%           % rows of X and Y belong to subject 1, the next three to subject
%           % 2, etc. They do not need to be ordered. Here there are 3
%	        % subject "blocks", one designated by 1s, another by 2s, etc.
%   n_subj = 4; % there are 3 unique subjects identified in subjid
%   fit_lme_options = {'FitMethod','REML','CovariancePattern','isotropic'}
%   [w,I,lme] = mlpcr(X,Y,'topLvl', {ones(length(Y),1),n_subj-1}, ...
%                    'withinSubj', {subjid,8}, 'fitlmeoptions', fit_lme_options);
%           
%   X_test; % n x d matrix of subject single trial maps
%
%   Y_pred = X_test*w{1} + X_test*x{2} + I{end}
%   %equivalently
%   Y_pred = X_test*w{3} + I{end};
%
%   % an improved prediction which leverages repeated measures for
%   % denoising and improved predictions
%   c = get_cntrng_mat(subjid) % a subject centering matrix
%   btNoise = c*X_test*w{1}; % w{1} is orthogonal to meaningful within 
%                            % subject variance by design, so any within
%                            % subject predictive variance of w{1} is noise
%                            % in an estimate of between subject variance.
%   Y_pred = X_test*w{3} + I{end} - btNoise;
%
%   Example 2:
%
%   PCR (same as
%   fmri_data/predict('algorithm_name','cv_pcr','numcomponents',dims):
%
%   [w,Intercept,lme] = mlpcr(train.X, train.Y,'fixed_effects', ...
%       {ones(length(train.Y),1), dims});
%
% USAGE NOTES ::
%
% lvlName and associated values must be specified in descending
% order from the topmost level to the bottommost level. This means the
% first level information to provide should be for the fixed effects.
%
% Different covariance structure can be used in the LME model. The full
% covariance may not in fact result in the best predictive performance, due
% to overfitting. Either way, isotropic covariance tends to produce good
% performance at low computational cost, and may also result in smoother 
% objective functions when evaluating how many PCA dimensions to retain. 
% This would make hyperparameter optimization quicker by requiring fewer
% samples to infer the overall function (e.g. using bayesian hyperparameter
% optimization). Diagonal covariance may also merit consideration. Full 
% covariance is typically computationally prohibitive.
%
% 'FitMethod' 'REML' is the default but 'ML' will run faster (~60%
% faster in one instance, probably highly variable, but this observation 
% gives a sense of potential magnitude of speed up).
%
% The level names are not actually used in any meaningful way, they're
%   just arbitrary labels, but should be distinct.
%
% Dependencies:
%   mlpca                   (required by mlpcr)
%
% Writen by Bogdan Petre (Dec 24, 2017)
% Updated for covariates and correct weighted average of level specific
% 	weights by Bogdan (Dec, 21, 2018)
% Corrected weight estimates now provided, with higher effects (e.g.
%   between effects) represented by vectors that are orthogonal to lower
%   levels (e.g. within effects). Alternative solution to total weight map
%   derived, equivalent to pinv([SCW][CM]') == pinv(pca_denoised_X) ==
%   sum(<level specific weight maps>). (May, 10, 2019)

function [weights, Intercept, lme, CM, SCW] = mlpcr(X,Y,varargin)
    %% parse arguments
    %sanity check on dimensions still needs to be implemented
    idx = 1;
    lmeopts = cell(0);
    depJobID = [];
    lmelvls = -1;
    verbose = false;
    fitintercept = {};
    mergeWeights = true;
    CM = {};
    SCW = {};
    vcomps = {};
    for i = 1:length(varargin)
        if ischar(varargin{i})
            switch varargin{i}
                case 'fitlmeoptions'
                    lmeopts = varargin{i+1};                
                    for j = 1:length(lmeopts)
                        if ischar(lmeopts{j})
                            switch lmeopts{j}
                                case 'CovariancePattern'
                                    if ischar(lmeopts{j+1})
                                        lmelvls = 1;
                                    else
                                        lmelvls = length(lmeopts{j+1});
                                    end
                            end
                        end
                    end
                case 'dependentJob'
                    depJobID = varargin{i+1};
                case 'verbose'
                    verbose = true;
                case 'mergeWeights' % this option is deprecated and no longer in use
                    if varargin{i+1} ~= true & varargin{i+1} ~= false
                        error('mergeWeights expected true/false.');
                    end
                    mergeWeights = varargin{i+1};
                otherwise 
                    dLabels{idx} = varargin{i}; % an arbitrary name to use as a label
                    levelOpts = varargin{i+1};
                    id{idx} = levelOpts{1}(:); % block identifying vector
                    n_d{idx} = levelOpts{2}; % number of dims to retain in mlpca
                    covs{idx} = [];
                    covs_Labels{idx} = [];
                    fitintercept{end+1} = true; % default
                    if length(levelOpts) > 2
                        aux_levelOpts = levelOpts(3:end);
                        for j = 1:length(aux_levelOpts)
                            if ischar(aux_levelOpts{j})
                                switch aux_levelOpts{j}
                                    case 'covs'
                                        if size(aux_levelOpts{i+1},1) == length(id{idx})
                                            covs{idx} = levelOpts{i+1};
                                        else
                                            if size(aux_levelOpts{i+1},2) == length(id{idx})
                                                warning(['Looks like covariates for ' dLabels{idx} ' are not oriented correctly. Rotating... Check to make sure this is right!']);
                                                covs{idx} = aux_levelOpts{i+1}';
                                            else
                                                error(['Did not understand 4th entry for in ' levelOpts{1} ' options.']);
                                            end
                                        end
                                    case 'nointercept'
                                        fitintercept{end} = false;
                                    case 'PC'
                                        CM{end+1} = aux_levelOpts{j+1};
                                    case 'CW'
                                        SCW{end+1} = aux_levelOpts{j+1};
                                    case 'varcomp'
                                        vcomps{end+1} = aux_levelOpts{j+1};
                                    otherwise
                                        error(['Did not understand entry for ' levelOpts{1}]);
                                end
                            end
                        end
                    end
                    % construct covariate labels for fitlme model invocation
                    if ~isempty(covs{idx})
                        this_covs = covs{idx};
                        these_labels = {};
                        for j = 1:size(this_covs,2)
                            these_labels = {these_labels{:}, [dLabels{idx}, '_cov' num2str(j)]};
                        end
                        covs_Labels{idx} = these_labels;
                    end
                    idx = idx+1;
            end
        end
    end
    n_lvls = length(dLabels);
    if lmelvls == -1
        lmelvls = n_lvls - 1;
    end
    
    %% sanity checks and preproc
    if ~isempty(CM)
        if length(CM) ~= n_lvls
            error(sprintf('Insufficient CWs supplied to level options. %d required, but only %d found.', n_lvls, length(CM)));
        else
            for i = 1:n_lvls
                CMd = [size(X,2),n_d{i}];
                if n_d{i} ~= 0
                    if size(CM{i},1) ~= CMd(1) || size(CM{i},2) < CMd(2)
                        error(sprintf('Dimensions of %s ''PC'' arg are %d x %d. Expected %d x d, d >= %d',...
                            dLabels{i},size(CM{i},1),size(CM{i},2),CMd(1),CMd(2)));
                    else
                        CM{i} = CM{i}(:,1:CMd(2));
                    end
                end
            end
        end
    end
    
    if ~isempty(SCW) 
        if length(SCW) ~= n_lvls
            error(sprintf('Insufficient PCs supplied to level options. %d required, but only %d found.', n_lvls, length(SCW)));
        else
            for i = 1:n_lvls
                SCWd = [size(X,1),n_d{i}];
                if n_d{i} ~=0
                    if size(SCW{i},1) ~= SCWd(1) || size(SCW{i},2) < SCWd(2)
                        error(sprintf('Dimensions of %s ''CW'' arg are %d x %d. Expected %d x d, d >= %d',...
                            dLabels{i},size(SCW{i},1),size(SCW{i},2),CMd(1),CMd(2)));
                    else
                        SCW{i} = SCW{i}(:,1:SCWd(2));
                    end
                end
            end
        end
    end
    
    if ~isempty(vcomps)
        if length(vcomps) ~= n_lvls
            error(sprintf('Insufficient varcomps supplied. %d levels were specified, but only %d varcomps were found',n_lvls,length(varcomps)));
        else
            for i = 1:length(vcomps)
                if size(vcomps{i},1) ~= size(X,1) || size(vcomps{i},2) ~= size(X,2)
                    error(sprintf('Varcomps must be %d x %d, but level %d varcomp was %d x %d.', size(X,1), size(X,2), i, size(vcomps{i},1), size(vcomps{i},2)));
                end
            end
        end
    end
    
    % construct MLPCA arguments from available information
    % Note this is used both by PCA and varcomps calls.
    pcaOpts = {};
    for i = 1:length(dLabels)
        pcaOpts = {pcaOpts{:}, dLabels{i}, id{i}, n_d{i}};
    end
    
    %% run MLPCA
    tstart = tic; 
    if isempty(SCW)
        [CM, SCW] = mlpca(X,pcaOpts{:});
        if verbose
            fprintf('Completed MLPCA in %.2fs\n',toc(tstart));
        end
    end
    
    lambda_cnt = zeros(length(SCW),1);
    for i = 1:length(SCW)
        lambda_cnt(i) = sum(var(SCW{i}) > 10^-9);
    end
    if fitintercept{1}
        df = size(X,1) - 1;
    else
        df = size(X,1);
    end
    if df < sum(lambda_cnt) + size(cell2mat(covs),2) + 1
        warning('Cannot fit requested %d total PCA dims and %d additional covariates to data with %d degrees of freedom.',sum(lambda_cnt),size(cell2mat(covs),2),df);
        del = sum(lambda_cnt) + size(cell2mat(covs),2) + 1 - df;
        warning('Dropping %d dims from bottom most level.',del);
        lambda_cnt(end) = lambda_cnt(end) - del;
    end
    for i = 1:n_lvls
        vars = min(n_d{i},lambda_cnt(i));
        if vars < n_d{i}
            warning('Requested %d dimensions for level %i but only %d returned by PCA. Using %d',n_d{i},i,lambda_cnt(i),lambda_cnt(i));
            n_d{i} = vars;
            CM{i} = CM{i}(:,1:n_d{i});
            SCW{i} = SCW{i}(:,1:n_d{i});
        end
    end
    
    %% orthogonalize within and between components, 
    % moving any shared variance to the deapest levels of the hierarchy 
    % (e.g. any between variance that is also varying within subject gets 
    % moved to within subject components)
    for i = 1:n_lvls % preprocess
        SCW{i} = double(SCW{i});
        CM{i} = double(CM{i});
    end
    % qr seems to orthogonalize eigenvectors iteratively and columnwise 
    % from left to right. On each iteration, a column is orthogonalized
    % relative to all columns left of it. So let's stack deepest
    % levels on the left, and most superficial levels on the right so that
    % qr preserves the deepest eigenvectors while preferentially modifying
    % more superficial layers.
    S = cell2mat(fliplr(SCW')); % merge all loadings
    C = cell2mat(fliplr(CM')); % merge all eigenvectors
    [q,r] = qr(C,0); % orthogonalize eigenvectors
    rotS = S*r'; % map loadings to new eigenvectors
    % reload loadings and orthogonal eigenvalues back to original CM and 
    % SCW matricies, overwriting old values with orthogonalized values.
    lvl_idx = [];
    for i = 1:n_lvls
        if n_d{i} > 0
            i_0 = 1+sum(cell2mat( n_d(1:(i-1)) ));
            i_f = sum(cell2mat( n_d(1:i) ));
            lvl_idx(i_0:i_f) = i*ones(n_d{i},1);
        end
    end
    lvl_idx = sort(lvl_idx,'descend');
    for i = 1:n_lvls
        if n_d{i} > 0
            SCW{i} = rotS(:,lvl_idx == i);
            CM{i} = q(:,lvl_idx == i);     
        end
    end 

    % MLPCA is a memory intensive step. It may not be possible to run
    % multiple simultaneously. However fitlme is a compute intensive step
    % and it is helpful to run multiple in parallel. The most efficient way
    % to run multiple MLPCRs then is to run MLPCA serially and fitlme in
    % parallel. Using spmd with job dependencies allows for this. The
    % following conditional will tell queued SPMD loops that this MLPCR's
    % MLPCA instance has completed, and that SPMD can initiate the next
    % MLPCR script. SPMD loops should have labRecieve() calls waiting.
    if ~isempty(depJobID)
        for i = 1:length(depJobID)
            labSend(1,depJobID);
        end
    end
    
    %% create table of data for fitlme
    col_lbls = cell(sum(cell2mat(n_d)),1); % this doesn't allocate enough cells if there are covariates
    idx = 1;
    for i = 1:n_lvls
        for j = 1:n_d{i}
            col_lbls{idx} = ['l', int2str(i), '_', int2str(j)];
            idx = idx+1;
        end
        for j = 1:length(covs_Labels{i})
            col_lbls{idx} = covs_Labels{i}{j};
            idx = idx+1;
        end
        col_lbls = col_lbls(1:(idx-1));
    end

    yOffset = mean(Y);
    painc = Y(:) - yOffset;

    design = [painc,id{:}];
    for i = 1:n_lvls
        design = [design,SCW{i}(:,1:n_d{i}),covs{i}];
    end
    try
        t = array2table(double(design),'VariableNames',{'painc',dLabels{:},col_lbls{:}});
    catch
        dLabels{:}
        col_labls{:};
        size(design)
        rethrow(e);
    end

    %% consruct model formula for fitlme
    % fixed effects
    if fitintercept{1}
        formula = 'painc ~ 1+';
    else
        formula = 'painc ~ -1+';
    end
    for i = 1:n_lvls
        for j = 1:n_d{i}
            formula = [formula, 'l', int2str(i), '_', int2str(j), '+'];
        end
        for j = 1:size(covs{i},2)
            formula = [formula, covs_Labels{i}{j}, '+'];
        end
    end
    
    % random effects
    for i = 2:(lmelvls+1) %for a given level of the model
        formula = [formula, '('];
        for j = i:(n_lvls) %iterate through this level and all lower levels
            for k = 1:n_d{j} % adding every component loading
                formula = [formula, 'l', int2str(j), '_', int2str(k), '+'];
            end
            for k = 1:size(covs{j},2)
                formula = [formula, covs_Labels{j}{k}, '+'];
            end
        end
        if fitintercept{i}
            formula = [formula, '1|', dLabels{i}, ')+']; 
        else
            % this I think is buggy, so trying an alternative
            % delete trailing '+'
            %del = strfind(formula,'+');
            %formula = [formula(1:del(end)-1), '|', dLabels{i}, ')+'];
            if endsWith(formula,'+')
                formula = [formula(1:end-1), '-1|', dLabels{i}, ')+'];
            else
                formula = [formula, '-1|', dLabels{i}, ')+'];
            end
        end
    end
    % delete trailing '+'
    del = strfind(formula,'+');
    formula = formula(1:del(end)-1);
    if verbose
        fprintf('%s\n',formula);
    end
    
    %% invoke fitlme    
    %lme = fitlme(t,'painc ~ bt1+...+btn+ wi1+...+win + (wi1+...+win+1|sid)');
    %Isotropic covariance for speed. Think of it as a regularization step if
    %it makes you uneasy. Diagonal covariance may also work but be slower.
    tic; 
    if n_lvls > 1
        try
            lme = fitlme(t,formula,lmeopts{:});
        catch %try both fitting methods to see if either works
            try 
                lme = fitlme(t,formula,lmeopts{:},'FitMethod','ML');
            catch
                lme = fitlme(t,formula,lmeopts{:},'FitMethod','REML');
            end
        end
    else
        lme = fitlm(t,formula,lmeopts{:});
    end
    if verbose
        fprintf('Fit regression model in %.2fs\n',toc);
    end

    %% project model back to voxel space
    coef = lme.Coefficients.Estimate;

    % Xx is used later
    [Xx, Xy] = size(X);
    
    % initialization is necessary in case between/within dimensions requested == 0
    Intercept = cell(1,n_lvls+1);
    weights = cell(1,n_lvls+1);
    for i = 1:n_lvls
        Intercept{i} = double(0);
        weights{i} = double(zeros(Xy,1));
    end
    B = cell(1,n_lvls);
    
    for i = 1:n_lvls
        if n_d{i} > 0
            % Index of first relevant coefficient:
            % first coef is (maybe) fixed effect intercept, so skip it, start with
            % second index + accounting for any prior levels
            i_0 = 1+fitintercept{1}+sum(cell2mat( n_d(1:(i-1)) ));
            % Index of last relevant coefficient:
            % first coef is intercept, then count up all levels including
            % this one
            i_f = fitintercept{1}+sum(cell2mat( n_d(1:i) ));
            
            B{i} = coef(i_0:i_f);
            weights{i} = double(CM{i}*B{i});
            Intercept{i} = 0 - double(mean(SCW{i}*B{i}));
            % Notes on Intercept:
            % should be 0 for top level, because we centered Y. Anything
            % above zero needs to be subtracted out, anything below zero
            % should be added in, so let's just compute where it falls and
            % subtract the dif.
            % for subsequent levels the Intercept variable just holds an
            % adjustment relative to the next level prior. So, 
            % level 1: Intercept{1}
            % level 2: Intercept{1} + Intercept{2}
            % level 3: sum(Intercept{1:3})
            % etc.
            % All should be zero, but just in case I've gotten the theory
            % wrong, we can compute them.
        end
    end
    
    weights{end} = sum(cell2mat(weights),2);
    Intercept{end} = yOffset;
    if verbose
        fprintf('Total time to generated MLPCR weight maps in %.2fs\n',toc(tstart));
    end
end
