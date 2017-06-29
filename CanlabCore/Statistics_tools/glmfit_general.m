function varargout = glmfit_general(Y, X, varargin)
% This function was designed to run a second-level GLS analysis on
% coefficients and variance estimates from a 1st-level model, but it may be
% used generally to implement weighted GLS with options to bootstrap or run
% a sign permutation test to get p-values.
%
% :Inputs:
%
%   **Y:**
%        data for n subjects on k variables
%
%        run tests on each column of y independently
%
%        if this is a 2nd-level analysis in a two-level model, Y is coefficients from 1st level (mu-hat)
%
%   **s2:**
%        estimates of variance (sigma-squared??) for k variables
%
%   **X:**
%        design matrix; include intercept (first) + predictors
%
% :Estimation method:
%  - 'unweighted' / 'weighted'
%  - 'robust' / 'nonrobust'
%  - 's2'
%
% NEW: followed by n-length cell array of var/cov mtx
% (xtxi*sigma2) estimates from previous level
%
% OLD: followed by n x k matrix of first-level variance estimates
%
% :Inference method:
%  - 't-test'
%  - 'bootstrap'
%  - 'signperm'
%
%
% :Optional Inputs:
%   **'names':**
%        names of each col. of y:  string followed by input
%
%   **'name'  : text string tag for analysis
%
%   **'verbose':**
%        print verbose output and tables
%
%   **'nresample':**
%        number of boot/sign perm samples to run; string followed by input
%
%   **'noint':**
%        do not require intercept to be first column of X.
%
% :Notes:
% Same as glmfit.m (2007a) without weights.'
%
% With weights, we estimate w and use Satterthwaite, whereas glmfit uses
% n - k for dfe (assumes correct weights are known.)
%
% :Outputs:
%
%   **stats.b_star:**
%        Weighted fixed effects estimate (if using 'weighted' option)
%
%   **stats.Y_star:**
%        Y_star';  % Empirical Bayes estimates; don't use for group
%        inference, but save;
%
%        if this is run in a context of a multi-level model, these are the individual subjects'
%        (first-level) Empirical Bayes slope estimates
%
% :Output to screen:
% Coeff is "beta", the effect magnitude
% STE is the standard error, a function of variance and df 
%   - T = beta / STE
%   - Z = beta / the STE for infinite df
%   - p
%
% :Examples:
% ::
%
%    Y = randn(100, 4); X = Y(:,1) + randn(100,1); X = [ones(size(X)) X];
%    stats = glmfit_general(Y, X, 'verbose');
%    stats
%    first_lev_var = rand(100, 4);
%
%    stats = glmfit_general(Y, X, 'weighted', first_lev_var, 'dfwithin', 20, ...
%    'verbose', 'names', {'DMPFC' 'ACC' 'VMPFC' 'SII'}, ...
%    'beta_names', {'Mean activity' 'Rating covariate'}, ...
%    'analysisname','My Sample ROI analysis');
%


% ..
%    Setup inputs, print info to screen if verbose
% ..

varargout = cell(1, nargout - 1);

%if ~exist('varargin', 'var') || isempty(varargin), varargin = {0}; end

[Y, X, names, analysisname, beta_names, ...
    robust_option, weight_option, inference_option, ...
    s2within, dfwithin, ...
    verbose, dosave, doplots, ...
    targetu, nresample, whpvals_for_boot, ...
    permsign] = ...
    setup_inputs(Y, X, varargin{:});

if verbose, totalt = clock; end


% -------------------------------------------------------------------------
% Set up efficient estimator of group mean; used in bootstrap, etc.
% MUST DEFINE these functions AFTER getting X
% -------------------------------------------------------------------------

% First level/Single level paths: set handle for mediation function to run
%if dorobust, mediationfun = @fast_robust_ab; else mediationfun = @fast_ols_ab; end

W = ones(size(Y)) ./ size(Y, 1);    % wts. sum to 1, b/c weighted mean function requires this

% First pass: W will always be all equal weights (all ones)
% Or, if unweighted, skip this step and proceed to main analysis
% -----------------------------------------------------------------

switch weight_option

    case 'unweighted'
        % do nothing; we already have equal weights  (IID case)
        btwn_var_est = [];

    case 'weighted'
        
        
        % Determine weights from 1st-level (prior level) variances and
        % variance estimates for glm at this level
        % ------------------------------------------------------------
        % FIRST PASS - unweighted
        %[b, s2between_ols] = GLScalc(Y, W, X, bforming_fcn);
        [b, s2between_ols] = scn_stats_helper_functions('gls', Y, W, X);


        % get weights function
        % redefine W
        
        % OLD: needs to be fixed...
        %[W, btwn_var_est] = get_weights_based_on_varcomponents(s2within, s2between_ols, dfwithin, verbose);
        
        % NEW: based on Raudenbush & Bryk (in progress...) s2within should
        % be series of V matrices in cell array
        
        if verbose, fprintf('R & B - type Empirical Bayes reweighting.\n'); end
        
        % Y is data, which are presumed to be estimates from a previous
        % level of analysis
        [Y_star, b, Vgam_hat, gam_t, btwn_var_est, W] = RB_empirical_bayes_params(Y', s2within);
        
        statsx.b_star = b';
        statsx.b_star_descrip = 'R & B reweighted fixed-effects estimates';
        statsx.Y_star = Y_star';  % individual subjects' (first-level) Empirical Bayes estimates; don't use for group inference, but save
        statsx.Y_star_descrip = 'Empirical Bayes estimates of data after re-weighting';
        
end
% Main estimation, with whatever weights we have
% Whether analysis is 'weighted' depends on contents of W
% (equal W for unweighted, variable W for weighted.)
% -----------------------------------------------------------------

% first get GLS betas and p-values, to determine weights and
% boot samples in bootstrap.
[b, s2between_ols, stats] = scn_stats_helper_functions('gls', Y, W, X);
stats.analysisname = analysisname;

if exist('statsx', 'var')
    FF = fieldnames(statsx);
    for i = 1:length(FF), stats.(FF{i}) = statsx.(FF{i}); end
end

switch inference_option
    case 't-test'
        %[b, s2between_ols, stats] = GLScalc(Y, W, X, bforming_fcn);
        % Nothing further to do.

    case 'bootstrap'
        % BOOTcalc
        stats = scn_stats_helper_functions('boot', Y, W, X, nresample, stats, whpvals_for_boot, targetu, verbose );
        stats.analysisname = [stats.analysisname '. Bootstrap inference'];

    case 'signperm'
        % SIGNcalc
        stats = scn_stats_helper_functions('signperm', Y, W, X, nresample, stats, permsign, verbose );
        stats.analysisname = [stats.analysisname '. Sign permutation inference for intercept'];

    otherwise
        error('Unknown inference_option.');

end

% Save output structure
stats.btwn_var_est = btwn_var_est;
stats.W = W;

varnames = {'analysisname', 'names', 'beta_names', 'robust_option', 'weight_option', 'inference_option', 's2within', 'dfwithin', ...
    'targetu', 'nresample', 'whpvals_for_boot', 'permsign', 'Y', 'X'};

stats.inputOptions = create_struct(varnames);

varargout{1} = stats;




if verbose

    fprintf('Total time: %3.2f s\n', etime(clock, totalt));

    % print output table
    % using stats structure

    scn_stats_helper_functions('print', stats)

end




% _________________________________________________________________________
%
%
%
% * Inline functions
%
%
%
%__________________________________________________________________________


    function newstruct = create_struct(varnames)

        newstruct = struct();

        for i = 1:length(varnames)
            eval(['newstruct.(varnames{i}) = ' varnames{i} ';']);
        end

    end



end  % end main function






















% _________________________________________________________________________
%
%
%
% * Sub-functions
%
%
%
%__________________________________________________________________________




% -------------------------------------------------------------------------
% Setup inputs, print info to screen if verbose
% -------------------------------------------------------------------------
function [Y, X, names, analysisname, beta_names, robust_option, weight_option, inference_option, ...
    s2within, dfwithin, ...
    verbose, dosave, doplots, ...
    targetu, nresample, whpvals_for_boot, ...
    permsign] = ...
    setup_inputs(Y, X, varargin)

% Initial compliance checks
% ----------------------------------------------------------------

% remove NaNs
[nout, Y, X] = nanremove(Y, X);

[n, k] = size(X);
nvars = size(Y, 2);
nobs_tmp = size(Y, 1);

% Check sizes
if nobs_tmp ~= n, error('Y and X must have same number of rows!'); end

if any(strcmp(varargin, 'noint'))
    % No intercept
else
    % Check intercept
    if ~(all(X(:, 1) == X(1)))
        error('First column of X must be an intercept column. (e.g., all ones)');
    end
end

beta_names = cell(1, k);
for i = 1:k
    beta_names{i} = sprintf('Predictor %02d', i);
end


% Defaults
% ----------------------------------------------------------------

% General Defaults
names = cell(1, nvars);     % variable names, columns of Y
for i = 1:nvars, names{i} = ['Y' num2str(i)]; end

analysisname = 'GLS analysis.';

% Estimation Defaults
robust_option = 'no';               % robust IRLS; 'no' 'yes'
weight_option = 'unweighted';       % 'weighted' 'unweighted'
s2within = [];
dfwithin = [];

% Inference defaults
inference_option = 't-test';        % 't-test' 'bootstrap' 'signperm'

% Display control defaults
verbose = 0;                % verbose output
dosave = 0;                 % save figures at end
doplots = 0;                % make plots
savefilename = 'glmfit_general_output.txt';

% Bootstrap defaults
targetu = .20;              % proportion contribution of boot procedure to p-value
nresample = 1000;         % initial bootstrap samples
whpvals_for_boot = 1:size(Y,2);   % indices of p-values, the min of which is used to determine boot samples needed
% lower p-values require more boot samples
% for the p-vals to be meaningful.

% Sign perm defaults
permsign = [];              % empty: setup new sign permutation indices
% if entered: keep same permutation matrix
% across repeated calls (much faster! but
% reduces accuracy in simulations!)

% Inputs
% ----------------------------------------------------------------
for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            % General Defaults
            case 'names', names = varargin{i+1};  varargin{i+1} = [];
            case 'analysisname', analysisname = varargin{i+1}; varargin{i+1} = [];
            case 'beta_names', beta_names = varargin{i+1}; varargin{i+1} = [];

                % Estimation Defaults
            case 'robust', robust_option = 'yes';
            case {'weight', 'weighted', 'var', 's2'}, weight_option = 'weighted'; s2within = varargin{i + 1};
            case 'unweighted'  % already got this, do nothing
                
            case {'dfwithin'}, dfwithin = varargin{i + 1};

                % Inference defaults
            case {'boot1', 'boot', 'bootstrap'}, inference_option = 'bootstrap';
            case {'sign perm', 'signperm', 'sign'}, inference_option = 'signperm';
            case {'t-test', 'ttest'}, inference_option = 't-test';


                % Display control defaults
            case 'plots', doplots = 1;
            case 'noplots', doplots = 0;
                
            case {'dosave', 'save', 'saveplots'}, dosave = 1;
            case 'nosave', dosave = 0;
            case 'verbose', verbose = 1;
            case 'noverbose', verbose = 0;
            case {'savefile', 'savefilename'}, savefilename = varargin{i + 1}; varargin{i+1} = [];

                % Bootstrap defaults
            case {'nresample' 'bootsamples', 'nperms'}, nresample = varargin{i+1};
            case {'pvals', 'whpvals_for_boot'}, whpvals_for_boot = varargin{i+1};
            case 'targetu', targetu = varargin{i+1}; %**** not used in output yet!

                % Sign perm defaults
            case {'permsign'}, permsign = varargin{i+1};

            case 'intercept'  % do nothing; because we pass this option in from calling function; avoid warning
            otherwise
                fprintf('Warning! Unknown input string option: %s', varargin{i});

        end
    end
end

% start diary, if appropriate
if dosave
    if verbose
        fprintf('Saving output in file %s\n', savefilename)
    end
    diary(savefilename);

end

% Check for conflicting input arguments
% ***IMPLEMENT LATER***

if strcmp(weight_option, 'weighted')
    if isempty(s2within) || isempty(dfwithin)
        error('You must enter ''s2'', s2within, ''dfwithin'', df within as input arguments to use weighted option.');
    end
end


if strcmp(robust_option, 'yes')
    warning('off', 'stats:statrobustfit:IterationLimit');
    if verbose, disp('Note: Turning off iteration limit warning for robustfit.'); end
end

if verbose
    fprintf('Analysis description: %s \n', analysisname);
    fprintf('GLS analysis\n\nObservations: %3.0f, Predictors: %3.0f\n', n, k);

    fprintf('Outcome names: ');
    for i = 1:nvars
        fprintf('%s\t',names{i});
    end
    fprintf(1,'\n');

    fprintf('Weighting option: %s\nInference option: %s\n\n', weight_option, inference_option);

    nms = {'No' 'Yes'};

    fprintf('Other Options:\n\tPlots: %s\n\tRobust: %s\n\tSave: %s\n\tOutput file: %s\n\t', ...
        nms{doplots+1}, robust_option, nms{dosave+1}, savefilename);
end


% Bootstrap / signperm specific options:
% ***implement later***

end






% -------------------------------------------------------------------------
% Get weights from combination of within-ss and btwn-ss variances
% -------------------------------------------------------------------------
% get_w
% naive estimate of weights, based on 1 / (s2within + s2between)
% runs OLScalc to get s2between

% -------------------------------------------------------------------------
% For multilevel model: Get weights
% 1st-level standard errors & variance components
% -------------------------------------------------------------------------
function [W, btwn_var_est] = get_weights_based_on_varcomponents(s2within, s2between, dfwithin, verbose)

if verbose, fprintf('\nEstimating variance components based on OLS variance estimates\n'); end


% ANOVA approach.  MLE estimate; not ReML.
% Use SST = SSB + SSW. Solve for SSB.
%
% dfwithin is within-subjects degrees of freedom, e.g., n. obs per
% subject - no. params est. per subject. e.g., N - k within each subject
%   % * must be col. vector
%
% n = number of subjects
% nvars, number of tests
% totaln, total observations; subjects x obs per subject


[n] = size(s2within, 1);        % obs. x variables

% make a col. vector
if max(size(dfwithin)) == 1
    dfwithin = ones(n, 1) .* dfwithin;

elseif isrow(dfwithin)
    dfwithin = dfwithin';

end

totaln = sum(dfwithin);                 % total n obs. within subjects summed over subjects      %%sum(n-1);

SST = totaln .* s2between;              % 1 x nvars, one est. for each test (column of Y)

SSW = dfwithin' * s2within;             % 1 x nvars, within-subjects sum sq. errors. dfwithin must be row vector

btwn_var_est = (SST - SSW) / ( n-1 );

btwn_var_est = max(0, btwn_var_est);  % restrict to positive values.


%stats2.std_beta = stats2.ste .* wistats.avg_sx ./ wistats.avg_sy;
%stats2.btwn_prop = stats2.std_beta ./ (stats2.std_beta + mean(wistats.ste));

W = 1  ./ ( s2within + repmat(btwn_var_est, n, 1) );

% normalize by sum so weights sum to 1, for consistency
% DOES NOT affect stats.  or scaling of betas.

W = W ./ repmat(sum(W), n, 1);

end

%






% % % %
% % % %
% % % % % -------------------------------------------------------------------------
% % % % % GLS estimation of betas
% % % % % -------------------------------------------------------------------------
% % % %
% % % % function b = glsfunction(Y, W, X)
% % % %
% % % %     [n, k] = size(X);
% % % %     nvars = size(Y, 2);
% % % %
% % % %     b = zeros(k, nvars);
% % % %
% % % %
% % % %     if k == 1 && all(X == X(1))
% % % %
% % % %         % intercept only; fast computation
% % % %         % Weighted mean function:
% % % %         % Very efficient when there are no predictors other than the intercept
% % % %         % Faster than looping, and faster than using mean() for equal weights
% % % %
% % % %         b = diag(W'*Y)';
% % % %
% % % %
% % % %     else
% % % %         % predictors; need full computation
% % % %         % Full GLS function, needed if there are predictors
% % % %
% % % %         for i = 1:nvars
% % % %
% % % %             Wi = diag(W(:, i));
% % % %
% % % %             b(:, i) = inv(X' * Wi * X) * X' * Wi * Y(:, i);
% % % %
% % % %         end
% % % %
% % % %     end
% % % %
% % % % end
% % % %
% % % %
% % % %
% % % % % -------------------------------------------------------------------------
% % % % % GLS estimation of betas, ols variances, and full stats if 3rd output
% % % % % is requested
% % % % % -------------------------------------------------------------------------
% % % %
% % % % function [b, s2between_ols, stats] = GLScalc(Y, W, X, bforming_fcn)
% % % %     % OLScalc
% % % %     % calculate group betas, group variance
% % % %     % Optional: full inference outputs (takes longer)
% % % %
% % % %
% % % %     b = bforming_fcn(Y, W);
% % % %
% % % %     % Optional additional computations, in case we want to return just b
% % % %     % and go really fast (for bootstrapping)
% % % %
% % % %     if nargout > 1
% % % %
% % % %         [n, k] = size(X);
% % % %         nvars = size(Y, 2);
% % % %
% % % %         % --------------------------------------
% % % %         % * Residuals
% % % %         % --------------------------------------
% % % %
% % % %         e = Y - X * b;         % residuals
% % % %
% % % %         %
% % % %         % OLS residual variance
% % % %         % If bootstrapping or permuting, use this to get weights
% % % %
% % % %         dfe_ols = (n - k) .* ones(1, nvars);
% % % %
% % % %
% % % %         s2between_ols = (1 ./ dfe_ols .* sum(e .^ 2));  % Estimates of Sigma^2 for each col. of Y
% % % %
% % % %     end
% % % %
% % % %     if nargout > 2
% % % %
% % % %         % --------------------------------------
% % % %         % * Degrees of freedom for each column of Y
% % % %         % --------------------------------------
% % % %         dfe = dfe_ols;
% % % %
% % % %
% % % %         % Get design-related component of STE (includes n)
% % % %         % replace dfe with Sattherwaite approx. if necessary
% % % %
% % % %         Xste = zeros(k, nvars);
% % % %
% % % %         isweighted = 0;
% % % %
% % % %         % ***bug here: recomputes s2between too many times...diff
% % % %         % numbers...
% % % %         for i = 1:nvars
% % % %             if all(W(:, i) == W(1, i))
% % % %                 % weights are equal
% % % %                 invxvx = inv(X' * X);
% % % %                 Xste(:, i) = diag(invxvx);                  % design-related contribution to variance, k x nvars
% % % %
% % % %                 % --------------------------------------
% % % %                 % * Residual variance
% % % %                 % --------------------------------------
% % % %                 %s2between(i) = diag(e' * e)' ./ dfe(i);               % var for each col of Y, 1 x nvars
% % % %                 s2between(i) = s2between_ols(i);
% % % %
% % % %             else
% % % %                 % all weights are not equal
% % % %                 isweighted = 1;
% % % %
% % % %                 Wi = diag(W(:, i));                         % Wi = V^-1, inverse of cov.matrix
% % % %
% % % %                 invxvx = inv(X' * Wi * X);
% % % %
% % % %                 Xste(:, i) = diag(invxvx);                  % design-related contribution to variance, k x nvars
% % % %                 R = Wi.^.5 * (eye(n) - X * invxvx * X' * Wi);        % Residual inducing matrix
% % % %
% % % %                 Q = R * inv(Wi);                            % Q = RV
% % % %                 dfe(i) = (trace(Q).^2)./trace(Q * Q);       % Satterthwaite approximation for degrees of freedom
% % % %
% % % %                 % --------------------------------------
% % % %                 % * Residual variance
% % % %                 % --------------------------------------
% % % %                 e = R * Y(:, i);                            % weighted residuals
% % % %                 s2between(i) = diag(e' * e)' ./ dfe(i);               % var for each col of Y, 1 x nvars
% % % %             end
% % % %
% % % %         end
% % % %
% % % %
% % % %
% % % %         % --------------------------------------
% % % %         % * Standard errors of coefficients
% % % %         % --------------------------------------
% % % %         sterrs =  ( Xste .* repmat(s2between, k, 1) ) .^ .5;
% % % %
% % % %         % -------------------------------------------------------------------------
% % % %         % Get statistic structure from OLS regression, including p-values and conf. intervals
% % % %         % -------------------------------------------------------------------------
% % % %
% % % %
% % % %         stats.mean = b(1, :);           % intercept;  mean response
% % % %         stats.mean_descrip = 'Intercept of each col. of Y; (mean response if predictors are centered)';
% % % %
% % % %         stats.b = b;
% % % %         stats.b_descrip = 'betas (regression coefficients), k predictors x nvars';
% % % %
% % % %         stats.var = s2between;
% % % %         stats.var_descrip = 'Residual variance of each col. of Y';
% % % %
% % % %         stats.ste = sterrs;
% % % %         stats.ste_descrip = 'Std. error of each beta for each col. of Y, k predictors x nvars';
% % % %
% % % %         stats.t = b ./ sterrs;
% % % %
% % % %         stats.dfe = dfe;
% % % %         stats.dfe_descrip = 'error DF for each col. of Y, Satterthwaite corrected if necessary;  1 x nvars';
% % % %
% % % %
% % % %
% % % %         stats.e = e;
% % % %         if ~isweighted
% % % %             stats.e_descrip = 'unweighted (OLS) residuals';
% % % %         else
% % % %             stats.e_descrip = 'weighted residuals (resid. from weighted GLS model.)';
% % % %         end
% % % %
% % % %         for i = 1:k
% % % %
% % % %             stats.p(i, :) = min(1, (2 .* (1 - tcdf(abs(stats.t(i, :)), stats.dfe))));
% % % %
% % % %             stats.p(i, :) = max(stats.p(i, :), eps);
% % % %
% % % %         end
% % % %
% % % %         stats.p_descrip = 'Two-tailed p-values';
% % % %
% % % %     end
% % % %
% % % %
% % % %
% % % %
% % % %
% % % %
% % % % end
% % % %

