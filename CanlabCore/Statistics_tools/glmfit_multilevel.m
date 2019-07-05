function stats = glmfit_multilevel(Y, X1, X2, varargin)
% :Usage:
% ::
%
%     stats = glmfit_multilevel(Y, X1, X2, varargin)
%
% Mixed-effects models differ in their assumptions and implementation
% details. glmfit_multilevel fits regressions for individual 2nd-level
% units (e.g., participants), and then uses a precision-weighted least
% squares approach to model group effects, treating 1st-level units as a
% random effect. This is appropriate when 2st-level units are participants
% and 1st-level units are observations (e.g., trials) within participants. 
% glmfit_multilevel was designed with this use case in mind.
% 
% In addition, glmfit_multilevel:
% - Requires enough 1st-level units to fit a separate model for each
% 2nd-level unit (participant). If this is not the case, other models
% (igls.m, LMER, etc.) are preferred.
%
% - glmfit_multilevel  is conservative in the sense that the degrees of 
% freedom in the group statistical test is always based on the number of 
% subjects - 2nd-level parameters. 
% The df is never higher than the sample size, which you would have with 
% mixed effects models that estimate the df from the data. This causes
% problems in many other packages, particularly when there are many 1st-level
% observations and they are uncorrelated, resulting in large and undesirable 
% estimated df. 
%
% - The correlations across 1st-level observations are not measured, and 
% 1st-level obs are assumed to be IID. This is valid when generalizing across 
% 2nd-level units, but may not be fully efficient (powerful) if 1st-level 
% units are correlated.
%
% :Inputs:
%
%   **Y:**
%        Is data in either:
%           -cell array, one cell per subject.
%            Column vector of subject outcome data in each cell.
%           -Matrix
%            One column per subject, with vector of subject outcome
%            data in that column
%
%   **X1 and X2:**
%        are first and 2nd level design matrices
%          - X1 in cell array, one cell per subject
%            design matrix for each subject in each cell.  
%            *columns must code for the same variable for all subjects*
%
%         - X2 in rect. matrix
%         - can be empty (intercept only)
%
% E.g., with one 2nd-level predictor:
% stats = glmfit_multilevel(Y, X1, X2, ...
% 'names', {'Int' 'Temp'}, 'beta_names', {'2nd-level Intercept (overall group effect)' '2nd-lev predictor: Group membership'});
%
% :Output:
%
%   **stats:**
%        is structure with results.
%        Intercept is always added as first column!
%        (do not add intercept to input predictors)
%
% See glmfit_general.m for varargin variable input options.
%
% :Examples:
% ::
%
%    len = 200; sub = 20;
%    x = zeros(len,sub);
%    x(11:20,:) = 2;                   % create signal
%    x(111:120,:) = 2;
%    c = normrnd(0.5,0.1,sub,1);       % slope between-subjects variations
%    d = normrnd(3,0.2,sub,1);         % intercept between-subjects variations
%    % Create y: Add between-subjects error (random effects) and measurement noise
%    % (within-subjects error)
%    for i=1:sub, y(:,i) = d(i) + c(i).*x(:,i) + normrnd(0,0.5,len,1);
%    end;
%
%    for i = 1:size(y, 2), YY{i} = y(:, i); end
%    for i = 1:size(y, 2), XX{i} = x(:, i); end
%
%    % one-sample t-test, weighted by inv of btwn + within vars
%    stats = glmfit_multilevel(YY, XX, [], 'verbose', 'weighted');
%
%    statsg = glmfit_multilevel(y, x, covti, 'names', {'L1 Intercept' 'L1 Slope'},...
%    'beta_names', {'Group Average', 'L2_Covt'});
%
% :Input Options:
%
% General Defaults
%  - case 'names', Names of first-level predictors, starting
%    with 'Intercept', in cell array
%
%  - case 'analysisname', analysisname = varargin{i+1}; varargin{i+1} = [];
%  - case 'beta_names', beta_names = Names of 2nd-level predictors, starting
%    with 'Intercept', in cell array
%
% Estimation Defaults
%  - case 'robust', robust_option = 'yes';
%  - case {'weight', 'weighted', 'var', 's2'}, weight_option = 'unweighted';
%
% Inference defaults
%  - case {'boot1', 'boot', 'bootstrap'}, inference_option = 'bootstrap';
%  - case {'sign perm', 'signperm', 'sign'}, inference_option = 'signperm';
%  - case {'t-test', 'ttest'}, inference_option = 't-test';
%
% Display control defaults
%  - case 'plots', doplots = 1; plotstr = 'plots';
%  - case 'noplots', doplots = 0; plotstr = 'noplots';
%
%  - case {'dosave', 'save', 'saveplots'}, dosave = 1; savestr = 'save';
%  - case 'verbose', verbose = 1; verbstr = 'verbose';
%  - case 'noverbose', verbose = 0; verbstr = 'noverbose';
% 
%  - case {'savefile', 'savefilename'}, savefilename = varargin{i + 1}; varargin{i+1} = [];
%
% Bootstrap defaults
%  - case 'nresample', nresample = varargin{i+1};
%  - case {'pvals', 'whpvals_for_boot'}, whpvals_for_boot = varargin{i+1};
%
% Sign perm defaults
%  - case {'permsign'}, permsign = varargin{i+1};
%



%    Programmer's notes:
%    9/2/09: Tor and Lauren: Edited to drop NaNs within-subject, and drop
%    subject only if there are too few observations to estimate.
% ..

  if ~iscell(Y)
    N = size(Y, 2);
    for i = 1:N
      YY{i} = Y(:, i); 
    end
    Y = YY;
    clear YY;
  end

  if ~iscell(X1)
    N2 = size(X1, 2);
    if N ~= N2, error('Sizes of X and Y do not match'); end
    for i = 1:N
      XX{i} = X1(:, i); 
    end
    X1 = XX;
    clear XX
  end
  
  N = length(Y);
  if N ~= length(X1)
    error('Enter one cell per subject for each of X and Y');
  end




 % first level: SETUP
 % -------------------------------------------------------------------
 % set up first-level X matrix (sample)
  if any(strcmp(varargin, 'noint')) % no-intercept version
    X1tmp = X1{1};
  else
    X1tmp = setup_X_matrix(X1{1}); % intercept first
  end

  k = size(X1tmp, 2);             % num predictors; assumed to be the same!!


 % Second level: SETUP
 % Need to remove 2nd-level units with NaN data at first level
 % -------------------------------------------------------------------
  wh_omit = false(1, N);
  for i = 1:N
     %if any(isnan(Y{i})) || any(isnan(X1{i}(:))), wh_omit(i) = 1; end
    
    can_be_nans = length(Y{i}) - k - 1;  % up to this many can be NaN, still leaving 1 degree of freedom
    if can_be_nans < 0, warning('Warning: you might be overparameterized!  Seems like you have more predictors than observations'); end
    if sum(isnan(Y{i}) | any(isnan(X1{i}), 2)) > can_be_nans
      wh_omit(i) = 1; 
    end
  end

  if any(wh_omit)
    if isempty(X2), X2 = ones(N, 1); end
    
    Y(wh_omit) = [];
    X1(wh_omit) = [];
    X2(wh_omit, :) = [];
    N = length(Y);
  end


  beta = zeros(k, N);
  sterr = zeros(k, N);
  t = zeros(k, N);
  p = zeros(k, N);
  dfe = zeros(1, N);
%phi = zeros(arorder, N);

 % first level: ESTIMATE
 % -------------------------------------------------------------------
  for i = 1:N
    [beta(:, i), sterr(:, i), t(:, i), p(:, i), dfe(:, i), phi(:,i), V{i}] = first_level_model(Y{i}, X1{i}, varargin{:});
       % V{i} is var/cov matrix (xtxi)*sigmasq
  end

  varnames = {'beta' 't' 'p' 'dfe' 'phi'};
  first_level = create_struct(varnames);
  first_level.ste = sterr;

 % second level: Finish SETUP and ESTIMATE
 % -------------------------------------------------------------------
 % set up second-level X matrix: intercept first
  X2 = setup_X_matrix(X2, beta(1,:)');

% set up second-level options
% names of outcomes become beta_names here b/c second level test on 1st
% level betas
  [beta_names1, analysisname, beta_names2, robust_option, weight_option, inference_option, ...
   verbose, dosave, doplots, ...
   verbstr, savestr, plotstr, ...
   targetu, nresample, whpvals_for_boot, ...
   permsign] = ...
  setup_inputs(beta', X2, varargin{:});
  switch weight_option
    case 'weighted'
  % Note: R & B-style : replaced sterr' with V

      stats = glmfit_general( ...
	      beta', X2, ...
	      'analysisname', analysisname, 'names', beta_names1, 'beta_names', beta_names2, ...
	      verbstr, savestr, plotstr, ...
	      weight_option, V, inference_option, 'dfwithin', dfe', ...
	      'targetu', targetu, 'nresample', nresample, ...
	      'whpvals_for_boot', whpvals_for_boot, 'permsign', permsign);

    case 'unweighted'

      stats = glmfit_general( ...
	      beta', X2, ...
	      'analysisname', analysisname, 'names', beta_names1, 'beta_names', beta_names2, ...
	      verbstr, savestr, plotstr, ...
	      weight_option, inference_option, ...
	      'targetu', targetu, 'nresample', nresample, ...
	      'whpvals_for_boot', whpvals_for_boot, 'permsign', permsign);

    otherwise
      error('Problem with weight_option. Please select either weighted or unweighted.')
  end

  stats.first_level = first_level;

  if doplots
    scn_stats_helper_functions('xyplot', X1, Y, weight_option, 'names', beta_names1, 'nostats');
    xlabel('X'); ylabel('Y');
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
      eval(['newstruct.' varnames{i} ' = ' varnames{i} ';']);
    end
  end

end %End of glmfit_multilevel 

function [b, sterr, t, p, dfe, phi, V] = first_level_model(y, X, varargin)

% defaults
% -------------------------------------------------------------------
verbose = 0;
verbstr = 'noverbose';
arorder = 0;                    % or Zero for no AR
interceptstr = 'intercept';

% optional inputs
% -------------------------------------------------------------------
for varg = 1:length(varargin)
    if ischar(varargin{varg})
        switch varargin{varg}

            % reserved keywords
            case 'verbose all', verbose = 1; verbstr = 'verbose';
            case 'verbose', % do nothing
            case {'ar', 'arorder'} , arorder = varargin{varg+1};
                %otherwise, disp(['Unknown input string option: ' varargin{varg}]);

            case 'noint', interceptstr = 'noint';
        end
    end
end

k = size(X, 2) + 1;

[whnan X y] = nanremove(X, y);

if isempty(X)
    % no data
    
    [b, t, p, sterr] = deal(NaN * ones(k, 1));
    dfe = deal(NaN);
    if arorder, phi = NaN * ones(arorder, 1); else phi = NaN; end
    V = [];
    
    return
    
end

% set up X matrix: intercept first
if ~strcmp(interceptstr, 'noint')
    X = setup_X_matrix(X, y);
end

if arorder
    % if we have missing observations or redundant columns, let's
    % regularize a bit so we can still estimate this, using a ridge prior
    % The degree of regularization is arbitrary.
    % V is used to estimate variance components and re-weight.
    if rank(X) < size(X, 2)
        disp('WARNING! RANK DEFICIENT.  THIS FUNCTION WILL RETURN AN ERROR.')
        X = [X; eye(size(X, 2)) ./ size(X, 1)];
        if ~strcmp(interceptstr, 'noint'), X(:, 1) = 1; end
        y = [y; ones(size(X, 2), 1) .* nanmean(y)];
    end
    
    [t, dfe, b, phi, sigma, sterr] = fit_gls(y, X, [], arorder);
    p = 2 * (1 - tcdf(abs(t), dfe));  % two-tailed

    V = inv(X' * X) * sigma .^ 2; % Var/Cov mtx, Precision^-1, used in weighted est. and empirical bayes
    
else
    % if we have missing observations or redundant columns, let's
    % regularize a bit so we can still estimate this, using a ridge prior
    % The degree of regularization is arbitrary.
    if rank(X) < size(X, 2)
        disp('WARNING! RANK DEFICIENT.  THIS FUNCTION WILL RETURN AN ERROR.')
        X = [X; eye(size(X, 2)) ./ size(X, 1)];
        if ~strcmp(interceptstr, 'noint'), X(:, 1) = 1; end
        y = [y; ones(size(X, 2), 1) .* nanmean(y)];
    end
    
    stats = glmfit_general(y, X, verbstr, interceptstr);
    t = stats.t; dfe = stats.dfe; b = stats.beta; phi = NaN; sterr = stats.ste; p = stats.p;

    V = inv(X' * X) * stats.var; % Var/Cov mtx, Precision^-1, used in weighted est. and empirical bayes
end

end

function X = setup_X_matrix(X, y)
  % set up X matrix: intercept first
  [n, k] = size(X);
  if n == 0
    n = size(y, 1);
  end
  equal_x = false(1, k);
  for i = 1:k
    if all(X(:, i) == X(1, i))
    % weights are equal
      equal_x(i) = 1;
    end
  end
  if any(equal_x)
    error('Warning: some columns of X have no variance.  Do not enter intercept in X; it will be added automatically as the first predictor.');
    X(:, equal_x) = [];
  end
  X = [ones(n, 1) X];
end

% -------------------------------------------------------------------------
% Setup inputs, print info to screen if verbose
% THIS IS FOR SECOND-LEVEL MODEL -- NOT FIRST
% -------------------------------------------------------------------------
function [names, analysisname, beta_names, robust_option, weight_option, inference_option, ...
    verbose, dosave, doplots, ...
    verbstr, savestr, plotstr, ...
    targetu, nresample, whpvals_for_boot, ...
    permsign] = ...
    setup_inputs(Y, X, varargin)

% Initial compliance checks
% ----------------------------------------------------------------
[n, k] = size(X);
nvars = size(Y, 2);
nobs_tmp = size(Y, 1);

% Check sizes
if nobs_tmp ~= n, error('Y and X must have same number of rows!'); end

% Check intercept
if ~(all(X(:, 1) == X(1)))
    error('First column of X must be an intercept column. (e.g., all ones)');
end

beta_names = cell(1, k);
for i = 1:k
    beta_names{i} = sprintf('2nd-level B%02d', i);
end


% Defaults
% ----------------------------------------------------------------

% General Defaults
names = cell(1, nvars);     % variable names, columns of Y
for i = 1:nvars, names{i} = ['1st-level B' num2str(i)]; end

analysisname = 'Second Level of Multilevel Model ';

% Estimation Defaults
robust_option = 'no';               % robust IRLS; 'no' 'yes'
weight_option = 'unweighted';       % 'weighted' 'unweighted'

% Inference defaults
inference_option = 't-test';        % 't-test' 'bootstrap' 'signperm'

% Display control defaults
verbose = 1;                % verbose output
dosave = 0;                 % save figures at end
doplots = 0;                % make plots
plotstr = 'noplots';
savestr = 'nosave';
verbstr = 'verbose';

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
            case {'weight', 'weighted', 'var', 's2'}, weight_option = 'weighted';

                % Inference defaults
            case {'boot1', 'boot', 'bootstrap'}, inference_option = 'bootstrap';
            case {'sign perm', 'signperm', 'sign'}, inference_option = 'signperm';
            case {'t-test', 'ttest'}, inference_option = 't-test';

                % Display control defaults
            case {'plot', 'plots'}, doplots = 1; plotstr = 'plots';
            case {'noplots', 'noplot'}, doplots = 0; plotstr = 'noplots';

            case {'dosave', 'save', 'saveplots'}, dosave = 1; savestr = 'save';
            case 'verbose', verbose = 1; verbstr = 'verbose';
            case 'noverbose', verbose = 0; verbstr = 'noverbose';

            case {'savefile', 'savefilename'}, savefilename = varargin{i + 1}; varargin{i+1} = [];

                % Bootstrap defaults
            case 'nresample', nresample = varargin{i+1};
            case {'pvals', 'whpvals_for_boot'}, whpvals_for_boot = varargin{i+1};


                % Sign perm defaults
            case {'permsign'}, permsign = varargin{i+1};

            case 'intercept'
                
            otherwise
                fprintf('Warning! Unknown input string option: %s', varargin{i});

        end
    end
end

% all the checking, etc. is done in glmfit_general


end
