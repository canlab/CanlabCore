function result = lasso_rocha(y, x, options)
% -------------------------------------------------------------------------
% function result = lasso(y, x, options)
% -------------------------------------------------------------------------
% PURPOSE:
% This is an implementation of the LARS algorithm for getting the LASSO 
% path in MATLAB.
%
% It can handle constant columns and linear dependent design matrices
% -------------------------------------------------------------------------
% INPUTS:
% y:     a n x 1 vector containing the values of the dependent variable
% x:     a n x p matrix containing the values of the predictors
% options: a structure containing options for the lasso
%          - max_steps (default = 'Inf'): number of steps until stopping
%          - show_progress (default = 0): boolean flagging whether progress should be output to screen (1) or not (0) 
% -------------------------------------------------------------------------
% OUTPUTS:
% A structure containing the following fields:
% nbetas:      the values of the coefficients for x and y normalized 
% betas:       the values of the coefficients in the original scales
% intercepts:  the values of the intercepts in the original scales
% xscale:      the constant by which each x was normalized
% yscale:      the normalizing constant for y
% ymean:       the mean of the y vector
% xmean:       the means of the columns of x
% npenalty:  the L1 norm of the coefficients along the path for the 
%              normalized Xs;
% penalty:   the L1 norm of the coefficients along the path for the
%              unnormalized Xs
% lambda:     values of the regularization parameter along the path
%              equals ||nx'*(ny-nx'*beta(\lambda)||_{\infty}
% nx:          a nxp design matrix with the normalized predictors
% ny:          a nx1 vector containing the values of the response variable
% sample_size: number of observations used to trace path
% dfs:         Zou and Hastie (2004)'s unbiased estimate of the
%              degrees of freedom along the path
% -------------------------------------------------------------------------
% Author: Guilherme V. Rocha and Peng Zhao
%         Department of Statistics
%         University of California, Berkeley
%         gvrocha@stat.berkeley.edu, gvrocha@gmail.com
% 2006/09
% -------------------------------------------------------------------------
%
% See also LASSO_CV, LASSO_COEFFICIENTS, BLASSOL2

% 0.0 Checking inputs:
% =========================================================================
if nargin<3
  options = [];
end;

% 0.1 Filling in options defaults:
% =========================================================================
if(~isfield(options, 'normalize_data'))
  options.normalize_data = 1;
end;

if(~isfield(options, 'max_steps'))
  options.max_steps = Inf;
end;

if(~isfield(options, 'show_progress'))
  options.show_progress = 0;
end;

if(~isfield(options, 'trace'))
  options.trace = 0;
end;

if(~isfield(options, 'trace'))
  options.trace = 0;
end;

% 0.1 Set options:
% =========================================================================
normalize_data  = options.normalize_data;
trace           = options.trace;
max_steps       = options.max_steps;
show_progress   = options.show_progress;

if nargin<2
  error('No predictors provided.');
  return;
end;

% 0.1 Starting up variables:
% =========================================================================
my_eps                         = 1e-8;                                     % numbers smaller than this will be considered zer
k                              = size(x, 2);                               % this is the number of regressors
n_obs                          = size(x, 1);                               % this is the number of observations
all_indices                    = (1:k)';

% 0.2 Normalizing both y and x:
% =========================================================================
if(normalize_data)
  xmeans = mean(x);
  xscale = sqrt(n_obs)*std(x);
  nx     = mrdivide((x-repmat(xmeans, n_obs, 1)), sparse(1:k,1:k,xscale,k,k,k));
  
  ymeans = mean(y);
  yscale = sqrt(n_obs)*std(y);
  ny     = full((y-repmat(ymeans, n_obs, 1))*sparse(1:1,1:1,1./yscale,1,1,1));
else
  xmeans = zeros(1, k);
  xscale = ones(1,k);
  nx     = x;
  
  ymeans = 0;
  yscale = 1;
  ny     = y;
end;

% 0.3 Starting up "control" variables:
% =========================================================================
constant_idxs      = find(xscale <= my_eps);
ignore_list        = constant_idxs;
variables_list     = sort(setdiff(all_indices, constant_idxs));            % list of "actual" variables (they DO vary)
nx(:, ignore_list) = zeros(n_obs, size(ignore_list, 2));                   % Set the constant columns to zero in normalized scale
actions            = -ignore_list;                                         % Registers the exclusion of the ignored variables in the action vector

if(show_progress)                                                                  % If show_progress was requested, list ignored variables to screen
  for i = 1:size(ignore_list,2)
    disp(['Adding ' num2str(ignore_list(i),3) ' to ignore list - constant']);
  end;
end;

% 0.4 Computes "correlations" - that is all we need:
% =========================================================================
xtr            = nx'*ny;                                                   % xty is a vector with the correlations between y and the explanatory variables
xtxA           = [];                                                       % xtxA stores the X'X(:,A) matrix - built incrememntally for speed

% 1 Initial point of the path:
% =========================================================================
betas          = zeros(k, 1);                                              % Path starts at zero
beta           = zeros(k, 1);                                              % Path starts at zero
lambdas        = norm(xtr,Inf);
drop           = 0;                                                        % No variable has been just dropped in the beginning
drop_ids       = [];                                                       % No variable has been just dropped in the beginning
model_size     = 0;                                                        % No variable has been select4ed in the beginning
A              = [];                                                       % A contains the set of active variables - it is empty at the start of the path
Ac             = sort(find(~ismember((1:k), union(A, ignore_list))));      % Ac contains the set of inactive variables: variables not ignored but not in the active list - variables in Ac can be included in model, variables in ignore_list can't
R              = [];                                                       % R is the Cholesky decomposition of the X'X(A, A) matrix
n_steps        = 0;
dlambda        = -Inf;
if(trace)
  xtrs = xtr';
end;

% 2.     normal equations are not satisfied (lambda=0)
%    AND n_steps is smaller than the maximal
%    AND lambda is still decreasing noticeably
% =========================================================================
while((lambdas(end)>my_eps)&(n_steps<max_steps)&(~(dlambda>-my_eps)))
  n_steps     = n_steps+1;
  lambda      = lambdas(end);
  new_action = [];
  % 2.a Add or drop variables
  % =======================================================================
  % 2.a.i Add variables
  % =======================================================================
  if(~drop)
    oldA = A;
    new_indices = Ac(find(abs(xtr(Ac))>=lambda-my_eps));                   % find all directions with maximum correlation and 
    for i = 1:size(new_indices, 2)                                         % tries adding each one of them individually (addition is only done after check for singularity)
      [rankok, newR] = updateR(nx, A, new_indices(i), R, my_eps);
      % 2.a.i.1  If variable being tried causes X'X(A,A) to become singular
      % ===================================================================
      if(~rankok)
        ignore_list = [ignore_list new_indices(i)];                        % the variable that causes singularity is added to ignore list
        Ac          = sort(find(~ismember((1:k), union(A, ignore_list)))); % the variable that causes singularity is excluded from possible choices
        if(show_progress)                                                          % if tracing was requested, ignored variable is reported to screen
          disp(['Adding ' num2str(new_indices(i),3) ' to ignore list - LD']);
        end;
        new_action = [new_action -new_indices(i)];                         % "dropping" of the variable is reported to action list
        
      % 2.a.i.2 If variable being tried does not cause X'X(A,A) to be singular
      % ===================================================================
      else
        R           = newR;                                                % X'X(A,A) is updated;
        A           = [A new_indices(i)];                                  % Active set is updated
        Ac          = sort(find(~ismember((1:k), union(A, ignore_list)))); % Inactive set is updated (remember that ignore list is not added to this set)
        A_sj        = sign(xtr(A));                                        % Signs of correlations of the active set are updated
        new_action  = [new_action new_indices(i)];                         % Addition is included in action list
        model_size  = model_size+1;                                        % Model size is updated
        if(show_progress)
          disp(['Adding variable ' num2str(new_indices(i),3) ...           % If tracing was requested, addition is reported to screen
                ' to model, model size = ' num2str(model_size,3) ', '...
                ' lambda = ' num2str(lambda, '%f')]);
        end;
      end;
    end;
    diffA                            = setdiff(A, oldA);                   % gets added terms
    new_xtxA                         = nx'*nx(:,diffA);                    % compute new terms in xtxA
    xtxA(:, end+1:end+length(diffA)) = new_xtxA;                           % update xtxA
    actions                          = [actions diffA];                    % add action to action list
        
  % 2.a.ii If a variable has just been dropped, 
  % =======================================================================
  else
    drop_indices        = find(abs(beta(A))<my_eps);
    keep_indices        = ~ismember(1:length(A),drop_indices);
    drop_ids            = A(drop_indices);
    A                   = A(keep_indices);                                 % updates active set
    Ac                  = sort(find(~ismember((1:k), union(A, ignore_list))));% updates inactive set
    A_sj                = A_sj(keep_indices);                              % updates correlation signs for active set
    R                   = downdateR(R, drop_indices);                      % updates X'X(A,A)
    xtxA(:,1:length(A)) = xtxA(:,keep_indices);
    xtxA                = xtxA(:,1:length(A));
    if(show_progress)                                                              % if show_progress was requested, show drop list to screen
      for i = 1:size(drop_ids,2)
        model_size = model_size - 1;
        disp(['Dropping variable ' num2str(drop_ids(i),3) ...
              ' from model, model size = ' num2str(model_size, 3) ', ' ...
              ' lambda = ' num2str(lambda, '%f')]);
      end;
    end;
    actions             = [actions -drop_ids];                             % add action to action list
    drop                = 0;                                               % resets drop flag
  end;
  
  % 2.b Computes direction in which estimates should move
  % =======================================================================
  dbeta         = R\(R'\A_sj);                                             % that is the direction that keeps correlation equal in size in active set

  % 2.c Computes how much to move in that direction
  % =======================================================================

  % 2.c.i First check when variable will be added to model
  % =======================================================================
  n1 = lambda - xtr(Ac);
  d1 = 1-xtxA(Ac,:)*dbeta;                                                    

  n2 = lambda + xtr(Ac);
  d2 = 1+xtxA(Ac,:)*dbeta;                                                    
  
  v1 = n1./d1;                                                             % v1 contains the size of the coefficient that makes each variable in the inactive set as correlated with the residuals as the active ones for positive movements
  v2 = n2./d2;                                                             % v2 contains the size of the coefficient that makes each variable in the inactive set as correlated with the residuals as the active ones for negative movements
  
  gamma_inc  = min([v1(v1>=my_eps);v2(v2>=my_eps)]);                       % gets the maximum distance that makes an inactive variable as correlated to residuals as an active variable
  if(isempty(gamma_inc));gamma_inc=Inf;end;
  
  % 2.c.ii Then check which gamma makes variable get dropped out of model 
  % =======================================================================
  gammas_exc = -(beta(A)./dbeta);
  gamma_exc  = min(gammas_exc(gammas_exc>my_eps));
  if(isempty(gamma_exc));gamma_exc=Inf;end;
  
  % 2.c.iii gamma is the minimal positive number that causes one of 2.c.i
  %         and ii above
  %         if no such positive number exists: run to terminal point
  % =======================================================================
  gamma = min(union(gamma_exc, gamma_inc));
  if(isempty(gamma)|gamma>lambda)                                          % If this happens, tries to go all the way to ``unregularized'' solution
    gamma = lambda;
  end;
  
  % 2.d Set drop flag if needed
  % =======================================================================
  if(gamma_exc<gamma_inc)                                                  % If one of the active variables will have its coefficient reaching zero before any of the inactive ones become as correlated to the residuals as the active set
    drop             = 1;                                                  % sets active set
  end;

  % 2.e update cofficients, correlations
  % =======================================================================
  update    = zeros(k,1);
  update(A) = gamma*dbeta;                                                 % calcuates update for the remaining variables in the active set
  beta      = betas(:,end) + update;                                       % updates beta
  xtr       = xtr - xtxA*update(A);                                        % updates correlation with residuals

  % 2.f store cofficients, correlations, 
  % =======================================================================
  betas     = [betas beta];                                                % add new step to list
  lambdas   = [lambdas norm(xtr, Inf)];
  dlambda   = diff(lambdas(end-1:end));
  if(trace)
    xtrs = [xtrs;xtr'];
  end;
end;

% 3. Storing the results in a structure
%==========================================================================
not_ignore_list                 = setdiff(all_indices, ignore_list);
variable_list                   = setdiff(all_indices, constant_idxs);

% Path:
result.nbeta                    = betas';
result.beta(:,constant_idxs)    = zeros(size(betas, 2), size(constant_idxs, 2)); 
result.beta(:,variable_list)    = yscale*betas(variable_list,:)'*diag(1./xscale(variable_list));
result.intercept                = ymeans-(xmeans*result.beta')';

% Indices:
result.npenalty                 = sum(abs(betas))';
result.penalty                  = sum(abs(result.beta'))';
result.lambda                   = lambdas';
result.df                       = sum(abs(result.beta')>my_eps)';

% Settings:
options.my_eps                  = my_eps;
result.options                  = options;

% Data:
result.ny                       = ny;
result.nx                       = nx;
result.yscale                   = yscale;
result.xscale                   = xscale;
result.ymean                    = ymeans;
result.xmean                    = xmeans;
result.sample_size              = n_obs;

% Debug:
if(trace)
  result.xtrs                     = xtrs;
end;

%==========================================================================
% End of lasso function
%==========================================================================

function [rankok, newR] = updateR(nx, A, new_indices, oldR, my_eps)
%==========================================================================
% function [rankok, newR] = updateR(xtx, xtR, oldR, my_eps)
%==========================================================================
% This function is used to update the Cholesky decomposition of the 
% X'X(A,A) matrix when a new variable is added to the model
% It checks whether the addition causes X'X(A,A) to become singular and
% returns a flag indicating that
%==========================================================================
xtR = nx(:,A)'*nx(:,new_indices);
xtx = nx(:,new_indices)'*nx(:,new_indices);

if(isempty(oldR))
  rankok = 1;
  newR   = chol(xtx);
  return;
end;
r   = (oldR')\xtR;
rpp = xtx - r'*r;
if(rpp <= my_eps)
  rankok = 0;
  newR = [];
else
  rankok = 1;
  newR   = [oldR r;zeros(1, size(oldR,2)) sqrt(rpp)];
end

function newR = downdateR(oldR, p)
%==========================================================================
% function newR = downdateR(oldR, p)
%==========================================================================
% This function is used to recompute the Cholesky decomposition of the
% X'X(A,A) matrix when a group of variables is deleted from the model
% The elements on the left and above the deleted variable are unchanged
%==========================================================================
  if(size(oldR,1)~=size(oldR,2))
    error('oldR must be square');
  end;
  num_delcols = max(size(p));
  for i = 1:num_delcols
    curr_p = p(i);
    k = size(oldR,2);
    R_11 = oldR(1:curr_p-1, 1:curr_p-1);
    R_13 = oldR(1:curr_p-1, curr_p+1:k);
    oldR_23 = oldR(curr_p, curr_p+1:k);
    oldR_33 = oldR(curr_p+1:k, curr_p+1:k);
    R_33    = chol(oldR_23'*oldR_23 + oldR_33'*oldR_33);
    oldR    = [R_11 R_13;zeros(k-curr_p,curr_p-1) R_33];
    p       = p-1; %one column was deleted from oldR so update of indices needed...
  end;
  newR = oldR;
  
