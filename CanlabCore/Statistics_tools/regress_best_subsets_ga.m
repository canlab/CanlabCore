function [wh_predictors, betas, b, stat] = regress_best_subsets_ga(X, Y)
% :Usage:
% GA-based best subsets regression
%
% ::
%
%     [wh_predictors, betas, b_subset, stat_subset] = regress_best_subsets_ga(X, Y)
%
% :Inputs:
%
%   **Y:**
%        is outcome data
%
%   **X:**
%        is predictor matrix
%
%   **wh_predictors:**
%        is the primary outcome -- it is vector of which predictors to include in the model
%
% the objective criterion is AIC

seeds = get_seeds(X, Y);

% Example: get_AIC(X, Y, seeds(:, 1))

% The fitness function is higher when AIC is lower
% The 1/AIC transformation would introduce a nonlinear weighting function;
% so would -log(AIC)
% log(1/AIC) is approximately linear.

%fitness_fcn = @(wh_preds)  1 ./ get_AIC(X, Y, wh_preds)
fitness_fcn = @(wh_preds)  log(1 ./ get_AIC(X, Y, wh_preds));

% for i = 1:size(seeds, 2)
%     a(i) = get_AIC(X, Y, seeds(:, i));
%     f(i) = fitness_fcn(seeds(:, i));
% end

% Run GA

example_vec = {seeds(:, 1)}; % this is 'inputs', what will be passed to the objective function (ofun) in the GA

for i = 1:size(seeds, 2), gaseeds{i} = seeds(:, i); end

[best_params,fit,beff,in,isconverged] = tor_ga(300, 50, example_vec, fitness_fcn, 'discrete', 'seeds', gaseeds, 'genconverge', 5);

wh_predictors = logical(best_params{1});

 Xs = X(:, wh_predictors);
 [b, dev, stat] = glmfit(Xs, Y); % subset plus intercept first
    
 betas = zeros(size(X, 2) + 1, 1);
 betas(logical([1; wh_predictors])) = b;  % add the intercept always
 
end


function aic = get_AIC(X, Y, wh_preds)

    X = X(:, logical(wh_preds));
    [nobs, npreds] = size(X);
    
    [b, dev, stat] = glmfit(X, Y);
    RSS = stat.resid' * stat.resid;

    % AIC = 2k - 2ln(L),  k = npreds, L = likelihood
    % AIC = 2k + n(ln(2piRSS/n) + 1)
    %
    % simplifying to relative values (we need only relative AIC, not
    % absolute:
    
    npreds = sum(wh_preds) + 1;  % include intercept
    aic = 2 * npreds + nobs * log(RSS);
    
end


function seeds = get_seeds(X, Y)
    
[nobs, npreds] = size(X);
z = false(npreds, 1);

% seed the GA with stepwise regression results and univariate results (top
% 1, 2, ... k)

[b, dev, stat] = glmfit(X, Y);
b = abs(b(2:end));
t = abs(stat.t(2:end));

% Set of candidate subsets with descending order of abs magnitude of betas
[bs, wh] = sort(b, 1, 'descend');
bseeds = false(npreds);

for i = 1:npreds
    my_subset = z; 
    my_subset(wh(1:i)) = true;
    bseeds(:, i) = my_subset;
end

% Set of candidate subsets with descending order of abs magnitude of
% t-values (last one is redundant with bseeds)
[ts, wh] = sort(t, 1, 'descend');
tseeds = false(npreds, npreds - 1);
for i = 1:npreds - 1
    my_subset = z; 
    my_subset(wh(1:i)) = true;
    tseeds(:, i) = my_subset;
end

% single-predictor seeds with each variable
eseeds = logical(eye(npreds));

% stepwise seeds
pentervals = [.05:.05:.95];
for i = 1:length(pentervals)
    [b, se, p, inmodel] = stepwisefit(X, Y, 'penter', pentervals(i), 'display', 'off');
    fseeds(:, i) = inmodel';
end

seeds = [bseeds tseeds eseeds fseeds];

% random seeds
numseeds = max(100, 200 - size(seeds, 2)); % at least 100 seeds
rseeds = rand(npreds, numseeds);
rseeds = rseeds > .5;

seeds = [seeds rseeds];

end

