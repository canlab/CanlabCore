function S = xval_SVM(varargin)
%  SVM for repeated-measures (within-person) classification, with repeated cross-val and nested cross-val options
%
% :Usage:
% ::
%
% S = xval_SVM(X, Y, id, varargin)
%
% Steps and features:
% -------------------------------------------------------------------------
% - Select holdout sets, keeping images from the same id together and stratifying on the outcome to be predicted
% - Fit overall model and get "sanity check" accuracy
% - Cross-validation with standard a priori hyperparameter choices
% - Plot scores and ROC curves
% - Nested cross-val with hyperparameter optimization [optional]
%   * Estimate model performance accounting for search across hyperparams
%   * Updates S.yfit, S.dist_from_hyperplane_xval, S.class_probability_xval
%   * Updates crossval_accuracy, Y_within_id, scores_within_id, scorediff, crossval_accuracy_within, classification_d
%   * Retrain to obtain estimate of best hyperparameters for future testing (update S.modeloptions)
% - Repeat cross-validation of optimized model with different random splits [optional]
% - Fit model on full dataset with final hyperparameters to get predictor weights (betas, b)
% - Bootstrap model parameter estimates (w) and get P-values, FDR-corrected significant features
% - Print output and text compatible with report-generation
%
% Notes:
% -------------------------------------------------------------------------
% - Uses Matlab's Stats/ML Toolbox SVM object, fitcsvm, and hyperparameter
% optimization. Other options, e.g., fitclinear, may be better for some
% situations (large num. variables)
%
% - Single-interval accuracy from cross-val and ROC_plot may differ because
% ROC_plot chooses a new score threshold that maximizes overall balanced
% accuracy. Forced-choice accuracy should be identical.
%
% Optimizing hyperparameters: 
% HyperparameterOptimizationOptions include defaults of:
% - 5-fold CV (not grouped by id) without repartitioning 
% - Bayesian optimization, using best estimates from smooth function
% - maximizes accuracy, so may need large(ish) samples to be effective.
%
% - works on linear classifiers only, but can explore nonlinear models
% with a small change in the code (now commented out). This could be added
% as an optional flag as well.
%
% - If optimizing hyperparameters AND repeating cross-validation, accuracy
% estimates will use nested cross-validation, and can take a long time to run
%
%
% Class probability estimates: 
% - crossval returns cross-validated yfit (class predictions), scores (dist_from_hyperplane_xval),
% class_probability_xval. Scores and class probabilties will diverge (may
% not be perfectly correlated) because the sigmoid scaling (Platt scaling) varies across folds.
% They will diverge more if there is no true signal.
%
% ..
%     Author and copyright information:
%
%     Copyright (C) 2020 Tor Wager
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
% ..
%
% :Inputs:
%
%   **X:**
%        obs x variables numeric matrix of predictors
%
%   **Y:**
%        obs x 1 numeric vector of outcomes, effects-coded (1, -1)
%
%   **id:**
%        obs x 1 numeric vector of grouping codes, e.g., for participants
%        all obs from each group will be included together in a training or test set
%
% :Optional Inputs:
%   **'doplot', [logical flag]:**
%        Create plots; default = true. 'noplot' to turn off.
%
%   **'doverbose', [logical flag]:**
%        Verbose output; default = true. 'noverbose' to turn off.
%
%   **'dooptimize', [logical flag]:**
%        Optimize hyperparameters; default = true. 'nooptimize' to turn off.
%
%   **'doprepeats', [num repeats]:**
%        Repeat cross-val with different partitions; default = 10.
%        Enter number of repeats, 'norepeats' to turn off.
%
%   **'modeloptions', [modeloptions cell]:**
%        Options for model structure and hyperparameters
%        Cell array of keyword-value pairs, as specified in fitcsvm (Matlab Stats/ML toolbox)
%        Default: {'KernelFunction', 'linear'}
%
% :Outputs:
%
%   **S:**
%        Structure with model output
%        of them (partially consistent with this function).
%                           Y: Actual (obs) outcome - should be 1, -1 for SVM analysis
%                        yfit: Predicted outcome - should be 1, -1 for SVM analysis
%                              Cross-validated, so can be used as series
%                              of predicted values based on brain
%                              measures (very useful!)
%                          id: Grouping variable for within-participant observations
%                      accfun: Function handle to get accuracy (single-interval)
%                           w: Model weights/betas, weights on input variables
%                      nfolds: Number of folds in holdout set
%                 cvpartition: Object with information about training/test sets
%                       teIdx: Testing IDs for each fold (holdout set)
%                       trIdx: Training IDs for each fold (holdout set)
%             dist_from_hyperplane_xval: Cross-validated distance perpendicular to class boundary
%                               higher = stronger prediction in favor of Class 1, lower = in favor of class -1.
%                               Useful! use as continuous measure of
%                               observation scores. Can calculate
%                               effect sizes from this, for example.
%       class_probability_xval: Cross-validated distance expressed as probabilities of Class 1, using Platt scaling
%            crossval_accuracy: Cross-validated accuracy (input options; no hyperparam opt)
%  classification_d_singleinterval: Classification effect size (input options; no hyperparam opt)
%crossval_accuracy_opt_hyperparams: Cross-validated accuracy with optimized hyper-parameters
%                  Y_within_id: Outcomes arranged by id, for within-person comparisons
%             scores_within_id: SVM scores arranged by id, for within-person comparisons
%                    scorediff: Within-person SVM scores arranged by id, for within-person comparisons
%     crossval_accuracy_within: Within-person cross-validated accuracy (input options; no hyperparam opt)
%      classification_d_within: Within-person classification effect size (input options; no hyperparam opt)
%                 S.boot_w_ste: Bootstrapped standard errors of model weights (feature-level)
%                S.boot_w_mean: Bootstrapped mean model weights (feature-level)
%                         S.wZ: Bootstrapped Z-scores of model weights (feature-level)
%                         S.wP: Bootstrapped P-values for individual model weights (feature-level)
%                 S.wP_fdr_thr: P-value threshold for FDR q < 0.05
%              S.boot_w_fdrsig: Logical vector for which weights are significant with FDR
%               S.w_thresh_fdr: Thresholded weights at FDR q < 0.05 based on bootstrapping (feature-level)
%                   S.SVMModel: ClasssificationSVM model object trained on full dataset, final chosen parameter set; for predicting in subsequent validation samples.
%                               Y-hat = (X/s)'* w + b
%                               ClassificationSVM objects store w, b, and s in the properties Beta, Bias, and KernelParameters.Scale, respectively.
%
% :Examples:
% ::
%
% Simulate some data with true signal, with two observations per person
% -------------------------------------------------------------------------
% n = 50; % participants
% k = 120; % features
% true_sig = [repmat(randn(1, k), n, 1); repmat(randn(1, k), n, 1)]; % First n are Class 1, second n are Class 2
% noise = 10 * randn(2 * n, k);         % Noise var >> true signal var
% X = true_sig + noise;                  % Predictors
% Y = [ones(n, 1); -ones(n, 1)];         % Outcome to classify, coded [1 -1]
% id = [(1:n)'; (1:n)'];                 % Grouping ID codes for participants
%
% S = xval_SVM(X, Y, id, 'nooptimize', 'norepeats', 'nobootstrap');    % Fastest for quick cross-validated performance
% S = xval_SVM(X, Y, id, 'nooptimize', 'norepeats', 'nboot', 100);     % Very quick test of bootstrapping with few samples (not for final models)
% S = xval_SVM(X, Y, id, 'nooptimize', 'dorepeats', 5, 'nobootstrap'); % Quick test of repeated cross-val, 5x
% S = xval_SVM(X, Y, id, 'nooptimize', 'nobootstrap');                 % Repeat cross-val only
% S = xval_SVM(X, Y, id, 'norepeats', 'nobootstrap');                  % Optimize with nested cross-val only
% S = xval_SVM(X, Y, id);                                              % Optimize, repeat with optimal model, bootstrap
%
% S = xval_SVM(X, Y, id, 'nooptimize', 'norepeats', 'nobootstrap', 'noverbose', 'noplot');  % Returns only the output S
%
% Train the same model without within-person observations (i.e., 1 observation per id):
% This does single-interval classification only
% 
% S = xval_SVM(X, Y, (1:length(Y))', 'norepeats', 'nobootstrap');
%
% :References:
%   See Mathworks functions
%
% :See also:
%   - fmri_data.predict, canlab_run_paired_SVM, other xval_ functions
%

% ..
%    Programmers' notes:
%    Created by Tor Wager, May 2020
%    ver 1.0 (no version number listed) - some bugs in optimization
%    ver 1.1 (Tor): fixed bugs, optimization changes and documentation, multiple small tweaks and some new plots
%             Updated to handle within-person or between-person obs only
%             ROC_plot should be at threshold 0 for single-interval
%             Get rid of plots for all folds in opt output
%             OptimizeHypParams returns point est. No need to rerun. Builds in 5-fold eval by default
%             Add ?Summary of nested cross-val accuracy with hyperparameter optimization?
%             Optim.: Allow control of verbosity and plots. No plots/verbose in nested eval.
%             Optim: use best estim, not best observed
%             Have a look at the data and add data plots. Sorted by cases
%   Future:
%             Enforce [1 -1] inputs
%             Troubleshoot prediction vs. class probability reversals
%                   % either scores or pscores appear to be inconsistent in terms of which column is class -1 and which is class 1. Check for instability in scores, and select pscores to be consistent with this on each fold.
%             Add permutation test
% 
% ..

%% ----------------------------------------------------------------------
% Version
% ----------------------------------------------------------------------

ver = 1.2;

%% ----------------------------------------------------------------------
% Parse inputs
% ----------------------------------------------------------------------
% Uses the inputParser object. Older schemes are below.

% Logical flags - parse specially because they are not name-value pairs
% Parser accepts only name/value pairs, so this code allows a workaround
% We must remove these inputs because they will cause inputparser to error
% ----------------------------------------------------------------------
wh = strcmp(varargin, 'noverbose'); if any(wh), doverbose1 = false; varargin(wh) = []; end
wh = strcmp(varargin, 'noplot'); if any(wh), doplot1 = false; varargin(wh) = []; end
wh = strcmp(varargin, 'noplots'); if any(wh), doplot1 = false; varargin(wh) = []; end
wh = strcmp(varargin, 'nooptimize'); if any(wh), dooptimize1 = false; varargin(wh) = []; end
wh = strcmp(varargin, 'norepeats'); if any(wh), dorepeats1 = 0; varargin(wh) = []; end
wh = strcmp(varargin, 'nobootstrap'); if any(wh), dobootstrap1 = 0; varargin(wh) = []; end

% The handling of parsing keyword-value pairs, checking attributes, and
% specifying default values is handled here, using Matlab's inputparser object:
ARGS = parse_inputs(varargin{:});

fn = fieldnames(ARGS);

for i = 1:length(fn)
    str = sprintf('%s = ARGS.(''%s'');', fn{i}, fn{i});
    eval(str)
end

% Replace logical flags because they will be overwritten by the parser
if exist('doverbose1', 'var'), doverbose = doverbose1; end
if exist('doplot1', 'var'), doplot = doplot1; end
if exist('dooptimize1', 'var'), dooptimize = dooptimize1; end
if exist('dorepeats1', 'var'), dorepeats = dorepeats1; end
if exist('dobootstrap1', 'var'), dobootstrap = dobootstrap1; end


%% Select holdout sets for outer loop
% - Keep images from the same id together
% - Stratify on the outcome to be predicted
% -------------------------------------------------------------------------

S = struct();   % Define structure for models and output; use variable names in fmri_data.predict when possible

S.Y = Y;
S.id = id;      % Subject grouping - not in fmri_data.predict

S.modeloptions = modeloptions;

S.accfun = @(Y, yfit) 100 .* nansum(Y == yfit) ./ sum(~isnan(Y));

if doverbose
    
    dashes = '---------------------------------------------------------------------------------';
    printhdr = @(str) fprintf('%s\n%s\n%s\n', dashes, str, dashes);
    
    printhdr(' ') %#ok<*UNRCH>
    printhdr(sprintf('xval_SVM ver %3.2f: Cross-validated Support Vector Machine classification', ver));
    printhdr(' ')
    
    disp(' ')
    disp(' ')
    printhdr('Selecting and analyzing cross-validation folds');
    
end


% Preliminary data plots
% -------------------------------------------------------------------------
if doplot
    prelim_data_plot(X, Y, id);
end

% Select train/test sets
% -------------------------------------------------------------------------

[S.trIdx, S.teIdx] = xval_stratified_holdout_leave_whole_subject_out(S.Y, S.id, 'doverbose', doverbose, 'doplot', doplot);
drawnow, snapnow

% Fit the overall a priori model: SVM with linear kernel
% -------------------------------------------------------------------------
if doverbose
    
    printhdr('Model training and cross-validation without optimization')
    fprintf('Training overall model. ')
    
    SVMModel = fitcsvm(X, S.Y, S.modeloptions{:}); %#ok<*NODEF>
    
    % non-cross-val accuracy: sanity check. should usually be 100% if training is overfit and n >> p
    label = predict(SVMModel, X);
    
    acc = S.accfun(S.Y, label);
    fprintf('Accuracy without cross-val: %3.0f%%\n', acc);
    
end

%% Cross-validation with standard a priori hyperparameter choices
% -------------------------------------------------------------------------
S.nfolds = length(S.trIdx);

% crossval returns cross-validated yfit (class predictions), scores (dist_from_hyperplane_xval),
% class_probability_xval. Scores and class probabilties will diverge (may
% not be perfectly correlated) because the sigmoid scaling (Platt scaling) varies across folds.
% They will diverge more if there is no true signal.

S = crossval(S, X, doverbose);

% Summarize accuracy
% -------------------------------------------------------------------------
% Given Y, yfit, cross-val scores (dist_from_hyperplane_xval), id
% Update crossval_accuracy,classification_d_singleinterval
% Y_within_id, scores_within_id, scorediff,
% crossval_accuracy_within, classification_d_within
S = update_accuracy_stats(S);

if doverbose
    
    fprintf('Accuracy (cross-validated): Single-interval: %3.0f%%, d = %3.2f\n', S.crossval_accuracy, S.classification_d_singleinterval);
    
    if S.mult_obs_within_person
        
        fprintf('                            Forced-choice within-person: %3.0f%% d = %3.2f\n', S.crossval_accuracy_within, S.classification_d_within);
        
    end
    
end

% Plot scores and ROC curves
% -------------------------------------------------------------------------
if doplot
    
    plot_scores_and_ROC(S);
    
    xval_plot_scores_vs_class_probability(S);
    
end

%% Nested cross-val with hyperparameter optimization
% -------------------------------------------------------------------------

if dooptimize
    
    if doverbose
        printhdr('Optimizing hyperparameters using nested cross-validation');
    end
    
    % Nested cross-validation (updates S.yfit, S.dist_from_hyperplane_xval,
    % S.class_probability_xval)
    % Estimate model performance accounting for search across hyperparams
    % ---------------------------------------------------------------------
    % Updates S.yfit, S.dist_from_hyperplane_xval, S.class_probability_xval
    S = crossval_nested(S, X, doverbose);
     
    % Given Y, yfit, cross-val scores (dist_from_hyperplane_xval), id
    % Update crossval_accuracy, Y_within_id, scores_within_id, scorediff, crossval_accuracy_within, classification_d
    S = update_accuracy_stats(S);
    
    % Retrain to obtain estimate of best hyperparameters for future testing (update S.modeloptions)
    % ---------------------------------------------------------------------
    % Modal parameters across folds is not
    % the best way to get the final hyperparameters. Instead, the model should
    % be re-fit using single cross-validation with all hparam options pick the
    % best relative model. Nested cross-val will estimate the accuracy of the
    % whole procedure, including the hparam search, but not return the final best
    % parameters.
    Mdl = fitcsvm(X, S.Y,...
        'OptimizeHyperparameters', {'BoxConstraint' 'Standardize'}, 'HyperparameterOptimizationOptions',...
        struct('AcquisitionFunctionName','expected-improvement-plus', 'Verbose', double(doverbose), 'ShowPlots', double(doplot))); %#ok<*IDISVAR>
    
    % Notes:
    % - Extract hyperparameter choices in a format that we can pass back in
    % to train or test a new model.
    % - We could use the best observed point, but it's preferred to use the
    % best estimated values from a smooth model fit to observed samples.
    
    best_modeloptions = convert_hyperparameter_choices_to_cell(Mdl, S.modeloptions);

    S.modeloptions = best_modeloptions;

    
    % Summarize accuracy
    % -------------------------------------------------------------------------
 
    if doverbose
        
        disp('Summary of nested cross-val accuracy with hyperparameter optimization')
        disp(' ')
        
        disp('Optimized hyperparameters for each fold')
        disp(S.hyperparams_by_fold)
        
        disp(' ')
        disp('Best options, stored in S.modeloptions:');
        
        disp(best_modeloptions)
        
        fprintf('Accuracy (nested cross-validation): Single-interval: %3.0f%%, d = %3.2f\n', S.crossval_accuracy, S.classification_d_singleinterval);
        
        if S.mult_obs_within_person
            
            fprintf('                            Forced-choice within-person: %3.0f%% d = %3.2f\n', S.crossval_accuracy_within, S.classification_d_within);
            
        end
        
    end
    
    % Plot scores and ROC curves
    % -------------------------------------------------------------------------
    if doplot
        
        plot_scores_and_ROC(S, 'Cross-validated SVM scores with hyperparameter optimization');
        
    end
    
    
end % do optimize

%% Repeat cross-validation of optimized model, if specified
% - if hyperparams optimized, then use "consensus params" for repeats
% - these are stored in S.modeloptions now, so will be used
% - If optimizing, re-do the first accuracy with the optimized
% hyperparameters
% -------------------------------------------------------------------------
if dorepeats > 1
    
    if doverbose
        disp(' ')
        fprintf('Repeating cross-val %d times. ', dorepeats)
    end
    
    % if dooptimize, startat = 1; else, startat = 2; end
    startat = 2;
    
    for i = startat:dorepeats     % do remainder of repeats without verbose.
        
        if doverbose, fprintf('%d ', i), end
        

        Sr = S;
        
        % Get new cvpartition (or outer loop if optimizing) 
        [Sr.trIdx, Sr.teIdx] = xval_stratified_holdout_leave_whole_subject_out(Sr.Y, Sr.id, 'doverbose', false, 'doplot', false);
        
        if dooptimize
            % If optimizing, repeat nested xval procedure
            % This can take a long time!
            
            Sr = crossval_nested(Sr, X, doverbose);
            
        else
            
            % Cross-validate
            Sr = crossval(Sr, X, false);
            
        end
        
        % Save accuracy
        % -------------------------------------------------------------------------
        % Given Y, yfit, cross-val scores (dist_from_hyperplane_xval), id
        % Update crossval_accuracy,classification_d_singleinterval
        % Y_within_id, scores_within_id, scorediff,
        % crossval_accuracy_within, classification_d_within
        
        Sr = update_accuracy_stats(Sr);

        S.crossval_accuracy(i) = Sr.crossval_accuracy; % Sr.accfun(Sr.Y, Sr.yfit);
        S.classification_d_singleinterval(i) = Sr.classification_d_singleinterval;
        
        if S.mult_obs_within_person
            
            S.crossval_accuracy_within(i) = Sr.crossval_accuracy_within; % Sr.accfun(Sr.Y, Sr.yfit);
            S.classification_d_within(i) = Sr.classification_d_within;

        end
        
    end
    
    if doverbose
        
        fprintf('Done!\n')
        
        fprintf('CV single-interval accuracy across %d reps, mean = %3.0f%%, std = %3.0f%%, min = %3.0f%%, max = %3.0f%%\n', dorepeats, ...
            mean(S.crossval_accuracy), std(S.crossval_accuracy), min(S.crossval_accuracy), max(S.crossval_accuracy));
        
        disp('Performance metrics saved in S.crossval_accuracy, classification_d_singleinterval')
        
        if S.mult_obs_within_person
            
            fprintf('\nCV forced_choice accuracy across %d reps, mean = %3.0f%%, std = %3.0f%%, min = %3.0f%%, max = %3.0f%%\n', dorepeats, ...
                mean(S.crossval_accuracy_within), std(S.crossval_accuracy_within), min(S.crossval_accuracy_within), max(S.crossval_accuracy_within));
            
             disp('Performance metrics saved in S.crossval_accuracy_within, classification_d_within')
            
        end
        
        
        disp(' ')
        
    end
    
end

%% Fit model on full dataset with selected options (final model if optimizing hyperparams) to get betas

S.SVMModel = fitcsvm(X, S.Y, S.modeloptions{:});
S.w = S.SVMModel.Beta;                                       % Model weights/betas

% Bootstrap betas
if dobootstrap
    
    if doverbose
        disp(' ')
        printhdr(sprintf('Bootstrapping model param estimates (b): %d samples', nboot));
    end
    
    rng('shuffle');
    
    S.boot_w = NaN .* zeros(length(S.w), nboot);
    
    for i = 1:nboot
        
        [Xb, Yb] = get_bootstrap_sample_grouped_by_id(S, X);
        
        SVMModel = fitcsvm(Xb, Yb, S.modeloptions{:});
        
        S.boot_w(:, i) = SVMModel.Beta;                                       % Model weights/betas
        
        
    end
    
    % Inference
    % (from fmri_data.predict)
    
    S.boot_w_ste = squeeze(nanstd(S.boot_w, 0, 2)); %1/20/16 add squeeze for multiclass case
    S.boot_w_mean = squeeze(nanmean(S.boot_w, 2)); %1/20/16 add squeeze for  multiclass case
    S.boot_w_ste(S.boot_w_ste == 0) = Inf;  % in case unstable regression returns all zeros
    
    S.wZ = S.boot_w_mean ./ S.boot_w_ste;  % tor changed from wmean; otherwise bootstrap variance in mean inc in error; Luke renamed to avoid confusion
    S.wP = 2 * (1 - normcdf(abs(S.wZ)));
    S.wP_fdr_thr = FDR(S.wP, .05);
    if isempty(S.wP_fdr_thr), S.wP_fdr_thr = -Inf; end
    
    S.boot_w_fdrsig = S.wP <= S.wP_fdr_thr;  % equals because can get exact vals in some cases...
    S.w_thresh_fdr = S.w;
    S.w_thresh_fdr(~S.boot_w_fdrsig) = 0;
    
    if doverbose
        
        disp('Summary of significant individual features (two-tailed)');
        
        Threshold = [.05 .01 .001 S.wP_fdr_thr]';
        
        for i = 1:length(Threshold)
            
            Num_Sig_Features(i, 1) = sum(S.wP <= Threshold(i));
            
        end
        
        t = table(Threshold, Num_Sig_Features);
        disp(t)
        
    end
end

%%
% permutations would go here...future project.
% consider exchangeability issue given within-person design...winkler and
% nichols...
% indx = randperm(size(X, 1));
% X = X(indx, :); Y = Y(indx); id = id(indx);
% S = xval_SVM(X, Y, id, 'nooptimize', 'norepeats');


end % main function


%%


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

% Subfunctions

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

function S = crossval(S, X, doverbose)

if doverbose, fprintf('..X-val, %d folds...', S.nfolds), end

for i = 1:S.nfolds
    
    if doverbose, fprintf('%d ', i), end
    
    % Fit to training data for this fold
    SVM_fold = fitcsvm(X(S.trIdx{i}, :), S.Y(S.trIdx{i}), S.modeloptions{:});
    
    % Apply holdout test set for this fold
    [label, score] = predict(SVM_fold, X(S.teIdx{i}, :));   % Get raw scores
    score = score(:, 2);                                    % Decision boundary is symmetrical
    
    S.dist_from_hyperplane_xval(S.teIdx{i}, 1) = score;     % Unscaled SVM scores
    
    S.yfit(S.teIdx{i}, 1) = label;                          % Predicted class, cross-val
    
    SVM_fold = fitPosterior(SVM_fold);                      % Works for fitcsvm, not fitclinear
    
    % Note:
    % either scores or pscores appear to be inconsistent in terms of which
    % column is class -1 and which is class 1. Check for instability in
    % scores, and select pscores to be consistent with this on each fold.
    if mean(label(score > 0)) < mean(label(score < 0))
        disp('WARNING!!! Scores for label 1 are < scores for label -1, scores are reversed!!! This should not happen. Check code/implementation.');
    end
    
    [~, pscore] = predict(SVM_fold, X(S.teIdx{i}, :));
    
    % which pscores correlate more positively with score? use this one.
     r = corr([score pscore]);
    [~, wh] = max(r(1, 2:3));
   
    pscore = pscore(:, wh);                                  % Platt scaling scores (class probability)
    
    S.class_probability_xval(S.teIdx{i}, 1) = pscore;
    
    
end

if doverbose, fprintf('Done!\n'), end

end % crossval


% -------------------------------------------------------------------------
% Nested cross-validation to optimize hyperparameters
% -------------------------------------------------------------------------


function S = crossval_nested(S, X, doverbose)

if doverbose, fprintf('..X-val, %d folds...', S.nfolds), end

for i = 1:S.nfolds
    
    if doverbose, fprintf('%d ', i), end
    
    % Optimize hyperparameters within this fold
    % -------------------------------------------
    
    % Option 1: Optimize all, including kernel (linear/nonlinear)
    
    %     Mdl = fitcsvm(X(S.trIdx{i}, :), S.Y(S.trIdx{i}),...
    %     'OptimizeHyperparameters','all','HyperparameterOptimizationOptions',...
    %     struct('AcquisitionFunctionName','expected-improvement-plus'));
    %
    %     S.bestp{i} = Mdl.HyperparameterOptimizationResults.XAtMinObjective; % best parameters - table
    
    % Option 2: Optimize specific hyperparameters, Slack param C
    % (BoxConstraint) and Standardize inputs, linear model only
    Mdl = fitcsvm(X(S.trIdx{i}, :), S.Y(S.trIdx{i}),...
        'OptimizeHyperparameters', {'BoxConstraint' 'Standardize'}, 'HyperparameterOptimizationOptions',...
        struct('AcquisitionFunctionName','expected-improvement-plus', 'Verbose', 0, 'ShowPlots', 0));
    
    % Notes: HyperparameterOptimizationOptions include defaults of:
    % - 5-fold CV (not grouped by id) without repartitioning 
    % - Bayesian optimization, using best estimates from smooth function
    
    % Add hyperparameter results to table
    % -------------------------------------
    % Use estimated rather than observed. This uses the smooth function estimated across observations

    opt_hyperparam_table = Mdl.HyperparameterOptimizationResults.XAtMinEstimatedObjective; 
    
    if i == 1
        S.hyperparams_by_fold = opt_hyperparam_table; % best parameters - table
        
    else
        S.hyperparams_by_fold(i, :) = opt_hyperparam_table;
    end
    
    fold_modeloptions = convert_hyperparameter_choices_to_cell(Mdl, S.modeloptions);
    
    S.fold_modeloptions{i} = fold_modeloptions;
    
    % *** REMOVE after checking
    % Re-train with final choices
    % -------------------------------------
%     fold_modeloptions = S.modeloptions;
%     
%     % Add optimal choices to the option set
%     for j = 1:size(opt_hyperparam_table, 2)
%         vname = opt_hyperparam_table.Properties.VariableNames{j};
%         
%         fold_modeloptions{end + 1} = vname;
%         paramval = opt_hyperparam_table.(vname)(1);
%         
%         % fix: categorical to text
%         if iscategorical(paramval)
%             paramval = char(string(paramval));
%             paramval = strcmp(paramval, 'true');
%             
%         end
%         
%         fold_modeloptions{end + 1} = paramval; %#ok<*AGROW>
%         
%     end
    
    % Fit to training data for this fold
    SVM_fold = fitcsvm(X(S.trIdx{i}, :), S.Y(S.trIdx{i}), fold_modeloptions{:});
    
    % Updates S.yfit, S. S.dist_from_hyperplane_xval, S.class_probability_xval
    
    % Apply holdout test set for this fold
    [label, score] = predict(SVM_fold, X(S.teIdx{i}, :));   % Get raw scores
    score = score(:, 2);                                    % Decision boundary is symmetrical
    
    S.dist_from_hyperplane_xval(S.teIdx{i}, 1) = score;     % Unscaled SVM scores
    
    S.yfit(S.teIdx{i}, 1) = label;                          % Predicted class, cross-val
    
    SVM_fold = fitPosterior(SVM_fold);                      % Works for fitcsvm, not fitclinear
    
    [~, pscore] = predict(SVM_fold, X(S.teIdx{i}, :));
    pscore = pscore(:, 2);                                  % Platt scaling scores (class probability)
    
    S.class_probability_xval(S.teIdx{i}, 1) = pscore;
    
end

if doverbose, fprintf('Done!\n'), end

end % crossval



% -------------------------------------------------------------------------
% Get scores
% -------------------------------------------------------------------------


function [scores_within_id, scorediff, d] = get_scores_within_id(S, varname)

myvar = S.(varname);
u = unique(S.id);

for i = 1:length(u)
    
    wh_id = S.id == u(i) & S.Y == -1;
    scores_within_id(i, 1) = nanmean(myvar(wh_id));
    
    wh_id = S.id == u(i) & S.Y == 1;
    scores_within_id(i, 2) = nanmean(myvar(wh_id));
    
end

scorediff = diff(scores_within_id')';
d = nanmean(scorediff) ./ nanstd(scorediff);

end %get scores


% -------------------------------------------------------------------------
% bootstrapping
% -------------------------------------------------------------------------


function [Xb, Yb] = get_bootstrap_sample_grouped_by_id(S, X)
% Given X, S.Y, S.id, return a bootstrap sample, keeping all obs from an id together

u = unique(S.id);
n = length(u);

wh = ceil(rand(n, 1) * n);  % Bootstrap sample of p ids, with replacement

[Xboot, Yboot] = deal(cell(n, 1));

for j = 1:n                 % add obs, allowing repeats
    
    wh_obs = ismember(S.id, u(wh(j)));
    
    Xboot{j} = X(wh_obs, :);
    Yboot{j} = S.Y(wh_obs);
    
end

Xb = cat(1, Xboot{:});
Yb = cat(1, Yboot{:});

end % function


% -------------------------------------------------------------------------
% optimized hyperparameter aggregation
% -------------------------------------------------------------------------

function  best_modeloptions = convert_hyperparameter_choices_to_cell(Mdl, other_modeloptions)
% Notes:
% - Extract hyperparameter choices in a format that we can pass back in
% to train or test a new model.
% - We could use the best observed point, but it's preferred to use the
% best estimated values from a smooth model fit to observed samples.

% opt_hyperparam_table = Mdl.HyperparameterOptimizationResults.XAtMinObjective;
opt_hyperparam_table = Mdl.HyperparameterOptimizationResults.XAtMinEstimatedObjective;

best_modeloptions = other_modeloptions;

% Add optimal choices to the option set
for j = 1:size(opt_hyperparam_table, 2)
    vname = opt_hyperparam_table.Properties.VariableNames{j};
    
    best_modeloptions{end + 1} = vname;
    paramval = opt_hyperparam_table.(vname)(1);
    
    % fix: categorical to text
    if iscategorical(paramval)
        paramval = char(string(paramval));
        paramval = strcmp(paramval, 'true');
        
    end
    
    best_modeloptions{end + 1} = paramval; %#ok<*AGROW>
    
end

end
    
% function best_modeloptions = get_modal_params_across_folds(S)
% 
% % This function is deprecated because modal parameters across folds is not
% % the best way to get the final hyperparameters. Instead, the model should
% % be re-fit using single cross-validation with all hparam options pick the 
% % best relative model. Nested cross-val will estimate the accuracy of the
% % whole procedure, including the hparam search, but not return the final best
% % parameters.
% 
% hbyfold = cat(1, S.fold_modeloptions{:});
% best_modeloptions = {};
% 
% for i = 1:size(hbyfold, 2)
%     
%     mydat = cat(1, hbyfold{:, i});
%     
%     if ischar(mydat) || iscategorical(mydat) || islogical(mydat)
%         
%         best_modeloptions{1, i} = mode(mydat);
%         
%     elseif isnumeric(mydat)
%         
%         best_modeloptions{1, i} = trimmean(mydat, 80);
%         
%     else
%         error('Unknown model option class! Extend this code.')
%         
%     end
%     
% end
% 
% end


% -------------------------------------------------------------------------
% accuracy
% -------------------------------------------------------------------------



function S = update_accuracy_stats(S)

S.crossval_accuracy = S.accfun(S.Y, S.yfit);

s1 = S.dist_from_hyperplane_xval(S.Y == 1);
s2 = S.dist_from_hyperplane_xval(S.Y == -1);

% d = difference / std pooled within, weighted by sample size
S.classification_d_singleinterval = (mean(s1) - mean(s2)) ./ sqrt( ( var(s1) .* length(s1) + var(s2) .* length(s2) ) ./ (length(s1) + length(s2)) );

% Do we have multiple obs within-person? If so return within-person stats
S.mult_obs_within_person = length(S.id) > length(unique(S.id));

if S.mult_obs_within_person
    
    varname = 'dist_from_hyperplane_xval';
    [scores_within_id, scorediff, d] = get_scores_within_id(S, varname);
    
    S.Y_within_id = get_scores_within_id(S, 'Y');
    
    S.scores_within_id = scores_within_id;
    S.scorediff = scorediff;
    S.crossval_accuracy_within = 100 * sum(scorediff > 0) ./ sum(~isnan(scorediff));
    S.classification_d_within = d;
    
end

end


% -------------------------------------------------------------------------
% plots
% -------------------------------------------------------------------------



function prelim_data_plot(X, Y, id)

create_figure('Data view', 1, 3);
imagesc(X); colorbar;
xlabel('Input features'); ylabel('Observation'); title('Predictor matrix (X)');
axis tight; set(gca, 'YDir', 'reverse');

subplot(1, 3, 2);
imagesc([Y scale(id)]);
axis tight; set(gca, 'YDir', 'reverse');
set(gca, 'XTick', [1 2], 'XTickLabel', {'Y' 'id'});
title('Outcome (Y) and id');

subplot(1, 3, 3);
r = corr(X');
[~, wh] = sort(Y);
imagesc(r(wh, wh));
n1 = sum(Y == -1);
n2 = sum(Y == 1);

h1 = drawbox(0, n1, 0, n1, 'k');
set(h1, 'FaceColor', 'none', 'EdgeColor', 'r', 'LineWidth', 2);

h1 = drawbox(n1, n2, n1, n2, 'k');
set(h1, 'FaceColor', 'none', 'EdgeColor', 'r', 'LineWidth', 2);

axis tight; set(gca, 'YDir', 'reverse');
% set(gca, 'XTick', [1 2], 'XTickLabel', {'Y' 'id'});
title('Inter-obs correlations sorted by Y [-1, 1]');

colorbar

drawnow, snapnow
end



function plot_scores_and_ROC(S, varargin)
% plot_scores_and_ROC(S, varargin) -> varargin = new title

create_figure('cross-val accuracy', 2, 2);
subplot(2, 2, 1); delete(gca); subplot(2, 2, 2); delete(gca);
axes('Position', [.13 .58 .77 .35]);
set(gca, 'FontSize', 16)
hold on

if S.mult_obs_within_person
    % Within-person scores
    
    n = size(S.scores_within_id, 1);
    x = [1:n; 1:n];
    plot(x(:, S.scorediff > 0), S.scores_within_id(S.scorediff > 0, :)', 'Color', [.3 .3 .3], 'LineWidth', 1);
    plot(x(:, S.scorediff < 0), S.scores_within_id(S.scorediff < 0, :)', 'r', 'LineWidth', 1);
    
    plot(x(1, :), S.scores_within_id(:, 1)', '^', 'Color', [.3 .5 1] / 2, 'MarkerFaceColor', [.3 .5 1],  'LineWidth', 1);
    plot(x(2, :), S.scores_within_id(:, 2)', 'v', 'Color', [1 .5 0] / 2, 'MarkerFaceColor', [1 .5 0],  'LineWidth', 1);
    
    xlabel(sprintf('ID, %3.0f%% single-interval acc, %3.0f%% forced-choice', S.crossval_accuracy, S.crossval_accuracy_within));
    
    disp('Black lines: Correct, Red lines: Errors. Red triangles: Scores for Class 1, Blue triangles: Scores for Class -1');
    
else
    
    % Scores, all between-person
    
    n = size(S.dist_from_hyperplane_xval, 1);
    x = 1:n;
    iscorrect = sign(S.dist_from_hyperplane_xval) == sign(S.Y);
    
    %plot(x(:, iscorrect), S.dist_from_hyperplane_xval(iscorrect, :)', 'o', 'MarkerSize', 10, 'Color', [.3 .3 .3], 'LineWidth', 1);
    plot(x(:, ~iscorrect), S.dist_from_hyperplane_xval(~iscorrect, :)', 'o', 'MarkerSize', 10, 'Color', [1 0 0], 'LineWidth', 1);
    
    h = plot_horizontal_line(0); set(h, 'LineStyle', '--');
    
    color1 = [.2 .8 .2];
    color2 = [.3 .3 1];
    
    plot(x(:, S.Y == 1), S.dist_from_hyperplane_xval(S.Y == 1, :)', 'v', 'MarkerSize', 6, 'Color', color1, 'MarkerFaceColor', color1, 'LineWidth', 1);
    plot(x(:, S.Y == -1), S.dist_from_hyperplane_xval(S.Y == -1, :)', '^', 'MarkerSize', 6, 'Color', color2, 'MarkerFaceColor', color2, 'LineWidth', 1);
    
    xlabel(sprintf('ID, %3.0f%% single-interval acc', S.crossval_accuracy));
    
    disp('Red cicles: Errors. Green triangles: Scores for Class 1, Blue triangles: Scores for Class -1');
    
    axis tight
    
end

ylabel('SVM Score');

title('Cross-validated SVM scores (no hyperparameter optimization)');
if ~isempty(varargin), title(varargin{1}); end


% ROC plot
subplot(2, 2, 3);
S.ROC_single_interval = roc_plot(S.dist_from_hyperplane_xval, logical(S.Y > 0), 'color', [.4 .4 .7], 'threshold', 0);
title('Single-interval ROC')
set(gca, 'FontSize', 16)

subplot(2, 2, 4);

if S.mult_obs_within_person
    
    % Paired forced-choice. Get complete cases - Remove NaNs id-wise
    outcomes = S.Y_within_id;
    scores = S.scores_within_id;
    [~, outcomes, scores] = nanremove(outcomes, scores);
    
    S.ROC_forced_choice = roc_plot(scores(:), logical(outcomes(:) > 0), 'color', [.4 .4 .7], 'twochoice');
    title('Forced-choice ROC')
    set(gca, 'FontSize', 16)
    
else
    % Scores, all between-person
    
    plot(S.Y(~iscorrect), S.dist_from_hyperplane_xval(~iscorrect, :)', 'o', 'MarkerSize', 10, 'Color', [1 0 0], 'LineWidth', 1);
    
    h = plot_horizontal_line(0); set(h, 'LineStyle', '--');
    
    color1 = [.2 .8 .2];
    color2 = [.3 .3 1];
    
    plot(ones(1, sum(S.Y == 1)), S.dist_from_hyperplane_xval(S.Y == 1, :)', 'v', 'MarkerSize', 6, 'Color', color1, 'MarkerFaceColor', color1, 'LineWidth', 1);
    plot(-ones(1, sum(S.Y == -1)), S.dist_from_hyperplane_xval(S.Y == -1, :)', '^', 'MarkerSize', 6, 'Color', color2, 'MarkerFaceColor', color2, 'LineWidth', 1);
    
    set(gca, 'XTick', [-1 1], 'XLim', [-1.5 1.5]);
    xlabel('True class');
    ylabel('Predicted class')
    
end

drawnow, snapnow

end



function xval_plot_scores_vs_class_probability(S)

% crossval returns cross-validated yfit (class predictions), scores (dist_from_hyperplane_xval),
% class_probability_xval. Scores and class probabilties will diverge (may
% not be perfectly correlated) because the sigmoid scaling (Platt scaling) varies across folds.
% They will diverge more if there is no true signal.

create_figure('xval scores vs. class probability estimates', 1, 2); 

plot(S.dist_from_hyperplane_xval, S.class_probability_xval, 'o');

for j = 1:length(S.teIdx)
    
    plot(S.dist_from_hyperplane_xval(S.teIdx{j}), S.class_probability_xval(S.teIdx{j}), 'o', 'MarkerFaceColor', rand(1, 3));

end

xlabel('Cross-validated SVM scores (colors are folds)');
ylabel('Cross-validated class prob estimates');

hh = plot_vertical_line(0); set(hh, 'LineStyle', '--');
hh = plot_horizontal_line(.5); set(hh, 'LineStyle', '--');

subplot(1, 2, 2)

lineh = plot(S.dist_from_hyperplane_xval(S.Y == 1), S.class_probability_xval(S.Y == 1), 'v', 'MarkerFaceColor', [.2 .8 .2]);
lineh = [lineh plot(S.dist_from_hyperplane_xval(S.Y == -1), S.class_probability_xval(S.Y == -1), '^', 'MarkerFaceColor', [.2 .2 1])];

xlabel('Cross-validated SVM scores (colors are folds)');
ylabel('Cross-validated class prob estimates');

hh = plot_vertical_line(0); set(hh, 'LineStyle', '--');
hh = plot_horizontal_line(.5); set(hh, 'LineStyle', '--');

legend(lineh, {'True Class = 1' 'True Class = -1'});

drawnow, snapnow

end




% -------------------------------------------------------------------------
% Inputs
% -------------------------------------------------------------------------



function ARGS = parse_inputs(varargin)

p = inputParser;


% Validation functions - customized for each type of input
% ----------------------------------------------------------------------

% valfcn_scalar = @(x) validateattributes(x, {'numeric'}, {'nonempty', 'scalar'});

valfcn_number = @(x) validateattributes(x, {'numeric'}, {'nonempty'}); % scalar or vector

valfcn_cell = @(x) validateattributes(x, {'cell'}, {'nonempty'}); % scalar or vector

% Validation: Region object, structure, or [x1 x2 x3] triplet
% valfcn_custom = @(x) isstruct(x) || isa(x, 'region') || (~isempty(x) && all(size(x) - [1 3] == 0) && all(isnumeric(x)));

valfcn_logical = @(x) validateattributes(x, {}, {'nonempty', 'scalar', '>=', 0, '<=', 1}); % could enter numeric 0,1 or logical

valfcn_effectscode = @(x) validateattributes(x, {'numeric'}, {'nonempty', '<=', 1, '>=', -1});

% Required inputs
% ----------------------------------------------------------------------
p.addRequired('X', valfcn_number);
p.addRequired('Y', valfcn_effectscode);
p.addRequired('id', valfcn_number);

% Optional inputs
% ----------------------------------------------------------------------
% Pattern: keyword, value, validation function handle

p.addParameter('doplot', true, valfcn_logical);
p.addParameter('doverbose', true, valfcn_logical);
p.addParameter('dooptimize', true, valfcn_logical);
p.addParameter('dorepeats', 10, valfcn_number);
p.addParameter('dobootstrap', true, valfcn_logical);
p.addParameter('nboot', 1000, valfcn_number);

p.addParameter('modeloptions', {'KernelFunction', 'linear'}, valfcn_cell);


% Parse inputs and distribute out to variable names in workspace
% ----------------------------------------------------------------------
% e.g., p.parse([30 1 0], [-40 0 10], 'bendpercent', .1);
p.parse(varargin{:});

ARGS = p.Results;

end % parse_inputs);
