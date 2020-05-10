function S = xval_SVM(varargin)
%  SVM for repeated-measures (within-person) classification, with repeated cross-val and nested cross-val options
%
% :Usage:
% ::
%
% S = xval_SVM(X, Y, id, varargin)
%
% Steps:
% - Select holdout sets, keeping images from the same id together and stratifying on the outcome to be predicted
% - Fit overall model and get predictor weights (betas, b)
% - Cross-validation with standard a priori hyperparameter choices
% - Plot scores and ROC curves
% - Nested cross-val with hyperparameter optimization [optional]
% - Repeat cross-validation of optimized model with different random splits [optional]
%
% Notes:
% - Uses Matlab's Stats/ML Toolbox SVM object, fitcsvm, and hyperparameter
% optimization. Other options, e.g., fitclinear, may be better for some
% situations (large num. variables)
%
% - Single-interval accuracy from cross-val and ROC_plot may differ because
% ROC_plot chooses a new score threshold that maximizes overall balanced
% accuracy. Forced-choice accuracy should be identical.
%
% - If optimizing hyperparameters AND repeating cross-validation, accuracy
% estimates will use optimized params, so are somewhat optimistically biased.
% if optimizing, apply final model to an independent dataset to test
% accuracy.
%
% - Not sure how to turn off output when optimizing hyperparameters
%
% - Optimizing hyperparameters works on linear classifiers only, but can explore nonlinear models
% with a small change in the code (now commented out). This could be added
% as an optional flag as well.
%
% - Optimizing hyperparameters maximizes accuracy, so may need large
% samples to be effective. 
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
%
%
% :Examples:
% ::
%
%    % give examples of code here
%    param1 = abc();
%    param2 = xyz();
%    [out1,out2] = func_call(param1, param2)
%
% :References:
%   See Mathworks functions
%
% :See also:
%   - fmri_data.predict, canlab_run_paired_SVM, other xval_ functions
%

% ..
%    Programmers' notes:
%    List dates and changes here, and author of changes
% ..

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
    printhdr('xval_SVM: Cross-validated Support Vector Machine classification');
    printhdr(' ')    
    
    disp(' ')
    disp(' ')
    printhdr('Selecting and analyzing cross-validation folds');
    
end

[S.trIdx, S.teIdx] = xval_stratified_holdout_leave_whole_subject_out(S.Y, S.id, 'doverbose', doverbose, 'doplot', doplot);
drawnow, snapnow

%% Fit the overall a priori model: SVM with linear kernel
% -------------------------------------------------------------------------
if doverbose
    
    printhdr('Model training and cross-validation without optimization')
    fprintf('Training overall model. ')

end

SVMModel = fitcsvm(X, S.Y, S.modeloptions{:});
S.b = SVMModel.Beta;                                       % Model weights/betas

% [SVMModel, S.ScoreParameters] = fitPosterior(SVMModel);    % Non-cross val params for score->class probability estimate

if doverbose
    
    % non-cross-val accuracy: sanity check. should usually be 100% if training is overfit and n >> p
    label = predict(SVMModel, X);

    acc = S.accfun(S.Y, label);
    fprintf('Accuracy without cross-val: %3.0f%%\n', acc);

end

%% Cross-validation with standard a priori hyperparameter choices
% -------------------------------------------------------------------------
S.nfolds = length(S.trIdx);

S = crossval(S, X, doverbose);

% Summarize accuracy
% -------------------------------------------------------------------------

S.crossval_accuracy = S.accfun(S.Y, S.yfit);

varname = 'dist_from_hyperplane_xval';
[scores_within_id, scorediff, d] = get_scores_within_id(S, varname);

S.Y_within_id = get_scores_within_id(S, 'Y');

S.scores_within_id = scores_within_id;
S.scorediff = scorediff;
S.crossval_accuracy_within = 100 * sum(scorediff > 0) ./ sum(~isnan(scorediff));
S.classification_d = d;

if doverbose
    
    fprintf('Accuracy (cross-validated): Single-interval: %3.0f%%, Forced-choice: %3.0f%%\nWithin-person d = %3.2f\n', S.crossval_accuracy, S.crossval_accuracy_within, S.classification_d)

end

%% Plot scores and ROC curves
% -------------------------------------------------------------------------
if doplot
    
    plot_scores_and_ROC(S);
    
end

%% Nested cross-val with hyperparameter optimization
% -------------------------------------------------------------------------

if dooptimize
    
    S = crossval_nested(S, X, doverbose);
    
    S.crossval_accuracy_opt_hyperparams(i) = S.accfun(S.Y, S.yfit);
    
    % Get modal consensus params (for future validation)
    hbyfold = cat(1, S.fold_modeloptions{:});
    best_modeloptions = {};
    
    for i = 1:size(hbyfold, 2)
        
        mydat = cat(1, hbyfold{:, i});
        
        if ischar(mydat) || iscategorical(mydat) || islogical(mydat)
            
            best_modeloptions{1, i} = mode(mydat);
            
        elseif isnumeric(mydat)
            
            best_modeloptions{1, i} = trimmean(mydat, 80);
            
        else
            error('Unknown model option class! Extend this code.')
            
        end
        
    end
    
    S.modeloptions = best_modeloptions;
    
    if doverbose
        
        disp('Optimized hyperparameters for each fold')
        disp(S.hyperparams_by_fold)
        
        disp(' ')
        disp('Best options, stored in S.modeloptions:');
        
        disp(best_modeloptions)
        
    end
    
end % do optimize

%% Repeat cross-validation of optimized model, if specified
% - if hyperparams optimized, then use "consensus params" for repeats
% - these are stored in S.modeloptions now, so will be used
% - If optimizing, re-do the first accuracy with the optimized
% hyperparameters
% -------------------------------------------------------------------------
if dorepeats > 1
    
    if doverbose, fprintf('Repeating cross-val %d times. ', dorepeats), end
    
    if dooptimize, startat = 1; else, startat = 2; end
    
    for i = startat:dorepeats     % do remainder of repeats without verbose.
        
        if doverbose, fprintf('%d ', i), end
        
        Sr = S;
        
        % Get new cvpartition
        [Sr.trIdx, Sr.teIdx] = xval_stratified_holdout_leave_whole_subject_out(Sr.Y, Sr.id, 'doverbose', false, 'doplots', false);
        
        % Cross-validate. If S.modeloptions has been replaced (in optimization above)
        Sr = crossval(Sr, X, false);
        
        % Save accuracy
        S.crossval_accuracy(i) = S.accfun(Sr.Y, Sr.yfit);
        
    end
    
    if doverbose
        fprintf('Done!\n')
        
        fprintf('CV accuracy across %d reps, mean = %3.0f%%, std = %3.0f%%, min = %3.0f, max = %3.0f\n', dorepeats, ...
            mean(S.crossval_accuracy), std(S.crossval_accuracy), min(S.crossval_accuracy), max(S.crossval_accuracy));
    end
    
end


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
    
    [~, pscore] = predict(SVM_fold, X(S.teIdx{i}, :));
    pscore = pscore(:, 2);                                  % Platt scaling scores (class probability)
    
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
        struct('AcquisitionFunctionName','expected-improvement-plus'));
    
    % Not sure how to turn off output when optimizing hyperparameters
    drawnow, snapnow; close, close
    
    
    % Add hyperparameter results to table
    % -------------------------------------
    opt_hyperparam_table = Mdl.HyperparameterOptimizationResults.XAtMinObjective;
    
    if i == 1
        S.hyperparams_by_fold = opt_hyperparam_table; % best parameters - table
        
    else
        S.hyperparams_by_fold(i, :) = opt_hyperparam_table;
    end
    
    % Re-train with final choices
    % -------------------------------------
    fold_modeloptions = S.modeloptions;
    
    % Add optimal choices to the option set
    for j = 1:size(opt_hyperparam_table, 2)
        vname = opt_hyperparam_table.Properties.VariableNames{j};
        
        fold_modeloptions{end + 1} = vname;
        paramval = opt_hyperparam_table.(vname)(1);
        
        % fix: categorical to text
        if iscategorical(paramval)
            paramval = char(string(paramval));
            paramval = strcmp(paramval, 'true');
            
        end
        
        fold_modeloptions{end + 1} = paramval; %#ok<*AGROW>
        
    end
    
    S.fold_modeloptions{i} = fold_modeloptions;
    
    % Fit to training data for this fold
    SVM_fold = fitcsvm(X(S.trIdx{i}, :), S.Y(S.trIdx{i}), fold_modeloptions{:});
    
    % CVSVMModel = fitcsvm(X, Y, 'CVpartition', c, ...
    %     'Learner', Mdl.ModelParameters.Learner, ...
    %     'lambda', Mdl.ModelParameters.Lambda, ...
    %     'Regularization', Mdl.ModelParameters.Regularization, ...
    %     'Prior', 'uniform');
    
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





function plot_scores_and_ROC(S)

create_figure('cross-val accuracy', 2, 2);
subplot(2, 2, 1); delete(gca); subplot(2, 2, 2); delete(gca);
axes('Position', [.13 .58 .77 .35]);
set(gca, 'FontSize', 16)
hold on

n = size(S.scores_within_id, 1);
x = [1:n; 1:n];
plot(x(:, S.scorediff > 0), S.scores_within_id(S.scorediff > 0, :)', 'Color', [.3 .3 .3], 'LineWidth', 1);
plot(x(:, S.scorediff < 0), S.scores_within_id(S.scorediff < 0, :)', 'r', 'LineWidth', 1);

plot(x(1, :), S.scores_within_id(:, 1)', '^', 'Color', [.3 .5 1] / 2, 'MarkerFaceColor', [.3 .5 1],  'LineWidth', 1);
plot(x(2, :), S.scores_within_id(:, 2)', 'v', 'Color', [1 .5 0] / 2, 'MarkerFaceColor', [1 .5 0],  'LineWidth', 1);

xlabel(sprintf('ID, %3.0f%% single-interval acc, %3.0f%% forced-choice', S.crossval_accuracy, S.crossval_accuracy_within));
ylabel('SVM Score');

title('Cross-validated SVM scores (no hyperparameter optimization');
disp('Black lines: Correct, Red lines: Errors. Red triangles: Scores for Class 1, Blue triangles: Scores for Class -1');

% ROC plot
subplot(2, 2, 3);
S.ROC_single_interval = roc_plot(S.dist_from_hyperplane_xval, logical(S.Y > 0), 'color', [.4 .4 .7]);
title('Single-interval ROC')
set(gca, 'FontSize', 16)

subplot(2, 2, 4);
% Paired forced-choice. Get complete cases - Remove NaNs id-wise
outcomes = S.Y_within_id;
scores = S.scores_within_id;
[~, outcomes, scores] = nanremove(outcomes, scores);

S.ROC_forced_choice = roc_plot(scores(:), logical(outcomes(:) > 0), 'color', [.4 .4 .7], 'twochoice');
title('Forced-choice ROC')
set(gca, 'FontSize', 16)

drawnow, snapnow

end




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
p.addParameter('dorepeats', true, valfcn_number);

p.addParameter('modeloptions', {'KernelFunction', 'linear'}, valfcn_cell);


% Parse inputs and distribute out to variable names in workspace
% ----------------------------------------------------------------------
% e.g., p.parse([30 1 0], [-40 0 10], 'bendpercent', .1);
p.parse(varargin{:});

ARGS = p.Results;

end % parse_inputs);
