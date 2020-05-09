
%%  SVM with repeated cross-val and nested cross-val options
%
% Select holdout sets for outer loop
% - Keep images from the same id together
% - Stratify on the outcome to be predicted
% -------------------------------------------------------------------------

doverbose = true;
dorepeats = 10;         % false, or number of repeats
optimizehyperparameters = true;

% Required Inputs: X, Y, id

S = struct();   % Define structure for models and output; use variable names in fmri_data.predict when possible

S.Y = Y;
S.id = id;      % Subject grouping - not in fmri_data.predict

S.modeloptions = {'KernelFunction', 'linear'};

S.accfun = @(Y, yfit) 100 .* nansum(Y == yfit) ./ sum(~isnan(Y));

[S.trIdx, S.teIdx, S.test_set_condf] = xval_stratified_holdout_leave_whole_subject_out(S.Y, S.id, 'doverbose', doverbose);

% Fit the overall a priori model: SVM with linear kernel
% -------------------------------------------------------------------------
if doverbose, fprintf('Training overall model. '), end

SVMModel = fitcsvm(X, S.Y, S.modeloptions{:});
S.b = SVMModel.Beta;                                       % Model weights/betas                           

[SVMModel, S.ScoreParameters] = fitPosterior(SVMModel);    % Non-cross val params for score->class probability estimate

% non-cross-val accuracy: sanity check. should usually be 100% if training is overfit and n >> p
[label, score] = predict(SVMModel, X);

acc = S.accfun(S.Y, label);
if doverbose, fprintf('Accuracy without cross-val: %3.0f%%\n', acc), end


% Cross-validation (outer loop)
% -------------------------------------------------------------------------
S.nfolds = length(S.trIdx);

S = crossval(S, X, doverbose);

% Summarize
S.crossval_accuracy = S.accfun(S.Y, S.yfit);
if doverbose, fprintf('Accuracy (cross-validated): %3.0f%%\n', S.crossval_accuracy), end
 
% Repeat, if specified
% -------------------------------------------------------------------------
if dorepeats > 1
    
    if doverbose, fprintf('Repeating cross-val %d times. ', dorepeats), end
    
    for i = 2:dorepeats     % do remainder of repeats without verbose.
        
        if doverbose, fprintf('%d ', i), end
        
        Sr = S;
        
        % Get new cvpartition
        [Sr.trIdx, Sr.teIdx] = xval_stratified_holdout_leave_whole_subject_out(Sr.Y, Sr.id, 'doverbose', false, 'doplots', false);

        % Cross-validate
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


%% Nested cross-val with hyperparameter optimization
% -------------------------------------------------------------------------

S = crossval_nested(S, X, doverbose);


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
    
    [SVM_fold,ScoreParameters] = fitPosterior(SVM_fold);    % Works for fitcsvm, not fitclinear
    
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


    % Add hyperparameter results to table
    % -------------------------------------
    opt_hyperparam_table = Mdl.HyperparameterOptimizationResults.XAtMinObjective;
    
    if i == 1
        S.hyperparams_by_fold = opt_hyperparam_table; % best parameters - table
        
    else
        S.hyperparams_by_fold(2, :) = opt_hyperparam_table;
    end
    
    % Re-train with final choices
    % -------------------------------------
    fold_modeloptions = S.modeloptions;
    
    % Add optimal choices to the option set
    for j = 1:size(opt_hyperparam_table, 2)
        vname = opt_hyperparam_table.Properties.VariableNames{j};
        
    	fold_modeloptions{end + 1} = vname;
        fold_modeloptions{end + 1} = opt_hyperparam_table.(vname)(1);

    end
    
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
    
    [SVM_fold,ScoreParameters] = fitPosterior(SVM_fold);    % Works for fitcsvm, not fitclinear
    
    [~, pscore] = predict(SVM_fold, X(S.teIdx{i}, :));
    pscore = pscore(:, 2);                                  % Platt scaling scores (class probability)
    
    S.class_probability_xval(S.teIdx{i}, 1) = pscore;
    
end

if doverbose, fprintf('Done!\n'), end

end % crossval
