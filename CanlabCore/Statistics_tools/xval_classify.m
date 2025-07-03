function S = xval_classify(X, labels, varargin)
% xval_classify  Cross-validated discriminant classification using fitcdiscr
%
% USAGE:
%   S = xval_classify(X, labels, varargin)
%
% Required inputs
% -------------------------------------------------------------------------
%   X                   [N×M] numeric matrix of observations × features
%   labels              [N×1] integer vector of class labels
%
% Optional inputs (Name-Value pairs)
% -------------------------------------------------------------------------
%   'id'                        [N×1] integer grouping variable to leave whole subject out (default = [])
%   'nFolds'                    positive integer, number of folds (default = 5)
%   'optimizeHyperparameters'   logical, perform nested CV hyperparameter optimization (default = false)
%   'hyperparameterOptions'     struct, options for fitcdiscr hyperparameter optimization (default = [] => use default auto settings)
%   'verbose'                   logical, display progress (default = true)
%   'doplot'                    logical, plot confusion matrix per fold (default = true)
%
% Outputs
% -------------------------------------------------------------------------
%   S           struct with fields:
%       Y                   [N×1] true labels
%       id                  [N×1] grouping variable (empty if not provided)
%       nfolds              positive integer, number of folds
%       trIdx               {1×nfolds} cell array of logical vectors for training set per fold
%       teIdx               {1×nfolds} cell array of logical vectors for test set per fold
%       models              {1×nfolds} cell array of fitcdiscr model objects
%       predictions         {1×nfolds} cell array of predicted labels per fold
%       trueLabels          {1×nfolds} cell array of true labels per fold
%       accuracy            [nfolds×1] classification accuracy per fold (percentage)
%
% DESCRIPTION:
%   Performs stratified k-fold cross-validated linear discriminant analysis (LDA)
%   using MATLAB's fitcdiscr. If 'id' is provided, uses
%   xval_stratified_holdout_leave_whole_subject_out to ensure all observations
%   sharing the same id are in the same test set. Otherwise, uses
%   stratified_holdout_set to assign balanced folds. Optionally performs nested
%   hyperparameter optimization within each training set if 'optimizeHyperparameters'
%   is true. Default optimization settings use:
%       fitcdiscr(..., 'OptimizeHyperparameters','auto', ...
%           'HyperparameterOptimizationOptions', struct('AcquisitionFunctionName','expected-improvement-plus'))
%
% The primary hyperparameters available for optimization in fitcdiscr are Delta, Gamma, and DiscrimType.  
% Delta controls a threshold below which linear coefficients are set to zero, effectively performing feature elimination  ￼ ￼
% Gamma determines the amount of regularization applied when estimating the covariance matrix, ranging from no regularization to full diagonalization
% DiscrimType specifies the covariance‐structure variant—linear, quadratic, diagonal, or pseudo­inverse—thus affecting both bias and computational behavior
%
% SEE ALSO:
%   fitcdiscr, xval_stratified_holdout_leave_whole_subject_out, stratified_holdout_set, confusionchart

% Author and copyright information:
%   By CANlab members and collaborators
%   GNU General Public License; see http://www.gnu.org/licenses/



    %% -------------------- Input parsing --------------------
    p = inputParser;
    p.FunctionName = mfilename;
    % Required
    addRequired(p, 'X', @(x) validateattributes(x, {'numeric'}, {'2d','nonsparse','real','nonempty'}));
    addRequired(p, 'labels', @(v) ...
        validateattributes(v, {'numeric'}, {'vector','integer','nrows',size(X,1)}));
    % Optional Name-Value parameters
    addParameter(p, 'id', [], @(v) ...
        validateattributes(v, {'numeric'}, {'vector','integer','nrows',size(X,1)}));
    addParameter(p, 'nFolds', 5, @(v) validateattributes(v, {'numeric'}, {'scalar','integer','positive'}));
    addParameter(p, 'optimizeHyperparameters', false, @(v) ...
        validateattributes(v, {'logical'}, {'scalar'}));
    addParameter(p, 'hyperparameterOptions', struct(), @(s) validateattributes(s, {'struct'}, {}));
    addParameter(p, 'verbose', true, @(v) validateattributes(v, {'logical'}, {'scalar'}));
    addParameter(p, 'doplot', true, @(v) validateattributes(v, {'logical'}, {'scalar'}));
    parse(p, X, labels, varargin{:});

    X = p.Results.X;
    labels = p.Results.labels(:);
    id = p.Results.id;
    nFolds = p.Results.nFolds;
    doOptimize = p.Results.optimizeHyperparameters;
    hpOptions = p.Results.hyperparameterOptions;
    verbose = p.Results.verbose;
    doplot = p.Results.doplot;

    N = size(X, 1);

    %% -------------------- Determine fold indices --------------------
    if ~isempty(id)
        % Use leave-whole-subject-out stratified folds
        [trIdxCell, teIdxCell] = xval_stratified_holdout_leave_whole_subject_out( ...
            labels, id, 'nfolds', nFolds, 'doverbose', verbose, 'doplot', doplot);
    else
        % Use stratified_holdout_set to get fold assignments
        foldIndices = stratified_holdout_set(labels, nFolds);
        trIdxCell = cell(1, nFolds);
        teIdxCell = cell(1, nFolds);
        for f = 1:nFolds
            teIdxCell{f} = foldIndices.teIdx{f};
            trIdxCell{f} = foldIndices.trIdx{f};
        end
    end

    %% -------------------- Initialize output structure --------------------
    S = struct();
    S.Y = labels;
    S.id = id;
    S.nfolds = nFolds;
    S.trIdx = trIdxCell;
    S.teIdx = teIdxCell;
    S.models = cell(1, nFolds);
    S.predictions = cell(1, nFolds);
    S.trueLabels = cell(1, nFolds);
    S.accuracy = zeros(nFolds, 1);

    %% -------------------- Cross-validation loop --------------------
    for f = 1:nFolds
        trainMask = logical(S.trIdx{f});
        testMask = logical(S.teIdx{f});
        Xtrain = X(trainMask, :);
        ytrain = labels(trainMask);
        Xtest  = X(testMask, :);
        ytest  = labels(testMask);

        if verbose
            fprintf('Fold %d/%d: Training on %d samples, Testing on %d samples...\n', ...
                f, nFolds, sum(trainMask), sum(testMask));
        end

        %% ---- Train fitcdiscr model ----
        if doOptimize
            % Determine optimization settings
            if isempty(fieldnames(hpOptions))
                % Default hyperparameter optimization
                optCell = { ...
                    'OptimizeHyperparameters', 'auto', ...
                    'HyperparameterOptimizationOptions', struct( ...
                        'AcquisitionFunctionName', 'expected-improvement-plus' ) ...
                };
            else
                % Convert user-provided struct to name-value cell
                optFields = fieldnames(hpOptions);
                optCell = cell(1, 2*numel(optFields));
                for k = 1:numel(optFields)
                    optCell{2*k-1} = optFields{k};
                    optCell{2*k}   = hpOptions.(optFields{k});
                end
            end
            model = fitcdiscr(Xtrain, ytrain, optCell{:});
        else
            model = fitcdiscr(Xtrain, ytrain);
        end
        S.models{f} = model;

        %% ---- Test on held-out fold ----
        ypred = predict(model, Xtest);
        S.predictions{f} = ypred;
        S.trueLabels{f} = ytest;
        S.accuracy(f) = 100 * mean(ypred == ytest);

        if verbose
            fprintf(' Fold %d Accuracy: %.2f%%\n', f, S.accuracy(f));
        end

        %% ---- Optional plotting ----
        % if doplot
        %     figure('Name', sprintf('Fold %d Confusion Matrix', f), 'NumberTitle', 'off');
        %     cm = confusionchart(ytest, ypred);
        %     cm.Title = sprintf('Fold %d Confusion Matrix (Acc = %.2f%%)', f, S.accuracy(f));
        %     cm.RowSummary = 'row-normalized';
        %     cm.ColumnSummary = 'column-normalized';
        % end 


    end  % fold



    %% -------------------- Post‐processing after all folds --------------------
% Aggregate predictions across folds into a single vector yfit
yfit = NaN(N,1);
for f = 1:S.nfolds
    teMask = logical(S.teIdx{f});
    yfit(teMask) = S.predictions{f};
end
S.yfit = yfit;

% Save true labels (already in S.Y) under a more concise name
S.y = S.Y;

% Compute overall accuracy (percentage of correctly classified observations)
S.overallAccuracy = 100 * mean(S.yfit == S.y);

% (Optional) Save any other relevant variables—these are already in S:
%   S.trIdx, S.teIdx, S.predictions, S.trueLabels, S.accuracy, S.models, etc.

if doplot

    % Aggregate and plot the overall confusion matrix
    figure('Name', 'Overall Confusion Matrix', 'NumberTitle', 'off');
    han = confusionchart(S.y, S.yfit);
    han.Title = sprintf('Overall Confusion Matrix (Acc = %.2f%%)', S.overallAccuracy);
    han.RowSummary = 'row-normalized';
    han.ColumnSummary = 'column-normalized';

    % Set font size and custom colors for diagonal vs. off-diagonal
    set(gca, 'FontSize', 16);
    set(han, 'DiagonalColor', [1 .7 .3], 'OffDiagonalColor', [.3 .3 .7]);

end


    %% -------------------- End of function --------------------
end