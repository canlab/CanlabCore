function [cverr, stats, optout] = predict(obj, varargin)
% Predict outcome (Y) from brain data and test cross-validated error rate for an fmri_data object
%
% :Usage:
% ::
%
%     [cverr, stats, optional_outputs] = predict(obj, varargin)
%
% :Features:
%   - flexible specification of algorithm by function name
%   - k-fold cross-validation, default = 5-fold, can enter custom fold membership
%   - folds are stratified on outcome
%   - choice of multiple error metrics (class loss, mse, etc.)
%   - by default, chooses error metric based on outcome type (classes vs. continuous-valued)
%   - returns all outputs for each fold returned by the algorithm in optout cell array variable
%   - bootstrapping of weights built in [optional keyword]
%   - select variable number of components (for pcr-based techniques)
%
% :Inputs:
% obj is mandatory, rest are optional
%   **obj:**
%        fmri_data or image_vector object, with fields .dat (data used to predict) and .Y (outcome)
%
% :Optional inputs: (with their default values)
%
%   **nfolds** = 5
%        number of folds
%
%   **nfolds** = [vector of integers]
%        can also input vector of integers for holdout set IDs
%
%   **error_type** = mcr
%        mcr, mse: misclassification rate or mean sq. error
%
%   **algorithm_name** = 'cv_pcr'
%        name of m-file defining training/test function
%
%   **useparallel** = 1
%        Use parallel processing, if available; follow by 1 for yes, 0 for no
%
%   **bootweights** = 0
%        bootstrap voxel weights; enter bootweights do bootstrapping of weight maps (based on all observations)
%
%   **savebootweights**
%        save bootstraped weights (useful for combining across multiple iterations of predict())
%
%   **bootsamples** = 100
%        number of bootstrap samples to use
%
%   **numcomponents** = xxx:
%        save first xxx components (for pca-based methods)
%
%   **nopcr**
%        for cv_lassopcr and cv_lassopcrmatlab: do not do pcr, use original variables
%
%   **lasso_num** = xxx
%        followed by number of components/vars to retain after shrinkage
%
%   **hvblock** = [h,v]
%        use hvblock cross-validation with a block size of 'h' (0 reduces to v-fold xval) and
%        number of test observations 'v' (0 reduces to h-block xval)
%
%   **rolling** = [h,v,g]
%        use rolling cross-validation with a block size of 'h' (0 reduces to v-fold xval) and
%        number of test observations 'v' (0 reduces to h-block xval), and a training size
%        of g * 2 surrounding hv
%
%   **verbose** = 1
%        Set to 0 to suppress output to command window
%
%   **platt_scaling**
%        calculate cross-validated platt scaling if using SVM.
%        Softmax parameters [A,B] are in other_output{3}
%
% :Algorithm choices:
%   You can input the name (as a string array) of any algorithm with the
%   appropriate inputs and outputs. i.e., this can either be one of the
%   built-in choices below, or the name of another m-file.
%   The format for algorithm functions is :
%   [yfit, other_outputs] = predfun(xtrain, ytrain, xtest, optional_inputs)
%   Each algorithm can take/interpret its own optional inputs.
%   For bootstrapping of weights, algorithms MUST RETURN 3 OUTPUTS
%   (programming 'feature')
%
%   To choose an algorithm, enter 'algorithm_name' followed by a text string
%   with a built-in algorithm name, or a function handle for a custom algorithm
%   Built-in algorithm choices include:
%
%   **cv_multregress:**
%        multiple regression
%
%   **cv_univregress:**
%        Average predictions from separate univariate regression of outcome on each feature
%
%   **cv_svr:**
%        Support vector regression with Spider package; requires spider
%
%   **cv_pcr:**
%        [default] Cross-validated principal components regression
%
%   **cv_mlpcr:**
%        Cross-validated multilevel principal components regression. See
%        'help mlpcr2' for full documentation. If run with default settings
%        returns the same result as cv_pcr (except when bootstrapping, 
%        see below), except with information pertaining to within and 
%        between predictive variance in optout. optout provides 8 outputs: 
%        total model, between model, within model, intercept (same for all 
%        models), between eigenvectors, between scores, within 
%        eigenvectors and within scores. Requires 'subjID' option followed 
%        by size(obj.dat,2) x 1 vector of block labels.
%        Optional: Concensus PCA, {'cpca', 1}. [Default]={'cpca, 0}.
%        Optional: Dimension selection, {'numcomponents', [bt, wi]}.
%                   [Default] = {'numcomponents',[Inf,Inf]} (df constrained)
%        Note: You probably want to bootstrap this manually if
%           bootstrapping. If bootstrapping using fmri_data/predict's
%           built in method you should note three things. First, You are
%           bootstrapping at the image level, not the block level (true
%           for all algorithms in fmri_data/predict, but carries special
%           implications for block aware algorithms). Second, higher order
%           within block PCA dimensions are HIGHLY unstable. If your
%           model's within block PCA dimension is not low you may have no
%           significant voxels. Third, if your full dataset is balanced you
%           may want to weight your bootstrap PCAs and regressions by using
%           the {'cpca',1} argument pair to compensate for imbalance in
%           bootstrap samples.
%
%   **cv_pls:**
%        Cross-validated partial least squares regression (only univariate
%        outcomes for now)
%
%   **cv_lassopcr:**
%        Cross-val LASSO-PCR; can enter 'lasso_num' followed by components to retain by shrinkage
%        NOTE: can enter 'EstimateParams' to use shrankage
%        lasso method based on the estimated optimal lambda
%        that minimizes the mean squared error (MSE) of nested
%        cross-validation models. Output of nested cv model is
%        saved in stats.other_output_cv{:,3}. Output includes
%        'Lambda' parameter and min MSE value.
%
%   **cv_lassopcrmatlab:**
%        Cross-val LASSO-PCR; can enter 'lasso_num' followed by components to retain by shrinkage
%        NOTE: this uses the matlab implementation of LASSO,
%        but can also run ridge or elastic net. Reduces to PCR
%        when no lasso_num is entered by default.  Use MSE for
%        predicting continuous data and MCR for classifying
%        binary data.
%        NOTE: You can input any optional inputs that lassoglm
%        takes.
%        Enter 'Alpha', (0,1] as optional inputs to
%        run ridge (Alpha approaches 0, but excluding 0), lasso (Alpha = 1), or elastic
%        net (Alpha between 0 and 1)
%        NOTE: Requires Matlab R2012a and higher.
%        NOTE: Optional input: 'EstimateParams' - this will
%        use grid search and nested cross validation to
%        estimate both Lambda and Alpha (independent of your specification 
%        to estimate Alpha or Lambda). Output is saved in
%        stats.other_output_cv{:,3}.  Output includes 'Alpha'
%        parameter which is the elastic net mixture value
%        between l1 and l2 regularization, 'Lambda' parameter,
%        which is amount of LASSO regularization/shrinkage, and
%        'errorMatrix', which is the amount of error for each
%        parameter combination. Use imagesc(obj.stats_other_output_cv{:,3}.errorMatrix)
%        to view matrix.  Min of this matrix is the best
%        fitting parameters.
%        NOTE: To estimate the optimal LASSO Lambda, use the optional 
%        argument 'cv' followed by the number of crosvalidation folds and 
%        the argument 'nfolds',1 (e.g. 'cv',5,'nfolds',1). This will utilize the internal 
%        cross-validataion of the lassoglm function to estimate Lambda.    
%
%   **cv_svm:**
%        Cross-val support vector machine using Spider package
%        NOTE: This is sensitive to scale of outputs! Use -1 , 1
%        NOTE: Optional inputs: Slack var parameter: 'C', 1 [default], 'C', 3 etc.
%        Distance from hyperplane saved in
%        stats.other_output_cv{:,2}.  Recommend using the reordered
%        cross-validated distance from hyperplane saved in stats.other_output{3}
%        stats.dist_from_hyperplane_xval =  cross-validated distance from hyperplane
%        stats.weight_obj = voxel (variable) weight object
%        e.g., orthviews(stats.weight_obj)
%        Intercept for calculating dist from hy is in stats.other_output_cv{:,3}
%        e.g., dist_hy = stats.weight_obj.dat' * obj.dat, where obj is a new set of test images
%        NOTE: To run nonlinear SVM using radial basis
%        function.  Add 'rbf' followed by size of sigma (e.g., 2).
%        NOTE: To estimate some of the parameters using
%        nested cross validation add 'EstimateParams' as optional input.
%        NOTE: To run multiclass SVM (i.e., one vs rest) add
%        'MultiClass' as optional input.  Important - Obj.Y must be a matrix (data x
%        class) with a column of 1 and -1 indicating each
%        class.  For example, if using 3 classes, then obj.Y
%        must have 3 columns.
%        NOTE: To run a balanced SVM where the number of cases for each class are unequal (i.e., one vs rest) add
%        'Balanced' as optional input, followed by a numerical value indicating the ridge amount (e.g., 0.01).
%
%   **cv_multilevel_glm:**
%        Runs glmfit_multilevel. Must pass in ''subjIDs'' followed by an array specifying which subject each trial belongs to
%        Subjects' trials must all be "adjacent", i.e., don't
%        put some of subject 1's trials at the beginning and
%        other trials at the end -- subjIDs does not handle
%        this case correctly. Also, 2ND LEVEL PREDICTORS NOT
%        CURRENTLY SUPPORTED.  code can be expanded to support this.
%        mean-centering X and/or Y will NOT impact the
%        predictor betas.  Note that it WILL impact the intercept
%        esimate as well as how much variance is explained
%        (pred_outcome_r).  Stratified CV partition not
%        supported either, pass in custom holdout set.
%
%
% :Outputs:
%
%   **Y:**
%        Copy of outcome data to be predicted
%
%   **algorithm_name:**
%        Name of algorithm; see options above
%
%   **function_call:**
%        String of the command evaluated to call the prediction function
%
%   **function_handle:**
%        Handle for the command evaluated to call the prediction function
%
%   **yfit:**
%        Predicted outcome data (cross-validated)
%
%   **err:**
%        Residuals/misclassification vector (cross-validated)
%
%   **error_type:**
%        Name of error metric used for cverr
%
%   **cverr:**
%        Cross-validated error
%
%   **nfolds:**
%        Number of folds in stratified cross-validation, or
%        vector of integers for membership in custom holdout set of each fold
%            - if k = 1, will estimate weights for full data object
%              and not crossvalidate (useful for bootstrapping)
%
%   **cvpartition:**
%        Cross-val partition object or structure with fold info
%
%   **teIdx:**
%        Cell array of logical vectors with test samples in each fold
%
%   **trIdx:**
%        Cell array of logical vectors with training samples in each fold
%
%   **other_output:**
%        Other outputs returned by the algorithm; number and nature depend on algo choice; e.g., beta weights, svr weights, etc.
%        For many algorithms, other_output{1} is a vector of
%        weights on variables (e.g., voxels)
%
%   **other_output_descrip:**
%        String description of other outputs
%
%   **other_output_cv:**
%        Other outputs for each cross-validation fold
%
%   **other_output_cv_descrip:**
%        Other output from algorithm - for each CV fold
%
%   **mse:**
%        For regression only; mean squared error
%
%   **rmse:**
%        For regression only; root mean squared error
%
%   **meanabserr:**
%        For regression only; mean absolute error
%
%   **pred_outcome_r:**
%        For regression only; prediction-outcome correlation
%
%   **WTS:**
%        bootstrapped weights on voxels
%
%   **weight_obj:**
%        for some algorithms, an fmri_data object with the predictive weights (from full sample)
%
% Here are the fields in the main output structure, with explanations:
%   struct with fields:
% 
%                             Y: Actual (obs) outcome - should be 1, -1 for SVM analysis
%                algorithm_name: Algorithm used (string, name)
%                 function_call: Actual function call string used (for tracking what was run)
%               function_handle: Actual function handle used (for tracking what was run)
%                          yfit: Predicted outcome - should be 1, -1 for SVM analysis
%                                Cross-validated, so can be used as series
%                                of predicted values based on brain
%                                measures (very useful!)
%                           err: Errors. For SVM, 0 = correct, 1 = error
%                    error_type: Type of error metric.  'mcr' = misclassification rate, for SVM
%                         cverr: Cross-validated mean error 
%                        nfolds: Holdout set type 'nfolds'
%                   cvpartition: Object with information about training/test sets
%                         teIdx: Testing IDs for each fold (holdout set)
%                         trIdx: Training IDs for each fold (holdout set)
%                  other_output: {[328798?1 double]  [39?1 double]  [4.7134]  [1?1 struct]}
%                                {1} : Weight map trained on all data
%                                {2} : ?
%                                {3} :Intercept trained on all data
%                                {4} : Training model specifications              
%          other_output_descrip: 'Other output from algorithm - trained on all data (these depend on algorithm)'
%               other_output_cv: Same as other output above, but for each holdout set separately
%                                Useful if you want to re-apply models from
%                                training subsets, re-creating
%                                cross-validated output
%       other_output_cv_descrip: 'Other output from algorithm - for each CV fold'
%                           phi: ??
%     dist_from_hyperplane_xval: Cross-validated distance perpendicular to class boundary
%                                   higher =
%                                   stronger prediction in favor of Class 1, lower = in favor of class 2.
%                                   Useful! use as continuous measure of
%                                   observation scores. Can calculate
%                                   effect sizes from this, for example. 
%                    weight_obj: statistic_image object with weight map.
%                    Use object methods to plot, etc. If we want to apply
%                    the weight map to new data, this is all we need. 
%                           WTS: More info, plus bootstrap statistics if
%                           bootstrapping was run.
%
% :Examples:
% ::
%
%    obj = fmri_data;
%    obj.dat = randn(30, 50); %    30 voxels, 50 images (observations)
%    obj.Y = obj.dat' * rand(30, 1) + randn(50, 1); %    toy Y, linear combo of X plus noise
%    [cverr, stats, regression_outputs] = predict(obj);
%
%    Simulated example with 100 observations, 1000 voxels, with bootstrapping
%    dat = fmri_data;
%    dat.Y = rand(100, 1);
%    dat.dat = repmat(dat.Y', 1000, 1) + 10*rand(1000, 100);
%    [err,stats] = predict(dat, 'bootweights', 'algorithm_name', 'cv_lassopcr');
%
%    [cverr, stats, regression_outputs] = predict(obj, 'nfolds', 3, 'error_type', 'meanabserr');
%    [cverr, stats, regression_outputs] = predict(obj, 'algorithm_name', 'cv_univregress', 'error_type', 'meanabserr');
%    [cverr, stats, optout] = predict(obj, 'algorithm_name', 'cv_lassopcr', 'lasso_num', 5, 'nfolds', 5, 'error_type', 'mse', 'bootweights');
%    [cverr, stats, optout] = predict(dat, 'algorithm_name', 'cv_svm', 'nfolds', 5, 'error_type', 'mse');
%    [cverr, stats, optout] = predict(dat, 'algorithm_name', 'cv_svm', 'rbf', 2, 'nfolds', 5, 'error_type', 'mse'); %SVM w/ radial basis function
%    [cverr, stats, optout] = predict(dat, 'algorithm_name', 'cv_svm', 'rbf', 2, 'EstimateParams', 'nfolds', 5, 'error_type', 'mse'); %SVM w/ radial basis function w/ parameters estimated using nested cross-valdiation
%    [cverr, stats, optout] = predict(dat, 'algorithm_name', 'cv_svm', 'nfolds', 5, 'MultiClass', 'error_type', 'mse');
%
%    Elastic net with first 10 components:
%    [cverr, stats, optout] = predict(dat_masked, 'algorithm_name', 'cv_lassopcrmatlab', 'nfolds', 5, 'error_type', 'mse', 'numcomponents', 10, 'Alpha', .5); stats.pred_outcome_r
%
%    Ridge with first 10 components:
%    [cverr, stats, optout] = predict(dat_masked, 'algorithm_name', 'cv_lassopcrmatlab', 'nfolds', 5, 'error_type', 'mse', 'numcomponents', 10, 'Alpha', 0.00001); stats.pred_outcome_r
%
%    Lasso with all components, but shrink to retain 2 components only:
%    [cverr, stats, optout] = predict(dat, 'algorithm_name', 'cv_lassopcrmatlab', 'nfolds', whfolds, 'nopcr', 'lasso_num', 2, 'Alpha', 1);
%    [cverr, stats, optout] = predict(dat, 'algorithm_name', 'cv_lassopcr', 'nfolds', whfolds, 'lasso_num', 2);
%
%    Lasso with the shrinkage methods based on the estimated optimal lambda that minimizes MSE of nested cross-validation models.
%    [cverr, stats, optout] = predict(dat, 'algorithm_name', 'cv_lassopcr', 'nfolds', whfolds, 'estimateparam');
%    [cverr, stats, optout] = predict(dat, 'algorithm_name', 'cv_lassopcr', 'nfolds', 5, 'estimateparam');
%
%    Lasso without doing PCR:
%    [cverr, stats, optout] = predict(dat, 'algorithm_name', 'cv_lassopcrmatlab', 'nfolds', whfolds, 'nopcr', 'lasso_num', 2, 'Alpha', 1);
%    [cverr, stats, optout] = predict(dat, 'algorithm_name', 'cv_lassopcr', 'nfolds', whfolds, 'lasso_num', 2, 'nopcr');
%    [cverr, stats, optout] = predict(dat, 'algorithm_name', 'cv_lassopcr', 'nfolds', 5, 'estimateparam', 'nopcr');
%
%    Lasso pcr using hvblock cross-validation on time-series, h = 3, v = 5;
%    [cverr, stats, optout] = predict(dat, 'algorithm_name', 'cv_lassopcr', 'hvblock',[3,5]);
%
%    Output display:
%    orthviews(stats.weight_obj)
%    line_plot_multisubject(stats.yfit, stats.Y, 'subjid', id_numbers);
%
% :See also:
%    predict_test_suite method for fmri_data, which runs predict with multiple
%    options and summarizes output.
%    xval_regression_multisubject, xval_lasso_brain
%
% ..
%    Original version: Copyright Tor Wager, Dec 2010


%    Programmers' notes:
%    Dec 2011: Changed to input double format to algorithms; some have
%    problems with single format
%
%    March 2012: Changed default method to cv_pcr, and updated input option
%    parsing to avoid silly mistakes in specifying algorithm.
%    Changed line 377 to if rank(X) < size(sc, 1). Tested: pinv and normal inv
%    give same results up to machine precision.
%
%    10/16/12: Luke Chang: Changed input to bootstrp to double because of problems with
%    single format (see Tor's edits on 11/10/12)
%
%    11/2/2012: tor wager: cv_svm output added to return distance from hyperplane in stats.other_output_cv{:,2}
%
%    11/7/12: Luke Chang: added cv_lassopcrmatlab to use new matlab lasso.
%    This uses coordinate descent methods compared to Least Angle Regression
%    in the Rocha Lasso algorithm.  This function can be used for prediction of continuous data or classification
%    of binary data.  For binary classification compare cverr to, Phi, or MSE in stats output. Still
%    working on a good way to quantify accuracy using the predicted
%    probabilities from logistic regression.  Also working on methods to
%    select optimal lambda for regularization.  Make sure you change the Rocha
%    lasso function name to 'lasso_rocha' to retain old functionality and prevent conflicts.
%
%    11/7/12: Luke Chang: added functionality to display xval iteration number
%    and elapsed time
%
%    11/8/12: Luke Chang: fixed calculation of mse. Previous version calculated sum of squared error
%
%    11/10/12: Tor Wager:
%              - forced double format at input processing stage
%              - added optional inputs to lassoglm in cv_lassopcrmatlab
%              - added functionality to use variable number of components
%              - forced remove_empty on object to avoid potential problems with empty variables
%              - added/improved documentation, cleaned up code structure for cv_lassopcrmatlab
%
%    11/28/12: Luke Chang:
%              -fixed bug with optional inputs in cv_lassopcrmatlab - lassoglm
%              didn't like empty cells in varargin
%              -fixed bug with replace_empty(obj) that was returning brain
%              weights that were a different size from the input***
%              -cleaned code structure for cv_lassopcrmatlab
%              -added functionality to estimate lambda and alpha for
%              lasspcrmatlab using grid search and nested cross validation
%
%    1/14/13: Tor Wager: SVM algorithm use updates
%              - fixed bug getting distance from hyperplane with svm
%              - fixed bug with input of slack var params ('C=x') in svm
%              - return cross-validated distance from hyperplane in special
%              output
%              - added stats.weight_obj; ...and intercept values in stats.other_output_cv{:,3}
%
%    2/21/13: Luke Chang:
%              - fixed bug with cv_lassopcr - now restimates betas selected
%              with lasso using OLS (see hastie, friedman, tibrishani pg 92)
%              - added new functionality to SVM, rbf kernel
%              - added new functionality to SVM can estimate slack parameter
%              (or sigma if using rbf) using nested xVal.  For some reason
%              this isn't working super well yet, could have something to do
%              with loss function.  Someone should look into this.
%
%    2/26/13: Luke Chang:
%              - added new functionality to cv_svm - can now perform
%              classification of multiple classes using one vs rest
%
%    3/5/13: Tor Wager:
%              - lassopcr : fixed bug in lasso_num to return correct num components
%              - added nopcr option to cv_lassopcr and cv_lassopcrmatlab
%              - cv_lassopcrmatlab: fixed bug: lasso_num not selecting correctly
%              was not re-fitting OLS solution after selecting components with lasso
%
%
%    3/6/13:  Tor Wager
%              - fixed minor bug introduced in last version with lassopcr/full rank
%              - added stats.weight_obj, an fmri_data object with the predictive weights (from full sample)
%
%    3/8/13: Tor Wager
%              - working on bootstrapping, not a full solution yet
%              - rocha lasso seems to be ok, but others not (cv_pcr, matlab
%              lasso) returning stable weights
%              - added code to use only non-redundant components in pcr, to
%              avoid warnings/instability during bootstrapping.
%
%    3/23/13: Wani Woo
%              - implemented 'EstimateParams' option in rocha lasso. With this
%                option, you can use the subset of coefficients that are non-zero
%                predictors based on the optimal lambda estimation. The optimal
%                lambda will be chosen based on the minimization of mean square
%                error (MSE) of netsed cross-validation.
%              - For this, cv_lassopcr now calls "lasso_cv.m" instead of "lasso_rocha.m"
%                , but lasso_cv calls lasso_rocha. lasso_cv gives the results of nested
%                cross-validation, which is using to select the optimal lambda,
%                in addition to all the outputs that lasso_rocha gives. For this
%                reason, I think it is beneficial to use lasso_cv instead of
%                lasso_rocha.
%              - In order to make this possible, I added a variable,
%                cv_assignment, into funhan. However, only cv_lassopcr is
%                actually using the variable.
%
%    6/10/13: Luke Chang
%              -Added balanced_ridge option to SVM.
%              -fixed some bugs with the cv_svm multiclass support
%
%    6/25/13: Luke Chang
%              -Added try/catch on lassopcr.  Will output nans if there is a
%              problem running PCA.  Should help with problems with svd
%              convergence when bootsrapping.
%              -turned off parallel by default for bootstrapping.  memory is
%              duplicated for each worker so will crash if not enough memory.
%              Also bootstrp seems to preallocate.
%
%    6/29/13: Luke Chang
%              -added option to save bootstrap weights 'savebootweights'.
%              This is useful if you want to aggregate bootstrap samples from
%              multiple iterations of predict() run at different times or
%              on different computers
%
%    7/2/13: Luke Chang
%              -added ability to not cross-validate by setting nfolds, k=1
%              -made a bunch of changes to bootstrapping to reduce memory
%              demands and fixed bugs to output weights.
%              -added nancorr to ignore nans when calculating correlation
%              -added rng 'shuffle' to ensure that bootstrapping will use
%              different inital seed.  VERY IMPORTANT for aggregating across
%              multiple boostrap sessions!
%
%    11/28/13: Luke Chang
%              -added ability to use hv block cross-validation, which is good
%              for timeseries data with stationary autocorrelation.
%              Use 'hvblock,[h,v]
%
%    12/16/13: Luke Chang
%              -added ability to use rolling block cross-validation with the 
%              ability to deal with autocorrelation, which is good
%              for timeseries data with stationary autocorrelation.
%              Use 'rolling,[h,v,g]
%    4/3/14: Luke Chang
%              -fixed bug with SVM, cross-validated distance from hyper plane
%              (hopefully, the distance from hyper plane is still correct)
%              -fixed bug with SVM 'nfolds',1
%
%    2/28/15: Luke Chang
%              -added platt scalling option for SVM
%
%    4/7/15: Wani Woo: SVM algorithm use updates
%              - fixed bug with input of slack var params ('C=x') in svm
%
%    5/7/15: Anjali and Wani: replaced princomp for PCA with SVD on transpose
%              training data for three algorithms, cv_pcr, cv_lassopcr,
%              cv_lassopcrmatlab. This reduces running time substantially.
%
%    9/4/2015: Tor: Created more functional statistic_image output in
%    weight_obj when bootstrapping.  Now p-values, etc. are included so you
%    can threshold.
% 
%    1/17/2016: Stephan: added optional stats-output from matlab lassoglm 
%                 to be returned in out.other_output{3}.stats
% 
%    2/6/2017: Stephan replaced 'pc(:,size(xtrain,1)) = [];' with 'pc(:,end) = [];
%    to accomodate predictor matrices with fewer features (voxels) than
%    data (trials/images)
%
%    6/2017: Tor and Phil Kragel: update parallel processing for
%    bootstrapping. Default is to use parallel.
%
%
%    4/22/2020: Phil Kragel, add PLS regression assuming univariate outcome
% ..

% ..
%    ---------------------------------------------------------------------
%    Defaults
%    ---------------------------------------------------------------------
% ..

nfolds = 5;
tsxval_hvblock = 0;
tsxval_rolling = 0;
%error_type = 'mcr';            % mcr, mse: misclassification rate or mean sq. error
algorithm_name = 'cv_pcr';      % name of m-file defining training/test function
useparallel = 'always';         % Use parallel processing, if available, for bootstrapping. Currently no parallel for xval. always = do it, anything else, don't use parallel
bootweights = 0;                % bootstrap voxel weights
verbose = 1;
doMultiClass = 0;               % option to run multiclass SVM

if length(unique(obj.Y)) == 2
    error_type = 'mcr';
else
    error_type = 'mse';
end

if size(obj.dat, 2) ~= size(obj.Y, 1)
    if ~strcmp('MultiClass',varargin) %Added by LC 2/26/13 for new svm multiclass functionality
        error('obj.dat must be [Predictors x observations] and obj.Y must be [observations x 1]');
    end
end

% Mandatory conditioning of data object
% ---------------------------------------------------------------------

%obj = remove_empty(obj);
%11/27/12: Luke Chang: This seems to be
%affecting the ability to use orthviews.  This should be either removed or
%a 'replace_empty(obj) needs to be added somewhere. % TOR: orthviews should
%be able to handle this by doing replace_empty there if needed

% force double to avoid various problems
obj.dat = double(obj.dat);
obj.Y = double(obj.Y);

% ---------------------------------------------------------------------
% Split inputs into those that control crossval and those that
% belong to the algorithm (can be different for different
% algorithms).
% ---------------------------------------------------------------------

predfun_inputs = {};
bootfun_inputs = {}; % for setting number of bootstrap samples
nfold = []; % to catch typos

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            % functional commands
            
            case {'noparallel'}
                useparallel = 'never';
                
            case {'nfolds', 'nfold' 'error_type', 'algorithm_name', 'useparallel', 'verbose'}
                str = [varargin{i} ' = varargin{i + 1};'];
                eval(str)
                varargin{i} = [];
                varargin{i + 1} = [];
                
            case {'cv_pcr', 'cv_multregress', 'cv_univregress', 'cv_svr', 'cv_lassopcr', 'cv_svm','cv_lassopcrmatlab' ,'cv_pls'}
                algorithm_name = varargin{i};
                
            case {'bootweights'}
                bootweights = 1; varargin{i} = [];
                
            case 'bootsamples'
                bootfun_inputs{end+1} = 'bootsamples';
                bootfun_inputs{end+1} = varargin{i+1};
                bootweights = 1;
                varargin{i} = [];
                varargin{i+1} = [];
                
            case 'savebootweights'
                bootfun_inputs{end+1} = 'savebootweights';
                bootweights = 1;
                varargin{i} = [];
           
            case 'MultiClass'
                doMultiClass = 1;
                predfun_inputs{end + 1} = varargin{i}; %PK - pass as input to cv_svm
            case 'rolling'
                tsxval_rolling = 1;
                inval = varargin{i + 1};
                h = inval(1);
                v = inval(2);
                g = inval(3);
                varargin{i} = [];
                varargin{i + 1} = [];
                
            case 'hvblock'
                tsxval_hvblock = 1;
                inval = varargin{i + 1};
                h = inval(1);
                v = inval(2);
                varargin{i} = [];
                varargin{i + 1} = [];
                
            otherwise
                % use all other inputs as optional arguments to specific
                % functions, to be interpreted by them
                predfun_inputs{end + 1} = varargin{i};
                
                if (i+1) <= length(varargin) && ~ischar(varargin{i + 1})
                    predfun_inputs{end + 1} = varargin{i + 1};
                end
                
        end
    end
end

if ~isempty(nfold) % catch a potentially problematic typo
    nfolds = nfold;
end

% Set up parallel processing input to bootstrap function
switch useparallel
    case 'always'
       bootfun_inputs{end+1} = 'parallel';
    otherwise
        % no parallel
end

% This is not used.
% opt = statset('crossval');
% opt.UseParallel = useparallel;

% ---------------------------------------------------------------------
% stratified partition, or custom holdout set
% ---------------------------------------------------------------------

% NOTE: COULD REPLACE ALL THIS WITH cvpart = stratified_holdout_set(Y, varargin)

if numel(nfolds) > 1 % a vector; assume integers for holdout sets
    % Custom holdout set
    fold_indicator = nfolds;
    u = unique(fold_indicator);
    nfolds = length(u);
    
    cvpart = struct('NumTestSets', nfolds);
    [trIdx, teIdx] = deal(cell(1, nfolds));
    
    for i = 1:length(u)
        teIdx{i} = fold_indicator == u(i);
        trIdx{i} = ~teIdx{i};
    end
    
elseif nfolds == 1 % special for 1 fold: use all obs for train and test; good for bootstrapping weights
    [trIdx, teIdx] = deal(cell(1, nfolds));
    trIdx={ones(size(obj.Y))};
    teIdx={ones(size(obj.Y))};
    cvpart = struct('NumTestSets', length(trIdx)); %LC: 4/3/14: added this as it was missing
        
elseif tsxval_hvblock == 1 %special case for timeseries CV using HVBlock : added LC: 11/28/13
    [trIdx, teIdx] = tscv(length(obj.Y), 'hvblock',[h,v]);
    cvpart = struct('NumTestSets', length(trIdx));
    
elseif tsxval_rolling == 1 %special case for timeseries CV using HVBlock : added LC: 12/16/13
    [trIdx, teIdx] = tscv(length(obj.Y), 'rolling',[h,v,g]);
    cvpart = struct('NumTestSets', length(trIdx));
    
else
    
    [trIdx, teIdx] = deal(cell(1, nfolds));
    
    % Classification: Stratified holdout set
    if doMultiClass
        
        %cvpart = cvpartition(obj.Y,'k',nfolds); %changed by LC 2/26/13 to account for multiclass svm
        cvpart = cvpartition(length(obj.Y),'k',nfolds);
        
        for i = 1:cvpart.NumTestSets
            
            trIdx{i} = cvpart.training(i);
            teIdx{i} = cvpart.test(i);
            
        end
        
    else
        % Regression: custom stratification
        % cvpartition object will not stratify continuous values
        % do our own
        
        cvpart = struct('NumTestSets', nfolds);
        
        [ys, indx] = sort(obj.Y);
        for k = 1:nfolds
            
            wh_holdout = indx(k:nfolds:end);
            if isempty(wh_holdout), error('Holdout set construction error: Dataset too small?'); end
            
            teIdx{k} = false(size(obj.Y));
            teIdx{k}(wh_holdout) = true;
            
            trIdx{k} = ~teIdx{k};
        end
    end
end

if verbose
    fprintf('Cross-validated prediction with algorithm %s, %3.0f folds\n', algorithm_name, nfolds)
end

% special for 1 fold: use all obs for train and test; good for
% bootstrapping weights
if nfolds == 1
    %trIdx{1} = teIdx{1}; %LC: this isn't actualy using all of the weights.
    % 3/23/13 Wani Woo added the following three lines to get cv_assignment. See Programmers' note for the detail.
    cv_assignment = double(teIdx{1});
else
    cv_assignment = sum((cat(2, teIdx{:}).*repmat((1:size(cat(2, teIdx{:}), 2)), size(cat(2, teIdx{:}), 1), 1))')';
end

% ---------------------------------------------------------------------
% build function handle, given custom inputs to function
% ---------------------------------------------------------------------

% 3/23/13 Wani Woo added cv_assignment in funstr to feed this variable to cv_lassopcr
funstr = ['@(xtrain, ytrain, xtest, cv_assignment) ' algorithm_name '(xtrain, ytrain, xtest, cv_assignment, predfun_inputs{:})'];

eval(['funhan = ' funstr ';'])


% ---------------------------------------------------------------------
% Cross-validated prediction and error
% ---------------------------------------------------------------------

% Use this loop instead of the crossval function
% it's not complicated and allows more control.
% yfit = crossval(funhan, obj.dat',obj.Y, 'Partition', cvpart, 'Options', opt);

yfit = zeros(size(obj.Y));

noptout = nargout(algorithm_name) - 1;
optout = cell(1, noptout);

t1 = clock;

% Fit on all data - apparent loss, optional outputs
[yfit_all, optout{:}] = funhan(obj.dat', obj.Y, obj.dat', cv_assignment);
% 3/23/13 Wani Woo added cv_assignment.

if verbose
    %fprintf(1,'\n_______________________________\n')
    [hour, minute, second] = sec2hms(etime(clock,t1));
    fprintf(1,'\nCompleted fit for all data in: %3.0f hours %3.0f min %2.0f secs \n',hour,minute,second);
end

% ---------------------------------------------------------------------
% Cross-validated prediction and error
% ---------------------------------------------------------------------

% nworkers = matlabpool('size');

% if nworkers
%     % slice
%     for i = 1:cvpart.NumTestSets
%         dat_tmp{i} = obj.dat(:, trIdx{i})';
%     end
%
%     [yfit_tmp, cvtmp] = deal(cell(1, cvpart.NumTestSets));
%
%     parfor i = 1:cvpart.NumTestSets
%
%         [yfit_tmp{i}, cvtmp{i}] = funhan(dat_tmp{i}, obj.Y(trIdx{i},:), obj.dat(:, teIdx{i})');
%
%     end
%
%     % reassemble
%     for i = 1:cvpart.NumTestSets
%         yfit(teIdx{i}) = yfit_tmp{i};
%         cv_optout{i, :} = cvtmp{i};
%     end
%     clear yfit_tmp dat_tmp cvtmp
%
% else

if nfolds ~= 1 %skip if using nfolds 1, added 7/2/13 by LC
    
    cv_optout = cell(cvpart.NumTestSets, noptout);
    
    for i = 1:cvpart.NumTestSets
        
        % 3/23/13 Wani Woo added the following three lines to get a proper cv_assignment variable for cross-validation.
        
        teIdx_cv = teIdx(1:end ~= i);
        cv_assignment = sum((cat(2, teIdx_cv{:}).*repmat((1:size(cat(2, teIdx_cv{:}), 2)), size(cat(2, teIdx_cv{:}), 1), 1))')';
        cv_assignment(cv_assignment == 0) = [];
        
        % 9/8/13 Yoni Ashar added to pass correct subjIDs to each fold
        if any(strcmp('subjIDs', predfun_inputs))
            
            subjIDs_index =  1 + find(strcmp('subjIDs', predfun_inputs));
            if ~exist('subjIDsFull', 'var') %first time through, save the full array
                subjIDsFull = predfun_inputs{subjIDs_index};
            end
            
            %overwrite subjIDs with appropriate subset of subjIDs
            predfun_inputs{subjIDs_index} = subjIDsFull(trIdx{i});
            eval(['funhan = ' funstr ';']) % "reload" predfun_inputs
        end
        
        t2 = clock;
        
        %    [yfit(teIdx{i},:), cv_optout{i, :}] = funhan(obj.dat(:, trIdx{i})', obj.Y(trIdx{i},:), obj.dat(:, teIdx{i})');
        [yfit(teIdx{i},:), cv_optout{i, :}] = funhan(obj.dat(:, trIdx{i})', obj.Y(trIdx{i},:), obj.dat(:, teIdx{i})', cv_assignment); %changed by lc on 2/26/13 to add svm multiclass support
        
        % 3/23/13 Wani Woo added cv_assignment.
        if verbose
            [hour, minute, second] = sec2hms(etime(clock,t2));
            fprintf(1,['Fold ' num2str(i) '/' num2str(cvpart.NumTestSets) ' done in: %3.0f hours %3.0f min %2.0f sec\n'],hour,minute,second);
        end
    end
    if verbose
        [hour, minute, second] = sec2hms(etime(clock,t1));
        fprintf(1,'\nTotal Elapsed Time = %3.0f hours %3.0f min %2.0f sec\n',hour, minute, second);
    end
else %added 7/213/ by lc for adding no cv option
    
    yfit = yfit_all;
    %     cvpart = [];
    cv_optout = [];
end


% ---------------------------------------------------------------------
% Get error
% ---------------------------------------------------------------------

switch error_type
    case {'mcr', 'class_loss'}
        %err = obj.Y ~= yfit;
        err = obj.Y ~= round(yfit); %10/7/12: Luke Chang: this will allow mcr to also be calculated using predicted probabilities
        
        %cverr = sum(err) ./ length(err);
        cverr = sum(err) ./ length(err);
        
        phi = corr(obj.Y, yfit); %10/7/12: Luke Chang: this will calculate phi correlation coefficient between two binary variables
        
    case {'mse' 'rmse', 'meanabserr','r'}
        err = obj.Y - yfit;
        
        %mse = mean(err' * err); %10/8/12: Luke Chang: I think this is only capturing sum of squared error
        mse = (err' * err)/length(err);  %This should be correct calculation of mean squared error
        rmse = sqrt(mse);
        meanabserr = nanmean(abs(err));  %if you are getting strange error here make sure yo are using matlab default nanmean (e.g., which nanmean)
        r = corrcoef(obj.Y, yfit, 'rows', 'pairwise');
        r = r(1, 2);
        
        eval(['cverr = ' error_type ';']);
        
    otherwise
        error('Illegal loss function type')
end


% ---------------------------------------------------------------------
% collect output
% ---------------------------------------------------------------------
stats = struct('Y', obj.Y, 'algorithm_name', algorithm_name, ...
    'function_call', funstr, 'function_handle', funhan, ...
    'yfit', yfit, 'err', err, 'error_type', error_type, 'cverr', cverr, ...
    'nfolds', 'nfolds', 'cvpartition', cvpart);

stats.teIdx = teIdx;
stats.trIdx = trIdx;

stats.other_output = optout;
stats.other_output_descrip = 'Other output from algorithm - trained on all data (these depend on algorithm)';
stats.other_output_cv = cv_optout;
stats.other_output_cv_descrip = 'Other output from algorithm - for each CV fold';

%LC: 4/3/14: Removed this as it is causing errors and seems redundant with line 732.
% %For SVM reorder distance from hyperplane: 12/6/12: Luke Chang: not very elegant solution
% if strcmp(algorithm_name,'cv_svm')
%     stats.other_output_cv{3}=ones(length(stats.Y),1);
%     for i = 1:cvpart.NumTestSets
%         stats.other_output{3}(find(teIdx{i}==1)) = stats.other_output_cv{i,2};
%     end
% end

switch error_type
    case {'mcr', 'class_loss'}
        
        stats.phi = phi;
        
        if strcmp(algorithm_name,'cv_lassopcrmatlab')
            stats.me = mean(abs(obj.Y-yfit)); %calculate average distance from probability and outcome - should be more sensitive measure of accuracy than mcr for probabilities
        end
        
    case {'mse' 'rmse', 'meanabserr'}
        
        stats.mse = mse;
        stats.rmse = rmse;
        stats.meanabserr = meanabserr;
        stats.pred_outcome_r = r;
end

% ---------------------------------------------------------------------
% Special output
% ---------------------------------------------------------------------
switch algorithm_name
    
    case 'cv_svm'
        % cross-validated distance from hyperplane
        if nfolds ~= 1
            dd = zeros(size(obj.Y));
            
            for i = 1:length(teIdx)
                dd(teIdx{i},:) = stats.other_output_cv{i, 2}; %added by lc 6/10/13 to accomodate multiclass
            end
            
            stats.dist_from_hyperplane_xval = dd;
        end
        
end

if length(stats.other_output{1}) == size(obj.dat, 1)
    % other_output{1} is right size for weight vector
    % assume it is and return weight_obj
    w = mean(obj);
    w.dat = stats.other_output{1};
    stats.weight_obj = w;
end


% ---------------------------------------------------------------------
% Bootstrap weights, if asked for
% ---------------------------------------------------------------------

if nfolds == 1
    %trIdx{1} = teIdx{1};
    % 3/23/13 Wani Woo added the following three lines to get cv_assignment. See Programmers' note for the detail.
    cv_assignment = double(teIdx{1});
else
    cv_assignment = sum((cat(2, teIdx{:}).*repmat((1:size(cat(2, teIdx{:}), 2)), size(cat(2, teIdx{:}), 1), 1))')';
end

if bootweights
    stats.WTS = boot_weights(funhan, obj, cv_assignment,doMultiClass, bootfun_inputs{:}); %1/20/16 add multiclass support
    fprintf('Returning weights and Z, p values in stats.WTS\n')
    
    % Replace weight_obj with better statistic_image version:
    
    stats.weight_obj = statistic_image('dat', stats.weight_obj.dat, ...
    'volInfo', stats.weight_obj.volInfo, ...
    'p', stats.WTS.wP', ...
    'ste', stats.WTS.wste', ...
    'dat_descrip', stats.function_call, ...
    'removed_voxels', stats.weight_obj.removed_voxels);

end


end % main function




% ---------------------------------------------------------------------
% ---------------------------------------------------------------------

% ALGORITHMS

% ---------------------------------------------------------------------
% ---------------------------------------------------------------------


function [yfit, b, intercept] = cv_multregress(xtrain, ytrain, xtest, cv_assignment, varargin)

b = glmfit(xtrain, ytrain); % Note: many times faster than inv(X'*X)X' for large k
yfit = b(1) + xtest * b(2:end);

intercept = b(1);
b = b(2:end);

end

% ----------------------------- algorithms -------------------------------

function [yfit, b, intercept] = cv_multilevel_glm(xtrain, ytrain, xtest, cv_assignment, varargin)

disp('2nd LEVEL PREDICTORS NOT CURRENTLY SUPPORTED')

i=find(strcmp('subjIDs', varargin));
if numel(i) ~= 1
    error('When using multilevel glm, must pass in ''subjIDs'' followed by an array specifying which subject every trial belongs to')
end
subjIDs = varargin{i+1};

[Xml, Yml] = deal(cell(size(unique(subjIDs))));

diffs = diff(subjIDs);

% build cell array of X and Y, using a varargin input identifying subjects.
count = 1;
for i=1:numel(ytrain)
    Xml{count}(end+1,:) = xtrain(i,:);
    Yml{count}(end+1,:) = ytrain(i);
    
    if i~=numel(ytrain) && diffs(i) >= 1, count = count+1; end % on to a new subject.  all of a subject's trials must be adjacent to each other.
end

% call ml glm
stats = glmfit_multilevel(Yml, Xml, [], varargin{:});

% test on xtest
intercept = stats.beta(1);
b = stats.beta(2:end)';
yfit = intercept + xtest * b;

end

% ----------------------------- algorithms -------------------------------

function [yfit, vox_weights, intercept] = cv_pls(xtrain, ytrain, xtest, cv_assignment, varargin)

% Choose number of components to save [optional]
wh = find(strcmp(varargin, 'numcomponents'));
if ~isempty(wh) && length(varargin) >= wh + 1
    
    numc = varargin{wh + 1};
    
else

    numc = [];
end

if ~isempty(numc)
[~,~,~,~,b]=plsregress(xtrain,ytrain,numc);
else %use max possible
[~,~,~,~,b]=plsregress(xtrain,ytrain);
end
vox_weights = b(2:end);

intercept = b(1);

yfit = intercept + xtest * vox_weights;

end


% ----------------------------- algorithms -------------------------------

function [yfit, vox_weights, intercept] = cv_pcr(xtrain, ytrain, xtest, cv_assignment, varargin)

[pc,~,~] = svd(scale(xtrain,1)', 'econ'); % replace princomp with SVD on transpose to reduce running time. 
pc(:,end) = [];                % remove the last component, which is close to zero.
                               % edit:replaced 'pc(:,size(xtrain,1)) = [];' with
                               % end to accomodate predictor matrices with
                               % fewer features (voxels) than trials. SG
                               % 2017/2/6                              
                               
% [pc, sc, eigval] = princomp(xtrain, 'econ');

% Choose number of components to save [optional]
wh = find(strcmp(varargin, 'numcomponents'));
if ~isempty(wh) && length(varargin) >= wh + 1
    
    numc = varargin{wh + 1};
    
    if numc > size(pc, 2)
        disp('WARNING!! Number of components requested is more than unique components in training data.');
        numc = size(pc, 2);
    end
    pc = pc(:, 1:numc);
end

sc = xtrain * pc;

if rank(sc) == size(sc,2)
    numcomps = rank(sc); 
elseif rank(sc) < size(sc,2)
    numcomps = rank(sc)-1;
end

% 3/8/13: TW:  edited to use numcomps, because sc is not always full rank during bootstrapping
X = [ones(size(sc, 1), 1) sc(:, 1:numcomps)];

if rank(X) <= size(sc, 1)
    b = pinv(X) * ytrain; % use pinv to stabilize; not full rank
    % this will only happen when bootstrapping or otherwise when there are
    % redundant rows
else
    b = inv(X'*X)*X'*ytrain;
end

% Programmers' notes: (tor)
% These all give the same answer for full-rank design (full component)
% X = [ones(size(sc, 1), 1) sc];
% b1 = inv(X'*X)*X'*ytrain;
% b2 = pinv(X)*ytrain;
% b3 = glmfit(sc, ytrain);
% tic, for i = 1:1000, b1 = inv(X'*X)*X'*ytrain; end, toc
% tic, for i = 1:1000, b2 = pinv(X)*ytrain; end, toc
% tic, for i = 1:1000, b3 = glmfit(sc, ytrain); end, toc
% b1 is 6 x faster than b2, which is 2 x faster than b3

vox_weights = pc(:, 1:numcomps) * b(2:end);

intercept = b(1);

yfit = intercept + xtest * vox_weights;

end

% ----------------------------- algorithms -------------------------------


function [yfit, vox_weights, intercept, out_cv] = cv_lassopcr(xtrain, ytrain, xtest, cv_assignment, varargin) % 3/23/13 Wani Woo added out_cv.

doSkip = 0; %skip predict and output null \\
dopcr = 1;
doNestedXval = 0; % 3/23/13 added by Wani Woo.

wh = find(strcmp(varargin, 'nopcr'));
if ~isempty(wh), dopcr = 0; end

wh = find(strcmpi(varargin, 'estimateparams')); % 3/23/13 added by Wani Woo.
if ~isempty(wh), doNestedXval = 1; end

if dopcr
    try %catch added 7/2/13 by LC
        [pc,~,~] = svd(scale(xtrain,1)', 'econ'); % replace princomp with SVD on transpose to reduce running time. 
        pc(:,end) = [];        % remove the last component, which is close to zero.
                               % edit:replaced 'pc(:,size(xtrain,1)) = [];' with
                               % end to accomodate predictor matrices with
                               % fewer features (voxels) than trials. SG
                               % 2017/2/6
        % [pc, sc, eigval] = princomp(xtrain, 'econ');

    catch exception
        sprintf(['ERROR! ' exception.message '\nSkipping lassoPCR'])
        doSkip = 1;  %skip predict and output null (useful for bootstrapping will ignore problems with convergence)
        sc = xtrain;
    end
    
    if ~doSkip %skip if problem with pca
        % Choose number of components to save [optional]
        wh = find(strcmp(varargin, 'numcomponents'));
        if ~isempty(wh) && length(varargin) >= wh + 1
            
            numc = varargin{wh + 1};
            
            if numc > size(pc, 2)
                disp('WARNING!! Number of components requested is more than unique components in training data.');
                numc = size(pc, 2);
            end
            pc = pc(:, 1:numc);
        end
        
        sc = xtrain * pc;
    end
    
else
    sc = xtrain;
end


% 3/14/13: tor changed because full rank was still returning unstable
% matrices (non-invertible) sometimes.
if rank(sc) == size(sc,2)
    numcomps = rank(sc); 
elseif rank(sc) < size(sc,2)
    numcomps = rank(sc)-1;
end

% numcomps = rank(sc, .001)-1;
%numcomps = rank(sc, mean(abs(sc(:))) ./ 10000);

%X = [ones(size(sc, 1), 1) sc(:, 1:numcomps)];

%11/7/12: LC:  Renamed to avoid conflict with matlab function
% 3/8/13: TW:  edited to use numcomps, because sc is not always full rank during bootstrapping
% 3/23/13: Wani: changed lasso_rocha to lasso_cv. Using options, now we can
%                actually reduce the features.

options.alignment = 'lambda';
options.cv_assignment = cv_assignment;

if ~doSkip
    
    % 3/23/13: Wani: add dopcr option.
    if dopcr
        if doNestedXval, out = lasso_cv(ytrain, sc(:, 1:numcomps), length(unique(cv_assignment)), options);
        else out = lasso_rocha(ytrain, sc(:, 1:numcomps)); end
    else
        if doNestedXval, out = lasso_cv(ytrain, sc, length(unique(cv_assignment)), options);
        else out = lasso_rocha(ytrain, sc); end
    end
    % out = lasso_rocha(ytrain, sc(:, 1:numcomps));
    
    n = size(out.beta, 1) - 1;
    wh = find(strcmp(varargin, 'lasso_num'));
    
    if ~isempty(wh) && length(varargin) >= wh + 1
        n = varargin{wh + 1};
        if n > size(out.beta, 1)
            error('You asked for more components/vars than there are in the dataset.')
        end
    elseif ~isempty(wh)
        error('Follow lasso_num input with number of vars to keep');
    end
    
    %output original lasso parameters - removed 2/20/13 by LC
    % b = [out.intercept(n); out.beta(n, :)'];
    % vox_weights = pc * b(2:end);
    % intercept = b(1);
    % yfit = intercept + xtest * vox_weights;
    
    %refit lasso selected betas using OLS (see hastie, Friedman, Tibrishani, pg92 - added 2/20/13 by LC
    % use n + 1 because 1st row has no non-zero betas
    % 3/23/13: Wani: Actually out.beta(n+1,:) at a default setting (i.e., last
    %                set of betas) is a full model. To use a reduced model, we
    %                need to first estimate the optimal lambda that minimizes
    %                the MSE of cross-validation models.
    
    if ~doNestedXval
        wh_beta = logical([out.intercept(n+1); out.beta(n+1, :)' ~= 0]);
        out_cv = [];
    else
        wh_beta = logical([out.CV.intercept; out.CV.beta' ~= 0]);
        out_cv.lambda = out.CV.lambda;
        out_cv.min_MSE = out.CV.min_MSE;
    end
    
    Xsubset = [ones(size(sc, 1), 1) sc(:, wh_beta(2:end))];
    betatmp = pinv(Xsubset) * ytrain;
    betas = zeros(size(wh_beta));
    betas(wh_beta) = betatmp;
    
    if dopcr
        vox_weights = pc(:, 1:numcomps) * betas(2:end);
    else
        vox_weights = betas(2:end);
    end
    
    intercept = betas(1);
    yfit = intercept + xtest * vox_weights;
else
    yfit = nan(size(ytrain));
    vox_weights = nan(size(xtest,2),1);
    intercept = nan;
    out_cv = [];
    
end

end

% ----------------------------- algorithms -------------------------------

function [yfit, vox_weights, intercept, out] = cv_lassopcrmatlab(xtrain, ytrain, xtest, cv_assignment, varargin)
% Notes:
% Added 11/7/12 by Luke Chang
% This will only run on Matlab 2012 or newer (i.e., 7.14.0) when the lasso function was implemented
% Lasso regularization %standardizes x & y by default
% This code can also be used to run ridge or elastic net regression by setting alpha to - 0 = ridge; 1 = lasso; anything else = elastic net
% This code can do prediction or classification.  Important to note that weights in classification are in log-odds
% Can estimate optimal lambda or alpha through nested cross validation
% Should add functionality to select lambda, alpha, or both.


% ---------------------------------------------------------------------
% Defaults
% ---------------------------------------------------------------------

alpha = 1; %indicates type of regularization to run for cv_lassopcrmatlab.  runs lasso by default
n = 1; %Set Lasso_num
chooseComponents = 0; %indicate whether to choose a custom number of components
chooseLassoNum = 0; %indicate whether to choose a custom number of lasso numbers
doNestedXval = 0; %use nested cross-validation to pick lambda
estAlpha = 0;
estLambda = 0;
dopcr = 1; % in  case we do not want to use PCR

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            case 'nopcr'
                dopcr = 0; % in  case we do not want to use PCR
                varargin{i} = [];
                
            case 'numcomponents'
                varargin{i} = [];
                numc = varargin{i + 1};
                if isempty(numc) || ~isnumeric(numc)
                    error('Follow numcomponents with number of components to keep');
                end
                varargin{i+1} = [];
                chooseComponents = 1;
            case 'lasso_num'
                varargin{i} = [];
                lassonum = varargin{i + 1};
                if isempty(lassonum) || ~isnumeric(lassonum)
                    error('Follow lasso_num input with number of vars to keep');
                end
                varargin{i+1} = [];
                chooseLassoNum = 1;
            case 'Alpha'
                varargin{i} = [];
                alpha = varargin{i+1};
                if isempty(alpha) || ~isnumeric(alpha) || alpha>1 || alpha<0
                    error('Follow Alpha input with elastic net mixing parameter');
                end
                varargin{i+1} = [];
            case 'EstimateParams'
                %Need to add functionality to selectively estimate alpha and/or lambda
                varargin{i} = [];
                params = varargin{i+1};
                varargin{i+1} = [];
                if max(strcmp(params,'Alpha') | strcmp(params,'alpha'))
                    estAlpha = 1;
                end
                if max(strcmp(params,'Lambda') | strcmp(params,'lambda'))
                    estLambda = 1;
                end
                doNestedXval = 1;
        end
    end
end

% check if using correct version of matlab
% ----------------------------------------------------------

if verLessThan('matlab', '7.14.0')
    error('''cv_lassopcrmatlab'' is only compatible with matlab 7.14.0 and newer, try using old version ''cv_lassocvr''');
end

% determine if running prediction (continous) or classification (binary)
% ----------------------------------------------------------

if length(unique(ytrain)) == 2 %run in classification mode
    fprintf(1, 'Classification mode (binomial) ');
    runmode = 'classification';
    distribution = 'binomial';
    varargin{end + 1} = 'Link';
    varargin{end + 1} = 'logit';
    
    
    %check if y is in correct format
    if length(unique(ytrain)) ~= 2  % tor changed this
        error('Classification mode requires binary y with values of 1 or 0');
    end
    
    % enforce logical in case of [1, -1] inputs
    ytrain = ytrain > 0;
elseif any(strcmp(lower(varargin),'poisson'))
    runmode = 'poisson regression';
    fprintf(1,'Regression mode (Poisson) ');
    distribution = 'poisson';
    varargin={};
else %Run in prediction mode
    runmode = 'regression';
    fprintf(1,'Regression mode (Gaussian) ');
    distribution = 'normal';
    
end

% Principal components
% ----------------------------------------------------------
if dopcr
    
    [pc,~,~] = svd(scale(xtrain,1)', 'econ'); % replace princomp with SVD on transpose to reduce running time. 
    pc(:,end) = [];                % remove the last component, which is close to zero
                                   % edit:replaced 'pc(:,size(xtrain,1)) = [];' with
                                   % end to accomodate predictor matrices with
                                   % fewer features (voxels) than trials. SG
                                   % 2017/2/6.
    % [pc, sc, eigval] = princomp(xtrain, 'econ');  %econ is much faster and should be used when n <= p and returns n-1 components.  additional components are meaningless

    % [optional] Choose number of components to save
    if chooseComponents
        if numc > size(pc, 2)
            disp('WARNING!! Number of components requested is more than unique components in training data.');
            numc = size(pc, 2);
        end
        pc = pc(:, 1:numc);
    end
    
    % 3/8/13: TW:  edited to use numcomps, because sc is not always full rank during bootstrapping
    numcomps = rank(pc);
    
    sc = xtrain * pc(:, 1:numcomps);  %this is slightly different from princomp output due to the automatic centering
    
else
    pc = xtrain;
    sc = xtrain;
    numcomps = size(sc, 2);
end

% Regression (penalized, normal or logistic)
% ----------------------------------------------------------


% lassoglm will accept a number of optional inputs, which you can put
% in as optional inputs to predict.m.  We need to eliminate invalid
% inputs earlier, however, or lassoglm will return an error.

varargin(cellfun(@isempty,varargin)) = [];  %clear empty cells to prevent errors in lasso

if ~doNestedXval
    
    [b, stats] = lassoglm(sc, ytrain, distribution, 'Alpha', alpha, varargin{:});
    
    % [Optional] Process inputs on desired penalization
    % ----------------------------------------------------------
    %This allows you to select a lambda performs identically to the original
    %rocha lassopcr (i.e., lambda decreases as numbers increase).  Eventually
    %as lambda approaches zero you will retain all original components.
    
    % n = stats.IndexMinDeviance; %this is the optimal lambda from nested cross validation
    
    if chooseLassoNum
        n_in_model = size(sc, 2) - sum(b == 0); % how many variables have non-zero b
        
        n = min(find(n_in_model == lassonum)); % lowest lambda (least shrinkage) with lassonum components retained
        
        % sometimes this may not exist, depending on sampling res of b (100 steps by default...)
        % if empty, find the next closest one by relaxing lambda a bit.
        indx = 1;
        while isempty(n)
            n = min(find(n_in_model == lassonum + indx));
            indx = indx + 1;
        end
        
        %n = 1 + size(pc, 2) - lassonum;  % # betas, including intercept
        if n < 1
            disp('WARNING!! Number of LASSO components retained is more than unique components in training data.');
            n = 1;
        end
    end
    
    % Collect outputs
    % ----------------------------------------------------------
    
    beta = b(:, n);
    
    % re-fit using OLS for non-zero components
    wh = beta ~= 0;
    
    if ~strcmp(runmode,'poisson regression')
    bb = pinv([ones(size(ytrain)) sc(:, wh)]) * ytrain;
    else
    bb = glmfit(sc(:, wh), ytrain,distribution); %pinv has problems, and isn't flexible maybe check binomial case..
    end
    beta(wh) = bb(2:end);
    
    if dopcr
        vox_weights = pc * beta;
    else
        vox_weights = beta;
    end
    
    intercept = bb(1);  %stats.Intercept(n);  % after OLS re-fitting
    out.alpha = stats.Alpha;
    out.lambda = stats.Lambda(n);
    out.stats  =  stats; % return all info from lassoglm
    
    switch runmode
        case 'classification'
            yfit = intercept + xtest * vox_weights;
            yfit = exp(yfit)./(1+exp(yfit)); %convert to probabilities
        case 'regression'
            yfit = intercept + xtest * vox_weights;

        case 'poisson regression'
            yfit = glmval([intercept; vox_weights], xtest, 'log');
        otherwise error('this should never happen.')
    end
    
else
    %%Use nested xvalidation and grid search to find optimal lambda and alpha parameters
    % NOTE: THIS COULD PROBABLY BE DONE USING LASSOGLM'S INTERNAL
    % CROSS-VALIDATION
    
    %Set Initial Parameters
    gridResolution = 20;
    nFold = 5;
    
    %Create Parameter search grid
    a = linspace(eps,1,gridResolution); %can't be 0
    l = linspace(eps,1,gridResolution);
    
    %set up nested xValidation - implement Tor's stratified code
    [trIdx, teIdx] = deal(cell(1, nFold));
    
    % Classification: Stratified holdout set
    if length(unique(ytrain)) < length(ytrain) / 2
        
        cvpart = cvpartition(ytrain,'k',nFold);
        
        for i = 1:cvpart.NumTestSets
            trIdx{i} = cvpart.training(i);
            teIdx{i} = cvpart.test(i);
        end
        
    else
        % Regression: custom stratification
        % cvpartition object will not stratify continuous values
        % do our own
        
        cvpart = struct('NumTestSets', nFold);
        
        [ys, indx] = sort(ytrain);
        
        for k = 1:nFold
            wh_holdout = indx(k:nFold:end);
            
            if isempty(wh_holdout), error('Holdout set construction error: Dataset too small?'); end
            
            teIdx{k} = false(size(ytrain));
            teIdx{k}(wh_holdout) = true;
            
            trIdx{k} = ~teIdx{k};
        end
    end
    beta = cell(length(a),length(l),nFold);
    stats = cell(length(a),length(l),nFold);
    yfit = cell(length(a),length(l));
    for j = 1:(length(a))
        for k = 1:(length(l))
            yfit{j,k} = ones(length(ytrain),1);
        end
    end
    
    %Run nested xValidation
    for i = 1:nFold
        for j = 1:length(a)
            for k = 1:length(l)
                [beta{j,k,i}, stats{j,k,i}] = lassoglm(sc(trIdx{i},:), ytrain(trIdx{i}), distribution, 'Alpha', a(j), 'Lambda', l(k), varargin{:});
                vox_weights = pc * beta{j,k,i};
                intercept = stats{j,k,i}.Intercept;
                switch runmode
                    case 'classification'
                        yfit{j,k}(teIdx{i}) = intercept + xtrain(teIdx{i},:) * vox_weights; %test fit on nested xval
                        yfit{j,k}(teIdx{i}) = exp(yfit{j,k}(teIdx{i}))./(1+exp(yfit{j,k}(teIdx{i}))); %convert to probabilities
                    case 'regression'
                        yfit{j,k}(teIdx{i}) = intercept + xtrain(teIdx{i},:) * vox_weights;
                    otherwise error('this should never happen.')
                end
            end
        end
    end
    
    %Create matrix of error to select optimal parameters
    err = zeros(length(a),length(l));
    for j = 1:length(a)
        for k = 1:length(l)
            switch runmode
                case 'classification'
                    %err(j,k) = sum(ytrain ~= round(yfit{j,k})) ./ length(ytrain); %classification error
                    err(j,k) = mean(abs(ytrain-yfit{j,k})); %calculate average distance from probability and outcome - should be more sensitive measure of accuracy than mcr for probabilities
                case 'regression'
                    %err(j,k) = ((ytrain - yfit{j,k})' * (ytrain - yfit{j,k}))/length(ytrain - yfit{j,k}); %average sum of squared error
                    err(j,k) = median(abs(ytrain - yfit{j,k})); %median absolute deviation - should be more robust to outliers
            end
        end
    end
    
    %find min of grid - if more than one select randomly from set.
    %Alternatives for multiple minima: 1) smooth err matrix; 2) select
    %minima that retains the most pc components.
    clear jj kk
    [jj,kk] = ind2sub(size(err), find(err == min(min(err))));
    if length(jj) > 1
        s = datasample(jj,1,'Replace',false);
        jj = jj(s);
        kk = kk(s);
    end
    
    %Rerun lasso with cross-validated parameters
    [b, stats] = lassoglm(sc, ytrain, distribution, 'Alpha', a(jj), 'Lambda', l(kk), varargin{:});
    
    % Collect outputs
    % ----------------------------------------------------------
    
    beta = b(:, n);
    vox_weights = pc(:, 1:numcomps) * beta;
    intercept = stats.Intercept(n);
    out.alpha = stats.Alpha;
    out.lambda = stats.Lambda;
    out.gridResolution = gridResolution;
    out.errorMatrix = err;
    
    switch runmode
        case 'classification'
            yfit = intercept + xtest * vox_weights;
            yfit = exp(yfit)./(1+exp(yfit)); %convert to probabilities
        case 'regression'
            yfit = intercept + xtest * vox_weights;
        otherwise error('this should never happen.')
    end
end

end % algorithm

% ----------------------------- algorithms -------------------------------


function [yfit, b, yfit_vox] = cv_univregress(xtrain, ytrain, xtest, cv_assignment, varargin)
% A simple strategy that treats voxels as independent and aggregates
% univariate predictions across all voxels.
% Could be a good strategy of independent, redundant info is contained in
% input variables and there is little covariance across measures.

[n, v] = size(xtrain);

b = zeros(2, v);
onesvec = ones(n, 1);

% Training
fprintf('Training...')
tic
for i = 1:v
    
    % X is brain, Y is outcome, for each voxel
    X = [onesvec xtrain(:, i)];
    
    %b = pinv(X) * Y
    b(:, i) = X \ ytrain;   % should be same, but faster
    
end
toc

% Testing
fprintf('Testing...')
nt = size(xtest, 1);
yfit_vox = zeros(nt, v);

for i = 1:nt
    yfit_vox(i, :) = b(1, :) + xtest(i, :) * b(2, :)';
end
fprintf('Done. \n');

% Integrate separate models into one
yfit = mean(yfit_vox, 2);

b = b';

end % function

% ----------------------------- algorithms -------------------------------


function [yfit, w, dummy, intercept] = cv_svr(xtrain, ytrain, xtest, cv_assignment, varargin)

dataobj = data('spider data', xtrain, ytrain);

% Define algorithm
%svrobj = svr({'C=1', 'optimizer="andre"', kernel({'rbf', 1})})

% slack parameter
slackstr = 'C=1';
wh = find(strcmp('C', varargin));
if ~isempty(wh), slackstr = ['C=' num2str(varargin{wh(1)+1})]; end

svrobj = svr({slackstr, 'optimizer="andre"'});

% Training
fprintf('Training...')
t1 = tic;
[res, svrobj] = train(svrobj, dataobj);
w = get_w(svrobj);
t2 = toc(t1);
fprintf('Done in %3.2f sec\n', t2)

% Testing
fprintf('Testing...')
res2 = test(svrobj, data('test data', xtest, []));
yfit = res2.X; % this is proportional to res.X*w (perfectly correlated), but scale is different
fprintf('\n');
% vec = [res.X res.Y dataobj.X * w' res2.X]
% corrcoef(vec)

w = w';
dummy = 0;
intercept = svrobj.b0;


end % function

% ----------------------------- algorithms -------------------------------
%This function has been replaced by a new one that can also run RBF - 2/21/13 LC

% function [yfit, w, dist_from_hyplane, b0] = cv_svm(xtrain, ytrain, xtest, varargin)
%
% dataobj = data('spider data', xtrain, ytrain);
%
% % Define algorithm
% %svrobj = svr({'C=1', 'optimizer="andre"', kernel({'rbf', 1})})
%
% % slack parameter
% slackstr = 'C=1';
% wh = find(strncmpi('C=', varargin, 2));
% if ~isempty(wh), slackstr = varargin{wh(1)}; end
%
% svrobj = svm({slackstr, 'optimizer="andre"'});
%
% % Training
% fprintf('Training...')
% t1 = tic;
% [res, svrobj] = train(svrobj, dataobj);
% w = get_w(svrobj);
% t2 = toc(t1);
% fprintf('Done in %3.2f sec\n', t2)
%
% % Testing
% fprintf('Testing...')
% res2 = test(svrobj, data('test data', xtest, []));
% yfit = res2.X; % this is proportional to res.X*w (perfectly correlated), but scale is different
% fprintf('\n');
% % vec = [res.X res.Y dataobj.X * w' res2.X]
% % corrcoef(vec)
%
% w = w';
% b0 = svrobj.b0;
%
% dist_from_hyplane = xtest * w + b0;
%
% %dummy = 0;
%
% end % function

% ----------------------------- algorithms -------------------------------
%NOTES:
%what parameters does rbf need?
%figure out how to estimate parameters, does spider have a xval gridsearch
%option?
%


function [yfit, w, dist_from_hyplane, intercept, svmpar] = cv_svm(xtrain, ytrain, xtest, cv_assignment, varargin)
%NOTES:
%-Estimate Parameters is working correctly, but doesn't seem like there is
%much variance in the loss function, which is resulting in the xVal
%procedure selecting the smallest possible values.  One thing to look into
%is the loss function being used.
%-Add multiclass option - make sure Y can be bigger than one column in xval
%loop

% ----------------------------------------------------------
% Defaults
% ----------------------------------------------------------

slackstr = 'C=1'; %slack parameter
dorbf = 0;
doNestedXval = 0; %nested cross validation
doMultiClass = 0; %multiclass classification
doBalanced = 0; %balanced using ridge
doPlatt = 0; %platt-scaling

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            case 'C' %slack parameter
                slackstr = ['C=' num2str(varargin{i+1})];
                varargin{i} = []; varargin{i+1} = [];
                % slack parameter
                % wh = find(strcmpi('C=', varargin));
                % if ~isempty(wh), slackstr = varargin{wh(1)}; end
            case 'rbf' %Use radial basis function as kernel
                varargin{i} = [];
                sig = varargin{i+1};
                varargin{i+1} = [];
                if isempty(sig) || ~isnumeric(sig)
                    error('Follow RBF input with sigma');
                end
                dorbf = 1;
            case 'EstimateParams' %Estimate parameters using nested x val
                varargin{i} = [];
                doNestedXval = 1;
                % train 3 svms with C=1,2,3 and validate with 3 fold cross validation
                %[r,a]=train(gridsel(param(svm,'C',[1,2,3]),{'score=cv;score.folds=3'}),gen(toy)) ;
            case 'MultiClass'
                varargin{i} = [];
                doMultiClass = 1;
            case 'Balanced'
                varargin{i} = [];
                ridgeAmt = varargin{i+1};
                varargin{i+1} = [];
                if isempty(ridgeAmt) || ~isnumeric(ridgeAmt)
                    error('Follow Balanced input with ridge amount (numeric)');
                end
                doBalanced = 1;
            case 'platt_scaling'
                doPlatt = 1;
                varargin{i} = [];
        end
    end
end

dataobj = data('spider data', xtrain, ytrain);

% ----------------------------------------------------------
% Define Algorithm
% ----------------------------------------------------------

if doNestedXval %Estimate parameters via xVal
    
    fprintf('Estimating Parameters via Nested Cross Validation...')
    t1 = tic;
    
    if ~dorbf %only estimate C parameter
        
        %Set Initial Parameters
        gridResolution = 10;
        nFold = 3;
        paramGrid = linspace(eps,5,gridResolution); %Create Parameter search grid
        
        if ~doMultiClass
            [r,a] = train(gridsel(param(svm('optimizer="andre"'),'C',paramGrid),{['score=cv;score.folds=' num2str(nFold)]}),dataobj);
        else
            [r,a] = train(gridsel(param(one_vs_rest(svm('optimizer="andre"')),'C',paramGrid),{['score=cv;score.folds=' num2str(nFold)]}),dataobj);
        end
        
        slackstr = ['C=' num2str(a.C)]; %slack parameter
        
    else %estimate C & Sigma parameter for RBF
        
        cGridRes = 10;
        sigGridRes = 10;
        nFold = 3;
        cGrid = linspace(eps,5,cGridRes); %Create Parameter search grid
        sigGrid = linspace(eps,10,sigGridRes); %Create Parameter search grid
        
        if ~doMultiClass
            [r,a] = train(gridsel(param(svm(kernel('rbf',sig)),{'C','kerparam'},{cGrid,sigGrid}),{['score=cv;score.folds=' num2str(nFold)]}),dataobj);
        else
            [r,a] = train(gridsel(param(one_vs_rest(svm(kernel('rbf',sig),'optimizer="andre"')),{'C','kerparam'},{cGrid,sigGrid}),{['score=cv;score.folds=' num2str(nFold)]}),dataobj);
        end
        
        slackstr = ['C=' num2str(a.C)]; %slack parameter
        sig = a.kerparam; %sigma parameter
        
    end
    
    t2 = toc(t1);
    fprintf('Done in %3.2f sec\n', t2)
    
end

%Create SVM objects with correct parameters (i.e., refit after finding optimal parameter)
if ~dorbf
    svmobj = svm({slackstr, 'optimizer="andre"'});
else
    svmobj = svm({slackstr, 'optimizer="andre"',kernel('rbf',sig)});
end

%Check if using multiclass svm, i.e., one vs rest
if doMultiClass
    svmobj = one_vs_rest(svmobj);
    display('Running in Multiclass Mode')
end

%Check if using balanced ridge
if doBalanced
    svmobj.balanced_ridge = ridgeAmt;
end

%Check if using Platt Scaling
if doPlatt
    platt_svmobj = platt(svmobj);
end

% Training
fprintf('Training...')
t1 = tic;
[res, svmobj] = train(svmobj, dataobj, loss);

% Testing
fprintf('Testing...')
res2 = test(svmobj, data('test data', xtest, []));
yfit = res2.X; % this is proportional to res.X*w (perfectly correlated), but scale is different
fprintf('\n');

t2 = toc(t1);
fprintf('Done in %3.2f sec\n', t2)

% vec = [res.X res.Y dataobj.X * w' res2.X]
% corrcoef(vec)

% ----------------------------------------------------------
% Collect Outputs
% ----------------------------------------------------------

if ~doMultiClass
    w = get_w(svmobj)';
    b0 = svmobj.b0;
    dist_from_hyplane = xtest * w + b0;
    svmpar.kernel = 'linear';
else
    for i = 1:size(res2.X,2)
        w(:,i) = get_w(svmobj{i})';
        b0(i) = svmobj{i}.b0;
        dist_from_hyplane(:,i) = xtest * w(:,i) + b0(i);
        
    end
end

intercept = b0; % 09/24/13 Wani added intercept
svmpar.C = slackstr;
if dorbf
    svmpar.sigma = sig;
    svmpar.kernel = 'rbf';
end

if doPlatt
    [platt_res, platt_svmobj] = train(platt_svmobj, dataobj, loss);
    platt_res2 = test(platt_svmobj, data('test data', xtest, []));
    yfit = platt_res2.X;
    svmpar.A = platt_svmobj.A;
    svmpar.B = platt_svmobj.B;
end

end % function



% ---------------------------------------------------------------------
% ---------------------------------------------------------------------

% BOOTSTRAPPING AND OTHER FUNCTIONS

% ---------------------------------------------------------------------
% ---------------------------------------------------------------------

function WTS = boot_weights(funhan, obj, cv_assignment,doMultiClass, varargin)

%default options
opt = statset('UseParallel', 'never');
doSaveWeights = 0;
bootsamples = 100;

WTS = struct;
rng 'shuffle' %pick random seed for randomization

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            case 'parallel'
                % doParallel = 1; %run in parallel
                varargin{i} = [];
                
                try
                    
                    pool = gcp;
                    opt = statset('UseParallel', 'always');
                    
                catch
                    
                    disp('Problem starting matlab pool - JAVA issue?');
                
                end
                
            case 'bootsamples'
                bootsamples = varargin{i + 1};
                varargin{i} = [];
                varargin{i + 1} = [];
                
            case 'savebootweights'
                doSaveWeights = 1; %Save bootstrap weights
                %weightname = varargin{i + 1};
                varargin{i} = [];
                %varargin{i+1} = [];
        end
    end
end


%let's try directly inputting these to save RAM -luke
%d = obj.dat';
%y = obj.Y;

% anonymous function to return weights, using user's input algorithm
% w = bootfunhan(d, y, cv_assignment);
bootfunhan = @(d, y, cv_assignment) boot_weights_fcn(d, y, cv_assignment, funhan);

%fprintf('Setting up bootstrap samples...');
%bootsam = setup_boot_samples(d, bootsamples);

fprintf('Bootstrapping weights, %3.0f samples...', bootsamples);
t0 = tic;
WTS.w = bootstrp(bootsamples, bootfunhan, obj.dat', obj.Y, cv_assignment, 'Options', opt);
fprintf('Done in %3.0f sec\n', toc(t0));

%[bootstat, bootsam] = bootstrp_havesamples(bootsam,bootfun,varargin)
%[bootstat, bootsam] = bootstrp_havesamples(bootsamples,bootfunhan, d, y);

%7/2/13: Luke changed to use nan mean and outputs mean of weights and to reduce
%duplicating variables

%1/20/16 Phil added code to deal with multiclass case
if doMultiClass
    tv=WTS.w;
    WTS.w=zeros(size(tv,1),size(obj.dat,1),size(cv_assignment,2));
    for i=1:size(tv,1)
        for j=1:size(cv_assignment,2)
             WTS.w(i,:,j)=tv(i,j:size(cv_assignment,2):end);
        end
    end%reshape because MATLAB's bootstrp makes a single row
end

% ToDO:
% This needs to be updated to a bias corrected bootstrap with kernel
% density estimation for sparse data.
WTS.wste = squeeze(nanstd(WTS.w)); %1/20/16 add squeeze for multiclass case
WTS.wmean = squeeze(nanmean(WTS.w)); %1/20/16 add squeeze for  multiclass case
WTS.wste(WTS.wste == 0) = Inf;  % in case unstable regression returns all zeros
WTS.wZ = WTS.wmean ./ WTS.wste;  % tor changed from wmean; otherwise bootstrap variance in mean inc in error; Luke renamed to avoid confusion
WTS.wP = 2 * (1 - normcdf(abs(WTS.wZ)));

if ~doSaveWeights %clear weights if they don't need to be saved (better for memory)
    WTS.w = [];
end


singlemn = bootfunhan(obj.dat', obj.Y, cv_assignment);

if ~doMultiClass %1/20/16 Phil added case for multiclass (average correlation across models)
    if min(size(singlemn)) > 1 % in case bootfunhan returns matrix; bootstrp vectorizes
        singlemn = singlemn(:)';
    end
    %r = corr(wmean, singlemn,'rows','pairwise');
    r=nancorr(WTS.wmean,singlemn); %7/2/13: LC ignore nans
else
    for i=1:size(singlemn,1)
        r(i)=nancorr(WTS.wmean(:,i),singlemn(i,:)');
    end
    r=mean(r);
end
% create_figure('bootweights', 1, 2); plot(wmean, mean(w), 'k.');
% xlabel('Weights - full sample');
% ylabel('Mean bootstrap weights');
% subplot(1, 2, 2); hist(wZ, 100);
% title('Z-scores of voxel weights')



fprintf('Correlation between bootstrapped mean weights and weights: %3.4f\n', r)
fprintf('Voxel weight Z-values range between %3.2f and %3.2f\n', min(WTS.wZ), max(WTS.wZ));

end

% ------------ bootstrapping and other functions -------------------------


function w = boot_weights_fcn(dat, Y, cv_assignment, funhan)
% runs algorithm on all data, returns weights; for bootstrapping

% works for cv_pcr, which has 2 outputs; diff setup here for diff algorithms
optout = cell(1, 2);

[yfit_all, optout{:}] = funhan(dat, Y, dat, cv_assignment);

w = optout{1}';
end

% ------------ bootstrapping and other functions -------------------------

function  bootsam = setup_boot_samples(x, bootsamples)
% bootsam = setup_boot_samples(x, bootsamples)

% initalize random number generator to new values; bootstrp uses this
% there are problems with randomization if not done!
%rand('twister', sum(100*clock))
rng('shuffle');

n = size(x, 1);

bootsam = ceil(n*rand(n, bootsamples)); % x value is the same for all boot samples

maxiter = 500;

% check and replace invalid ones with no variance in predictor
% could do this later for mediator(s)/outcome, but may slow
% down/pose other problems?
wh_bad = ~any(diff(x(bootsam))); %| ~any(diff(m(bootsam)))
icount = 1;

while any(wh_bad)
    bootsam(:, wh_bad) = [];
    %ix = size(bootsam, 2) + 1;
    ntoadd = bootsamples - size(bootsam, 2);
    
    bootsam = [bootsam ceil(n*rand(n, ntoadd))]; % fill out rest of samples
    
    wh_bad = ~any(diff(x(bootsam)));
    icount = icount + 1;
    
    if icount > maxiter
        warning('PREDICT:cannotGetBootSamples', 'Cannot get bootstrap samples; max iterations exceeded.  Sample size too small with categorical predictors?');
        bootsam(:, wh_bad) = [];
        break
    end
end

end

% ------------ bootstrapping and other functions -------------------------

function [hour, minute, second] = sec2hms(sec)
%SEC2HMS  Convert seconds to hours, minutes and seconds.
%
%   [HOUR, MINUTE, SECOND] = SEC2HMS(SEC) converts the number of seconds in
%   SEC into hours, minutes and seconds.

hour   = fix(sec/3600);      % get number of hours
sec    = sec - 3600*hour;    % remove the hours
minute = fix(sec/60);        % get number of minutes
sec    = sec - 60*minute;    % remove the minutes
second = sec;
end

% ------------ bootstrapping and other functions -------------------------

function rho = nancorr(a,b)
%calculate pearson correlation between a & b for nonnan values
keep = logical(~isnan(a) .* ~isnan(b));
rho = (mean((a(keep)-mean(a(keep))) .* (b(keep)-mean(b(keep)))))/(std(a(keep))*std(b(keep)));
end
