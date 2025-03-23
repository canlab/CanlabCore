classdef predictive_model
% Predictive model object for training and testing brain data
%   Example: Fit a linear discriminant to the Fisher iris data, and look at
%   a confusion matrix comparing the true and predicted Species values:
%       t = readtable('fisheriris.csv','format','%f%f%f%f%C');
%       d = fitcdiscr(t,'Species')
%       confusionmat(t.Species,predict(d,t))
%
%   See also fitcdiscr,
%   classreg.learning.classif.CompactClassificationDiscriminant.
%
% Methods:
% train
% test 
% crossval
% set_training_test
% report
% plot  % predictions vs. outcomes
% montage
% select_features
% bootstrap  (and reproducibility)
% error_analysis (misclassified images, within/between error)
%
% This is a stub - the start of an incomplete function. For now, use predict.m to
% train models, which returns stats structures. 
%
% Properties:
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


properties
    
    % Provenance
    training_data_description string = 'Description of training dataset';
    training_data_references  string; % Publication(s)
    model_description string = 'Description of model';
    model_references string          % Publication(s)
    
    % Before fitting
    Y (:, :) double                 % Outcome [198×1 double]
    id (:, 1) double                % Participant (or grouping) ID for crossval, bootstrapping, permutation
    model_type string               % Linear multivariate regression, SVM, etc.
    algorithm_name                  %: compatible with predict( )  'cv_lassopcr'
    function_call                   % : '@(xtrain, ytrain, xtest, cv_assignment) cv_lassopcr(xtrain, ytrain, xtest, cv_assignment, predfun_inputs{:})'
    function_handle                 %  [function_handle]
    error_type string               % : 'mse'
    nfolds = 5                      % string 'nfolds' or number of folds
    cvpartition struct              % [1×1 struct]
    teIdx (1, :) cell               % Cell: 1 x number of folds: {1×5 cell}
    trIdx (1, :) cell      
    model_options cell              % Inputs passed in during training: {'param_name', value, 'name', value}
    
    % After fitting
    
    % data and residuals
    yfit (:, 1) double 
    err (:, 1) logical 
    cverr (1, 1) double
    mse (1, 1) double
    rmse (1, 1) double
    meanabserr (1, 1) double
    pred_outcome_r (1, 1) double
    
    % Model parameters
    weight_obj fmri_data                    % statistic_image? if bootstrapped
    model_intercept (1, 1) double
    model_encoding_obj fmri_data            % fmri_data object with voxels x test participants (groups) 
    cv_weight_obj fmri_data
    cv_intercept (1, :)
                                    %other_output: {[223707×1 double]  [38.4627]  []}
                                    %other_output_descrip: 'Other output from algorithm - trained on all data (these depend on algorithm)'
                                    %     other_output_cv: {5×3 cell}
                                    %     other_output_cv_descrip: 'Other output from algorithm - for each CV fold'
    
    % bootstrapped
    % ****
    
end

methods
    
    
end


end
