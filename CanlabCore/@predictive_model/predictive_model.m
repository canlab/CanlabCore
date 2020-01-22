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
% error_analysis
%
% This is a stub - the start of an incomplete function. For now, use predict.m to
% train models, which returns stats structures. 

properties
    
    % Provenance
    training_data_description 
    training_data_references  % Publication(s)
    model_description
    model_references          % Publication(s)
    
    % Before fitting
    Y                  % Outcome [198×1 double]
    model_type         % Linear multivariate regression, SVM, etc.
    algorithm_name     %: 'cv_lassopcr'
    function_call      % : '@(xtrain, ytrain, xtest, cv_assignment) cv_lassopcr(xtrain, ytrain, xtest, cv_assignment, predfun_inputs{:})'
    function_handle: [function_handle]
    error_type: 'mse'
    nfolds: 'nfolds'
    cvpartition: [1×1 struct]
    teIdx: {1×5 cell}
    trIdx: {1×5 cell}
    model_options    % {'param_name', value, 'name', value}
    % After fitting
    
    % data and residuals
    yfit  %
    err: [198×1 double]
    cverr: 1.6876e+03
    mse: 1.6876e+03
    rmse: 41.0808
    meanabserr: 34.0429
    pred_outcome_r: 0.6412
    
    % Model parameters
    weight_obj: [1×1 fmri_data]
    model_intercept
    other_output: {[223707×1 double]  [38.4627]  []}
    other_output_descrip: 'Other output from algorithm - trained on all data (these depend on algorithm)'
    other_output_cv: {5×3 cell}
    other_output_cv_descrip: 'Other output from algorithm - for each CV fold'
    
    % bootstrapped
    % ****
    
end

methods
    
    
end


end
