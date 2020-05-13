% function [dim loss] = mlpcr_out_of_id_mse(kfolds,id,X,Y,...)
%
% Helper function for bayesopt optimization of MLPCR hyperparams. Wrapper
% for mlpcr_cv_pred that generates CV partitions on each call without
% fragmenting groups defined by "id" across folds. Returns mse computed by
% mlpcr_cv_pred.
%
% Input ::
%
%   kfolds          - kfolds to use when evaluating loss function
%
%   id              - n x 1 vector of labels elements that share the same
%                     label will be in either used for testing or training 
%                     sets but not both. e.g. if these are subject labels,
%                     then mlpcr_cv_pred will estimate out of subject
%                     generalization accuracy using kfold cross validation
%                     across subjects.
%
%   X               - n x p matrix of predictors
%
%   Y               - n x 1 vector of outcomes
%
%   optionalArgs    - any additional arguments will be passed on to mlpcr
%
% Output::
%
%   mse             - mse estimate for mlpcr_cv_pred
%

function [mse, pred] = mlpcr_out_of_id_mse(kfolds,id,X,Y,varargin)    
    % define CV fold test instances, without fragmenting grps defined by id
    % labels
    uniq_grp = unique(id);
    cc = cvpartition(length(uniq_grp),'Kfold',kfolds);
    C = zeros(length(Y),kfolds);
    for i = 1:kfolds
        test_set = uniq_grp(cc.test(i));

        C(:,i) = ismember(id,test_set);
    end

    % eval expected best point
    [pred,STATS] = mlpcr_cv_pred(C, X, Y, varargin{:});
    mse = STATS.mse;
end