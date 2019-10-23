% [weights, Intercept, lme, hp] = mlpcr_full(kfold,grps,bayesopt_opts,mlpcr_options)
%
% Performs MLPCR with hyperparameters estimated using bayesian 
% optimization of mean squared error (MSE). See mlpcr help for additional 
% details on usage and what to pass as mlpcr_options arguments.
%
% Hyperparameters can be specified explicitly or as type 
% optimizableVariable. If specified as optimizableVariable they will be
% optimized and the best result will be returned. Explicitly specified
% hyperparameters will simply be passed on to mlpcr as is.
%
% Note: mlpcr_full() relies on cross validation (CV) to estimate
% hyperparameters. Do not use on data where CV estimates are invalid, for
% instance on resampled datasets.
%
% Input ::
%
%   kfolds          - number of folds to use when evaluating objective
%                       function.
%
%   grps            - Optimization relies on kfold cross validation. Grps
%                       is an n x 1 vector of labels (e.g. subjects) which
%                       are not fragmented across cv folds. Elements with
%                       the same label will either be in a test fold or a
%                       training fold, but not both. Thus, in a
%                       multisubject study providing subject labels will 
%                       result in hyperparameter optimization for 
%                       out-of-subject performance.
%
%   bayesopts_opts  - Options to pass to bayesopt() call. "help bayesopts"
%                       for details. If no acquisition function is
%                       specified expected-improvement-plus is used by
%                       default. Pass empty cell array to use defaults.
%                     Note: if running with {'UseParallel',true} invoke
%                       parpool manually and then run multithreadWorkers()
%                       before invoking mlpcr_full for best performance.
%
%   mlpcr_options   - Options to pass to mlpcr. Substitute
%                       optimizableVariable objects for any hyperparameter
%                       you want optmized.
%
% Output ::
%
%   weights         - Cell array of weights returned by mlpcr (using 
%                       optimal hyperparameters). See "help mlpcr" for 
%                       details.
%
%   Intercept       - Cell array of intercepts returned by mlpcr (using 
%                       optimal hyperparameters. See "help mlpcr" for
%                       details.
%
%   lme             - LinearMixedModel object returned by fitlme() when
%                       fitting PCA components to data after optimizing
%                       hyperparameters.
%
%   hp              - BayesianOptimization object returned by bayesopt(). 
%                       Useful for inspecting optimization results, 
%                       resuming optimization for additional iterations,
%                       or getting estimate of objective function 
%            			at optimal hyperparameters.
%                     Note: MSE estimate for final model is stored in
%                       hp.MinObjective.
%
% Dependencies:
%   R2016b machine learning toolbox (required by this script for bayesopt)
%   mlpcr_out_of_id_mse     (required by this script)
%   mlpcr_cv_pred           (required by mlpcr_out_of_id_mse)
%   mlpcr                   (required by mlpcr_cv_pred)
%   mlpca                   (required by mlpcr)
%   get_nested_var_comps    (required by mlpcr_cv_pred)
%   get_cntrng_mats         (required by mlpca and get_nested_var_comps)
%
% Writen by Bogdan Petre (Feb 20, 2018)

function [weights, Intercept, lme, hp] = mlpcr_full(kfold,grps,bayesopt_opts,X,Y,varargin)
    weights = {};
    Intercept = {};
    lme = {};

    optVars = [];
    mlpcr_arg = {X,Y};
    for i = 1:length(varargin)
        [new_arg, new_optVars] = extractOptVars(varargin{i});
        mlpcr_arg = [mlpcr_arg, new_arg];
        optVars = [optVars, new_optVars];
    end
    
    % construct obj fxn invocation
    execstr = construct_obj_fxn(1, 'mlpcr_arg', mlpcr_arg{:});
    eval(['objfxn = @(x1)(mlpcr_out_of_id_mse(kfold,grps,', execstr, '));']);
    
    % default acquisition function
    AcFxn = 'expected-improvement-plus';
    % If acquisition function was specified, use that instead, and remove
    % from bayesopt_opt (will be passed explicitly)
    for i = 1:length(bayesopt_opts)
        if ischar(bayesopt_opts)
            switch bayesopt_opts{i}
                case 'AcquisitionFunctionName'
                    AcFxn = bayesopt_opts{i+1};
                    bayesopt_opts = bayesopt_opts(1:i-1,i+3:end);
            end
        end
    end
    hp = bayesopt(objfxn,optVars,'AcquisitionFunctionName',AcFxn,bayesopt_opts{:});
    
    execstr = strrep(execstr,'(x1(','(hp.XAtMinEstimatedObjective(');
    
    eval(['[weights,Intercept,lme] = mlpcr(', execstr, ');']);
end
   
% Parses arguments from mlpcr_full
% fixedArgs will be identical to varargin, except that optimizableVariable
% types will be replaced with nans, and the optimizableVariable will be
% copied to optVars.
function [fixedArgs, optVars] = extractOptVars(varargin)
    optVars = [];
    fixedArgs = {};
    for i = 1:length(varargin)
        % nans are used as flags for optimizableVariables internally, so make sure there aren't any preexisting.
        if isnumeric(varargin{i})
            if isnan(varargin{i}) 
                error('Found ''nan'' in argument list. This is not supported.');
            end
        end
        switch class(varargin{i})
            case 'optimizableVariable'
                optVars = [optVars, varargin{i}];
                fixedArgs = [fixedArgs, {nan}];
            case 'cell'
                [subFixedArgs, subOptVars] = extractOptVars(varargin{i}{:});
                optVars = [optVars, subOptVars];
                fixedArgs = [fixedArgs, {subFixedArgs}];
            otherwise
                fixedArgs = [fixedArgs, varargin(i)];
        end
    end
end

% ov_idx - optimization variable index
% fv_idx - fixed variable index
function [execstr, ov_idx] = construct_obj_fxn(ov_idx, name, varargin)
    execstr = [];
    for i = 1:length(varargin)
        if ~isa(varargin{i},'cell')
            if isnan(varargin{i})
                execstr = [execstr, 'table2array(x1(1,' int2str(ov_idx) ')),'];
                ov_idx = ov_idx + 1;
            else
                execstr = [execstr, name, '{', int2str(i), '},'];
            end
        else
            % Note: ov_idx's carry through
            [new_execstr, ov_idx] = construct_obj_fxn(ov_idx, [name, '{' int2str(i) '}'], varargin{i}{:});
            
            execstr = [execstr, '{', new_execstr, '},']; 
        end
    end
    execstr = execstr(1:end-1); % drop trailing comma
end

