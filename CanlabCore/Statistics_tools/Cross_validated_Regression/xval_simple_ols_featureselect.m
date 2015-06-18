function pred_value = xval_simple_ols_featureselect(X, Y, pthreshold, holdout_method)
% Check: a very simple leave-one-out cross-validated regression
% pred_value = xval_simple_ols_loo(X, Y, p-value selection for univariate feature selection, holdout_method)
%
% pred_value: cross-validated predictions of outcome data
% X: n x variables matrix of predictors
% Y: n x 1 vector of outcomes
%
% Tor Wager, June 2010
%
% Go to any LASSO output directory and run this:
% 
% maskInfo = iimg_read_img(fullfile(pwd, 'mask.img'), 2);
% dat{1} =  iimg_get_data(maskInfo, anticimages);
% pred_value = xval_simple_ols_loo(X, Y)


[N, k] = size(X);
pred_value = zeros(N, 1);

holdout_set = nested_select_holdout_set;
 
create_figure('test', 1, 2);

fprintf('Fold: %03d', 0);
for wh_fold = 1:length(holdout_set) 
    
    fprintf('\b\b\b%3.0f ', wh_fold);
    
    % select training data
    Xi = X; 
    Yi = Y; 
    
    Yi(holdout_set{wh_fold}) = [];
    Xi(holdout_set{wh_fold}, :) = [];              % leave out the missing observation(s)

    Xtest = X(holdout_set{wh_fold}, :);
    Ytest = Y(holdout_set{wh_fold}, :);          % only used later when we test prediction

    ntrain = size(Xi, 1);
    nholdout = size(Xtest, 1);
        
    % select features based on univariate correlations
    [r, p, Tstat] = correlation_fast_series(Xi, Yi);
    wh_features = p <= pthreshold;
    nfeatures(wh_fold) = sum(wh_features);
    if nfeatures(wh_fold) == 0
        disp('Warning: no features pass threshold');
        wh_features = p <= prctile(p, 10);
    end

    Xi = Xi(:, wh_features);
    
    Xi = [ones(ntrain, 1) Xi]; % add intercept
    
    % Make prediction
    b = pinv(Xi) * Yi;
    pred_value(holdout_set{wh_fold}, 1) = [ones(nholdout, 1) Xtest(:, wh_features)] * b;
 
    create_figure('test', 1, 2, 1); subplot(1, 2, 1); plot(Xi*b, Yi, 'ko'); 
    plot(pred_value(holdout_set{wh_fold}, 1), Ytest, 'ro', 'MarkerFaceColor', 'r');
    subplot(1, 2, 2); 
   title('Black circles = training set; Red = holdout obs');
    drawnow;
    
end

cm = colormap(jet(N));
figure; hold on; 
for i = 1:N
    plot(pred_value(i), Y(i), 'ko', 'MarkerFaceColor', cm(i, :));
end
xlabel('Predicted outcome (xval)'); ylabel('Outcome');
title('Color = order in data series');




function holdout_set = nested_select_holdout_set
        % purpose: return holdout_set variable

        nobs_s = length(Y);

        switch lower(holdout_method)
            case 'loo'
                holdout_set = cell(1, nobs_s);
                for i = 1:nobs_s, holdout_set{i} = i; end

            case 'l4o_covbalanced'
                disp('Selecting holdout sets: Balancing on covariates and outcome, and also trying to ensure that each obs is selected equally often.');
                holdout_proportion = 4 ./ nobs_s; % prop for leave-4-out
                nfolds = 40;
                if isempty(cov_val)
                    wh_holdout = xval_select_holdout_set(Y, [], nfolds, holdout_proportion, verbose);
                else
                    wh_holdout = xval_select_holdout_set(Y, cov_val, nfolds, holdout_proportion, verbose);
                end
                holdout_set = cell(1, nfolds);
                for k = 1:nfolds
                    holdout_set{k} = wh_holdout(:, k);
                end

            case 'categorical_covs'
                
                holdout_set = xval_select_holdout_set_categoricalcovs(cov_val);
                
            case 'balanced4'
                nfolds = 4;
                holdout_set = cell(1, nfolds);
                [ys, indx] = sort(Y);
                for k = 1:nfolds

                    holdout_set{k} = indx(k:nfolds:end);
                    if isempty(holdout_set{k}), error('Holdout set construction error: Dataset too small?'); end

                end

            otherwise error('Unknown holdout method.  See help.');
        end
end
    
end % main function
        
