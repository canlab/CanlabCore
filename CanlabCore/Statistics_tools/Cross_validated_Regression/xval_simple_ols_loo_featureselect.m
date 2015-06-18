function pred_value = xval_simple_ols_loo(X, Y, pthreshold)
% Check: a very simple leave-one-out cross-validated regression
% pred_value = xval_simple_ols_loo(X, Y, p-value selection for univariate feature selection)
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

 
create_figure('test', 1, 2);

fprintf('Fold: %03d', 0);
for i = 1:N
    
    fprintf('\b\b\b%3.0f ', i);
    
    % select training data
    Xi = X; Xi(i, :) = []; 
    Yi = Y; Yi(i) = [];
    
    % select features based on univariate correlations
    [r, p, Tstat] = correlation_fast_series(Xi, Yi);
    wh_features = p <= pthreshold;
    nfeatures(i) = sum(wh_features);
    if nfeatures(i) == 0
        disp('Warning: no features pass threshold');
        wh_features = p <= prctile(p, 10);
    else
        Xi = Xi(:, wh_features);
    end
    
    Xi = [ones(N-1, 1) Xi]; % add intercept
    
    % Make prediction
    b = pinv(Xi) * Yi;
    pred_value(i, 1) = [1 X(i, wh_features)] * b;
 
    create_figure('test', 1, 2, 1); subplot(1, 2, 1); plot(Xi*b, Yi, 'ko'); 
    plot(pred_value(i, 1), Y(i), 'ro', 'MarkerFaceColor', 'r');
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


% figure; hold on; 
% for i = 1:N
%     plot(intercept_vals(i), pred_value(i), 'ko', 'MarkerFaceColor', cm(i, :));
% end

% 
% figure; hold on; 
%     plot(pred_value(whorder(1:24)), Y(whorder(1:24)), 'bo', 'MarkerFaceColor', 'b');
% plot(pred_value(whorder(25:end)), Y(whorder(25:end)), 'ro', 'MarkerFaceColor', 'r');
