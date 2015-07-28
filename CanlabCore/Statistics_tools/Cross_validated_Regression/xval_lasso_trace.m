function STATS = xval_lasso_trace(Y, data, whfolds, names)
% STATS = xval_lasso_trace(Y, data, whfolds, names)
%
% y = outcome
% data = obs x variables predictor matrix
% whfolds = integer vector of test set membership for each obs, obs x 1
% names = cell array of names for each predictor variable
%
% Tor Wager, July 2011

colors = scn_standard_colors(size(data, 2));

lassopath = '/Users/tor/Documents/matlab_code_external/machine_learning/lasso_rocha/lasso';
addpath(lassopath)

%% Basic cross-validated prediction

% STATS = xval_regression_multisubject('ols', {Y}, {data}, 'holdout_method', 'balanced4');
%
% STATS = xval_regression_multisubject('ols', {Y}, {data}, 'holdout_method', whfolds);

STATS = xval_regression_multisubject('ols', {Y}, {data}, 'holdout_method', whfolds);

fprintf('r = %3.2f\tr^2 = %3.2f\tpe = %3.4f\n', STATS.r_each_subject, STATS.r_squared, STATS.pred_err);

%% Now leave one var out at a time

for i = 1:size(data, 2)
    
    data2 = data;
    data2(:, i) = [];
    
    STATS = xval_regression_multisubject('ols', {Y}, {data2}, 'holdout_method', whfolds, 'noverbose');
    
    xval_pe(i, 1) = STATS.pred_err;
    
end

%% Print Leave out one var output

disp('Leave-one-variable out prediction error, sorted from highest to lowest beta weight.')

STATS = xval_regression_multisubject('ols', {Y}, {data}, 'holdout_method', whfolds, 'noverbose');

[b, wh] = sort(abs(STATS.subjbetas{1}(2:end)), 'descend');
names2 = names(wh);
xval_pe2 = xval_pe(wh);

fprintf('\n%s\t%s\t%s\n', 'Name', 'xvalBeta', 'LeaveOutVarPE');

for i = 1:length(wh)
    
    fprintf('%s\t%3.2f\t%3.4f\n', names2{i}, b(i), xval_pe2(i));
    
end

create_figure('all vars', 2, 1);

w = STATS.mean_vox_weights;  % STATS.subjbetas{1}(2:end)
hh = barh(w);
set(hh, 'FaceColor', [.5 .5 .5]);
set(gca, 'YDir', 'Reverse', 'YTick', 1:length(names), 'YTickLabel', names, 'YLim', [0 length(names)+1], 'FontSize', 24)
xlabel('Beta weights');
title('OLS');

% prediction-outcome correlation with all vars
r_all = STATS.r_each_subject;

%% LASSO

out = lasso_rocha(Y, data);

[k, v] = size(out.beta); % steps, vars

subplot(2, 1, 2);
w = out.beta(end, :)';
hh = barh(w);
set(hh, 'FaceColor', [.5 .5 .5]);
set(gca, 'YDir', 'Reverse', 'YTick', 1:length(names), 'YTickLabel', names, 'YLim', [0 length(names)+1], 'FontSize', 24)
xlabel('Beta weights');
title('LASSO - NO PENALTY');

%% LASSO trace plot 

create_figure('Lassopath', 3, 1); 

for i = 1:size(out.beta, 2)
    plot(0:(k-1), out.beta(:, i), 'o-', 'LineWidth', 3, 'color', colors{i});
end

legend(names)
ylabel('Regression coefficient');
drawnow

xval_pe = NaN;

for i = 2:k % omit first step with all zero betas

    STATS = xval_regression_multisubject('lasso', {Y}, {data}, 'holdout_method', whfolds, 'noverbose', 'lassopath', out.lambda(i));
    
    % xval_regression_multisubject will fit the OLS model after selecting variables. 
    xval_pe(i, 1) = STATS.pred_err;

    r(i, 1) = STATS.r_each_subject; % prediction-outcome correlation
    
    % for BIC
    v(i-1, 1) = STATS.var_full;
    f(i-1, 1) = i - 1;            % number of free params
    n(i-1, 1) = length(Y);        % this should be independent; or else estimate df?...
    
    % alternatives:  all produces same results, 3/5/13
%     dat = fmri_data;
%     dat.Y = Y; dat.dat = data';
%     [cverr, stats, optout] = predict(dat, 'algorithm_name', 'cv_lassopcr', 'nfolds', whfolds, 'lasso_num', i-1, 'nopcr');
%     [cverr, stats, optout] = predict(dat, 'algorithm_name', 'cv_lassopcrmatlab', 'nfolds', whfolds, 'lasso_num', i-1, 'Alpha', 1, 'nopcr');
end

subplot(3, 1, 2);
plot(0:(k-1), xval_pe, 'ko-', 'LineWidth', 3);

hh = plot_horizontal_line(STATS.pred_err_null);
set(hh, 'LineStyle', ':');
xlabel('Number of variables in model')
ylabel('Error');
title('Cross-validated prediction error');

subplot(3, 1, 3);
plot(0:(k-1), r.^2, 'ko-', 'LineWidth', 3);

hh = plot_horizontal_line(r_all.^2);
set(hh, 'LineStyle', ':');
xlabel('Number of variables in model')
ylabel('r^2');
title('Variance explained');

% BIC: This assumes observations are independent...

bic = n.*log(v) + f.*log(n);
[mm, i] = min(bic);

print_matrix([(0:k-1)' r out.lambda [nan; bic]], {'Vars in model' 'pred-Y r' 'lasso shrinkage lambda', 'bic'});

% Do one more time for output
STATS = xval_regression_multisubject('lasso', {Y}, {data}, 'holdout_method', whfolds, 'noverbose',  'lassopath', out.lambda(i+1));




end % function
