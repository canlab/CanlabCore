function [CVMdl, statstable, profile_data, mprop, perm_acc_stats] = knn_classify_cv_by_profile(X, y, analysis_name, condition_labels)

% Parameters

nperms = 100;
nneighbors = 10;
dashes = '----------------------------------------------';

% Z-score each image so that it is mean 0 and std 1 across regions
% This way, we are testing whether the distributions across regions differ.
profile_data = scale(X')';

Mdl = fitcknn(profile_data, y, 'NumNeighbors', nneighbors); % better for correlated variables, as we have
% 'Distance', 'mahalanobis', does not work as well.

% Train cross-validated model
CVMdl = crossval(Mdl, 'KFold', 5);

% Predict class - CV
[predLabel_cv, class_scores] = kfoldPredict(CVMdl);

cvacc = 100 * (1 - kfoldLoss(CVMdl));

[m,aprime,corr,far,missclass, mprop, statstable] = confusion_matrix(y, predLabel_cv);

% Display output
% ---------------------------------------------------------------
disp(analysis_name)
disp(dashes)

% Permutation: Test chance classification and get CI
% ---------------------------------------------------------------
cvacc_p = NaN * zeros(nperms, 1);
fprintf('Permuting %d times: %5.0f', nperms, 0);

for j = 1:nperms
    
    fprintf('\b\b\b\b\b%5.0f', j);
    
    y_p = y(randperm(length(y)));
    Mdl_p = fitcknn(profile_data, y_p, 'NumNeighbors', nneighbors);
    CVMdl_p = crossval(Mdl_p, 'KFold', 5);
    cvacc_p(j, 1) = 100 * (1 - kfoldLoss(CVMdl_p));
    
end

fprintf('\n')
mean_perm_acc = mean(cvacc_p);
ci_95_perm_acc = prctile(cvacc_p, [2.5 97.5]);
acc_perm_p = sum(cvacc_p >= cvacc) ./ nperms; % P-value for overall acc, one-tailed

if acc_perm_p == 0 % then use normal approx.
    acc_perm_p = 1 - normcdf((cvacc - mean(cvacc_p)) ./ std(cvacc_p));
end

perm_acc_stats.names = {'cv_acc' 'mean_perm_acc' 'ci_95_perm_acc' 'acc_perm_p'};
perm_acc_stats.values = [cvacc mean_perm_acc ci_95_perm_acc acc_perm_p];

% Display output
% ---------------------------------------------------------------
disp(analysis_name)
disp(dashes)

fprintf('5-fold CV Accuracy\nOverall: %3.3f%%\tChance (perm) is: %3.3f%%, CI: [%3.3f %3.3f] P = %3.6f\n', cvacc, mean_perm_acc, ci_95_perm_acc(1), ci_95_perm_acc(2), acc_perm_p);
statstable.Properties.RowNames = condition_labels;
disp(statstable)

fprintf('\n')

disp('5-fold CV confusion matrix. Rows: true, Columns: Predicted')
print_matrix(mprop, condition_labels, condition_labels);

disp(dashes)
fprintf('\n')

end
