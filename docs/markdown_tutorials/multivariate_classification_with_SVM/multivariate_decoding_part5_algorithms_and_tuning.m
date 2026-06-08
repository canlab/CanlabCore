%[text] # Multivariate decoding — Part 5: algorithms, tuning, and inference
%[text] Part 5 of the multivariate-decoding series. One dataset, many estimators: binary and multiclass (ECOC) classification and regression (SVR/lasso/ridge), compared over identical folds, with grid\_search tuning and stability\_selection inference. See the matching .md for the full narrative, figures, and a "when to use which" guide.
%%
%[text] ## 1. Setup
hw = load_image_set('DPSP_hotwarm', 'noverbose');
rf = load_image_set('DPSP_rejectorfriend', 'noverbose');
X = double(hw.dat'); Y = hw.Y; id = hw.metadata_table.subj_id;
%%
%[text] ## 2. Binary classification — compare algorithms on identical folds
%[text] Cross-validate several classifiers over the same splitter so the comparison is fold-matched. (predict\_test\_suite does this for you in one call.)
cv = cv_splitter.stratified_group_kfold(5);
for a = {'svm','linear_svm','logistic','lda'}
    pm = crossval(predictive_model('algorithm', a{1}, 'task', 'classification'), ...
                  X, Y, 'groups', id, 'cv', cv, 'scoring', 'balanced_accuracy');
    fprintf('%-12s cv bal-acc = %.3f\n', a{1}, pm.error_metrics.balanced_accuracy.value);
end
%%
%[text] ## 3. Multiclass classification with ECOC
%[text] For >2 classes, ecoc reduces the K-class problem to binary SVMs and decodes the votes. Stack the two DPSP tasks into a 4-class problem (Hot/Warm/Rejecter/Friend).
Xall  = [double(hw.dat'); double(rf.dat')];
Ycls  = [(hw.Y==1)*1 + (hw.Y==-1)*2; (rf.Y==1)*3 + (rf.Y==-1)*4];
idall = [hw.metadata_table.subj_id; rf.metadata_table.subj_id];
pm_ecoc = crossval(predictive_model('algorithm','ecoc','task','classification'), Xall, Ycls, ...
                   'groups', idall, 'cv', cv_splitter.group_kfold(5), 'scoring', 'accuracy');
ecoc_accuracy = 100 * pm_ecoc.error_metrics.accuracy.value
confusionchart(pm_ecoc);
%%
%[text] ## 4. Regression — compare SVR / lasso / ridge
%[text] DPSP ships binary labels; here we treat the signed label as a continuous target purely to compare regression algorithms on the same signal (in a real study, swap in a continuous outcome). r asks "does the prediction track the outcome", R^2 "on the right scale".
for a = {'svr','linear_svr','lasso','ridge'}
    pm = crossval(predictive_model('algorithm', a{1}, 'task', 'regression'), ...
                  X, Y, 'groups', id, 'cv', cv_splitter.group_kfold(5));
    fprintf('%-12s r = %.3f  R^2 = %.3f\n', a{1}, ...
        predictive_model.metric_value(pm.error_metrics, 'prediction_outcome_r'), ...
        pm.error_metrics.r2.value);
end
%[text] Gaussian-process regression (fitrgp) cannot take ~200k voxels — reduce first with a PCA step in an @pipeline (which refits the PCA per fold): pipeline({{'pca','k',30}}, predictive\_model('algorithm','gp','task','regression')).
%%
%[text] ## 5. Hyperparameter tuning with grid_search (nested CV for the honest estimate)
%[text] grid\_search cross-validates a parameter grid and refits at the winner. The score at the chosen point is optimistically biased; nest the search (outer folds + inner grid\_search per fold) for an unbiased estimate.
g.BoxConstraint = [0.01 0.1 1 10];
pm_gs = predictive_model('algorithm','svm','task','classification','cv',cv_splitter.stratified_group_kfold(5));
pm_gs = grid_search(pm_gs, X, Y, g, 'groups', id, 'verbose', false);
best_C = pm_gs.diagnostics.grid_search.best_params{2}
%%
%[text] ## 6. High-dimensional inference with stability_selection
%[text] For wide regularised models the bootstrap z/p collapses; stability selection counts how often each voxel lands in the top-k by |weight| across resamples — the robust "where is the signal" map.
pm_ss = stability_selection(predictive_model('algorithm','linear_svm','task','classification'), ...
            X, Y, 'nboot', 100, 'k', 2000, 'threshold', 0.6, 'groups', id, 'verbose', false);
n_stable = pm_ss.diagnostics.stability_selection.n_stable
freq = hw; freq.dat = pm_ss.diagnostics.stability_selection.selection_freq;
create_figure('stability'); axis off; montage(freq);

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline"}
%---
