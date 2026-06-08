%[text] # Multivariate decoding — Part 4: cross-classification
%[text] Part 4 of the multivariate-decoding series (Parts 1–3 build up single-task decoding; Part 5 covers algorithms & tuning). See the matching .md for the full narrative and figures.
%[text] Cross-classification tests representational overlap: train a classifier on one task, freeze it, apply it to another. Here we ask the Woo et al. (2014) question — does a brain pattern that separates **physical pain** (Hot vs Warm) also separate **social rejection** (Rejecter vs Friend)? Above-chance off-diagonal accuracy means a shared representation; how much below the within-system diagonal tells us how separable they are.
%%
%[text] ## 1. Load both tasks (same voxel space)
hw = load_image_set('DPSP_hotwarm', 'noverbose');        % Hot (+1) vs Warm (-1)
rf = load_image_set('DPSP_rejectorfriend', 'noverbose'); % Rejecter (+1) vs Friend (-1)
Xhw = double(hw.dat'); Yhw = hw.Y; idhw = hw.metadata_table.subj_id;
Xrf = double(rf.dat'); Yrf = rf.Y; idrf = rf.metadata_table.subj_id;
assert(size(hw.dat,1) == size(rf.dat,1), 'datasets must share voxel space');
%%
%[text] ## 2. Within-system baselines (cross-validated)
%[text] crossval auto-computes the within-subject forced-choice accuracy when groups are supplied — the diagonal of the generalization matrix.
cv = cv_splitter.stratified_group_kfold(5);
pm_hw = crossval(predictive_model('algorithm','svm','task','classification'), Xhw, Yhw, 'groups', idhw, 'cv', cv);
pm_rf = crossval(predictive_model('algorithm','svm','task','classification'), Xrf, Yrf, 'groups', idrf, 'cv', cv);
within_pain   = pm_hw.error_metrics.crossval_accuracy_within.value
within_reject = pm_rf.error_metrics.crossval_accuracy_within.value
%%
%[text] ## 3. Cross-classification: train one, freeze, apply to the other
%[text] Because the two datasets share voxel space and the test images were never in training, a plain full-sample fit is the right training step (no CV needed across datasets).
pm_pain = fit(predictive_model('algorithm','svm','task','classification'), Xhw, Yhw);  % pain pattern
pm_rej  = fit(predictive_model('algorithm','svm','task','classification'), Xrf, Yrf);  % rejection pattern
[~, score_rf] = predict(pm_pain, Xrf);   % pain pattern applied to rejection
[~, score_hw] = predict(pm_rej,  Xhw);   % rejection pattern applied to pain
%%
%[text] ## 4. Score with forced choice, NOT raw accuracy
%[text] The SVM bias term does not transfer across datasets, so raw accuracy at the native threshold misleads — only the score *ranking* transfers. Score with within-subject forced choice (each person their own control).
fc = @(scores, Y, id) local_forced_choice(scores, Y, id);
cross_pain_to_rej = fc(score_rf(:,end), Yrf, idrf)
cross_rej_to_pain = fc(score_hw(:,end), Yhw, idhw)
%%
%[text] ## 5. The 2x2 generalization matrix
%[text] Rows = training task, columns = testing task, cells = within-subject forced-choice accuracy. Diagonal = within-system (CV); off-diagonal = cross-applied frozen model. Off-diagonal well above 50% but below the diagonal is the *shared-but-separable* signature.
G = [within_pain, cross_pain_to_rej; cross_rej_to_pain, within_reject];
T = array2table(G, 'RowNames', {'trainPain','trainReject'}, 'VariableNames', {'testPain','testReject'})
%%
%[text] For the common case, xval\_cross\_classify packages this whole flow and returns a predictive\_model with the results under pm.cross\_classify.

function acc = local_forced_choice(scores, Y, id)
% Fraction of subjects whose +1 image scores above their -1 image.
    [~, ~, g] = unique(id(:), 'stable');
    u = unique(g); hit = nan(numel(u), 1);
    for i = 1:numel(u)
        sp = scores(g==u(i) & Y== 1);
        sn = scores(g==u(i) & Y==-1);
        if ~isempty(sp) && ~isempty(sn), hit(i) = mean(sp) > mean(sn); end
    end
    acc = 100 * mean(hit, 'omitnan');
end

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline"}
%---
