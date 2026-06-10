%[text] # Multivariate decoding — Part 4: cross-classification
%[text] Part 4 of the multivariate-decoding series (Parts 1–3 build up single-task decoding; Part 5 covers algorithms & tuning). See the matching .md for the full narrative and figures.
%[text] Cross-classification tests representational overlap: train a classifier on one task and apply it to another. Here we ask the Woo et al. (2014) question — does a brain pattern that separates **physical pain** (Hot vs Warm) also separate **social rejection** (Rejecter vs Friend)? Above-chance off-diagonal accuracy means a shared representation; how far below the within-system diagonal tells us how separable they are. Because Hot/Warm and Rejecter/Friend are the **same subjects**, we use a **shared cross-validation fold structure** so the within- and cross-condition numbers are directly comparable.
%%
%[text] ## 1. Load both tasks (same voxel space)
hw = load_image_set('DPSP_hotwarm', 'noverbose');        % Hot (+1) vs Warm (-1)
rf = load_image_set('DPSP_rejectorfriend', 'noverbose'); % Rejecter (+1) vs Friend (-1)
Xhw = double(hw.dat'); Yhw = hw.Y; idhw = hw.metadata_table.subj_id;
Xrf = double(rf.dat'); Yrf = rf.Y; idrf = rf.metadata_table.subj_id;
assert(size(hw.dat,1) == size(rf.dat,1), 'datasets must share voxel space');
%%
%[text] ## 2. Within-system baselines (quick, separate CV per task)
%[text] crossval auto-computes within-subject forced-choice accuracy when groups are supplied. These quick baselines use a separate CV per task; §3 re-derives the within numbers under a **shared** fold structure so every cell of the matrix is computed the same way.
cv = cv_splitter.stratified_group_kfold(5);
pm_hw = crossval(predictive_model('algorithm','svm','task','classification'), Xhw, Yhw, 'groups', idhw, 'cv', cv);
pm_rf = crossval(predictive_model('algorithm','svm','task','classification'), Xrf, Yrf, 'groups', idrf, 'cv', cv);
fprintf('within (separate CV): pain %.1f%%, reject %.1f%%\n', ...
    pm_hw.error_metrics.crossval_accuracy_within.value, ...
    pm_rf.error_metrics.crossval_accuracy_within.value);
%%
%[text] ## 3. Cross-classification with shared cross-validation folds
%[text] A full-sample pain model applied to all rejection images isn't *leakage* (the test images were never trained on), but it's a **bad comparison**: it is trained on all 59 subjects vs the diagonal's ~47, and — same subjects across tasks — every test subject's identity was seen in training. Instead we share one fold structure: for each fold, train on condition A's training subjects, then apply that fold's model to the held-out subjects' condition-A data (within) AND their condition-B data (cross). Every subject is scored once, by a model that saw none of their data in either task. Same folds, same training sizes, same held-out subjects → comparable cells.
[~, ~, ghw] = unique(idhw, 'stable');
sp = cv.split(Xhw, Yhw, ghw);                          % folds defined on subjects
n = size(Xhw, 1);
sc_hw_within = nan(n,1); sc_pain2rej = nan(n,1);       % pain-trained scores
sc_rf_within = nan(n,1); sc_rej2pain = nan(n,1);       % reject-trained scores
for k = 1:numel(sp)
    tr = sp(k).trIdx; te = sp(k).teIdx;
    test_subj = unique(idhw(te));  tr_subj = unique(idhw(tr));
    rf_te = ismember(idrf, test_subj);                 % held-out subjects' rejection rows
    rf_tr = ismember(idrf, tr_subj);                   % training subjects' rejection rows
    % train on pain -> score held-out pain (within) + rejection (cross)
    mp = fit(predictive_model('algorithm','svm','task','classification'), Xhw(tr,:), Yhw(tr));
    [~, s] = predict(mp, Xhw(te,:));     sc_hw_within(te)   = s(:,end);
    [~, s] = predict(mp, Xrf(rf_te,:));  sc_pain2rej(rf_te) = s(:,end);
    % train on rejection (same fold's training subjects) -> rejection (within) + pain (cross)
    mr = fit(predictive_model('algorithm','svm','task','classification'), Xrf(rf_tr,:), Yrf(rf_tr));
    [~, s] = predict(mr, Xrf(rf_te,:));  sc_rf_within(rf_te) = s(:,end);
    [~, s] = predict(mr, Xhw(te,:));     sc_rej2pain(te)     = s(:,end);
end
%%
%[text] ## A note on the bias term (intercept)
%[text] |predict(pm, Xnew)| applies the model's **intercept** as well as its weights — in-sample, on a held-out test set, and in every CV fold (each fold uses its own training intercept). So the scores are full decision values |w·x + b|. The intercept is calibrated to the **training** task's score distribution, which is exactly why raw accuracy doesn't transfer (below). If you ever apply a pattern **by hand** (|Xnew\*w|, a dot-product / |apply_mask| / |image_similarity|), note that |pm.weights.w| (and the |statistic_image| from |weight_map_object|) is the **slope only, no intercept** — correct for a weight *map*, but your score will miss the offset. The intercept lives in |pm.ml_model.Bias| (SVM / fitclinear) or |pm.ml_model.intercept| (PCR / lassoPCR). Easiest: just call |predict(pm, Xnew)|.
%%
%[text] ## 4. Score with forced choice, NOT raw accuracy
%[text] The bias is calibrated to the training task, so the decision boundary lands in the wrong place on the testing task and raw accuracy misleads — only the score *ranking* transfers. Score with within-subject forced choice (each person their own control). The helper returns the per-subject hits too, for the confidence intervals in §5b.
within_pain       = forced_choice(sc_hw_within, Yhw, idhw)
cross_pain_to_rej = forced_choice(sc_pain2rej,  Yrf, idrf)
cross_rej_to_pain = forced_choice(sc_rej2pain,  Yhw, idhw)
within_reject     = forced_choice(sc_rf_within, Yrf, idrf)
%%
%[text] ## 5. The 2x2 generalization matrix
%[text] Rows = training task, columns = testing task, cells = within-subject forced-choice accuracy — all cross-validated under the shared folds, so directly comparable. Off-diagonal above 50% but below the diagonal is the *shared-but-separable* signature.
G = [within_pain, cross_pain_to_rej; cross_rej_to_pain, within_reject];
T = array2table(G, 'RowNames', {'trainPain','trainReject'}, 'VariableNames', {'testPain','testReject'})
%%
%[text] ## 5b. Inference on the generalization error
%[text] Forced-choice accuracy is a **binomial** proportion — each of the n=59 subjects is one correct/incorrect trial — so an exact (Clopper–Pearson) 95% CI comes from the per-subject hits via |binofit|. The error bars let you read each cell against chance and against each other.
labels = {'Pain\rightarrowPain','Pain\rightarrowReject','Reject\rightarrowPain','Reject\rightarrowReject'};
[~, h1] = forced_choice(sc_hw_within, Yhw, idhw);
[~, h2] = forced_choice(sc_pain2rej,  Yrf, idrf);
[~, h3] = forced_choice(sc_rej2pain,  Yhw, idhw);
[~, h4] = forced_choice(sc_rf_within, Yrf, idrf);
H = {h1,h2,h3,h4}; vals = nan(1,4); lo = nan(1,4); hi = nan(1,4);
for j = 1:4
    h = H{j}; h = h(~isnan(h));
    [p, ci] = binofit(sum(h), numel(h));
    vals(j) = 100*p; lo(j) = 100*ci(1); hi(j) = 100*ci(2);
end
create_figure('generalization error');
bar(1:4, vals, 0.6, 'FaceColor', [.45 .6 .8]); hold on
errorbar(1:4, vals, vals-lo, hi-vals, 'k', 'linestyle','none', 'LineWidth',1.5, 'CapSize',10);
yline(50, '--', 'chance', 'LineWidth',1.5, 'Color',[.5 .5 .5]);
set(gca, 'XTick',1:4, 'XTickLabel', labels, 'FontSize',12); ylabel('forced-choice accuracy (%)'); ylim([40 100]);
title('Generalization accuracy \pm 95% binomial CI (n=59)');
%[text] An alternative is **repeated cross-validation** (|cv\_splitter.repeated\_kfold(5, n)|, or loop a different rng) and take mean ± SD across repeats — that additionally captures fold-assignment variability. For a clean paired design the binomial CI is usually sufficient.
%%
%[text] ## 6. The one-call legacy wrapper
%[text] |xval\_cross\_classify| packages a cross-classification flow and returns a predictive\_model with results under |pm.cross\_classify|. The step-by-step version above exposes the shared-fold scheme so you can swap algorithm, scorer, or scoring rule.

function [acc, hits] = forced_choice(scores, Y, id)
% Fraction of subjects whose +1 image scores above their -1 image.
    [~, ~, g] = unique(id(:), 'stable');
    u = unique(g); hits = nan(numel(u), 1);
    for i = 1:numel(u)
        sp = scores(g==u(i) & Y== 1);
        sn = scores(g==u(i) & Y==-1);
        if ~isempty(sp) && ~isempty(sn), hits(i) = mean(sp) > mean(sn); end
    end
    acc = 100 * mean(hits, 'omitnan');
end

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline"}
%---
