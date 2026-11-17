%% Multi-level LME modeling of RSA structure
% Reproduces the 08072024 Run-Level RDM Analysis workflow: model pairwise
% RSM similarity as a function of same-condition / same-bodysite / same-
% session, with subject as a random effect.

dat = make_synthetic_rsa_data('n_sub', 9, 'n_ses', 3);

%% 1. Omnibus fixed-effects: which factors structure the representation?
[mdl_lm, tbl_lm] = dat.rsa_lm('predictors', {'sub','sesno','condition','bodysite'}, ...
                              'subject_var', 'sub', 'pair_scope', 'all');
T_partial = rsa_partial_r2(mdl_lm, tbl_lm);
fprintf('\nPartial R^2 per factor (omnibus fixed-effects):\n');
disp(T_partial)

%% 2. Per-subject fits -> 2nd-level inference
T_by_sub = dat.rsa_lm_by_subject('predictors', {'condition','bodysite','sesno'}, ...
                                 'subject_var', 'sub', 'verbose', false);
sc = T_by_sub.beta(strcmp(T_by_sub.term, 'SameCondition'));
[~, pcond, ~, st] = ttest(sc);
fprintf('SameCondition beta across subjects: mean=%.3f, t(%d)=%.2f, p=%.4g\n', ...
    mean(sc), st.df, st.tstat, pcond);

%% 3. Random-effects LME (the principled multi-level model)
mdl = dat.rsa_lme('Y ~ SameCondition + SameBodysite + SameSesno + (SameCondition | Sub)', ...
                  'predictors', {'condition','bodysite','sesno'}, 'subject_var', 'sub');
disp(mdl.Coefficients)

%% 4. Variance decomposition + per-subject BLUPs
icc = rsa_lme_icc(mdl);
fprintf('\nVariance components:\n'); disp(icc.summary)
blups = rsa_lme_blups(mdl);
fprintf('Per-subject SameCondition BLUPs:\n');
disp(blups(strcmp(blups.Term, 'SameCondition'), :))

%% 5. Nested model comparison (AIC/BIC/LRT ladder)
seq = rsa_model_sequence('Y ~ 1 + (1|Sub)');
seq = seq.add_term('SameCondition');
seq = seq.add_term('SameBodysite');
seq = seq.add_term('SameSesno');
[T_models, best] = dat.rsa_compare_models(seq.formulas, ...
    'predictors', {'condition','bodysite','sesno'}, 'subject_var', 'sub', ...
    'select_by', 'aic', 'verbose', false);
disp(T_models)
fprintf('Best model by AIC: %s\n', seq.formulas{best});
