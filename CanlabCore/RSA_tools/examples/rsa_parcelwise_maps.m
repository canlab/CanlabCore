%% Parcelwise RSA -> brain maps
% Reproduces the Sun et al. Figs 7G / S8 / S9 pipeline: per-parcel RSA
% inference projected onto the brain as statistic_image maps.
%
% This example needs the data to be in a real atlas space. We build a small
% synthetic dataset directly in a 10-parcel atlas subset so it runs quickly.

atlas = select_atlas_subset(load_atlas('canlab2024'), 1:10);
codes = double(atlas.dat);

% Plant condition/bodysite structure in two parcels only (3 and 7)
rng(11);
n_sub = 6; conditions = {'hot','warm','imag'}; bodysites = {'face','arm','leg','chest'};
n_cond = 3; n_bs = 4; n_vox = size(atlas.dat, 1);
sig_mask = ismember(codes, [3 7]);
cond_pat = randn(n_cond, n_vox); bs_pat = randn(n_bs, n_vox);
X = []; sub_v = {}; cond_v = {}; bs_v = {}; idx = 0;
for s = 1:n_sub
    for c = 1:n_cond
        for b = 1:n_bs
            for rep = 1:2
                idx = idx + 1; base = 0.3*randn(n_vox, 1);
                base(sig_mask) = base(sig_mask) + 0.8*cond_pat(c, sig_mask)' + 0.5*bs_pat(b, sig_mask)';
                X(:, idx) = base; %#ok<AGROW>
                sub_v{idx,1} = sprintf('sub-%02d', s); %#ok<AGROW>
                cond_v{idx,1} = conditions{c}; bs_v{idx,1} = bodysites{b}; %#ok<AGROW>
            end
        end
    end
end
dat = fmri_data(atlas, 'noverbose');
dat.dat = X;
dat.metadata_table = table(sub_v, cond_v, bs_v, 'VariableNames', {'sub','condition','bodysite'});

%% Contrast path: within-hot per parcel, FDR across parcels
spec = { 'within_hot', 'hot', [];  'hot_vs_warm', 'hot', 'warm' };
results = dat.rsa_parcelwise('atlas', atlas, ...
    'group_by', {'condition','bodysite'}, 'subject_var', 'sub', ...
    'metric', 'spearman', 'contrasts', spec, 'correction', 'fdr', 'tail', 'right');

disp(results.contrast_table)

% The within_hot map is a statistic_image -- threshold and montage it
% results.maps.within_hot.threshold(0.05, 'unc').montage;
fprintf('within_hot map: %d sig voxels at FDR\n', ...
    sum(results.maps.within_hot.sig));

%% LME path: per-parcel mixed model, map each fixed-effect term
results_lme = dat.rsa_parcelwise('atlas', atlas, ...
    'group_by', {'condition','bodysite'}, 'subject_var', 'sub', ...
    'lme', 'Y ~ SameCondition + SameBodysite + (1|Sub)', ...
    'predictors', {'condition','bodysite'});

disp(results_lme.lme_table)
% results_lme.maps.SameCondition.montage;

%% Searchlight RSA (restrict to a mask for speed)
% Build a small mask and run a model-RDM searchlight
mask = fmri_data(atlas, 'noverbose'); mask.dat = double(ismember(codes, [1 2 3 7]));
dat_masked = apply_mask(dat, mask);
sl = dat_masked.searchlight_rsa('condition', 'radius', 3, ...
    'group_by', {'condition','bodysite'}, 'subject_var', 'sub', ...
    'metric', 'spearman', 'permutations', 100);
fprintf('Searchlight condition map built (%d centers).\n', sl.n_spheres);
% sl.maps.condition.montage;
