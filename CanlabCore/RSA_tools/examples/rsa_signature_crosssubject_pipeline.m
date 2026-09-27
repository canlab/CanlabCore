%% Cross-subject bodysite SVM-signature RSA pipeline (corrected / defensible)
% Consolidates the WASABI AcceptMap / DistractMap "Bodysite Pain SVM Signature
% Visualization Reports" signature-RSA workflow into the RSA_tools API, with
% the statistical fixes documented in RSA_tools/SIGNATURE_RSA_AUDIT.md.
%
% SCIENTIFIC SCOPE -- read this first.
%   These maps are SVM *weight* maps trained to discriminate the within-site
%   contrast (e.g. hot vs warm). Their cross-map similarity measures whether two
%   decoders point the same way -- NOT whether bodysite identity is represented
%   somatotopically. For the somatotopy question, prefer data-level RSA on the
%   evoked patterns (compute_rsm), or pass Haufe / structure-coefficient maps
%   into the builder instead of raw weights. This script answers the narrower,
%   well-posed question: "are the learned signatures consistent across subjects,
%   and more so for the matched bodysite than for mismatched bodysites?"

%% 1. Load the per-subject x bodysite SVM outputs for each ROI ------------
base = "\\dartfs-hpc\rc\lab\C\CANlab\labdata\projects\WASABI\WASABI_N_of_Few\analysis\WASABI-NofFew_BodyMap\results\SPM_FAST\PFR_painpathways";
tmp1 = load(fullfile(base, "pf_roi_svm_s1.mat"),    'svm_stats_cell');
tmp2 = load(fullfile(base, "pf_roi_svm_dpIns.mat"), 'svm_stats_cell');
svm_s1 = tmp1.svm_stats_cell;
svm_dp = tmp2.svm_stats_cell;

subjects = {'sub-SID000002','sub-SID000743','sub-SID001567','sub-SID001641', ...
            'sub-SID001651','sub-SID001684','sub-SID001804','sub-SID002328'};
bodysite_names = {'leftface','rightface','leftarm','rightarm', ...
                  'leftleg','rightleg','chest','abdomen'};

%% 2. Build the RSA matrices ---------------------------------------------
% Correlation is the PRIMARY (scale-free) metric. Dot product is magnitude-driven
% and kept only as a secondary diagnostic.
RSA_s1 = build_crosssubject_signature_rsa(svm_s1, subjects, bodysite_names, 'metric','correlation');
RSA_dp = build_crosssubject_signature_rsa(svm_dp, subjects, bodysite_names, 'metric','correlation');

% Optional: use Haufe / structure-coefficient maps instead of raw weights for a
% more interpretable representational comparison (recommended secondary analysis):
%   ext = @(c) structure_coefficients(c{1}.weight_obj, c{1}.dat_used);  % needs training data
%   RSA_dp_haufe = build_crosssubject_signature_rsa(svm_dp, subjects, bodysite_names, 'map_extractor', ext);

%% 3. Visualize the geometry ---------------------------------------------
plot_signature_rsa_matrix(RSA_s1, 'title', 'S1');
plot_signature_rsa_matrix(RSA_dp, 'title', 'dpIns');

%% 4. Same-site vs different-site (subject-level inference) ---------------
% Subject-level test (df = nSubjects-1), NOT ttest2 on non-independent pairs.
o_s1 = subject_level_crosssubject_effect(RSA_s1, 'tail','right');
o_dp = subject_level_crosssubject_effect(RSA_dp, 'tail','right');

% Permutation test (shuffles site labels within subject) -- distribution-free.
perm_s1 = permutation_test_site_specificity(RSA_s1, 'nperm', 5000);
perm_dp = permutation_test_site_specificity(RSA_dp, 'nperm', 5000);

plot_same_vs_different_site({RSA_s1, RSA_dp}, 'names', {'S1','dpIns'});

%% 5. dpIns vs S1 cross-subject consistency (subject-level, paired) -------
% One value per subject per ROI, then paired t across the 8 subjects.
sub_dp = mean(o_dp.sitewise_same_z, 2, 'omitnan');   % z-space per subject
sub_s1 = mean(o_s1.sitewise_same_z, 2, 'omitnan');
[~, p_roi, ci_roi, st_roi] = ttest(sub_dp, sub_s1, 'tail','right');
fprintf('\ndpIns > S1 cross-subject consistency: paired t(%d)=%.2f, p=%.4g, 95%%CI[%.3f %.3f]\n', ...
    st_roi.df, st_roi.tstat, p_roi, ci_roi(1), ci_roi(2));

%% 6. Sitewise cross-subject consistency by bodysite ---------------------
plot_sitewise_crosssubject_similarity({RSA_s1, RSA_dp}, 'names', {'S1','dpIns'});

%% 7. (Optional) Brain-wide parcelwise specificity -----------------------
% Requires one svm_stats_cell per parcel: parcel_svm_cells{pidx} = subject x site.
% atlas = load_atlas('canlab2024');
% out = signature_specificity_parcelwise(parcel_svm_cells, subjects, ...
%           bodysite_names, atlas, 'nperm', 2000);
% disp(out.table)
% montage(out.consistency_map);          % cross-subject consistency (r)
% montage(out.specificity_fdr);          % same-diff specificity, FDR q<.05

fprintf('\nSignature cross-subject RSA pipeline complete.\n');
