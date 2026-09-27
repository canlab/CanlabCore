function out = signature_specificity_parcelwise(parcel_svm_cells, subjects, bodysite_names, atlas_obj, varargin)
% signature_specificity_parcelwise  Brain-wide cross-subject bodysite specificity.
%
% For each parcel of a parcellation, builds the cross-subject signature RSA from
% that parcel's per-subject x bodysite maps, computes cross-subject consistency
% and same>different specificity with a defensible subject-level test (and
% optional permutation p), then projects the parcel-level summaries onto the
% brain and FDR-corrects the significance map.
%
% :Usage:
% ::
%   out = signature_specificity_parcelwise(parcel_svm_cells, subjects, ...
%             bodysite_names, atlas_obj)
%   out = signature_specificity_parcelwise(..., 'nperm', 2000, 'metric','correlation')
%
% :Inputs:
%   **parcel_svm_cells:** nParcels x 1 cell. parcel_svm_cells{pidx} is itself a
%        subject x bodysite svm_stats_cell for that parcel (same layout that
%        build_crosssubject_signature_rsa consumes).
%   **subjects:**         cellstr of subject IDs.
%   **bodysite_names:**   cellstr of bodysite names.
%   **atlas_obj:**        atlas object whose parcels correspond, in label order,
%        to parcel_svm_cells (numel(parcel_svm_cells) == num_regions(atlas_obj)).
%
% :Optional Inputs:
%   **'metric':**     'correlation' (default) | 'cosine' | 'dotproduct'.
%   **'nperm':**      permutations per parcel (default 0 = skip permutation,
%                     use subject-level t p-values; set e.g. 2000 for perm).
%   **'parcel_names':** cellstr labels for the table/maps; default atlas.labels.
%   **'doverbose':**  progress (default true).
%
% :Outputs:
%   **out:** struct with
%      .table              nParcels-row table: parcel, n_ok, consistency_r,
%                          specificity_z, t, df, p, p_fdr
%      .consistency_map    statistic_image: cross-subject same-site consistency (r)
%      .specificity_map    statistic_image: same-diff specificity (z), p from test
%      .specificity_fdr    thresholded specificity_map (FDR q<.05, right tail)
%      .RSA                nParcels x 1 cell of per-parcel RSA structs (for reuse)
%      .params
%
% :Examples:
% ::
%   atlas = load_atlas('canlab2024');
%   % parcel_svm_cells{p} produced by training one SVM set per parcel.
%   out = signature_specificity_parcelwise(parcel_svm_cells, subjects, ...
%             bodysite_names, atlas, 'nperm', 2000);
%   montage(out.specificity_fdr);
%   disp(out.table)
%
% :See also: build_crosssubject_signature_rsa, subject_level_crosssubject_effect,
%            permutation_test_site_specificity, assign_vals_to_atlas
%
% ..
%    2026 Michael Sun. GPL v3.
% ..

p = inputParser;
p.addParameter('metric', 'correlation', @(x) ischar(x) || isstring(x));
p.addParameter('nperm', 0, @(x) isnumeric(x) && isscalar(x) && x >= 0);
p.addParameter('parcel_names', [], @(x) isempty(x) || iscell(x));
p.addParameter('doverbose', true, @(x) islogical(x) || isnumeric(x));
p.parse(varargin{:});
metric    = lower(char(p.Results.metric));
nperm     = round(p.Results.nperm);
doverbose = logical(p.Results.doverbose);

nParcels = numel(parcel_svm_cells);
if isempty(p.Results.parcel_names)
    parcel_names = atlas_obj.labels(:);
else
    parcel_names = p.Results.parcel_names(:);
end
if numel(parcel_names) ~= nParcels
    error('signature_specificity_parcelwise:nameMismatch', ...
        '%d parcels but %d parcel names.', nParcels, numel(parcel_names));
end

consistency_r = nan(nParcels,1);
specificity_z = nan(nParcels,1);
tval = nan(nParcels,1);
dfval = nan(nParcels,1);
pval = nan(nParcels,1);
n_ok = zeros(nParcels,1);
RSAc = cell(nParcels,1);

for pidx = 1:nParcels
    if doverbose && mod(pidx, 25) == 0
        fprintf('  parcel %d/%d\n', pidx, nParcels);
    end
    try
        RSA = build_crosssubject_signature_rsa(parcel_svm_cells{pidx}, subjects, ...
            bodysite_names, 'metric', metric, 'doverbose', false);
        RSAc{pidx} = RSA;

        o = subject_level_crosssubject_effect(RSA, 'doverbose', false);
        consistency_r(pidx) = mean(o.same_r, 'omitnan');
        specificity_z(pidx) = o.mean_spec_z;
        tval(pidx) = o.t;
        dfval(pidx) = o.df;
        n_ok(pidx) = sum(isfinite(o.spec_z));

        if nperm > 0
            pr = permutation_test_site_specificity(RSA, 'nperm', nperm, ...
                'tail', 'right', 'doverbose', false);
            pval(pidx) = pr.p;
        else
            pval(pidx) = o.p;   % subject-level t p (right tail)
        end
    catch ME
        if doverbose
            fprintf('  parcel %d (%s) skipped: %s\n', pidx, parcel_names{pidx}, ME.message);
        end
    end
end

% ---------- FDR across parcels ----------
p_fdr = nan(nParcels,1);
good  = isfinite(pval);
if any(good)
    p_fdr(good) = local_fdr_bh(pval(good));
end

% ---------- Table ----------
out = struct();
out.table = table(parcel_names, n_ok, consistency_r, specificity_z, tval, dfval, pval, p_fdr, ...
    'VariableNames', {'parcel','n_subjects','consistency_r','specificity_z','t','df','p','p_fdr'});

% ---------- Brain maps ----------
% Consistency: descriptive r per parcel (fmri_data is fine, no inference).
out.consistency_map = assign_vals_to_atlas(atlas_obj, [], consistency_r, ...
    'output_type', 'fmri_data', 'dat_descrip', 'cross-subject same-site consistency (r)');

% Specificity: statistic_image carrying the per-parcel p so thresholding works.
out.specificity_map = assign_vals_to_atlas(atlas_obj, [], specificity_z, ...
    'p_vals', pval, 'output_type', 'statistic_image', ...
    'dat_descrip', 'same-diff bodysite specificity (z)');

% FDR-thresholded version: threshold on the BH-corrected p directly.
% (Use function syntax — statistic_image has a `threshold` property AND method.)
spec_fdr = out.specificity_map;
spec_fdr.p = p_fdr;                 % swap in BH-corrected p
spec_fdr.p(~isfinite(spec_fdr.p)) = 1;
out.specificity_fdr = threshold(spec_fdr, 0.05, 'unc');   % q<.05 (p already FDR-corrected)

out.RSA = RSAc;
out.params = struct('metric', metric, 'nperm', nperm);

if doverbose
    nsig = sum(out.table.p_fdr < 0.05);
    fprintf('signature_specificity_parcelwise: %d/%d parcels with FDR q<.05 specificity.\n', ...
        nsig, nParcels);
end

end


function pcorr = local_fdr_bh(pvec)
% Benjamini-Hochberg FDR-corrected p-values (monotone, capped at 1).
pvec = pvec(:);
n = numel(pvec);
[ps, ord] = sort(pvec);
ranks = (1:n)';
pc = ps .* n ./ ranks;
% enforce monotonicity from the largest p downward
for i = n-1:-1:1
    pc(i) = min(pc(i), pc(i+1));
end
pc = min(pc, 1);
pcorr = nan(n,1);
pcorr(ord) = pc;
end
