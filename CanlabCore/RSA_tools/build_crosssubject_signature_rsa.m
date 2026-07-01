function RSA = build_crosssubject_signature_rsa(svm_stats_cell, subjects, bodysite_names, varargin)
% build_crosssubject_signature_rsa  Cross-subject x cross-bodysite signature RSA matrix.
%
% Builds the full (nSub*nSites) x (nSub*nSites) similarity matrix among learned
% per-subject x bodysite signature maps (e.g. SVM weight maps), with bookkeeping
% that lets every downstream test exclude same-subject and self comparisons.
%
% WHAT THIS COMPUTES (and what it does NOT)
%   The maps are *backward-model* SVM weight maps that discriminate the within-
%   site contrast (e.g. hot vs warm). Their cross-map similarity therefore
%   measures "do two decoders point the same way", NOT "is bodysite identity
%   represented somatotopically". For the somatotopy question, run RSA on the
%   *evoked response patterns* (compute_rsm) or pass Haufe / structure-coefficient
%   maps here instead of raw weights. See SIGNATURE_RSA_AUDIT.md.
%
% DESIGN CHOICES THAT MAKE THIS DEFENSIBLE
%   - Default metric is 'correlation' (scale-free, the primary RSA-like measure).
%     'cosine' and 'dotproduct' are magnitude-sensitive and provided as secondary
%     diagnostics only.
%   - All maps are pulled into one voxels x conditions matrix on a COMMON voxel
%     support, so the matrix is exactly symmetric and every cell uses the same
%     voxels (unlike apply_mask(a,b,...) whose included-voxel set depends on
%     argument order when supports differ).
%   - Self/diagonal entries are filled but flagged; same-subject pairs are flagged.
%     Downstream inference must use RSA.cross_subject_mask, never the diagonal.
%
% :Usage:
% ::
%   RSA = build_crosssubject_signature_rsa(svm_stats_cell, subjects, bodysite_names)
%   RSA = build_crosssubject_signature_rsa(..., 'metric','correlation')
%   RSA = build_crosssubject_signature_rsa(..., 'map_extractor', @(c) c{1}.weight_obj)
%
% :Inputs:
%   **svm_stats_cell:**  nSub x 1 cell; svm_stats_cell{s} is nSites x 1 cell.
%        By default each leaf svm_stats_cell{s}{b} is unwrapped via the
%        map_extractor to an fmri_data / image_vector map.
%   **subjects:**        cellstr of subject IDs (length nSub).
%   **bodysite_names:**  cellstr of bodysite names in matching order (length nSites).
%
% :Optional Inputs:
%   **'metric':**        'correlation' (default) | 'cosine' | 'dotproduct'.
%   **'map_extractor':** function handle mapping svm_stats_cell{s}{b} to a map
%                        object. Default @(c) c{1}.weight_obj (the obligatory
%                        single-element cell wrap, then the full-data weight map).
%   **'doverbose':**     print progress / QC (default true).
%
% :Outputs:
%   **RSA:** struct with fields
%        .M                  nTotal x nTotal similarity matrix (r-space if 'correlation')
%        .subject_idx        nTotal x 1 subject index per row/col
%        .site_idx           nTotal x 1 bodysite index per row/col
%        .labels             nTotal x 1 'subj | site' labels
%        .nSub, .nSites, .metric, .subjects, .bodysite_names
%        .cross_subject_mask nTotal x nTotal logical, true where subjects differ
%        .same_site_mask     nTotal x nTotal logical, true where sites match
%        .nvox_common        number of voxels in the common support actually used
%        .qc                 struct of diagnostics (per-map voxel counts, dropped maps)
%
% :Examples:
% ::
%   subjects = {'sub-SID000002','sub-SID000743'};            %#ok<NASGU>
%   bodysites = {'leftface','rightface'};                    %#ok<NASGU>
%   % RSA = build_crosssubject_signature_rsa(svm_stats_cell_dpIns, subjects, bodysites);
%   % plot_same_vs_different_site(RSA);
%
% :See also: get_crosssubject_site_effect, subject_level_crosssubject_effect,
%            permutation_test_site_specificity, collapse_rsa_to_bodysite_matrix
%
% ..
%    2026 Michael Sun. GPL v3.
% ..

% ---------- INPUT PARSER ----------
% Keyword-pair options. 'plot'/'verbose' analogues per CanlabCore convention.
p = inputParser;
p.addRequired('svm_stats_cell', @iscell);
p.addRequired('subjects',       @(x) iscell(x) || isstring(x));
p.addRequired('bodysite_names', @(x) iscell(x) || isstring(x));
p.addParameter('metric', 'correlation', @(x) ischar(x) || isstring(x));
p.addParameter('map_extractor', @(c) c{1}.weight_obj, @(x) isa(x,'function_handle'));
p.addParameter('doverbose', true, @(x) islogical(x) || isnumeric(x));
p.parse(svm_stats_cell, subjects, bodysite_names, varargin{:});

metric        = lower(char(p.Results.metric));
map_extractor = p.Results.map_extractor;
doverbose     = logical(p.Results.doverbose);
subjects      = cellstr(subjects);
bodysite_names= cellstr(bodysite_names);

validmetrics = {'correlation','cosine','dotproduct'};
if ~ismember(metric, validmetrics)
    error('build_crosssubject_signature_rsa:metric', ...
        'metric must be one of: %s', strjoin(validmetrics, ', '));
end

nSub   = numel(subjects);
nSites = numel(bodysite_names);
nTotal = nSub * nSites;

% ---------- INDEX BOOKKEEPING ----------
subject_idx = nan(nTotal,1);
site_idx    = nan(nTotal,1);
labels      = cell(nTotal,1);
k = 0;
for s = 1:nSub
    for b = 1:nSites
        k = k + 1;
        subject_idx(k) = s;
        site_idx(k)    = b;
        labels{k}      = sprintf('%s | %s', subjects{s}, bodysite_names{b});
    end
end

% ---------- EXTRACT MAPS INTO A COMMON-SPACE MATRIX ----------
% W is voxels x nTotal. We require every map to share the same in-mask voxel
% count (true for same-ROI signatures). If they differ, we error with guidance
% rather than silently comparing mismatched spaces.
W = [];
nvox_per_map = nan(nTotal,1);
dropped = false(nTotal,1);

for c = 1:nTotal
    s = subject_idx(c); b = site_idx(c);
    try
        mapobj = map_extractor(svm_stats_cell{s}{b});
    catch ME
        warning('build_crosssubject_signature_rsa:extract', ...
            'Could not extract map for subject %d site %d (%s). Column set to NaN.', ...
            s, b, ME.message);
        dropped(c) = true;
        continue
    end

    % Pull a full, aligned voxel vector. replace_empty re-expands image_vector
    % objects to original space so columns are voxel-aligned even if some maps
    % had removed voxels. Plain structs with a .dat field are used as-is.
    if isobject(mapobj) && ismethod(mapobj, 'replace_empty')
        mapobj = replace_empty(mapobj);
    end
    v = double(mapobj.dat(:));

    if isempty(W)
        W = nan(numel(v), nTotal);
    elseif numel(v) ~= size(W,1)
        error('build_crosssubject_signature_rsa:space', ...
            ['Map for subject %d site %d has %d voxels but expected %d. ', ...
             'All maps must share the same voxel space. Resample first ', ...
             '(e.g. resample_space) or build one RSA per ROI.'], ...
            s, b, numel(v), size(W,1));
    end
    W(:,c) = v;
    nvox_per_map(c) = sum(isfinite(v) & v~=0);
end

% ---------- COMMON VALID-VOXEL SUPPORT ----------
% Keep voxels finite in every retained column and not identically zero across
% all of them. This guarantees a symmetric matrix computed on shared voxels.
keepcol = ~dropped;
valid_rows = all(isfinite(W(:,keepcol)), 2) & any(W(:,keepcol) ~= 0, 2);
Wv = W(valid_rows, :);
nvox_common = sum(valid_rows);

if nvox_common < 3
    error('build_crosssubject_signature_rsa:support', ...
        'Common valid voxel support is only %d voxels — cannot compute similarity.', nvox_common);
end
if doverbose
    fprintf('build_crosssubject_signature_rsa: %d conditions, %d common voxels, metric=%s\n', ...
        nTotal, nvox_common, metric);
    if any(dropped)
        fprintf('  WARNING: %d map(s) could not be extracted and are NaN in M.\n', sum(dropped));
    end
end

% ---------- SIMILARITY ----------
M = nan(nTotal, nTotal);
switch metric
    case 'correlation'
        % Pearson correlation among columns over the common support. Symmetric.
        M(keepcol, keepcol) = corr(Wv(:,keepcol));
    case 'cosine'
        Wn = Wv(:,keepcol);
        Wn = Wn ./ vecnorm(Wn);
        M(keepcol, keepcol) = Wn' * Wn;
    case 'dotproduct'
        M(keepcol, keepcol) = Wv(:,keepcol)' * Wv(:,keepcol);
end
% Enforce exact numerical symmetry (guards against round-off).
M = (M + M') / 2;

% ---------- MASKS FOR DOWNSTREAM INFERENCE ----------
cross_subject_mask = subject_idx ~= subject_idx';   % true where subjects differ
same_site_mask     = site_idx   == site_idx';        % true where sites match

% ---------- PACK ----------
RSA = struct();
RSA.M                  = M;
RSA.subject_idx        = subject_idx;
RSA.site_idx           = site_idx;
RSA.labels             = labels;
RSA.nSub               = nSub;
RSA.nSites             = nSites;
RSA.metric             = metric;
RSA.subjects           = subjects;
RSA.bodysite_names     = bodysite_names;
RSA.cross_subject_mask = cross_subject_mask;
RSA.same_site_mask     = same_site_mask;
RSA.nvox_common        = nvox_common;
RSA.qc                 = struct('nvox_per_map', nvox_per_map, 'dropped', dropped);

end
