function [site_similarity_r, site_similarity_z] = get_sitewise_crosssubject_similarity(RSA)
% get_sitewise_crosssubject_similarity  Per-bodysite cross-subject similarity pools.
%
% For each bodysite b, collects every CROSS-SUBJECT pair of that same bodysite
% (subject i site b vs subject j site b, i ~= j). Same-subject/self excluded.
%
% DESCRIPTIVE ONLY — the pools are non-independent pairwise cells (see
% get_crosssubject_site_effect). For inference across sites use the subject-level
% per-site means returned by subject_level_crosssubject_effect.
%
% :Usage:
% ::
%   [site_r, site_z] = get_sitewise_crosssubject_similarity(RSA)
%
% :Inputs:
%   **RSA:** struct from build_crosssubject_signature_rsa.
%
% :Outputs:
%   **site_similarity_r:** nSites x 1 cell; each cell is a column vector of
%                          cross-subject same-site similarities (r-space).
%   **site_similarity_z:** Fisher-z version ('correlation' metric only; else
%                          empty cells).
%
% :See also: get_crosssubject_site_effect, plot_sitewise_crosssubject_similarity
%
% ..
%    2026 Michael Sun. GPL v3.
% ..

nSites = RSA.nSites;
M  = RSA.M;
ut = triu(true(size(M)), 1);
xs = RSA.cross_subject_mask & ut;

site_similarity_r = cell(nSites,1);
site_similarity_z = cell(nSites,1);
iscorr = strcmp(RSA.metric, 'correlation');

for b = 1:nSites
    site_b   = RSA.site_idx == b;
    bb       = site_b & site_b';                 % both ends are site b
    sel      = xs & bb;
    vals     = M(sel);
    vals     = vals(isfinite(vals));
    site_similarity_r{b} = vals;
    if iscorr
        site_similarity_z{b} = rsa_fisher_z(vals);
    else
        site_similarity_z{b} = [];
    end
end

if ~iscorr
    warning('get_sitewise_crosssubject_similarity:metric', ...
        'Fisher-z only valid for correlation metric (got ''%s''); z cells empty.', RSA.metric);
end

end
