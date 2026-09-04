function [G, G_z] = collapse_rsa_to_bodysite_matrix(RSA, exclude_within_subject)
% collapse_rsa_to_bodysite_matrix  Collapse subject x site RSA to a site x site matrix.
%
% Averages the full (nSub*nSites)^2 signature RSA into an nSites x nSites
% bodysite-by-bodysite matrix. By default (recommended) within-subject pairs are
% excluded, so the result is a purely CROSS-SUBJECT representational similarity
% matrix: G(i,j) = mean similarity of site-i maps to site-j maps across different
% subjects. The diagonal G(i,i) is then the cross-subject same-site consistency
% (NOT trivial self-correlation, because i and j come from different subjects).
%
% Averaging is done in Fisher-z space for the correlation metric (then mapped
% back to r), and in raw space otherwise.
%
% :Usage:
% ::
%   [G, G_z] = collapse_rsa_to_bodysite_matrix(RSA)
%   G = collapse_rsa_to_bodysite_matrix(RSA, false)   % include within-subject
%
% :Inputs:
%   **RSA:** struct from build_crosssubject_signature_rsa.
%   **exclude_within_subject:** logical, default true.
%
% :Outputs:
%   **G:**   nSites x nSites mean similarity matrix in r-space (or raw for
%            non-correlation metrics).
%   **G_z:** the matrix in z-space (correlation metric only; else []).
%
% :See also: build_crosssubject_signature_rsa, get_sitewise_crosssubject_similarity
%
% ..
%    2026 Michael Sun. GPL v3.
% ..

if nargin < 2 || isempty(exclude_within_subject)
    exclude_within_subject = true;
end

nSites = RSA.nSites;
M      = RSA.M;
iscorr = strcmp(RSA.metric, 'correlation');

% Base pair mask: exclude self (diagonal). Optionally exclude same-subject.
base = ~eye(size(M));
if exclude_within_subject
    base = base & RSA.cross_subject_mask;
end

G   = nan(nSites, nSites);
G_z = nan(nSites, nSites);

for i = 1:nSites
    for j = 1:nSites
        sel  = base & (RSA.site_idx == i) & (RSA.site_idx' == j);
        vals = M(sel);
        vals = vals(isfinite(vals));
        if isempty(vals), continue; end
        if iscorr
            zz       = rsa_fisher_z(vals);
            G_z(i,j) = mean(zz, 'omitnan');
            G(i,j)   = tanh(G_z(i,j));      % report mean in r-space
        else
            G(i,j)   = mean(vals, 'omitnan');
        end
    end
end

if ~iscorr, G_z = []; end

end
