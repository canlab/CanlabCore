function plot_signature_rsa_matrix(RSA, varargin)
% plot_signature_rsa_matrix  Full subject x site RSA matrix and collapsed site matrix.
%
% Draws (left) the full (nSub*nSites)^2 signature RSA with subject block
% boundaries, and (right) the cross-subject-only bodysite x bodysite collapsed
% matrix. For the correlation metric the color limits are [-1 1]; otherwise auto.
%
% :Usage:
% ::
%   plot_signature_rsa_matrix(RSA)
%   plot_signature_rsa_matrix(RSA, 'title', 'dpIns')
%
% :Inputs:
%   **RSA:** struct from build_crosssubject_signature_rsa.
%
% :Optional Inputs:
%   **'title':** string prefix for panel titles (default RSA.metric).
%
% :See also: add_subject_boundaries, collapse_rsa_to_bodysite_matrix
%
% ..
%    2026 Michael Sun. GPL v3.
% ..

p = inputParser;
p.addParameter('title', RSA.metric, @(x) ischar(x) || isstring(x));
p.parse(varargin{:});
ttl = char(p.Results.title);

iscorr = strcmp(RSA.metric, 'correlation');
if iscorr, clim = [-1 1]; else, clim = []; end

% Use the lab diverging colormap if available, else default.
try cm = colormap_tor([0 0 1],[1 0 0],[1 1 1]); catch, cm = parula; end

figure('Color','w','Position',[100 100 1200 480]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

% --- full matrix ---
ax1 = nexttile;
if isempty(clim), imagesc(ax1, RSA.M); else, imagesc(ax1, RSA.M, clim); end
axis(ax1,'square'); colorbar(ax1); colormap(ax1, cm);
add_subject_boundaries(ax1, RSA.nSub, RSA.nSites);
title(ax1, sprintf('%s: full subject\\times site RSA', ttl));
xlabel(ax1, 'Subject \times bodysite'); ylabel(ax1, 'Subject \times bodysite');

% --- collapsed (cross-subject only) ---
G = collapse_rsa_to_bodysite_matrix(RSA, true);
ax2 = nexttile;
if isempty(clim), imagesc(ax2, G); else, imagesc(ax2, G, clim); end
axis(ax2,'square'); colorbar(ax2); colormap(ax2, cm);
title(ax2, sprintf('%s: cross-subject bodysite\\times bodysite', ttl));
xticks(ax2, 1:RSA.nSites); yticks(ax2, 1:RSA.nSites);
xticklabels(ax2, RSA.bodysite_names); yticklabels(ax2, RSA.bodysite_names);
xtickangle(ax2, 45);

end
