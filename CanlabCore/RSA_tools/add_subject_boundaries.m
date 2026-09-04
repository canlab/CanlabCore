function add_subject_boundaries(ax, nSub, nSites)
% add_subject_boundaries  Draw subject block separators on a signature RSA image.
%
% Overlays black grid lines at the boundaries between subjects on an imagesc plot
% of a (nSub*nSites) x (nSub*nSites) RSA matrix whose rows/cols are ordered
% subject-major, bodysite-minor.
%
% :Usage:
% ::
%   imagesc(RSA.M, [-1 1]); ax = gca;
%   add_subject_boundaries(ax, RSA.nSub, RSA.nSites);
%
% ..
%    2026 Michael Sun. GPL v3.
% ..

hold(ax, 'on');
bounds = (1:nSub-1) * nSites + 0.5;
for bnd = bounds
    xline(ax, bnd, 'k-', 'LineWidth', 1);
    yline(ax, bnd, 'k-', 'LineWidth', 1);
end
hold(ax, 'off');

end
