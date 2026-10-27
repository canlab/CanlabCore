function out = plot_sitewise_crosssubject_similarity(RSA_list, varargin)
% plot_sitewise_crosssubject_similarity  Per-bodysite cross-subject consistency.
%
% Plots cross-subject same-bodysite similarity broken down by bodysite, for one
% or more ROIs. Points are subject-level means (each subject's mean cross-subject
% same-site similarity for that bodysite); error bars are SEM across subjects.
%
% :Usage:
% ::
%   plot_sitewise_crosssubject_similarity(RSA)
%   plot_sitewise_crosssubject_similarity({RSA_s1, RSA_dp}, 'names', {'S1','dpIns'})
%
% :Inputs:
%   **RSA_list:** a single RSA struct or cell array of RSA structs (shared
%                 bodysite ordering / names).
%
% :Optional Inputs:
%   **'names':** cellstr ROI labels.
%   **'axes':**  target axes (default new figure).
%
% :Outputs:
%   **out:** struct with per-ROI nSites x 1 subject-level means (r-space) and SEMs.
%
% :See also: get_sitewise_crosssubject_similarity, subject_level_crosssubject_effect
%
% ..
%    2026 Michael Sun. GPL v3.
% ..

if ~iscell(RSA_list), RSA_list = {RSA_list}; end
nROI = numel(RSA_list);

p = inputParser;
p.addParameter('names', arrayfun(@(i) sprintf('ROI%d', i), 1:nROI, 'UniformOutput', false), @iscell);
p.addParameter('axes', [], @(x) isempty(x) || ishandle(x));
p.parse(varargin{:});
names = p.Results.names;
ax    = p.Results.axes;

if isempty(ax), figure('Color','w'); ax = gca; end
hold(ax, 'on');

bodysite_names = RSA_list{1}.bodysite_names;
nSites = RSA_list{1}.nSites;
iscorr = strcmp(RSA_list{1}.metric, 'correlation');

means = nan(nSites, nROI);
sems  = nan(nSites, nROI);

for i = 1:nROI
    o = subject_level_crosssubject_effect(RSA_list{i}, 'doverbose', false);
    % sitewise_same_z is nSub x nSites (z-space for correlation).
    sw = o.sitewise_same_z;
    if iscorr
        sw_r = tanh(sw);          % per-subject, per-site, back to r
    else
        sw_r = sw;
    end
    means(:,i) = mean(sw_r, 1, 'omitnan')';
    sems(:,i)  = (std(sw_r, 0, 1, 'omitnan') ./ sqrt(sum(isfinite(sw_r), 1)))';
    errorbar(ax, 1:nSites, means(:,i), sems(:,i), '-o', 'LineWidth', 2);
end

xticks(ax, 1:nSites);
xticklabels(ax, bodysite_names);
xtickangle(ax, 45);
ylabel(ax, sprintf('Cross-subject similarity (r)\n[%s metric]', RSA_list{1}.metric));
title(ax, 'Across-subject consistency of bodysite representations (subject-level \pm SEM)');
yh = yline(ax, 0, 'k--'); yh.HandleVisibility = 'off';   % keep out of legend
legend(ax, names, 'Location', 'best');
box(ax, 'off');
hold(ax, 'off');

out = struct('names', {names}, 'bodysite_names', {bodysite_names}, ...
    'means', means, 'sems', sems);

end
