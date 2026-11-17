function out = plot_same_vs_different_site(RSA_list, varargin)
% plot_same_vs_different_site  Same-site vs different-site cross-subject similarity bars.
%
% Plots cross-subject same-bodysite vs different-bodysite signature similarity
% for one or more ROIs/RSA structs side by side. Bars are subject-level means
% and error bars are SEM ACROSS SUBJECTS (the honest unit of replication), not
% across non-independent pairwise cells.
%
% :Usage:
% ::
%   plot_same_vs_different_site(RSA)                          % single ROI
%   plot_same_vs_different_site({RSA_s1, RSA_dp}, 'names', {'S1','dpIns'})
%
% :Inputs:
%   **RSA_list:** a single RSA struct or a cell array of RSA structs.
%
% :Optional Inputs:
%   **'names':** cellstr of ROI labels (default 'ROI1','ROI2',...).
%   **'axes':**  target axes handle (default: new figure).
%
% :Outputs:
%   **out:** struct with per-ROI subject-level means (same_r, diff_r), SEMs, and
%            the subject-level test (t, p, df) from subject_level_crosssubject_effect.
%
% :See also: subject_level_crosssubject_effect, plot_sitewise_crosssubject_similarity
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

if isempty(ax)
    figure('Color','w'); ax = gca;
end
hold(ax, 'on');

same_m = nan(nROI,1); diff_m = nan(nROI,1);
same_se = nan(nROI,1); diff_se = nan(nROI,1);
stats = cell(nROI,1);

for i = 1:nROI
    o = subject_level_crosssubject_effect(RSA_list{i}, 'doverbose', false);
    stats{i} = o;
    same_m(i)  = mean(o.same_r, 'omitnan');
    diff_m(i)  = mean(o.diff_r, 'omitnan');
    same_se(i) = std(o.same_r, 'omitnan') / sqrt(sum(isfinite(o.same_r)));
    diff_se(i) = std(o.diff_r, 'omitnan') / sqrt(sum(isfinite(o.diff_r)));
end

% Grouped bars: per ROI a pair (same, diff).
xc = (1:nROI) * 3;          % cluster centers
w  = 0.6;
for i = 1:nROI
    bar(ax, xc(i)-0.5, same_m(i), w);
    bar(ax, xc(i)+0.5, diff_m(i), w);
    errorbar(ax, xc(i)-0.5, same_m(i), same_se(i), 'k.', 'LineWidth', 1.2);
    errorbar(ax, xc(i)+0.5, diff_m(i), diff_se(i), 'k.', 'LineWidth', 1.2);
    % significance annotation from subject-level test
    star = '';
    if stats{i}.p < .001, star = '***'; elseif stats{i}.p < .01, star = '**'; ...
    elseif stats{i}.p < .05, star = '*'; end
    yt = max(same_m(i)+same_se(i), diff_m(i)+diff_se(i));
    text(ax, xc(i), yt*1.1, sprintf('%s p=%.3g', star, stats{i}.p), ...
        'HorizontalAlignment','center', 'FontSize', 9);
end

xticks(ax, xc);
xticklabels(ax, names);
ylabel(ax, sprintf('Cross-subject similarity (r)\n[%s metric]', RSA_list{1}.metric));
title(ax, 'Same-bodysite vs different-bodysite (subject-level \pm SEM)');
yh = yline(ax, 0, 'k--'); yh.HandleVisibility = 'off';   % keep out of legend
legend(ax, {'Same site','Different site'}, 'Location', 'best');
box(ax, 'off');
hold(ax, 'off');

out = struct('names', {names}, 'same_m', same_m, 'diff_m', diff_m, ...
    'same_se', same_se, 'diff_se', diff_se, 'stats', {stats});

end
