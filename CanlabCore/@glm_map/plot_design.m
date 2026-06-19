function plot_design(obj, varargin)
% Plot the design matrix and per-regressor VIFs for a glm_map object.
%
% Shows the design matrix X as an image (columns = regressors) and, when
% available, a bar plot of the variance inflation factor for each regressor.
%
% :Usage:
% ::
%
%     plot_design(obj)
%
% :Inputs:
%
%   **obj:**
%        A glm_map object with a design available (obj.X non-empty).
%
% :Outputs:
%
%        A figure with the design matrix and (if computed) VIF bars.
%
% :Examples:
% ::
%
%     g = glm_map('X', [ones(30,1) zscore((1:30)')], 'level', 2);
%     g = diagnostics(g, 'noverbose');
%     plot_design(g);
%
% :See also:
%   - diagnostics, fmri_glm_design_matrix.plot
%
% ..
%    Programmers' notes:
%    2026 - Initial implementation.
% ..

X = obj.X;
if isempty(X)
    error('glm_map:NoDesign', 'No design matrix available (obj.X is empty).');
end

rn = obj.regressor_names;
if isstruct(obj.diagnostics) && isfield(obj.diagnostics, 'Variance_inflation_factors')
    vifvals = obj.diagnostics.Variance_inflation_factors;
else
    vifvals = [];
end
has_vif = ~isempty(vifvals);

create_figure('glm_map design', 1, 1 + has_vif);

% --- Design matrix image ---
subplot(1, 1 + has_vif, 1);
imagesc(X);
colormap(gray);
colorbar;
xlabel('Regressor');
ylabel('Image / observation');
title('Design matrix X');
if ~isempty(rn) && numel(rn) == size(X, 2)
    set(gca, 'XTick', 1:size(X, 2), 'XTickLabel', rn, 'XTickLabelRotation', 45);
end

% --- VIF bars ---
if has_vif
    subplot(1, 2, 2);
    bar(vifvals);
    ylabel('Variance inflation factor');
    xlabel('Regressor');
    title('VIF per regressor');
    if ~isempty(rn) && numel(rn) == numel(vifvals)
        set(gca, 'XTick', 1:numel(vifvals), 'XTickLabel', rn, 'XTickLabelRotation', 45);
    end
    hold on;
    yl = get(gca, 'YLim');
    plot(get(gca, 'XLim'), [4 4], 'r--');   % conventional VIF = 4 reference
    set(gca, 'YLim', [0 max(yl(2), 4.5)]);
end

end % plot_design
