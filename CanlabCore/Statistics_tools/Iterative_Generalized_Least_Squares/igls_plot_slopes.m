% -------------------------------------------------------------------------
% Plot individual slopes of regressions, using conf. interval for x range
% -------------------------------------------------------------------------
function igls_plot_slopes(out, X, varargin)

    % data
    % -------------------------------------------------------------------
    N = out.sub; % subjects

    betas = out.beta_indiv([2 1], :);  % format for compatibility with mediation.m 
    % slope is first row, intercept is 2nd
    
    groupbeta = out.beta([2 1]);
    
    w = ones(N, 1);
    if ~isempty(varargin)
        % we have stats structure with (maybe) weights
        if isfield(out, 'w')
            w = out.w;
            w = w .* N;
        end
    end

    ncols = 1;
    
    % plot
    % -------------------------------------------------------------------
    fh = create_figure('Slope_Plot');
  
    subplot(1, ncols, 1);
    linehan = plot_individual_slopes(betas, X, w(:,1));
    set(linehan, 'Color', [.7 .7 .7]);
    
    linehan_group = plot_individual_slopes(groupbeta, X(:), 3);
    set(linehan_group, 'LineWidth', 4);
    
    %xlabel(vnames{1}), ylabel(vnames{3});
    title('Slopes');

% %     subplot(1, ncols, 2);
% %     plot_individual_slopes(bbetas, M, w(:,2));
% %     xlabel(vnames{3}), ylabel(vnames{2});
% %     title('b: M->Y controlling X');
% % 
% %     subplot(1, ncols, 3);
% %     plot_individual_slopes(cpbetas, X, w(:,3));
% %     xlabel(vnames{1}), ylabel(vnames{2});
% %     title('c'': X->Y controlling M');
% % 
% %     subplot(1, ncols, 4);
% %     plot_individual_slopes(cbetas, X, w(:,4));
% %     xlabel(vnames{1}), ylabel(vnames{2});
% %     title('c: X->Y');
end


function linehan = plot_individual_slopes(betas, X, w)

    % sort so that we plot from lowest to highest weights
    [w, sorti] = sort(w);
    betas = betas(:,sorti);
    if iscell(X), X = X(sorti); else X = X(:,sorti); end

    % get colors based on weights
    N = size(betas, 2);
    minwt = .4;     % make sure all lines are visible.
    w = rescale_range(w, [minwt 1]);

    % line widths: median split, top half gets 1, bottom gets 1/2.
    linew = (w(:,1) >= median(w(:,1))) +.5;

    w = 1 - w;  % for colors, 0 is black
    colors = repmat([1 1 1], N, 1) .* repmat(w, 1, 3);
    colors(colors < 0) = 0; % include to remove rounding error

    for i = 1:N
        if iscell(X), x = X{i}; else x = X(:,i); end

        % 95% conf. interval for x
        [nanvec, x] = nanremove(x);
        mx = mean(x);
        s = std(x) * tinv(.975, length(x)-1);
        x = [mx - s mx + s];
        y = betas(2, i) + betas(1, i) * x;
        linehan(i) = plot(x, y, '-', 'Color', colors(i,:), 'LineWidth', linew(i));
    end

    drawnow
end

function rx = rescale_range(x, y)
    % re-scale x to range of y
    m = range(y)./range(x);

    if isinf(m)
        % no range/do not rescale
        rx = x;
    else
        b = y(1) - m * x(1);
        rx = m*x + b;
    end
end
