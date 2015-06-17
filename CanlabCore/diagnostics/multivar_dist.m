% [ds, S, p] = multivar_dist(X)
%
% multivariate normality checking
% and diagnostic plots
%
% given matrix X with cases = rows, cols = variables
%
% Outputs:
% ----------------------------------------------------
% ds: is matrix of squared distances, case numbers, and
%   expected chi2 values (in columns in this order)
%   rows are cases
%   *NOTE: Sorted in order of ascending distance!
% S: estimated covariance matrix
% mv_distance: squared distances in original order of rows
% p: p-values in original order of rows
%
% by Tor Wager

function [ds, S, p] = multivar_dist(X)
    % -----------------------------------------------------
    % * determine multivariate standard deviation matrix S
    % -----------------------------------------------------

    % center
    Xs = X - repmat(mean(X), size(X, 1), 1);

    % covariance matrix S
    S = (Xs' * Xs) ./ (size(Xs, 1)-1);

    % -----------------------------------------------------
    % * get squared distance
    % -----------------------------------------------------

    % squared distance, Johnson & Wichern p. 201
    % (X-mean(X))' * inv(S) * (X - mean(X)) generalized to matrices
    d = Xs * inv(S) * Xs';
    d = diag(d);            % what do the off-diagonals signify?

    % -----------------------------------------------------
    % * compare with chi2 distribution
    % -----------------------------------------------------

    % calculate chi2 threshold
    % chi2 value compares the number of points within the ellipsoid
    % contour to those outside; roughly alpha of the squared distances
    % should be within the ellipsoid (for rough general test of normality).
    % outliers will have very high chi2 values

    t = chi2inv(.5, size(X, 2));  % for general test
    p50 = 100 * (sum(d > t) ./ length(d));
    fprintf(1, 'Expected 50%% of points within 50%% normal ellipsoid, found %3.2f%%\n', p50);

    t = chi2inv(.95, size(X, 2));  % for general test
    p95 = sum(d > t);
    fprintf(1, 'Expected %3.2f outside 95%% ellipsoid, found %3.0f\n', .05*length(d), p95);

    % -----------------------------------------------------
    % * get case numbers and sort by distance
    % -----------------------------------------------------
    d(:,2) = (1:length(d))';
    ds = sortrows(d, 1);


    % -----------------------------------------------------
    % * get chi2 quantiles for qd plot
    % -----------------------------------------------------
    q = (((1:size(ds, 1)) - .5) ./ size(ds, 1))';
    q = chi2inv(q, size(X, 2));
    ds(:,3) = q;

    % -----------------------------------------------------
    % * get distance^2 in original order and p-values
    % -----------------------------------------------------
    ds = sortrows(ds, 2);
    p = 1 - chi2cdf(ds(:,1), size(X, 2));


    % -----------------------------------------------------
    % * plot the results in a figure
    % -----------------------------------------------------
    figure('color', 'w');
    subplot(1, 2, 1); hold on; grid on
    plot([1:size(d, 1); 1:size(d, 1)], [zeros(size(d, 1), 1) d(:,1)]', 'b', 'LineWidth', 1.5);
    xlabel('Case number');
    ylabel('Squared stat. distance from origin');
    wh = (d(:,1) > t); d2 = d; d2(:,1) = d2(:,1) .* wh;  % zero out the non-"significant" chi2 cases
    plot([1:size(d2, 1); 1:size(d2, 1)], [zeros(size(d2, 1), 1) d2(:,1)]', 'r', 'LineWidth', 1.5);
    title('d^2, red cases outside 95% normal ellipsoid');
    plot([0 size(d, 1)], [t t], 'k');
    set(gca, 'YLim', [0 max(t+.5, max(d(:,1)))]);

    subplot(1, 2, 2); hold on; grid on
    plot(ds(:,3), ds(:,1), 'MarkerSize', 0.01, 'Color', 'w');
    for i = 1:size(ds, 1)
        text(ds(i,3), ds(i, 1), num2str(ds(i, 2)), 'Color', 'k');
    end
    xlabel('Expected chi2 value'), ylabel('Squared distance');
    plot([0 max([ds(:,3);ds(:,1)])], [0 max([ds(:,3);ds(:,1)])], 'k', 'LineWidth', 1.5);
    title('Line with slope = 1 is normal');
end

