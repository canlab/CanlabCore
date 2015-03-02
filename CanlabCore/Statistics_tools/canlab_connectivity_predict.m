function OUT = canlab_connectivity_predict(dat, subject_grouping, varargin)
% Connectivity and multivariate pattern-based prediction for multi-subject timeseries data
% Currently runs predictive algorithm(s) on pairwise correlations among regions
%
% Usage:
% -------------------------------------------------------------------------
% OUT = canlab_connectivity_predict(dat, subject_grouping, ['outcome', outcome_dat])
%
% Features:
% - Within-subject correlation matrices and 'random effects' statistics
% - [optional] Prediction with LASSO-PCR/SVR/SVM of outcomes from pairwise connectivity
% - Time-lagged cross-correlation options
% - Graph theoretic measures
% - [optional] Prediction with LASSO-PCR/SVR/SVM of outcomes from graph measures
% - Can easily be extended to handle partial regression/correlation coefficients
%
% Author and copyright information:
% -------------------------------------------------------------------------
%     Copyright (C) 2015  Tor Wager
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
% Inputs:
% -------------------------------------------------------------------------
% dat           concatenated data matrix of time points (t) within subjects x variables (e.g., ROIs)
%               [subj x time] x [variables]
%
% subject_grouping  [subj x time]-length integer vector of which
%               observations belong to which subjects, e.g., [1 1 1 ... 2 2 2 ... 3 3 3 ...]'
%
% Optional:
% 'outcome'     followed by outcome data for multivariate prediction.  connectivity
%               values and graph metrics are used to predict outcome data.
%
% 'shift_by'    Followed by integer value for max number of time points to shift
% 'partialr'    Use partial correlation instead of raw correlation
%
% Outputs:
% -------------------------------------------------------------------------
% OUT           A structure containing subject correlation matrices, the mean
%               matrix, and raw and FDR-thresholded group matrix
%               Also contains matrices with [subjects x variables] pairwise
%               correlation elements and graph metrics
%
% Examples:
% -------------------------------------------------------------------------
%
% % Use partial correlations:
% OUT = canlab_connectivity_predict(dat, subject_grouping, 'partialr');
%
% Omit graph met
% OUT = canlab_connectivity_predict(dat, subject_grouping, 'outcome', y, 'nograph');
%
% See also:
% parcel_cl, parcel_cl_nmds_plots, canlab_force_directed_graph,
% canlab_connectivity_preproc

% Programmers' notes:
% Created 2/5/15 by tor wager
%
% TO-DOS:  Thresholds and sig matrix should probably be individualized
%          - plotting methods if you input Clusters/regions

% -------------------------------------------------------------------------
% DEFAULTS AND INPUTS
% -------------------------------------------------------------------------

% Defaults
% -----------------------------------

docluster = 1;
dograph = 1;
doplot = 1;
y = [];             % this is for the outcome data, if any

spath = which('use_spider.m');
if isempty(spath)
    disp('Warning: spider toolbox not found on path; prediction may break')
end

spath = which('all_shortest_paths.m');
if isempty(spath)
    disp('Warning: matlab_bgl toolbox not found on path; Will not run graph measures')
end

% optional inputs with default values
% -----------------------------------

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            
            case {'partialr', 'shift_by'}  % do nothing; these are passed in to xcorr_multisubject
            case 'noplot', doplot = 0;
            case 'nograph', dograph = 0;
            case 'outcome', y = varargin{i+1}; varargin{i+1} = [];
                
            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end


% Cross-correlation stats
% -----------------------------------
s = unique(subject_grouping)';
n = length(s);

% get time points for each subject
for i = s
    t(i) = sum(subject_grouping == s(i));
end

% convert to cells for compatibility with xcorr_multisubject
C = mat2cell(dat, t, size(dat, 2));

% stats on pairwise associations
OUT = xcorr_multisubject(C, varargin{:});

plot_individual_subjects(OUT)

% Image the correlation matrix
%-----------------------------------
newcm = colormap_tor([0 0 .7], [1 .5 0], [1 1 1]);

ylim = OUT.stats.mean(:);
ylim = prctile(abs(ylim), 95);

create_figure('associations', 1, 3);
imagesc(OUT.stats.mean, [-ylim ylim]);
colorbar
set(gca, 'YDir', 'Reverse')
axis tight
colormap(newcm)
title('Group average association matrix');
xlabel('Variable No.');

subplot(1, 3, 2)
tthr = OUT.stats.t .* double(OUT.stats.fdrsig ~= 0);
tthr(isinf(tthr)) = 0; %max(tthr(~isinf(tthr)));
ylim = tthr(:);
ylim = prctile(abs(ylim), 95);
imagesc(tthr, [-ylim ylim]);
colorbar
set(gca, 'YDir', 'Reverse')
axis tight
title('t-values for significant associations');
xlabel('(FDR q < .05)');

% do clusters, if requested
%-----------------------------------

if docluster
    % Cluster variables based on average intercorrelation
    parcelindx = clusterdata(OUT.stats.mean,'linkage','ward','savememory','on', 'maxclust', 12);
    
    subplot(1, 3, 3)
    plot_sorted_correlation_matrix(OUT.stats.mean, parcelindx);
    
    subplot(1, 3, 3)
    ylim = OUT.stats.mean(:);
    ylim = prctile(abs(ylim), 95);
    set(gca, 'Clim', [-ylim ylim]);
    
    OUT.parcelindx = parcelindx;
    OUT.parcelcolors = scn_standard_colors(length(unique(parcelindx)'));

end

% get data for prediction
%-----------------------------------

k = size(OUT.pairwise_assoc{1}, 1);
for i = 1:n
    
    OUT.connectdata(i, :) = squareform(OUT.pairwise_assoc{1}(:, :, i) - eye(k));
    
end

% within- and between-clusters

% graph-theoretic measures
if dograph
    [OUT.betweenness_centrality, OUT.shortestpath, OUT.weighted_degree] = get_graph_metrics(OUT.pairwise_assoc{1}, OUT.stats.fdrsig);
    % bc : betweenness-centrality for [subjects x variables]
    % shortestpath : shortest path for [subjects x variable pairs]
    % wdeg : weighted degree (average correlation) for [subjects x variables]
    
    plot_graph_metrics(OUT);
    
end


% predict, if we have outcome
%-----------------------------------
if ~isempty(y)
    if dograph
        create_figure('prediction', 2, 2);
    else
        create_figure('prediction');
    end
    
    pdat = fmri_data;
    pdat.dat = OUT.connectdata';
    pdat.Y = y;
    [cverr, stats, optout] = predict(pdat, 'algorithm_name', 'cv_lassopcr', 'nfolds', 5, 'error_type', 'mse');
    OUT.PREDICT.pairwise_association = stats;
    
    plot(stats.yfit, stats.Y, 'ko', 'MarkerFaceColor', [.7 .3 0]);
    refline
    title('Prediction from pairwise associations');
    xlabel('Predicted'); ylabel('Observed');
    
    if dograph
        
        pdat.dat = OUT.betweenness_centrality';
        [cverr, stats, optout] = predict(pdat, 'algorithm_name', 'cv_lassopcr', 'nfolds', 5, 'error_type', 'mse');
        OUT.PREDICT.betweenness_centrality = stats;
        subplot(2, 2, 2);
        plot(stats.yfit, stats.Y, 'ko', 'MarkerFaceColor', [.7 .3 0]);
        refline
        title('Prediction from betweenness-centrality');
        xlabel('Predicted'); ylabel('Observed');
        
        pdat.dat = OUT.shortestpath';
        pdat.dat(isinf(pdat.dat)) = max(pdat.dat(~isinf(pdat.dat(:)))) + 1;
        [cverr, stats, optout] = predict(pdat, 'algorithm_name', 'cv_lassopcr', 'nfolds', 5, 'error_type', 'mse');
        OUT.PREDICT.shortestpath = stats;
        subplot(2, 2, 3);
        plot(stats.yfit, stats.Y, 'ko', 'MarkerFaceColor', [.7 .3 0]);
        refline
        title('Prediction from shortest path');
        xlabel('Predicted'); ylabel('Observed');
        
        pdat.dat = OUT.weighted_degree';
        [cverr, stats, optout] = predict(pdat, 'algorithm_name', 'cv_lassopcr', 'nfolds', 5, 'error_type', 'mse');
        OUT.PREDICT.weighted_degree = stats;
        subplot(2, 2, 4);
        plot(stats.yfit, stats.Y, 'ko', 'MarkerFaceColor', [.7 .3 0]);
        refline
        title('Prediction from weighted degree');
        xlabel('Predicted'); ylabel('Observed');
        
    end
    
end


end % function


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%
% Sub-functions
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

function [bc shortestpath wdeg] = get_graph_metrics(C3d, mask)
% [bc D deg] = get_graph_metrics(C3d, mask)
%
% given C3d 3-d matrix of associations and mask matrix for significant
% associations, calculate some graph metrics for each subject
%
% bc : betweenness-centrality for [subjects x variables]
% shortestpath : shortest path for [subjects x variable pairs]
% wdeg : weighted degree (average correlation) for [subjects x variables]

n = size(C3d, 3);
k = size(C3d, 2);

[bc wdeg] = deal(zeros(n, k));
shortestpath = zeros(n, k*(k-1)/2);

fprintf('Getting graph metrics...   ');

for i = 1:n
    
    fprintf('\b\b\b%03d', i);
    
    C = C3d(:, :, i);
    
    C = C .* double(mask ~= 0);
    
    % Enforce format for graph
    C = (C + C') ./ 2;
    C = sparse(C);
    
    % Graph Stats
    % -----------------------------------------------
    bc(i, :) = betweenness_centrality(abs(C))';  % abs if ignoring neg connections
    
    % shortest paths: used as estimate of connectivity
    D = get_path_length(C);
    shortestpath(i, :) = squareform(D);
    
    wdeg(i, :) = full(sum(C));  %full(sum(C ~= 0));
    
end

fprintf('...done.\n');

end

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

function D = get_path_length(r)
% compute matrix of path lengths for n x n matrix of regions, D (connectivity)
%
% r is a correlation matrix
% rset is a cell array of length k, for k sets, with vectors describing the
% indices of members of each set.
%
% e.g., r = region_r{i};

r(isnan(r)) = 0;
r = (r' + r) ./ 2;  % enforce symmetry, just in case


% shortest paths: used to calculate connectivity
rtmp = sparse(r);
rtmp(rtmp < 0) = 0;

% Note: Uses Floyd-Warshall method if 10% non-zero elements or more
D = all_shortest_paths(rtmp);

end % function

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


function plot_sorted_correlation_matrix(r, parcelindx)

nclasses = unique(parcelindx)';

basecolors = scn_standard_colors(length(nclasses));

[classsort, wh] = sort(parcelindx);

rsorted = r(wh, wh);

% put spaces in to mark off classes
%-----------------------------------
%wh_spaces = [0; diff(parcelindx)];
rnew = [];
pinew = [];
for i = nclasses
    rnew = [rnew; rsorted(classsort == i, :)];
    rnew(end + 1, :) = 0;
    
    pinew = [pinew; classsort(classsort == i)];
    pinew(end + 1, 1) = 0;
    
end

rnew2 = [];
for i = nclasses
    rnew2 = [rnew2 rnew(:, classsort == i)];
    rnew2(:, end + 1) = 0;
end
r = rnew2;

% Image the correlation matrix
%-----------------------------------
imagesc(r); %colorbar
set(gca,'YDir', 'Reverse')
axis tight
title('Correlations sorted by cluster');
axis off

% Color bars for class ID
%-----------------------------------
axpos = get(gca, 'Position');
axh = axes('Position', [axpos(1)-.03 axpos(2) .02 axpos(4)]);
set(axh, 'YDir', 'Reverse', 'YLim', [1 length(pinew)]);
hold on;
for i = nclasses
    
    yy = find(pinew == i);
    plot(.7 * ones(size(yy)), yy, '-', 'Color', basecolors{i}, 'LineWidth', 10);
    
    text(-4, round(median(yy)), num2str(i), 'FontSize', 18, 'Color', 'k');
    
end
axis off
drawnow

end



function plot_graph_metrics(OUT)


[h, p] = ttest(OUT.connectdata);
thr = FDR(p, .05);
sig = p < thr;
fprintf('FDR q < .05 sig. pairwise connections: %3.0f out of %3.0f connections.\n', sum(sig), length(sig));

create_figure('graph_metrics', 2, 2);

x = mean(OUT.connectdata);
[xsort, wh] = sort(x(sig), 'descend');

%[h,bins] = hist(x(sig), min(10, round(x(sig)/100)));
plot(xsort, 'LineWidth', 3);
title('Avg. association for significant associations only')


subplot(2, 2, 2);

x = mean(OUT.shortestpath);
[xsort, wh] = sort(x, 'descend');
plot(xsort, 'LineWidth', 3);
title('Avg. shortest path connecting all pairs')


subplot(2, 2, 3);

x = mean(OUT.betweenness_centrality);
[xsort, wh] = sort(x, 'descend');
lineplot_columns(OUT.betweenness_centrality(:, wh));

title('Betweenness-centrality')

subplot(2, 2, 4);

x = mean(OUT.weighted_degree);
[xsort, wh] = sort(x, 'descend');
lineplot_columns(OUT.weighted_degree(:, wh));

title('Weighted degree')


end % function



% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


function plot_individual_subjects(OUT)

m = OUT.pairwise_assoc{1};
m2 = m(:, :, 1);

nc = max(1, floor(size(m2, 2) ./ 5));
nc = zeros(size(m2, 2), nc);

for i = 2:size(m, 3)
m2 = [m2 nc m(:, :, i)];
end

mm = abs(m(:));
clim = prctile(mm(mm ~= 1 & mm ~= 0), 98);

create_figure('individual subjects'); 
imagesc(m2, [-clim clim]); 
colorbar; set(gca, 'YDir', 'Reverse');
axis tight
axis off

%cm = colormap_tor([0 0 1], [1 1 0], [1 1 1]);
cm = colormap_tor([0 0 .7], [1 .5 0], [1 1 1]);
colormap(cm)

end


%     disp('Clustering variables based on full time-series data');
%     disp('This can take a long time...');
%
%     % standardize within-subject data to avoid issues of scale
%     disp('Standardizing subject data');
%     for i = 1:n
%         C{i} = zscore(C{i});
%         cdat = cat(1, C{:});
%     end
%
%     % compress variables if necessary, avoid sparse space
%     disp('Compressing variable space');
%
%     % if many observations, we need to compress observations to make
%     % computationally tractable.
%
%     c = clusterdata(cdat,'linkage','ward','savememory','on', 'maxclust', 4);
