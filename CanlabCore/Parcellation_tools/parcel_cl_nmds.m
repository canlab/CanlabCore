function [parcel_cl_avgs, NMDS, class_clusters] = parcel_cl_nmds(parcel_cl_avgs)
% :Usage:
% ::
%
%     [parcel_cl_avgs, NMDS, class_clusters] = parcel_cl_nmds(parcel_cl_avgs)
%
% :Examples:
% ::
%
%    load Parcellation_info/parcellation.mat
%    [parcel_cl_avgs, NMDS, class_clusters] = parcel_cl_nmds(parcel_cl_avgs)
%
%    parcel_cl_nmds_plots(parcel_cl_avgs, NMDS, 'save')
%    parcel_cl_nmds_plots(parcel_cl_avgs, NMDS, 'save', 'savedir', 'Parcellation_info')
%
% :Complete methods example:

%% 1. Dimension reduction
% 
% Clustering of multivariate data is most stable when the data is not sparse, i.e., the dimensionality is low relative to the number of observations. To limit the dimensionality of the data, a spatio-temporal dimension-reduction step is first performed on the [n x v x N] data matrix of  AUC data for n trials x v voxels x N participants (here, n = 48 trials [usually], v = 17,112 voxels , and N = 26 participants).  A temporal data reduction is first performed to identify components with correlated AUC trial time series within each participant, followed by a spatial reduction to identify components with correlated spatial patterns across subjects.  First, the [n x v] matrix of AUC data for each participant was subjected to PCA, using the [v x v] correlation matrices.  Based on the scree plots across subjects, we saved the first 7 eigenvectors (spatial maps).  These eigenvectors explained 74 +- 2.5% (st. dev. across subjects) of the variance in the full dataset. These eigenvectors were scaled by their variances (eigenvalues) and concatenated across subjects to form an [v x N*7] matrix of eigenvectors.   This matrix was subjected to second (across participant) PCA step to identify components with similar spatial maps across participants.  We retained 12 eigenvectors (maps) based on the scree plot, which explained 65% of the variance across individuals. Component scores in this space were used for clustering.  This is a data reduction step, and the results are not expected to depend strongly on the number of eigenvectors retained at either step, as long as most of the variance in the data is explained.
% 
% 2. Parcellation
% 
% Hierarchical agglomerative clustering with average linkage was used to group voxels into parcels--sets of contiguous voxels with similar profiles--in the [v x 12] matrix of component maps.  The goal of parcellation was to reduce the space from voxels to parcels (regions) for non-metric multidimensional scaling (NMDS)-based clustering of regions, so that NMDS is computationally tractable. Voxels whose trial AUC time series did not correlate with that of other voxels in the same parcel at p < .001 (in a random effects analysis across participants) were pruned from each parcel. Out of 140 parcels total, 127 parcels with more than 3 voxels after pruning were retained for subsequent analysis.  
% 
% 3. NMDS and clustering
% 
% Parcels were now treated as the unit of analysis, and the AUC trial time series data were extracted from each parcel for each subject.  Values in the  [n trials x parcels x N subjects]  data matrix were z-scored within participant to remove inter-subject differences in scaling, then concatenated into an [n*N x parcels] data matrix.  Correlations among parcels were converted into a [parcels x parcels] matrix of distances using the formula distance = (1 - r) / 2. The NMDS stress plot (which operates on ranked distances and is therefore more robust to outliers and does not require a strictly Euclidean distance space) was examined and 7 dimensions (which explained 79% of the variance in distances) were retained for the final cluster analysis.  
% Hierarchical agglomerative clustering with average linkage was used to cluster the parcels in this space into interconnected networks with correlated AUC trial time series. To choose the number of clusters in the final solution (k), for every possible choice of clusters between  k = 2 and 20, we compared the cluster solution to the average and standard deviation of 1,000 clustering iterations with permuted parcel time series.  This yielded a Z-score ([actual solution - mean permuted solution] /  standard deviation of permuted solution for each value of k. The best solution was k = 13, with a value of Z = 5.57 compared to the null-hypothesis single-cluster solution (p < .0001).
%
% ..
%    Documentation not complete. please update me.
%    Tor Wager, Oct 2008
% ..


    disp('Getting networks of parcels') % Get networks of these parcels

    clear data
    N = length(parcel_cl_avgs(1).timeseries);
    for i = 1:length(parcel_cl_avgs), for j = 1:N, data{j}(:,i) = parcel_cl_avgs(i).timeseries{j}; end, end

    % Correlate like this:
    % ---------------------------------------
    % NMDS = xcorr_multisubject(data);
    % NMDS.stats.D = (1 - NMDS.stats.mean) ./ 2;

    % OR like this, if we want to be closer to data and don't need inferences across ss on correlations:
    % ---------------------------------------
    data = [];
    for i = 1:length(parcel_cl_avgs)

        % zscore within, or at least center to avoid individual diffs in
        % baseline driving correlations
        for s = 1:N
            parcel_cl_avgs(i).timeseries{s} = scale(parcel_cl_avgs(i).timeseries{s});
        end

        data(:, i) = cat(1, parcel_cl_avgs(i).timeseries{:});
    end

    desired_alpha = .01;
    fprintf(['Saving FDR-corrected connections at q < %3.3f corrected\n' desired_alpha]);
    
    NMDS.parcel_data = data;
    [NMDS.stats.c, NMDS.stats.pvals] = corrcoef(data);
    [NMDS.stats.fdrthr, NMDS.stats.fdrsig] = fdr_correct_pvals(NMDS.stats.pvals, NMDS.stats.c, desired_alpha);

    NMDS.stats.D = (1 - NMDS.stats.c) ./ 2;

    % Pick number of dimensions
    maxclasses = round(length(parcel_cl_avgs) ./ 10);
    maxclasses = min(maxclasses, 20);
    maxclasses = max(maxclasses, 5);

    disp('Choosing number of dimensions')
    [NMDS.stats_mds.GroupSpace,NMDS.stats_mds.obs,NMDS.stats_mds.implied_dissim] = shepardplot(NMDS.stats.D, maxclasses);

    % Clustering with permutation test to get number of clusters

    disp(['Clustering regions in k-D space: 2 to ' num2str(maxclasses) ' classes.'])

    NMDS.stats_mds = nmdsfig_tools('cluster_solution',NMDS.stats_mds, NMDS.stats_mds.GroupSpace, 2:maxclasses, 1000, []);


    % Set Colors
    basecolors = {[1 0 0] [0 1 0] [0 0 1] [1 1 0] [0 1 1] [1 0 1] ...
        [1 .5 0] [.5 1 0] [.5 0 1] [1 0 .5] [0 1 .5] [0 .5 1]};

    NMDS.basecolors = basecolors;
       
%     % Make figure
%     disp('Visualizing results')
%     create_figure('nmdsfig')

%     nmdsfig(NMDS.stats_mds.GroupSpace,'classes',NMDS.stats_mds.ClusterSolution.classes, ...
%         'names', [],'sig',NMDS.stats.fdrsig, 'fill', 'colors', basecolors);

    % disp('Saving NMDS structure with networks in parcellation.mat')
    % save(fullfile(mysavedir, 'parcellation.mat'), '-append', 'class_clusters', 'parcel*')

    % re-define class clusters
    clear class_clusters
    for i = 1:max(NMDS.stats_mds.ClusterSolution.classes)
        wh = find(NMDS.stats_mds.ClusterSolution.classes == i);
        class_clusters{i} = parcel_cl_avgs(wh);

        % refine class (network) membership
        [parcel_cl_avgs(wh).from_class] = deal(i);
    end

    % disp('Saving NMDS structure and final class clusters in parcellation.mat')
    % save(fullfile(mysavedir, 'parcellation.mat'), '-append', 'NMDS', 'class_clusters')


end


function [pthr,sig] = fdr_correct_pvals(p, r, desired_alpha)

    psq = p; psq(find(eye(size(p,1)))) = 0;
    psq = squareform(psq);
    pthr = FDR(p, desired_alpha);
    if isempty(pthr), pthr = 0; end

    sig = sign(r) .* (p < pthr);

    sig(isnan(sig)) = 0;
end




