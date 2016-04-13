function parcel_clusters(clpos_data, clneg_data)
% :Usage:
% ::
%
%     parcel_clusters(clpos_data, clneg_data)
%
% No outputs. Saves all output in separate directory.
%
% First, get eigenvectors for each subject.
%
% We're interested obtaining PARCELS of voxels that tend to co-activate, 
% or have the same activation profile.
%
% We can find these by using clustering algorithms to group voxels with
% similar profiles.
%
% Because we have a many voxel x many voxel covariance matrix for each
% subject (lots of data!), it's important to reduce the dimensionality of
% the problem and peform clustering on a REDUCED_DIMENSIONAL space. 
%
% We use PCA to do this.  Instead of clustering activation profiles (e.g.,
% time-courses) directly, we cluster eigenvector loadings for each voxel on
% a reduced set of components that explains most of the variance in the
% data.
%
% Similar voxels will have similar loadings across the set of 
% eigenvectors.  e.g., two voxels may load high on components [1 3 and 5],
% and low on components [10 and 13].  If they have the same pattern of
% loadings, they should be considered part of the same CLASS.  Groups
% of voxels that are contiguous in space and are members of the same CLASS 
% are called parcels.
%
% Images and outputs are saved in their own subdirectory called
% Parcellation_info
%
% The main outputs are: 
% parcel_cl             % parcels, one cell per subject, one parcel per
%                       element within cells. same format as clpos_data
%
% parcel_cl_avgs        % parcels, one parcel per element within cells. 
%                         same format as clpos_data2
%
% parcel_cl_avgs(x).timeseries contains one cell per subject, with data
% averaged across voxels within that parcel for that subject
% This kind of output is useful, because you can input it directly into
% other mediation analyses.
% [paths, stats2] = mediation(SETUP.data.X, SETUP.data.Y, parcel_cl_avgs(1).timeseries, 'plots', 'verbose', 'names', {'Hi-Low Cue' 'Pain Report' 'Parcel'}, 'boot');
% cluster_orthviews(parcel_cl_avgs(1), {[0 1 0]}, 'add');
% ::
%
%     cd('/Volumes/SCNAlpha/Data_and_Tools/SpeechTask/analysis/wb_multisubject_correl_HR_corrected/mediation_Xprepvsb_Mbrain_Yhr')
%     load cl_b_fdr05_002_01_k3_1_1_prune
%
% then run
%
% There are 2 dimension-reduction steps:
%   1. within-subjects
%   2. is on eigenvectors concatenated across subjects

initial_eigval_limit = 15;      % Number of eigenvalues to save initially for each subject
                                % Need this to reduce computational burden
       
n_eigs = 7;                     % Number of eigenvectors to use per subject in the clustering algorithm.  

% We form a matrix of [voxels x subjects*eigenvectors] (called group_eigenvectors)
% For example, 7 eigenvectors x 10 subjects means clustering will be done
% in a 70 dimensional space.
% This matrix is subjected to a second step of PCA data reduction.  A
% smaller number of scores are saved, and these scores represent sets of
% eigenvector loadings in the original 70 dimensional space.


n_eigs_across = 12;             % Number of dimensions to save out of the original subjects*eigenvectors set
                                % Clustering is done on component scores.
                                % These are canonical eigenvector patterns
                                % that capture regular variations across
                                % subjects. Think of them as compressed
                                % eigenvectors.
                                
% More n_eigs_across  means that clustering will be done in a more complex space.
% Increasing this will tend to split the voxels into smaller parcels, and
% decreasing it will tend to lump them into larger parcels.

n_class_range = [2:20];       % test solutions from this many to this many CLASSES

nclasses = 20;

% The final number of classes you want to use in clustering.
% Increasing this will tend to split the voxels into smaller parcels, and
% decreasing it will tend to lump them into larger parcels.

do_nclasses_search = 0;         % search over n_class_range?

doprune = 1;                    % prune voxels and parcels that are too small or whose voxels don't inter-correlate

mysavedir = 'Parcellation_info';
if ~exist(mysavedir, 'dir'), mkdir(mysavedir), end


create_figure('eigenvalues');

nsubjects = length(clpos_data);
fprintf('Subjects: %3.0f\n', nsubjects);

if ~isempty(clneg_data)
    test_subj_dat = [cat(2, clpos_data{1}(:).all_data) cat(2, clneg_data{1}(:).all_data)];
else
    test_subj_dat = [cat(2, clpos_data{1}(:).all_data)];
end

nvox = size(test_subj_dat, 2);
fprintf('Voxels in mask area: %3.0f\n', nvox);


group_eigenvalues = zeros(nsubjects, initial_eigval_limit);

group_scores = cell(1, nsubjects);


fprintf('PCA: Subject ')

%Clustering of multivariate data is most stable when the data is not sparse, 
% i.e., the dimensionality is low relative to the number of observations. 
% To limit the dimensionality of the data, a spatio-temporal dimension-reduction 
% step is first performed on the [n x v x N] data matrix of  AUC data for 
% n trials x v voxels x N participants.  A temporal data reduction is first performed 
% to identify components with correlated AUC trial time series within each participant, 
% followed by a spatial reduction to identify components with correlated spatial patterns 
% across subjects.  First, the [n x v] matrix of AUC data for each participant was 
% subjected to PCA, using the [v x v] correlation matrices.  Based on the scree plots 
% across subjects, we saved the first [n_eigs] eigenvectors.  
% These eigenvectors explained xxx +- xxx% (st. dev. across subjects) of the variance 
% in the full dataset. These eigenvectors were scaled by their variances (eigenvalues) 
% and concatenated across subjects to form an [v x N*7] matrix of eigenvectors, 
% where N=27 subjects.   This matrix was subjected to another (spatial) PCA step to 
% identify components with similar spatial maps across participants.  
% We retained [n_eigs_across] eigenvectors based on the scree plot, which explained xx% of the 
% variance across individuals. This is a data reduction step, and the results are not 
% expected to depend strongly on the number of eigenvectors retained at either step, 
% as long as most of the variance in the data is explained.

% Temporal reduction
% -------------------------------------------------------------------------
for s = 1:nsubjects

    fprintf('%3.0f', s)

    if ~isempty(clneg_data)
        subj_dat = [cat(2, clpos_data{s}(:).all_data) cat(2, clneg_data{s}(:).all_data)];
    else
        subj_dat = [cat(2, clpos_data{s}(:).all_data)];
    end

    nanvec = any(isnan(subj_dat)) | all(subj_dat == 0);
    wh_bad = find(nanvec);
    if any(wh_bad)
        fprintf('Warning! Subject %3.0f s has missing data (0 or NaN) for these voxels: ', s)
        fprintf('%3.0f ', wh_bad);
        fprintf('\n')

        subj_dat(:, wh_bad) = [];
    end

    % to do PCA on correlation matrix rather than cov
    subj_dat = zscore(subj_dat);

    %[U, eigenvalues, eigenvectors] = svd(subj_dat, 'econ'); % almost, but
    %scaling isn't right, so just use princomp, which does it all

    [eigenvectors, score, eigenvalues] = princomp(subj_dat, 'econ');

    clear subj_dat

    % insert bad voxels back in
    if any(wh_bad)
        for i = 1:size(eigenvectors, 2)
        ev(:, i) = naninsert(nanvec, eigenvectors(:, i));
        end
        
        ev(isnan(ev)) = 0;
        eigenvectors = ev;
    end

    % The first [initial_eigval_limit] 
    group_eigenvalues(s, :) = eigenvalues(1:initial_eigval_limit)';

    % we want to scale the eigenvectors by their variances (eigenvalues),
    % so that components that account for more variation in the data are weighted more heavily.
    eigenvectors = eigenvectors(:, 1:initial_eigval_limit) * diag(eigenvalues(1:initial_eigval_limit));

    group_eigenvectors{s} = eigenvectors;

    plot(group_eigenvalues(s, :), 'ko-');
    drawnow
    
    % Calculate variance explained for first [initial_eigval_limit]
    % eigenvalues.
    ev = cumsum(group_eigenvalues');
    ev = ev ./ repmat(sum(group_eigenvalues', 1), size(ev, 1), 1);
    evn = ev(n_eigs, :);  % explained variance for these eigs; approximate as it is only proportion of first n eigs
    fprintf('Temporal reduction: Eigenvectors explain %3.2f%% +- %3.2f%% of variance.', 100*mean(evn), 100*std(evn));
    
end

clear eigenvalues eigenvectors


%%

%n_eigs = input('Enter number of eigenvectors to save: ');

plot_vertical_line(n_eigs);
drawnow

if ~exist(mysavedir, 'dir'), mkdir(mysavedir); end

saveas(gcf, fullfile(mysavedir, 'Eigenvalues'), 'png');


for i = 1:nsubjects
    group_eigenvectors{s}(:, n_eigs + 1 : end) = [];
end

group_eigenvectors = cat(2, group_eigenvectors{:});

%%

%input('Enter range of clusters: ');
%niter = 100;

names = [];

% c = nmdsfig_tools('cluster_solution', [], group_eigenvectors, n_class_range, niter, names);

%% 2nd Dimension redution step : to get scores in lower-dim space for
% clustering

% Spatial reduction
% -------------------------------------------------------------------------

[across_eigenvectors, across_scores, across_eigenvalues] = princomp(group_eigenvectors, 'econ');

create_figure('group_eigenvalues', 1, 3);
plot(across_eigenvalues(1:initial_eigval_limit), 'ko-');
title('Across subjects eigenvalues');

% Calculate variance explained for first [n_eigs_across] eigenvalues.
ev = cumsum(across_eigenvalues);
ev = ev ./ sum(across_eigenvalues);
evn = ev(n_eigs_across, :);  % explained variance for these eigs; approximate as it is only proportion of first n eigs
fprintf('Spatial reduction: Eigenvectors explain %3.2f%% of variance.', 100*evn);

%n_eigs_across = input('Enter number of eigenvectors to save: ');

scores_to_cluster = across_scores(:, 1:n_eigs_across);

plot_vertical_line(n_eigs_across);

subplot(1, 3, 2)
imagesc(scores_to_cluster);

drawnow
%% CLUSTER voxels

classes = [];
sil_vals = {};
mean_sil = [];

if do_nclasses_search

    fprintf('Clustering : ')

    % problem with silhouette is that it really finds natural break-point in
    % data; so looks good for 2 clusters with pos/neg groups in data usually.

    clear s

    for i = n_class_range

        fprintf('%3.0f ', i);

        classes = clusterdata(scores_to_cluster, 'linkage', 'average', 'maxclust', i);
        sil_vals{i} = silhouette(scores_to_cluster, classes);
        mean_sil(i) = mean(sil_vals{i});

    end

    fprintf('\n');

    subplot(1, 3, 3)
    plot(n_class_range, mean_sil(n_class_range), 'ko-', 'MarkerFaceColor', [.2 .6 1], 'LineWidth', 2);

    saveas(gcf, fullfile(mysavedir, 'Group_eigs_and_clustering'), 'png');

end


cl = [clpos_data clneg_data];

disp(['Saving data file: ' mysavedir filesep 'parcellation.mat'])
save(fullfile(mysavedir, 'parcellation'), 'cl', 'n*', '*eig*', '*score*', 'classes', '*sil*')


%% % Get parcels of contiguous regions
% Save averages over voxels for each subject, within each region

fprintf('Getting parcels and associated data: Clustering with %3.0f classes\n', nclasses)

classes = clusterdata(scores_to_cluster, 'linkage', 'average', 'maxclust', nclasses);

clear parcel_cl
disp('Getting contiguous regions: These become parcels');

fprintf('Subject ');
for s = 1:nsubjects

    fprintf('%3.0f', s)

    if ~isempty(clneg_data)
        cl = [clpos_data{s} clneg_data{s}];
    else
        cl = clpos_data{s};
    end

    CLU = clusters2CLU(cl);

    for i = 1:max(classes)

        my_cl = CLU;
        my_cl.XYZmm = my_cl.XYZmm(:, classes == i);
        my_cl.XYZ = my_cl.XYZ(:, classes == i);
        my_cl.Z = my_cl.Z(:, classes == i);

        my_cl.all_data = my_cl.all_data(:, classes == i);

        class_clusters{i} = tor_extract_rois([], my_cl, my_cl);

        [class_clusters{i}.from_class] = deal(i);
    end

    parcel_cl{s} = cat(2, class_clusters{:});

end
fprintf('\n')

% Prune parcels here
if doprune
    [parcel_cl, meanpval, meancor] = prune_parcels(parcel_cl);
    
    disp('Saving meanpval and meancor for pruned parcels in parcellation.mat')
    save(fullfile(mysavedir, 'parcellation.mat'), '-append', 'meanpval', 'meancor')
end


parcel_cl_avgs = parcel_cl{1};

% another convenient format
for i = 1:length(parcel_cl_avgs)
    parcel_cl_avgs(i).all_data = cell(1, nsubjects);
    parcel_cl_avgs(i).timeseries = cell(1, nsubjects);
    for s = 1:nsubjects
        parcel_cl_avgs(i).all_data{s} = parcel_cl{s}(i).all_data;
        parcel_cl_avgs(i).timeseries{s} = parcel_cl{s}(i).timeseries;
    end
end
    
fprintf('\n')

disp('Saving parcel_cl and parcel_cl_avgs in parcellation.mat')
save(fullfile(mysavedir, 'parcellation.mat'), '-append', 'class_clusters', 'parcel*')

%volInfo = iimg_read_img('mask.img', 2);

%% RE-do networks on these parcels
% re-define class clusters
% refine class (network) membership
% ---------------------------
[parcel_cl_avgs, NMDS, class_clusters] = parcel_cl_nmds(parcel_cl_avgs);

disp('Saving NMDS structure and final class clusters in parcellation.mat')
save(fullfile(mysavedir, 'parcellation.mat'), '-append', 'NMDS', 'class_clusters')

%% Plots
% ---------------------------
% Plot: data panel
% Orthviews of parcels
% Montages of parcels
parcel_cl_nmds_plots(parcel_cl_avgs, NMDS, 'save', 'savedir', 'Parcellation_info')

end



%% Sub-functions


function [parcel_cl, meanpval, meancor] = prune_parcels(parcel_cl)

    nsubjects = length(parcel_cl);
    
    %% remove parcels with too few voxels
    vcutoff = 3;
    p_cutoff = .002;

    whomit = false(size(parcel_cl{1}));
    nvox = cat(1, parcel_cl{1}.numVox);
    whomit(nvox < vcutoff) = 1;
    fprintf('Eliminated %3.0f parcels smaller than %3.0f voxels\n', sum(whomit), vcutoff)
    fprintf('Keeping %3.0f parcels\n', sum(~whomit))
    for s = 1:nsubjects
        parcel_cl{s}(whomit) = [];
    end

    nparcels = length(parcel_cl{s});
    [meancor, meanpval, vox_removed] = deal(zeros(nparcels, 1));
    omit_parcels = false(nparcels, 1);

    %
    % Get data across subjects for one parcel
    for p = 1:nparcels

        fprintf('Parcel %3.0f ', p);

        % Get data across subjects for one parcel (p)
        % ---------------------------------------------
        dat = cell(nsubjects, 1);
        for s = 1:nsubjects
            dat{s} = parcel_cl{s}(p).all_data;
        end
        %dat = cat(1, dat{:});

        % Correlate, and get stats
        % ---------------------------------------------
        %[c, pvals] = corrcoef(dat);
        clear c
        for i = 1:nsubjects
            c(:, :, i) = corrcoef(dat{i});
        end

        % Note: p vals will not be correct if there are NaNs!
        
        mc = nanmean(c, 3);
        se = nanstd(c, 0, 3) ./ sqrt(nsubjects);

        mc = mc .* (1 - eye(size(mc)));
        mc = squareform(mc);

        se = se .* (1 - eye(size(se)));
        se = squareform(se);

        t = mc ./ se;
        pvals = 2 .* (1 - tcdf(abs(t), nsubjects - 1));

        meanpval(p) = mean(pvals);
        meancor(p) = mean(mc);

        pvals = squareform(pvals);

        % Omit the bad (unrelated) voxels
        % ---------------------------------------------
        pv = mean(pvals, 2);

        wh = pv > p_cutoff;

        vox_removed(p) = sum(wh);

        if sum(wh) > length(wh) - vcutoff + 1  
            % we have zero or 1 valid voxels left; omit the whole parcel
            omit_parcels(p) = 1;
        elseif any(wh)
            % omit bad voxels

            for s = 1:nsubjects
                parcel_cl{s}(p).XYZmm(:, wh) = [];
                parcel_cl{s}(p).XYZ(:, wh) = [];
                parcel_cl{s}(p).all_data(:, wh) = [];

                parcel_cl{s}(p).timeseries = nanmean(parcel_cl{s}(p).all_data, 2);
                parcel_cl{s}(p).numVox = sum(~wh);

            end

            pvals(wh, :) = [];
            pvals(:, wh) = [];
            mc = squareform(mc);
            mc(wh, :) = [];
            mc(:, wh) = [];

            meanpval(p) = mean(squareform(pvals));
            meancor(p) = mean(squareform(mc));
        end



        fprintf(' mean inter-voxel corr: %3.2f, mean p = %3.4f, vox excluded = %3.0f', meancor(p), meanpval(p), vox_removed(p))
        if omit_parcels(p), fprintf(' OMITTED'); end
        fprintf('\n')

    end

    % remove whole bad parcels
    for s = 1:nsubjects
        parcel_cl{s}(omit_parcels) = [];
    end

    meanpval(omit_parcels) = [];
    meancor(omit_parcels) = [];
    
    % remove parcels that are now too small
    whomit = false(size(parcel_cl{1}));
    nvox = cat(1, parcel_cl{1}.numVox);
    whomit(nvox < vcutoff) = 1;
    fprintf('Eliminated %3.0f parcels smaller than %3.0f voxels\n', sum(whomit), vcutoff)
    fprintf('Keeping %3.0f parcels\n', sum(~whomit))
    for s = 1:nsubjects
        parcel_cl{s}(whomit) = [];
    end

    nvox = cat(1, parcel_cl{1}.numVox);
    create_figure('Pruned parcels', 2, 1);
    plot(nvox, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', [.5 .5 1]);
    title('Number of voxels in each parcel');
    xlabel('Parcel index number')
    plot_horizontal_line(vcutoff)

    subplot(2, 1, 2)
    plot(meancor, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', [.5 .5 1]);
    title('Mean within-subject correlation value with other voxels in parcel')
    drawnow

end

