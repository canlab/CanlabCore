function [clusters, subcl] = cluster_barplot(P, clusters, varargin)
% :Usage:
% ::
%
%    [clusters, subcl] = cluster_barplot(P, clusters, varargin)
%
% this function not only plots, but separates subclusters using pca / cluster_princomp.m
% based on pattern across all conditions, covariance (not correlation),
%
% :Inputs:
%
%   **P:**
%        cell array of strings containing image file names (data extracted from these)
%
%   **clusters:**
%
%   **[opt] - 'subclusters':**
%        to get sub-clustering based on pca and clustering of voxels
%
%   **[opt] - cell array:**
%        of strings with condition names
%
%   **[opt] - 'split':**
%        value is 1: to split into 2 plots (first half and last half of P)
%
%        value is 2: to plot individual subjects over bars
%
%   **[opt] - 'center':**
%        center parameter values in plot (subtract row means)
%        this gives closer to correct "within subjects" error bars
%        and may be used when the overall parameter values have no meaning
%
%
%   **[opt] - 'covs':**
%        followed by between-subject covariates (e.g., behavioral regressors)
%         plots remove these before plotting means and std. errors
%
%   **[opt] - 'max':**
%        to make plots based on max z-values within region for each dataset P
%        not compatible with 'split' and 'center' (ignores these commands)
%        right now, special for inhib - see also inhib2_cluster_barplot (good function)
%
%   **[opt] - 'indiv':**
%        to threshold based on individual t-statistics
%
%        FOLLOW with cell array of t-images -- usually, there will be one cell, with images
%        for each subject in rows, to define voxels for each ss.
%
%        BUT Tnames can be the same length as
%        contrast images, one t-img per subject per contrast, if
%        desired.
%
% :Outputs:
%
%   **clusters:**
%        clusters struture, with BARPLOT substructure added
%
%   **subcl:**
%        substructure contains data extracted and image file names
%
% This program uses XYZmm millimeter coordinates in clusters to find voxels
% So clusters and data files may have different dimensions.
%
% :Examples:
% ::
%
%    cluster_barplot(EXPT.SNPM.P(4:6), clusters(2:3))
%    cluster_barplot(EXPT.SNPM.P(7:12), clusters(2:3), {'ObjE' 'AttE' 'InteractE' 'ObjI' 'AttI' 'InteractI'}, 'split')
%    [clusters, subclusters] = cluster_barplot(EXPT.SNPM.P(17:24), clusters, 'subclusters', 'split')
%    RS2_8vs2_placeboCP = cluster_barplot(EXPT.SNPM.P([8 10 12 14
%    16]), RS2meta, 'indiv', T);
%
% :See Also: mask2clusters.m, a simpler version that just extracts clusters from a mask file.
%
% ..
%    by Tor Wager, 5/15/03
%    modified 3/19/04 by Tor to add individual subject plots
% ..

    % ..
    %    Set up optional inputs
    % ..
    splitat = 0;
    dosubcl = 0;
    docenter = 0;
    covs = [];
    domax = 0;
    doindiv = 0;

    for j = 1:length(P), cnames{j} = num2str(j); end

    for i = 1:length(varargin),
        if iscell(varargin{i}), cnames = varargin{i};
        elseif strcmp('split', varargin{i}), if varargin{i+1}==2, splitat=2; else splitat=round(length(P)./2); end
        elseif strcmp('subclusters', varargin{i}), dosubcl = 1; disp('Finding sub-clusters')
        elseif strcmp('center', varargin{i}), docenter = 1;
        elseif strcmp('covs', varargin{i}), covs = varargin{i+1};
        elseif strcmp('max', varargin{i}), domax = 1;
        elseif strcmp('indiv', varargin{i}), doindiv = 1;
            Tnames = varargin{i+1}; varargin{i+1} = 'saved.';
        end
    end

    if doindiv
        % replicate names, in case only one set is entered
        if ~iscell(Tnames)
            Tnamestmp = Tnames;
            Tnames={};
            Tnames{1}=Tnamestmp;
        end
        while length(Tnames) < length(cnames)
            Tnames{end+1} = Tnames{1};
        end
    end

    %if ~isempty(strmatch('split', varargin)), splitat=round(length(P)./2);, else splitat=length(P);, end

    % -------------------------------------------------
    % * load image files
    % -------------------------------------------------
    fprintf(1, 'Loading images.')
    CLU = clusters2CLU(clusters);
    dat = cell(1, length(clusters));
    alldat = cell(1, length(clusters));


    V = spm_vol(P{1}(1, :));
    if any(any(V.mat - CLU.M))
        disp('CLUSTERS and IMAGE DIMENSIONS DIFFERENT!  Resizing clusters.'),
        %CLU = transform_coordinates(CLU, V.mat);
        V.M = V.mat;
        for i = 1:length(clusters),
            clusters(i).XYZ = mm2voxel(clusters(i).XYZmm, V)'; clusters(i).M=V.mat;
            clusters(i).voxSize = diag(V.M(1:3, 1:3))';
            clusters(i).XYZmm = voxel2mm(clusters(i).XYZ, clusters(i).M);
            clusters(i).Z = ones(1, size(clusters(i).XYZ, 2));
        end
        CLU = clusters2CLU(clusters);
        for i = 1:length(clusters)
            clusters(i).u=CLU.u;
            clusters(i).VOX=CLU.VOX;
        end
        doslow = 1;
    else
        doslow = 0;
    end

    % check if any contiguous clusters we lose
    %tmpCL = tor_extract_rois([], CLU, CLU);
    %if length(tmpCL) < length(clusters), doslow = 1;, end
    %clear tmpCL

    for i = 1:length(P)

        if ~doslow
            % fast way: if all voxel sizes match
            cl{i} = tor_extract_rois(P{i}, clusters);  %CLU, CLU);
        end

        % concatenated data for each cluster
        for j = 1:length(clusters)

            if doslow
                % slow way: to keep clusters separate that otherwise touch
                cl{i}(j) = tor_extract_rois(P{i}, clusters(j), CLU);
            end

            % alldat{cluster} is data x voxels, data is subj within conditions strung
            % together
            alldat{j} = [alldat{j}; cl{i}(j).all_data];

            if length(cl{i}(j).timeseries) < size(dat{j}, 1),
                cl{i}(j).timeseries = [cl{i}(j).timeseries; NaN*zeros(size(dat{j}, 1)-length(cl{i}(j).timeseries), 1)];
            end

            dat{j} = [dat{j} cl{i}(j).timeseries];
        end

        %V{i} = spm_vol(P{i});
        %vols{i} = spm_read_vols(V{i});
    end

    %V{1}(1).M = V{1}(1).mat;
    %imgdims = size(vols{1});
    %mask = zeros(imgdims(1:3));


    % -------------------------------------------------
    % * subclusters
    % -------------------------------------------------
    fprintf(1, 'Plotting clusters.')
    for i = 1:length(clusters)

        if dosubcl && clusters(i).numVox > 4
            % -------------------------------------------------
            % * subclusters
            % -------------------------------------------------

            [subcl{i}] = getsubcl(clusters(i), alldat{i}, length(P), size(P{1}, 1), cnames, splitat, i, docenter, covs);  % does subclusters within it
            if length(subcl{i}) > 1,
                subcluster_montage(subcl{i});    % now plot all subcluster locations
                saveas(gcf, ['cl' num2str(i) '_subcl_montage'], 'fig')
                saveas(gcf, ['cl' num2str(i) '_subcl_montage'], 'tif')
            end
        else

            % no subclusters

            ste = nanstd(dat{i}) ./ sqrt(size(dat{i}, 1));

            if isfield(clusters, 'shorttitle')
                mytit = clusters(i).shorttitle;
            else
                mytit = clusters(i).title;
            end

            str = sprintf('%s %2.0f (%3.0f, %3.0f, %3.0f), %3.0f voxels', ...
                mytit, i, clusters(i).mm_center(1), ...
                clusters(i).mm_center(2), clusters(i).mm_center(3), clusters(i).numVox);

            if domax
                % set up
                len = size(clusters(1).BARPLOT.data, 1);
                tmpcl = []; tmpcl2 = [];
                for j = 1:round(size(alldat{1}, 1) ./ len)
                    tmpcl{j}(i) = clusters(i); tmpcl{j}(i).all_data = alldat{i}(1+(j-1)*len:j*len, :);
                    tmpcl{j}(i).timeseries = mean(tmpcl{j}(i).all_data, 2);

                    [ste, tmp] = ste(tmpcl{j}(i).all_data);
                    clusters(i).Z(j, :) = spm_t2z(tmp, size(tmpcl{j}(i).all_data, 1)-1);

                    wh = find(tmp == max(tmp)); wh = wh(1);
                    clusters(i).BARPLOT.grp_peakdat(j, :) = tmpcl{j}(i).all_data(:, wh);

                    %tmpcl2{j}(i) = cluster_ttest(tmpcl{j}(i), covs);  % add covariate here if necessary
                end
                % needs checking and/or debugging
                %inhib2_barplot(clusters(i), tmpcl2{1}(i), tmpcl2{2}(i), tmpcl2{3}(i), 'common');
                keyboard
            end



            f = makefigure(dat{i}, str, cnames, splitat, docenter, covs);
            try
                saveas(gcf, ['cl' num2str(i) '_bar'], 'fig')
                saveas(gcf, ['cl' num2str(i) '_bar'], 'tif')
            catch
                disp('error saving figure.');
            end

            clusters(i).BARPLOT.files = P;
            clusters(i).BARPLOT.data = dat{i};
            clusters(i).BARPLOT.all_data = alldat{i};

            if dosubcl
                disp('Fewer than 5 voxels; skipping subclustering.')

                % make all field names match.
                % won't work if this is the first one!
                subcl{i} = subcl{1}(1);
                N = fieldnames(subcl{i});
                for NN = 1:length(N)
                    eval(['subcl{i}(1).' N{NN} ' = [];'])
                end

                N = fieldnames(clusters(i));
                for NN = 1:length(N)
                    eval(['subcl{i}(1).' N{NN} ' = clusters(i).' N{NN}])
                end
            end
            close all
        end
        clusters(i).covs = covs;
    end


    fprintf(1, '\n')

    % extract spatial peak data and locations
    for i = 1:length(cl)
        cltmp = extract_ind_peak([], cl{i});
        for j = 1:length(cltmp)
            clusters(j).BARPLOT.peakdata{i} = cltmp(j).peakdata;
            clusters(j).BARPLOT.peakXYZ{i} = cltmp(j).peakXYZ;
            clusters(j).BARPLOT.peakXYZmm{i} = cltmp(j).peakXYZmm;
        end
    end


    % extract individual significant region center of mass
    if doindiv
        disp('getting individual significance regions')
        indpeak = [];

        for c = 1:size(Tnames, 2)    % for each contrast tested
            for i = 1:size(Tnames{1}, 1) % for each subject
                % pass in clusters with all data saved for contrast c, Tnames
                % for contrast c; timeseries saves data from this

                [tmp] = cluster_tmask(cl{c}, Tnames{c}(i, :), i, alldat);
                cl{c} = tmp;
            end

            % put stuff where it belongs
            for j = 1:length(clusters)
                clusters(j).BARPLOT.indiv_data(:, c) = cl{c}(j).timeseries;
                clusters(j).BARPLOT.indivXYZ{c} = cl{c}(j).INDIV.center;
                clusters(j).BARPLOT.indivXYZmm{c} = cl{c}(j).INDIV.mm_center;
                clusters(j).BARPLOT.indiv_sig{c} = cl{c}(j).INDIV.sigt;
                clusters(j).BARPLOT.indivXYZall{c} = cl{c}(j).INDIV.XYZ;
            end
        end
    end
end



% -----------------------------------------------------------------------------
% * sub-functions
% -----------------------------------------------------------------------------


function [subclusters] = getsubcl(clusters, adat, numimglists, numimgs, varargin)
    % always pass in only one cluster in clusters!! (vec of length 1)
    % if > 1 output, subclusters is created and plotted here
    % if subclusters, pass in cnames and splitat

    if length(varargin) > 0, cnames = varargin{1}; splitat = varargin{2};
        clnum=varargin{3}; docenter = varargin{4};
    end
    if length(varargin) > 4, covs = varargin{5}; else covs = []; end

    % don't scale if this is across subjects, because voxels
    % with more variation across ss SHOULD drive the pca
    %adat{j} = scale(adat{j});   % scale so high-noise voxels
    %                            % don't drive the pca

    % separate subclusters
    subclusters.BARPLOT = [];   % add for consistency of field names
    clusters.all_data = adat;

    [clusters, subclusters] = cluster_princomp(clusters, []);   %  , 0, 0, 1);

    subclusters = subclusters{1};

    % re-format averages within sub-regions for plotting

    for j = 1:length(subclusters)
        vind = 1;
        if size(subclusters(j).all_data, 2) > 1,
            tmp = nanmean(subclusters(j).all_data')';
        else
            tmp = subclusters(j).all_data;
        end

        for k = 1:numimglists
            subclusters(j).BARPLOT.dat(:, vind) = ...
                tmp((k-1)*numimgs+1:k*numimgs);
            vind = vind + 1;
        end

        % now plot it
        ste = nanstd(subclusters(j).BARPLOT.dat) ./ sqrt(size(subclusters(j).BARPLOT.dat, 1));
        str = sprintf('Cl %3.0f (%4.0f vox) subcl at (%3.0f, %3.0f, %3.0f), %3.0f voxels', ...
            clnum, clusters.numVox, subclusters(j).mm_center(1), ...
            subclusters(j).mm_center(2), subclusters(j).mm_center(3), subclusters(j).numVox);

        f = makefigure(subclusters(j).BARPLOT.dat, str, cnames, splitat, docenter, covs);
        saveas(gcf, ['cl' num2str(clnum) '_subc' num2str(j) '_bar'], 'fig')
        saveas(gcf, ['cl' num2str(clnum) '_subc' num2str(j) '_bar'], 'tif')
    end
end



function f = makefigure(dat, varargin)

    if length(varargin) > 1, cnames=varargin{2}; else cnames = []; end
    if length(varargin) > 2, splitat=varargin{3}; else splitat = 0; end
    if length(varargin) > 3, docenter=varargin{4}; else docenter = 0; end
    if length(varargin) > 4, covs = varargin{5}; covs(:, end+1) = 1; else covs = []; end
    if covs == 1, covs = []; end    % if empty, make it empty!

    dat1 = dat;     % save original data for ind subj plot

    if ~isempty(covs)
        covs(:, 1:end-1) = scale(covs(:, 1:end-1));
        for i = 1:size(dat, 2),
            b = pinv(covs) * dat(:, i);
            r = dat(:, i) - covs * b + b(end);
            dat(:, i) = r;
        end
    end

    if docenter
        disp('Centering rows of data for plotting.')
        dat = dat - repmat(nanmean(dat')', 1, size(dat, 2));

    end

    doprintt = 1;

    if doprintt
        % print t values
        for i = 1:size(dat, 2)
            [h, p, ci, stats] = ttest(dat(:, i));
            nums{i} = sprintf('%3.2f', stats.tstat);
        end
    end

    ste = nanstd(dat) ./ sqrt(size(dat, 1));
    dat = nanmean(dat);
    mysum = sign(dat) .* (ste + abs(dat) + .15 .* abs(dat));

    if splitat > 2
        f = figure('Color', 'w');

        h(1) = subplot(1, 2, 1); set(gca, 'FontSize', 16); hold on; grid on;
        bar(dat(1:splitat)); tor_bar_steplot(dat(1:splitat), ste(1:splitat), {'b'});
        set(gca, 'XTick', 1:splitat)
        xlabel('Conditions', 'FontSize', 18), ylabel('fMRI Signal', 'FontSize', 18)
        if length(varargin) > 1, set(gca, 'XTickLabel', varargin{2}(1:splitat)), end
        if length(varargin) > 0, title(varargin{1}, 'FontSize', 20), end

        for i = 1:splitat, text(i-.5, mysum(i), nums{i}, 'Color', 'k', 'FontWeight', 'b', 'FontSize', 14); end


        h(2) = subplot(1, 2, 2); set(gca, 'FontSize', 16); hold on; grid on;
        bar(dat(splitat+1:end)); tor_bar_steplot(dat(splitat+1:end), ste(splitat+1:end), {'b'});
        set(gca, 'XTick', 1:length(dat)-splitat)
        xlabel('Conditions', 'FontSize', 18)
        if length(varargin) > 1, set(gca, 'XTickLabel', varargin{2}(splitat+1:end)), end

        for i = 1:length(dat)-splitat, text(i, mysum(splitat+i), nums{splitat+i}, 'Color', 'k', 'FontWeight', 'b', 'FontSize', 14); end

        equalize_axes(h);
        set(gcf, 'Position', [195         308        1259         674]), drawnow

    elseif splitat == 2

        % ------------------------
        % this plot does individual subject estimates
        % ------------------------
        barplot_columns(dat1)



        % add reg lines, if covs and sig
        for i = 1:size(dat1, 2)
            [B, dev, stat]=glmfit(covs(:, 1), dat1(:, i)); tmp = corrcoef(covs(:, 1), dat1(:, i));
            fprintf(1, 'Cond. %3.0f: Bo: t = %3.2f, se = %3.4f, p = %3.4f  B1: r = %3.2f, t = %3.2f, se = %3.4f, p = %3.4f\n', ...
                i, stat.t(1), stat.se(1), stat.p(1), tmp(1, 2), stat.t(2), stat.se(2), stat.p(2));

        end

        if ~isempty(cnames), set(gca, 'XTickLabel', cnames), end

    else

        % ------------------------
        % standard bar plot
        % ------------------------

        f = figure('Color', 'w'); set(gca, 'FontSize', 16); %hold on; grid on;
        h = bar(dat); set(h, 'FaceColor', [.7 .7 .7])
        tor_bar_steplot(dat, ste, {'k'});

        doinhib2 = 0;
        if doinhib2
            % special insert for inhib2 (triple inhibition)
            set(gca, 'XTickLabel', {'GNG' 'Flanker' 'SRC' 'Saccade'})
        end

        dointext = 0;
        if dointext
            % special insert for intext
            cla
            xp = dat; h = bar([xp(1:3);xp(4:6)], 'grouped'); cm = [0 1 0;1 0 0; 0 0 1];colormap(cm)
            xe = ste;
            tor_bar_steplot(xp(1:3), xe(1:3), {'k'}, .55, .225)
            tor_bar_steplot(xp(4:6), xe(4:6), {'k'}, 1.55, .225)
            set(gca, 'FontSize', 16, 'XTickLabel', {'External' 'Internal'});legend(h, {'Object switching' 'Attribute Switching' 'Interaction'})
            ylabel('BOLD Contrast')
        end

        %set(gca, 'XTick', 1:size(dat, 2))
        %xlabel('Conditions', 'FontSize', 18), ylabel('fMRI Signal', 'FontSize', 18)
        %if length(varargin) > 1, set(gca, 'XTickLabel', varargin{2}), end
        if length(varargin) > 0, title(varargin{1}, 'FontSize', 20), end

        set(gcf, 'Position', [464   283   930   827]);
        drawnow
    end
end





