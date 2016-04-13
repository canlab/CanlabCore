function parcel_cl_nmds_plots(parcel_cl_avgs, NMDS, varargin)
% Plot: data panel
%
% case 'save', dosave = 1;
% case {'savedir', 'mysavedir'}, mysavedir = varargin{i+1};
% case {'figs', 'dotfigs'}, savedotfigs = 1;
%
% Doc not complete yet.  Please update me!
%
% See parcel_clusters.m
% See parcel_cl_nmds.m
%
% load Parcellation_info/parcellation.mat
%
% :Examples:
% ::
%
%    parcel_cl_nmds_plots(parcel_cl_avgs, NMDS, 'save')
%    parcel_cl_nmds_plots(parcel_cl_avgs, NMDS, 'save', 'savedir', 'Parcellation_info')
%    parcel_cl_nmds_plots(parcel_cl_avgs, NMDS, 'save', 'savedir', 'Parcellation_info_tor_mask_try1', 'savedotfig')


    % ..
    %    SETUP
    % ..

    dosave = 0;
    mysavedir = [];
    savedotfigs = 0;

    for i = 1:length(varargin)
        if ischar(varargin{i})
            switch varargin{i}

                % functional commands
                case 'save', dosave = 1;

                case {'savedir', 'mysavedir'}, mysavedir = varargin{i+1};

                case {'figs', 'dotfigs'}, savedotfigs = 1;

                otherwise, warning(['Unknown input string option:' varargin{i}]);
            end
        end
    end

    N = length(parcel_cl_avgs(1).timeseries);

    if ~isfield(NMDS, 'basecolors')
        NMDS.basecolors = {[1 0 0] [0 1 0] [0 0 1] [1 1 0] [0 1 1] [1 0 1] ...
            [1 .5 0] [.5 1 0] [.5 0 1] [1 0 .5] [0 1 .5] [0 .5 1]};
    end

    basecolors = NMDS.basecolors;

    % NMDS plot
    % --------------------------------------------------
    % Make figure
    disp('Visualizing results')
    create_figure('nmdsfig');

    nmdsfig(NMDS.stats_mds.GroupSpace,'classes',NMDS.stats_mds.ClusterSolution.classes, ...
        'names', [],'sig',NMDS.stats.fdrsig, 'fill', 'colors', basecolors);
    % --------------------------------------------------

    % --------------------------------------------------
    % Get class clusters (convenient format)
    % --------------------------------------------------

    % re-define class clusters
    clear class_clusters
    for i = 1:max(NMDS.stats_mds.ClusterSolution.classes)
        wh = find(NMDS.stats_mds.ClusterSolution.classes == i);
        class_clusters{i} = parcel_cl_avgs(wh);

        % refine class (network) membership
        [parcel_cl_avgs(wh).from_class] = deal(i);
    end


    nclasses = length(class_clusters); % redefine based on results

    % --------------------------------------------------
    % get index of which parcels (clusters) are in which classes
    % --------------------------------------------------

    parcelindx = cat(1, parcel_cl_avgs.from_class);
    [parcelindx, ivals] = sort(parcelindx);
    cl = parcel_cl_avgs(ivals);

    % --------------------------------------------------
    % Get data
    % --------------------------------------------------

    data = [];
    for i = 1:length(cl)

        % zscore within, or at least center to avoid individual diffs in
        % baseline driving correlations
        for s = 1:N
            cl(i).timeseries{s} = scale(cl(i).timeseries{s});
        end

        data(:, i) = cat(1, cl(i).timeseries{:});
    end

    for s = 1:N-1, len(s) = length(cl(1).timeseries{s}); end
    len = cumsum([1 len]);

    % --------------------------------------------------
    % Panel TS figure
    % --------------------------------------------------

    create_figure('Panel_timeseries', nclasses, 1);
    for i = 1:nclasses
        subplot(nclasses, 1, i);
        hold on;
        mycolor = basecolors{i};

        [lh, fh] = plot_error(data(:, parcelindx == i)');
        mymean = mean(data(:, parcelindx == i)');

        set(lh, 'Color', mycolor);
        set(fh, 'FaceColor', mycolor);
        set(gca, 'YLim', [min(mymean) max(mymean)]);
        axis tight
        hold on
        plot_vertical_line(len);

        if i == 1
            for j = 1:length(len), text(len(j), max(mymean), ['Subj' num2str(j)]); end
        end

        axis off
        drawnow
    end

    % --------------------------------------------------
    % Correlation figure
    % --------------------------------------------------


    % get within-parcel correlations
    %-----------------------------------
    for i = 1:length(cl)
        dd = cat(1, cl(i).all_data{:});
        rr = corrcoef(dd);
        rr = rr .* (1 - eye(size(rr)));
        within_corr(i) = mean(squareform(rr));
    end

    create_figure('Correlations');
    [r, p] = corrcoef(data);

    % put withins on diagonal
    r = r .* (1 - eye(size(r)));
    r = r + diag(within_corr);
    %r(p > .05) = NaN;

    % put spaces in to mark off classes
    %-----------------------------------
    %wh_spaces = [0; diff(parcelindx)];
    rnew = [];
    pinew = [];
    for i = 1:nclasses
        rnew = [rnew; r(parcelindx == i, :)];
        rnew(end + 1, :) = NaN;

        pinew = [pinew; parcelindx(parcelindx == i)];
        pinew(end+1) = NaN;

    end

    rnew2 = [];
    for i = 1:nclasses
        rnew2 = [rnew2 rnew(:, parcelindx == i)];
        rnew2(:, end + 1) = NaN;
    end
    r = rnew2;

    % Image the correlation matrix
    %-----------------------------------
    imagesc(r); colorbar
    set(gca,'YDir', 'Reverse')
    axis tight
    newcm = colormap_tor([0 0 .7], [1 1 0], [1 .5 0]);
    colormap(newcm)
    title('Correlations among parcels');
    xlabel('Parcel index number');

    % Color bars for class ID
    %-----------------------------------
    axpos = get(gca, 'Position');
    axh = axes('Position', [.05 axpos(2) .04 axpos(4)]);
    set(axh, 'YDir', 'Reverse', 'YLim', [1 length(pinew)]);
    hold on;
    for i = 1:nclasses

        yy = find(pinew == i);
        plot(.7 * ones(size(yy)), yy, '-', 'Color', basecolors{i}, 'LineWidth', 10);

        text(-2, round(median(yy)), num2str(i), 'FontSize', 18, 'Color', 'k');

    end
    axis off
    drawnow


    % --------------------------------------------------
    % Orthviews of parcels
    % --------------------------------------------------

    nclasses = length(class_clusters); % redefine based on results


    while length(basecolors) < nclasses
        basecolors{end+1} = basecolors{end - 11} + basecolors{end - 6} ./ 2;
    end

    parcelindx = cat(1, parcel_cl_avgs.from_class);

    colors = nmdsfig_tools('class_colors',parcelindx,basecolors);

    addstr = 'noadd';

    for i = 1:length(parcel_cl_avgs)

        cluster_orthviews(parcel_cl_avgs(i), colors(i), addstr, 'solid');

        addstr = 'add';
    end

    % --------------------------------------------------
    % Montages of parcels
    % --------------------------------------------------

    disp(['Saving montages in ' mysavedir])

    for i = 1:length(class_clusters)

        mycolor = basecolors(i);

        montage_clusters([], class_clusters{i}, mycolor);
        drawnow

        savename = fullfile(mysavedir, ['Montage_class' num2str(i)]);
        saveas(gcf, savename,'png')

    end

    if dosave
        save_figures(1, savedotfigs, mysavedir);
    end
    
    % --------------------------------------------------
    % Glass brains
    % --------------------------------------------------
    for i = 1:length(class_clusters)
        create_figure('glass'); 
        
        mysig = NMDS.stats.fdrsig;
        mysig = mysig(:, parcelindx == i);
        mysig = mysig(parcelindx == i, :);
        
        cluster_nmdsfig_glassbrain( ...
            class_clusters{i}, ...
            ones(length(class_clusters{i}), 1), ...
            NMDS.basecolors(i), ...
            mysig, [], 'blobs');

        if dosave
           saveas(gcf, fullfile(mysavedir, ['glass_network_' num2str(i)]), 'png');
           if savedotfigs
               saveas(gcf, fullfile(mysavedir, ['glass_network_' num2str(i)]), 'fig');
           end
        end

    end



end



% -------------------------------------------------------------------------
% Save all figures
% -------------------------------------------------------------------------
function save_figures(verbose, savedotfigs, mysavedir)

    if nargin == 0, verbose = 1; end
    if nargin == 1, savedotfigs = 0; end
    if nargin == 2 || isempty(mysavedir), mysavedir = '.'; end

    fignames = {'nmdsfig' 'Panel_timeseries' 'Correlations'};

    for i = 1:length(fignames)
        name = fignames{i};

        h = findobj('Tag', name);
        if ~isempty(h)

            if length(h) > 1
                if verbose, disp(['Warning: More than one plot with tag ' name '. Saving highest-numbered figure.']); end
                h = max(h);
            end

            figure(h)
            scn_export_papersetup(500);

            if savedotfigs
                if verbose
                    fprintf('Saving: %s%s\n', fullfile(mysavedir, name), '.fig');
                end

                saveas(h, fullfile(mysavedir, name), 'fig');
            end

            if verbose
                fprintf('Saving: %s%s\n', fullfile(mysavedir, name), '.png');
            end
            saveas(h, fullfile(mysavedir, name), 'png');
        end
    end

end
