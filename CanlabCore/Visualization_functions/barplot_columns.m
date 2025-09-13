function [handles, dat, xdat, statstable, order_idx, pairwise_results] = barplot_columns(dat, varargin)
% :Usage:
% ::
%
%    [graphics handles, adjusted data, x-data, statstable, order_idx] = barplot_columns(dat, [other optional arguments])
%
% This function makes a barplot of columns of data, with standard error
% bars.  Optional arguments include removing continuous covariates before plotting,
% robust (IRLS) estimation of means and correlations with covariates, and
% within-subject error bars based on the subject x condition interaction
% (overall), which is not quite the standard error contrasts of interest,
% but is the standard error for a 1-way repeated measures ANOVA.
%
% Very flexible control with many options for estimation and display:
% - Plot type options (line, bar, violin / distribution)
% - Point plotting, individual lines, and outlier options
% - Robustness options (IRLS)
% - Error bar options (within-subject, 95% CIs)
% - Display axis, title, widths, and color control
% - Covariate options (option to add)
%
% E.g.:  (and see options below)
% - default is violin plot (distribution) plus bars
% - optional: lines, not bars
% - optional:  plots circles around points at z >= 1.96
% - optional: add case numbers or labels instead of points
%
% Input formats (dat input argument):
% - Matrix or cell vector (good for unequal numbers of observations)
% - if dat is a cell array, each entry becomes one "bar".  Useful if n
% observations is different for each column.
% - table object can be entered. This function will use Variable Names in table.
%
% :Examples: Just plot means and SE
% ::
%
%    h = barplot_columns(tmp,'Cluster 1',[],1);
%
% :Optional arguments:
%
%   1.
%
%   2-k. String Arguments - in any order
%
%   Covariates
%        - 'covs' : Followed by matrix of covariates to remove by regression
%                   Note: covs are mean-centered
%        - 'wh_reg' : Which regressor to leave in and sort data points by
%                     OR enter 0 to remove all covariates before plotting.
%
%   Figure control
%        - 'nofig' : do not make figure
%        - 'skipallplots': skip making figure and all plots; table only
%
%   Plot type options
%        - 'line' : Make line plot instead of bar plot
%        - 'violin': add violin plot to each bar, with data points
%        - 'noviolin': suppress default violin plot behavior
%        - 'nobars' : suppress bars (e.g., for violin plots only)
%
%   Point plotting, individual lines, and outlier options
%        - 'dolines' : plot lines showing individual effects
%        - 'noind' : do not plot individual scores
%        - 'plotout': circle potential outliers at z>1.96 in red
%        - 'removeout': Remove outliers at z>1.96 from statistical computations. Not compatible with 'plotout' since that data will be removed from plotting.
%        - 'number' : plot case numbers instead of points
%        - 'MarkerSize' : followed by marker size
%        - 'MarkerAlpha' : followed by marker transparency value (alpha)
%        - 'stars', 'dostars' : plot stars for significance above each column (default)
%           Sig. levels are coded as * = p<.05, ** = p<.01, *** = q<.05 FDR across variables tested in this plot
%        - 'nostars' : do not plot stars
%        - 'custom_p' : Custom cell-array to determine p-values if manually generating plots. Don't abuse this.
%        - 'pairwisetest' : Conducts pairwise statistical comparisons between all columns.
%                   Performs uncorrected t-tests (either paired or unpaired depending on data shape):
%                   - Paired t-tests are used when both columns have the same number of non-NaN elements.
%                   - Unpaired t-tests are used when column lengths differ or missing data prevents alignment.
%                   Effect sizes (Cohen's d), t-statistics, and p-values are printed to the console and
%                   returned in a table (`pairwise_results`) in the base workspace.
%                   Significant results (p < .05) are visualized on the plot using bracket connectors
%                   with stars (*, **, ***) placed above the bar pairs.
%
%   Robustness options
%        - 'dorob' : do robust IRLS means and correlations
%
%   Error bar options
%        - 'within' : within-subjects standard errors, followed by contrast
%                   matrix
%        - '95CI'   : error bars are 95% CI instead of SE
%        - 'custom_se': followed by a singleton vector of errors to use,
%               one per bar
%
%   Display axis, title, widths, and color control
%        - 'order': 'ascend' or 'descend' by the mean value
%        - 'x' : followed by x-axis values for bars
%        - 'noxlim' : Suppress automatic setting of x-limit (e.g., when adding to existing plot)
%        - 'color' : followed by color for bars (text: 'r' or [r g b]) OR
%               cell array with names of colors cell for each line/bar
%        - 'title' : followed by title for figure (or empty for no title)
%        - 'width' : followed by bar width
%
% :Examples:
% ::
%
%    barplot_columns(ctmp, 'title', 'RT effects by Switch Type',overall_sw,'nofig','dorob')
%
% Standard Errors ARE NOT Adjusted for covariate, right now.
%
% Example: within-subjects std. errors
% ::
%
%    barplot_columns(dat, 'title', 'Means', [], 'nofig', 'within', c);
%
% The example below uses color, width, and xposition arguments to make a grouped
% ::
%
%    barplot showing effects for two groups:
%    exp_dat = EXPT.error_rates(EXPT.group==1,:);
%    control_dat = EXPT.error_rates(EXPT.group==-1,:);
%    barplot_columns(exp_dat, 'title', 'Error rates', [], 'nofig', 'noind', 'color', 'r','width', .4);
%    barplot_columns(control_dat, 'title', 'Error rates', [], 'nofig', 'noind', 'color', 'b','width', .4, 'x', (1:9)+.5);
%    set(gca, 'XLim', [0 10], 'XTick', 1:9)
%
%    barplot_columns(nps_by_study, 'title', 'NPS by study', [], 'doind', 'colors', mycolors, 'nofig');
%
%    create_figure('example_plots', 1, 4);
%
%    Y{:,1} = rand(20,1);
%    Y{:,2} = rand(100,1);
%
%    [h, L, MX, MED, bw, F, U] = violinplot(Y,'facecolor',[1 .5 0; 0 .5 1],'edgecolor',[1 .5 0; 0 .5 1].*.75,'mc', [1 .5 0].*.5, 'x', [1 3], 'medc', []);
%    title('Violinplot.m', 'FontSize', 16)
%
%    subplot(1, 4, 2)
%    barplot_columns(Y, 'nofig')
%    title('barplot\_columns.m default', 'FontSize', 16)
%
%    subplot(1, 4, 3)
%    barplot_columns(Y, 'nofig', 'violin', 'colors', {[1 .5 0] [0 .5 1]})
%    title('barplot\_columns.m colored', 'FontSize', 16)
%
%    subplot(1, 4, 4)
%
%    Y{:,1} = randn(50,1) + 5;
%    Y{:,2} = Y{1} + .3 * randn(50,1) + 3;
%
%    barplot_columns(Y, 'nofig', 'noviolin', 'colors', {[1 .5 0] [0 .5 1]}, 'dolines')
%    title('barplot\_columns.m parallel coords', 'FontSize', 16)
%
% Covariate(s): Remove regressors (covs) "group":
% barplot_columns(mydata, figtitle, 'colors', DAT.colors, 'dolines', 'nofig', 'names', DAT.conditions, 'covs', group, 'wh_reg', 0);
%
% Covariate(s): Leave in regressor (cov) # 1 "group" and sort points by its values:
% barplot_columns(mydata, figtitle, 'colors', DAT.colors, 'dolines', 'nofig', 'names', DAT.conditions, 'covs', group, 'wh_reg', 1);
%
% See also: lineplot_columns, barplot_colored, line_plot_multisubject, violinplot
%
%   Edited 06/12/2024 Added 'order' option. - Michael Sun, PhD
%   Edited 08/22/2024 Added 'custom_p' option. - Michael Sun, PhD
%   Edited 06/10/2025 Added 'pairwisetest' option. - Michael Sun, PhD

% ..
%    Defaults
% ..

dashes = '---------------------------------------------';
dofig = 1;
doind = 1;
removeout = 0;
plotout = 0;
dorob = 0;
xdat = [];
dolines = 0;
dobars = 1;
dowithin = 0;
custom_se = false;
custom_p = false;
donumber = 0;
dojitter = 1; % jitter is for numbers only
mycolor = [.8 .8 .8];
myalpha = 1;
barwidth = .8;
dolineplot = 0;
do95CI = 0;
do95within = 0;
order = 'none';
dopairwise = false;
pairwise_results = [];
order_idx = [];

nanwarningflag = 1;
doviolin = 1;
mytitle = [];
covs = [];
doxlim = 1;
names = {};
wh_reg = 1; % regressor of interest - 0 for "no regressor of interest", remove all
mymarkersize = 6;
dostars = true;
handles = [];
doprinttable = 1;
skipallplots = false;


% ----------------------------------------------------
% > handle table input - save names
% ----------------------------------------------------
if isa(dat, 'table')
    names = format_strings_for_legend(dat.Properties.VariableNames);
    dat = table2array(dat);
end

% ----------------------------------------------------
% > handle cell input - concatenate and pad with NaN
% ----------------------------------------------------

dat = enforce_padded_matrix_form(dat);

xvals = 1:size(dat, 2);

% ----------------------------------------------------
% > Set up additional input arguments
% ----------------------------------------------------

if length(varargin) > 0
    for i = 1:length(varargin)

        if iscell(varargin{i}), continue, end

        % Figure control
        if strcmp(varargin{i},'nofig') || strcmp(varargin{i},'nofigure'), dofig = 0;  end

        if strcmp(varargin{i},'skipallplots'), skipallplots = true;  end


        % Plot type options
        if strcmp(varargin{i},'line'), dolineplot = 1;  end
        if strcmp(varargin{i},'nobars'), dobars = 0;  end
        if strcmp(varargin{i}, 'violin')
            doviolin = 1;
            doind = 0;  % ind is already done in violin
        end

        if strcmp(varargin{i}, 'noviolin') || strcmp(varargin{i}, 'noviolins')
            doviolin = 0;
        end

        % Point plotting, individual lines, and outlier options
        if strcmp(varargin{i},'noind'), doind = 0;  end
        if strcmp(varargin{i},'removeout'), removeout = 1;  end
        if strcmp(varargin{i},'plotout'), plotout = 1;  end
        if strcmp(varargin{i},'dolines'), dolines = 1;  end
        if strcmp(varargin{i},'number'), donumber = 1;  end
        if strcmp(varargin{i},'MarkerSize') || strcmp(varargin{i},'markersize')
            mymarkersize = varargin{i + 1}; varargin{i + 1} = [];
        end
        if strcmpi(varargin{i},'MarkerAlpha') || strcmpi(varargin{i},'MarkerFaceAlpha')
            myalpha = varargin{i + 1}; varargin{i + 1} = [];
        end

        if strcmp(varargin{i}, 'stars') || strcmp(varargin{i}, 'dostars'), dostars = true; end
        if strcmp(varargin{i}, 'nostars'), dostars = false; end
        if strcmp(varargin{i}, 'custom_p'), PValues = varargin{i+1}; custom_p = true; end
        if strcmpi(varargin{i}, 'pairwisetest'), dopairwise = true; end



        % Robustness options
        if strcmp(varargin{i},'dorob'), dorob = 1;  end

        % Error bar options
        if strcmp(varargin{i}, 'within'), dowithin = 1; end %cons = varargin{i + 1}; end
        if strcmp(varargin{i}, '95within'), dowithin = 1; do95within = 1; end
        if strcmp(varargin{i}, 'custom_se'), within_ste = varargin{i+1}; custom_se = true; end
        if strcmp(varargin{i}, '95CI'), do95CI = 1; end

        % Display axis, title, widths, and color control
        if strcmp(varargin{i}, 'x'), xvals = varargin{i + 1}; end
        if strcmp(varargin{i}, 'color') || strcmp(varargin{i}, 'colors')
            mycolor = varargin{i + 1};
            varargin{i + 1} = [];
        end
        if strcmp(varargin{i}, 'width'), barwidth = varargin{i + 1}; end
        if strcmp(varargin{i}, 'title'), mytitle = varargin{i + 1}; end

        if strcmp(varargin{i}, 'noxlim'), doxlim = 0;  end

        % Labels
        if strcmp(varargin{i}, 'names'), names = varargin{i + 1}; varargin{i + 1} = []; end

        % Order
        if strcmp(varargin{i}, 'order'), order = varargin{i + 1}; varargin{i + 1} = []; end

        % Covariate options
        if strcmp(varargin{i}, 'covs')
            covs = varargin{i + 1};
            if ~isempty(covs), covs = scale(covs,1); end
        end

        % verbose options
        if strcmpi(varargin{i}, 'notable'), doprinttable = 0; end


        if strcmp(varargin{i}, 'wh_reg')
            wh_reg = varargin{i + 1};
        end

    end % for
end % varargin


% ----------------------------------------------------
% > Set Order
% ----------------------------------------------------
if strcmp(order, 'ascend') || strcmp(order, 'descend')
    [~, order_idx] = sort(nanmean(dat, 1), order);
    dat = dat(:, order_idx);

    names=names(order_idx);
    if iscell(mycolor) & numel(mycolor)==numel(order_idx)
        mycolor=mycolor(order_idx);
    end
    if custom_se == true
        within_ste=within_ste(order_idx);
    end
end

% ----------------------------------------------------
% > Build design matrix X for controlling for covariates
% ----------------------------------------------------

dat = double(dat);

% find all-NaN columns, replace with zeros
wh = find(all(isnan(dat), 1));
dat(:, wh) = 0;

% get final design matrix, intercept is last column
[nn, ny] = size(dat);
k = size(covs, 2);

% add intercept as last column
if ~isempty(covs)
    % Remove intercept if we added one manually
    wh_oldintercept = find(all(diff(covs) < eps));
    covs(:, wh_oldintercept) = [];
end

X = [covs ones(nn,1)];
wh_intercept = k + 1;

% ----------------------------------------------------
% > Get means and standard error of means
% With robust option, if specified, and removing
% covariates, if there are any.
% ----------------------------------------------------

Std_Error = [];

% key vars are :
% Mean_Value, Std_Error, T, P (for stars), Cohens_d (for table)

for i = 1:ny


    if removeout
        % Remove outliers
        if plotout
            warning('Cannot use plotout if removeout is set.');
        end
    
        z = (dat(:,i) - nanmean(dat(:,i))) ./ nanstd(dat(:,i));
        wh = find(abs(z) >= 1.96);
        
        dat(wh,i)=NaN;

        
    end




    if doprinttable
        if ~isempty(names) && length(names) >= i
            fprintf(1,'Col %3.0f: %s\t', i, names{i});
        else
            fprintf(1,'Column %3.0f:\t', i);
        end
    end

    % ----------------------------------------------------
    % > Get [robust or non-robust] mean and standard error
    %   Return y, data from column, with nans in
    % ----------------------------------------------------

    % remove nans from this column - dat and design (X) including covs if entered
    tmpy = dat(:, i);
    tmpx = X;

    [wasnan, tmpx, tmpy] = nanremove(tmpx, tmpy);
    all_nan_index{i} = wasnan;

    n(i) = size(tmpy,1);

    % get mean and standard error of intercept (robust or OLS)
    % y is adjusted for all non-intercept covs
    % stats has weights, stats.w, which are all 1 for OLS

    [x, newy, r, p, Std_Error(i, 1), Mean_Value(i, 1), stats] = partialcor(tmpx, tmpy, wh_intercept, 0, dorob);
    T(i, 1) = stats.t(wh_intercept);
    P(i, 1) = stats.p(wh_intercept);

    % very low P-values: do not use exactly zero, because of problems
    P(i, 1) = max(P(i, 1), 10 * eps);

    Cohens_d(i, 1) = stats.t(wh_intercept) ./ (size(x, 1) .^ .5);  % mean(tmpy) ./ std(tmpy), but adjusts for covs

    % Convert to 95% CI, if requested
    if do95CI
        Std_Error(i) = Std_Error(i) * tinv(.975, nn-1);
    end

    y(:, i) = naninsert(wasnan, newy);

    %%%not needed y(:,i) = y(:,i) + Mean_Value(i);   % add mean
    myweights(:,i) = naninsert(wasnan, stats.w);

    % If we have a covariate of interest specified by wh_reg, re-calculate
    % adjusted y to leave in this covariate. This will be used to sort and
    % plot y values.
    % wh_reg should be 0 for all covs removed, or a number to select which
    % regressor to leave in and sort by
    % ----------------------------------------------------
    % > Use partialcor to remove covariates if requested
    %   Return y, adjusted y-values
    % ---------------------------------------------------wid-

    if ~isempty(covs) && wh_reg

        % if we have covs, leave in cov. of interest (cov1)
        % y is adjusted for all non-intercept covs
        [x, y(:, i)] = partialcor(tmpx, tmpy, wh_reg, 0, dorob);

    end

    if dopairwise
        fprintf('\n%s\nPairwise t-tests (uncorrected)\n%s\n', dashes, dashes);
        combs = nchoosek(1:ny, 2);
        pairwise_p = NaN(size(combs, 1), 1);
        pairwise_t = NaN(size(combs, 1), 1);
        pairwise_d = NaN(size(combs, 1), 1);
        pairwise_labels = strings(size(combs, 1), 1);
    
        for i = 1:size(combs, 1)
            c1 = combs(i, 1);
            c2 = combs(i, 2);
            y1 = dat(:, c1);
            y2 = dat(:, c2);
        
            % Remove NaNs
            [~, y1c, y2c] = nanremove(y1, y2);
        
            % Determine whether to do paired or unpaired test
            is_paired = isequal(length(y1c), length(y2c));
        
            if is_paired
                [~, p, ~, st] = ttest(y1c, y2c);
                d = mean(y1c - y2c) / std(y1c - y2c);  % Cohen's d for paired
            else
                [~, p, ~, st] = ttest2(y1, y2);  % use full original vectors w/ nans
                pooled_sd = sqrt(((std(y1, 'omitnan')^2 + std(y2, 'omitnan')^2) / 2));
                d = (nanmean(y1) - nanmean(y2)) / pooled_sd;  % Cohen's d for unpaired
            end
        
            pairwise_p(i) = p;
            pairwise_t(i) = st.tstat;
            pairwise_d(i) = d;
        
            label1 = sprintf('Col %d', c1);
            label2 = sprintf('Col %d', c2);
            if ~isempty(names)
                label1 = names{c1};
                label2 = names{c2};
            end
            pairwise_labels(i) = sprintf('%s vs. %s', label1, label2);
        
            if is_paired
                test_type = 'paired';
            else
                test_type = 'unpaired';
            end
            
            fprintf('%-30s t = % 5.2f, p = % .4f, d = % 5.2f (%s)\n', ...
                pairwise_labels(i), st.tstat, p, d, test_type);
        end
    
        pairwise_results = table(pairwise_labels, pairwise_t, pairwise_p, pairwise_d, ...
            'VariableNames', {'Comparison', 'T', 'P', 'Cohens_d'});
    
        % Optionally: assign this to base workspace or attach to statstable
        assignin('base', 'pairwise_results', pairwise_results);
    end

end

if custom_p
    if iscell(PValues)
        P=[PValues{:}]';
    else
        P=PValues';
    end
    % very low P-values: do not use exactly zero, because of problems
    P(find(P < 10 * eps)) = 10 * eps;
end


dat = y;  % adjusted data, for plot

if dowithin && ~custom_se

    within_ste = barplot_get_within_ste(dat);
    if do95within
        within_ste = within_ste * tinv(.975, nn-1);
    end

    Std_Error = repmat(within_ste, 1, size(dat, 2));

elseif custom_se

    Std_Error = within_ste;

end

% ----------------------------------------------------
% > Print Table
% ----------------------------------------------------
if doprinttable
    fprintf(1, '\n%s\nTests of column means against zero\n%s\n', dashes, dashes);
end
Name = names';
if ~iscolumn(Name), Name = Name'; end

if ~iscolumn(Std_Error), Std_Error = Std_Error'; end % can happen for within-subject err

if isempty(Name), for i = 1:length(T), Name{i, 1} = sprintf('Col %3.0f', i); end, end

statstable = table(Name, Mean_Value, Std_Error, T, P, Cohens_d);
if do95CI || do95within, statstable.Properties.VariableNames{3} = 'half_CI_95perc'; end
if doprinttable
    disp(statstable)
end

% ----------------------------------------------------
% > Make figure
% ----------------------------------------------------
if skipallplots
    return
end

if dofig
    handles.fig_han = create_figure('barplot');
else
    handles.fig_han = get(gcf);

end

handles.axis_han = gca;
set(handles.axis_han, 'FontSize',18);
hold on

% ----------------------------------------------------
% > BARPLOT (or line plot)
% ----------------------------------------------------

if dolineplot

    if iscell(mycolor)
       for i = 1:length(xvals)
        handles.line_han{i} = plot(xvals, Mean_Value, 'o-', 'Color', mycolor{i}, 'MarkerFaceColor', mycolor{i}, 'MarkerSize', 8);
        handles.errorbar_han{i} = errorbar(xvals, Mean_Value, Std_Error, Std_Error);
        set(handles.errorbar_han{i}, 'LineWidth', 2, 'Color', mycolor{i});
       end


    else

        handles.line_han = plot(xvals, Mean_Value, 'o-', 'Color', mycolor, 'MarkerFaceColor', mycolor, 'MarkerSize', 8);
        handles.errorbar_han = errorbar(xvals, Mean_Value, Std_Error, Std_Error);
        set(handles.errorbar_han, 'LineWidth', 2, 'Color', mycolor);
    end

elseif dobars

    handles.bar_han1 = bar(xvals, Mean_Value, barwidth);

    if iscell(mycolor)
        % each bar a different color

        for i = 1:length(xvals)
            handles.bar_han{i} = bar(xvals(i), Mean_Value(i), barwidth, 'FaceColor', mycolor{i});

            handles.errorbar_han{i} = errorbar(xvals(i), Mean_Value(i), Std_Error(i), 'Color', mycolor{i} ./ 2, 'LineWidth', 3);

        end

    else
        %all bars the same color

        set(handles.bar_han1, 'FaceColor', mycolor); %,'LineWidth',2)

        for i = 1:length(xvals)
            handles.errorbar_han{i} = errorbar(xvals(i), Mean_Value(i), Std_Error(i), 'Color', mycolor ./ 2, 'LineWidth', 3);
        end

    end

end

% add violin if entered
if doviolin

    if iscell(mycolor)
        mycolor = cat(1, mycolor{:});
    elseif size(mycolor, 1) < ny
        mycolor = repmat(mycolor, ny, 1);
    end

    Y = enforce_cell_array(y);

    % Do not plot indiv points here; we will do later if requested. So use 'noind' option.
    try
        violinplot(Y, 'noind', 'facecolor', mycolor, 'edgecolor', mycolor.*.75, 'mc', mycolor.*.5, 'x', xvals, 'medc', [], varargin{:});
        legend off
    catch
        warning('Default violinplot cannot be plotted.')
    end

end

if doxlim
    set(gca,'XLim', [min(xvals)-.5 max(xvals) + .5]);
end

set(gca, 'XTick', xvals, 'XTickLabel', xvals)
xlabel('Condition'), ylabel('Outcome value')
title(mytitle, 'FontSize', 24);

if ~isempty(names)
    if verLessThan('matlab','8.4')
        set(gca, 'XTickLabel', names);
    else % not sure if this works with 8.4, but doesn't work with 8.3. Update conditional if needed
        set(gca, 'XTickLabel', names, 'XTickLabelRotation', 45);
    end
end

if doind

    if iscell(mycolor)
        mycolor = cat(1, mycolor{:});
    elseif size(mycolor, 1) < ny
        mycolor = repmat(mycolor, ny, 1);
    end

    % ----------------------------------------------------
    % > Plot individual points and lines
    % ----------------------------------------------------

    % these are cells - get violin points in the same way that violinplot does
    xvalues = get_violin_points(xvals, y);

    for i = 1:ny
        xvalues{i} = naninsert(all_nan_index{i}, xvalues{i});
    end

    % ----------------------------------------------------
    % Plot parallel coordinate lines across columns, if requested

    if dolines

        % this works because points are matched, equal numbers in each column...
        xvalues_for_lines = cat(2, xvalues{:});

        for j = 1:size(dat, 1) % j is row (observation)

            handles.parallel_line_han{j} = plot(xvalues_for_lines(j, :), dat(j, :), 'k', 'LineWidth', .5, 'Color', [.7 .7 .7]);

        end
    end

    % ----------------------------------------------------
    % > Sort individual scores by covariate if entered, for point plot
    % ----------------------------------------------------
    % sort individual scores by covariate, if covs
    if ~isempty(covs) && wh_reg

        [sortedcov,indx] = sort(covs(:,wh_reg));
        dat = dat(indx,:);
        sortedw = myweights(indx,:);

        x = (0:(nn-1)) ./ (nn+5) - (.5 - 1/(nn-5));
    else
        x = zeros(1,nn) - .1;
        sortedw = myweights;
    end

    % NOTE: X IS NOT USED ANYMORE, LEGACY FOR FUTURE CLEAN-UP...
    x = [x; x];

    hold on

    % ----------------------------------------------------
    % Plot individual points, if requested

    [handles.point_han1, handles.text_han, handles.point_han] = deal(cell(size(dat, 1), ny));

    for i = 1:ny % i is column

        % marker
        mym = 'o';

        % get color
        % ----------------------------------------------------
        if iscell(mycolor) && ~ischar(mycolor)
            mycolcolor = mycolor{i} ./ 2;

        elseif iscell(mycolor) && ischar(mycolor{1})
            mycolcolor = [.2 .2 .2];

        else
            mycolcolor = mycolor(i, :) ./ 2;
        end

        % if matrix
        if ismatrix(mycolcolor) && size(mycolcolor, 1) == ny
            mycolcolor = mycolcolor(i, :);
        end

        % Plot points
        % ----------------------------------------------------

        % 11/2023, tor: fast method, avoiding slow plot of different weights on
        % points if all weights are the same
        if ~donumber && all(stats.wts == 1)
            % simple column
            myc = mycolcolor;
            wh = ~isnan(dat(:, i));


            myhan = plot(xvalues{i}(wh), dat(wh, i), mym, 'MarkerSize', mymarkersize, 'Color', mycolcolor ./2 , 'LineWidth', 1, 'MarkerFaceColor', myc);
            handles.point_han{1, i} = myhan;

        else
            % point-by-point
            for j = 1:size(dat, 1) % j is row (observation)
                % Plot points for each row, for this outcome (column)
                % If first outcome column, plot lines if requested

                % Color by weight
                % Uses original color if weights are 1, shades towards white for weights < 1

                mywt = sortedw(j,i);
                myc = mywt * mycolcolor + (1-mywt) * ([1 1 1] - mycolcolor);  % (myc .* sortedw(j,i));

                pointsexist = ~(any(isnan(x(:, j))) || isnan(dat(j, i)));

                if donumber && pointsexist
                    % Text labels with numbers

                    handles.point_han1{j, i} = plot(xvalues{i}(j), dat(j, i), 'w.'); % to set axis scale appropriately

                    handles.text_han{j, i} = text(xvalues{i}(j), dat(j, i), num2str(j), 'Color', [0 0 0], 'FontSize', 14);

                elseif pointsexist

                    if verLessThan('matlab','8.4')
                        %handles.point_han{j, i} = scatter(xvalues{i}(j), dat(j, i), mymarkersize, mycolcolor ./2 , mym, 'LineWidth', 1, 'MarkerFaceColor', myc);
                        handles.point_han{j, i} = plot(xvalues{i}(j), dat(j, i), mym, 'MarkerSize', mymarkersize, 'Color', mycolcolor ./2 , 'LineWidth', 1, 'MarkerFaceColor', myc);

                    else % doesn't work with 8.3. Not sure about 8.4. Update conditional if needed

                        % Note: tor: changed to plot, scatter is very slow. Eliminate alpha.
                        handles.point_han{j, i} = plot(xvalues{i}(j), dat(j, i), mymarkersize, mycolcolor ./2 , mym, 'LineWidth', 1, 'MarkerFaceColor', myc);

                        if myalpha > .9
                            %handles.point_han{j, i} = plot(xvalues{i}(j), dat(j, i), mymarkersize, mycolcolor ./2 , mym, 'LineWidth', 1, 'MarkerFaceColor', myc);
                            handles.point_han{j, i} = plot(xvalues{i}(j), dat(j, i), mym, 'MarkerSize', mymarkersize, 'Color', mycolcolor ./2 , 'LineWidth', 1, 'MarkerFaceColor', myc);
                        else
                            %                         This can be very slow with scatter.m sometimes...
                            handles.point_han{j, i} = scatter(xvalues{i}(j), dat(j, i), mymarkersize, mycolcolor ./2 , mym, 'LineWidth', 1, 'MarkerFaceColor', myc,'MarkerFaceAlpha',myalpha,'MarkerEdgeAlpha',myalpha);
                        end

                    end


                end % if points exist

            end % j data points

        end % simple column vs point-by-point

        if plotout
            % Plot outliers

            z = (dat(:,i) - mean(dat(:,i))) ./ std(dat(:,i));
            wh = find(abs(z) >= 1.96);

            if ~isempty(wh)
                handles.outlier_han{i} = plot(x(:,wh) + i,[dat(wh,i) dat(wh,i)]', 'ro', 'MarkerSize', mymarkersize .* 2, 'LineWidth', 2);
            end
        end

    end

end % plot individuals


% Add stars
% ----------------------------------------------------
if dostars

    handles.star_handles = star_plot(P, xvals, mymarkersize);

end



% --- Add pairwise comparison brackets and stars ---
if dopairwise
    % Get current Y limits
    yl = ylim();
    
    % Use actual plot data (mean + error bars) for upper reference
    y_max_plot = max(Mean_Value(:) + Std_Error(:), [], 'omitnan');
    
    % Start brackets a bit above plot
    y_start = y_max_plot + 0.05 * range(yl);
    y_spacing = 0.05 * range(yl);
    
    bracket_level = 1;
    top_star_y = y_start;  % track highest Y value used
    
    for i = 1:size(combs, 1)
        if pairwise_p(i) < 0.05
            col1 = combs(i, 1);
            col2 = combs(i, 2);
            x1 = xvals(col1);
            x2 = xvals(col2);
            x_middle = mean([x1, x2]);
    
            % Compute Y coordinate for this bracket
            y_bracket = y_start + (bracket_level - 1) * y_spacing;
            y_star = y_bracket + y_spacing + 0.01 * range(yl);
            top_star_y = max(top_star_y, y_star);
    
            % Draw bracket
            bracket_height = y_spacing * 0.6;  % shorter verticals
            top_y = y_bracket + bracket_height;
            
            % Draw bracket manually: left vertical, top horizontal, right vertical
            plot([x1 x1], [y_bracket, top_y], 'k', 'LineWidth', 1.5);  % left vertical
            plot([x1 x2], [top_y, top_y],   'k', 'LineWidth', 1.5);    % top line
            plot([x2 x2], [top_y, y_bracket], 'k', 'LineWidth', 1.5);  % right vertical
    
            % Significance stars
            if pairwise_p(i) < 0.001
                star_text = '***';
            elseif pairwise_p(i) < 0.01
                star_text = '**';
            else
                star_text = '*';
            end
    
            text(x_middle, y_star, star_text, ...
                 'HorizontalAlignment', 'center', ...
                 'FontSize', 16, 'FontWeight', 'bold');
    
            bracket_level = bracket_level + 1;
        end
    end
    
    % Now safely extend ylim to include topmost star
    ylim([yl(1), top_star_y + 0.05 * range(yl)]);
end


end % main function



function linehandles = plot_violin_points(x, Y, lc, fc)
% x = vector of x positions for each "column"
% Y = cell array of input data, one cell per "column"
% U, F = outputs from ksdensity, normalized, or [] to recalculate
% lc = len(Y) x 3 vector of point fill colors
% fc = len(Y) x 3 vector of point line colors
%
% added by Tor Wager, Sept 2015
% lc is line color, will be used as fill
% fc is fill color, will be used for lines, in a strange twist of fate
% designed to increase contrast

% Enforce cell, no NaNs
% ------------------------------------------------
Y = enforce_cell_array(Y);

xvalues = get_violin_points(x, Y);

linehandles = [];

for i = 1:size(Y, 2)

    myfillcolor = lc(i, :); % line color for this plot
    mylinecolor = fc(i, :); % line color for this plot

    myY = Y{i};     % data points

    % set point size
    pointsize = 1000 ./ length(myY);
    pointsize(pointsize < 1) = 1;
    pointsize(pointsize > 12) = 12;

    linehandles{i} =  plot(xvalues{i}, myY, 'o', 'Color', mylinecolor, 'MarkerSize', pointsize, 'MarkerFaceColor', myfillcolor);

end % column

end % function




function Y = enforce_cell_array(Y)

k = size(Y, 2);

if ~iscell(Y)

    for i = 1:k

        Ytmp{i} = Y(:, i);

        Ytmp{i}(isnan(Ytmp{i})) = [];

    end

    Y = Ytmp;

end
end % function




function xvalues = get_violin_points(x, Y)
% x = vector of x positions for each "column"
% Y = cell array of input data, one cell per "column"
%
% added by Tor Wager, Sept 2015

nbins = 10;

k = size(Y, 2);

xvalues = cell(1, k);

% Enforce cell, no NaNs
% ------------------------------------------------
Y = enforce_cell_array(Y);


% calculate  density values
% ------------------------------------------------
for i = 1:k

    [f, u, bb]=ksdensity(Y{i});

    f=f/max(f)*0.3; %normalize
    F(:,i) = f;
    U(:,i) = u;

end

% get x positions
% ------------------------------------------------

for i = 1:k

    myx = x(i);     % x-value for this bar in plot

    myU = U(:, i);  % x-values of ksdensity output
    myF = F(:, i);  % y-values (density) of ksdensity output

    myY = Y{i};     % data points
    mybins = linspace(min(myY), max(myY), nbins);

    % starting and ending values
    st = [-Inf mybins(1:end-1)];
    en = [mybins];

    for j = 1:nbins
        % define points within a bin or 'slab'

        whpoints = myY > st(j) & myY <= en(j);

        if sum(whpoints) == 0, continue, end

        whu = myU > st(j) & myU <= en(j);

        mylimit(j) = nanmean(myF(whu));  % average density for this 'slab' of points
        %if isnan(mylimit(j)), mylimit(j) = 0; end % can happen...

        % this will be the limit on x-values

        % re-use nbins here for convenience - now it means bins ACROSS the
        % graph within the violin shape, though
        % interleave points on either side of the midline
        % make the first point the actual midline value

        my_xvals = linspace(myx - mylimit(j), myx + mylimit(j), sum(whpoints))';

        if length(my_xvals) == 1, my_xvals = myx;  end

        %         if mod(length(my_xvals), 2)
        %             % odd number, use the midline point
        %             my_xvals = [myx; my_xvals];  % (1:end-1)];
        %         end

        % build coordinates:
        % xlocs is left to right, ylocs is y-values to plot

        ylocs = myY(whpoints);
        %xlocs = repmat(my_xvals, ceil(length(ylocs) ./ length(my_xvals)), 1);

        xlocs = my_xvals(1:length(ylocs));

        % save in original point list
        xvalues{i}(whpoints, 1) = xlocs;

        %linehandles{i, j} =  plot(xlocs, ylocs, 'o', 'Color', mylinecolor, 'MarkerSize', pointsize, 'MarkerFaceColor', myfillcolor);

    end % slab

end % column

end % function


function dat = enforce_padded_matrix_form(dat)

%handle input of different lengths -- passed in as cell array
if iscell(dat)
    maxlen = 0;
    for i = 1:length(dat)
        if length(dat{i}) > maxlen
            maxlen = length(dat{i});

        end
    end

    dat2 = repmat(NaN, maxlen, length(dat));
    for i=1:length(dat)
        mylen = length(dat{i});

        dat2(1:mylen,i) = dat{i};

        if mylen < maxlen, nanwarningflag = 0; end % turn off warning in this case
    end

    dat = dat2;
end

end


function star_han = star_plot(P, xvals, mymarkersize)

star_han = [];
%starwid = .06;
starsize = round(mymarkersize);

% Get y position for all stars
my_ylim = get(gca, 'YLim');
yval = my_ylim(2) - .05 * range(my_ylim);

for i = 1:length(P)

    if P(i) < FDR(P,.05), mystr = '***';
    elseif P(i) < .01, mystr = '**';
    elseif P(i) < .05, mystr = '*';
    else mystr = ''; xadj = 0;
    end

    star_han(i) = text(xvals(i), yval, mystr, 'FontSize', starsize, 'FontWeight', 'b');

    % Center text - Adjust for width of text string
    star_han(i) = center_text(star_han(i));


end % loop through columns of plot

end % star_plot function


function text_handle = center_text(text_handle)
% Center text - Adjust for width of text string
% Given text handle, subtract 1/2 the width of text string from x and y positions
% Borrowed/adapted from John Barber, http://www.mathworks.com/matlabcentral/fileexchange/30671-calcticks

textExt = get(text_handle, 'Extent');
textHeight = textExt(4);
textWidth = textExt(3);

% If using a proportional font, shrink text width by a fudge factor to
% account for kerning.
% ax = gca; -kragel changed to be compatible with older versions of matlab

if ~strcmpi(get(gca,'FontName'),'FixedWidth')
    textWidth = textWidth*0.8;
end

mypos = get(text_handle, 'Position');
mypos(1) = mypos(1) - textWidth ./ 2;
%mypos(2) = mypos(2) + textHeight ./ 2; % from bottom - adjust up
set(text_handle, 'Position', mypos);

end

