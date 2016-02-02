function [hout,dat,xdat, h] = barplot_columns(dat,varargin)
% :Usage:
% ::
%
%    [axishandle,adjusted data,x-data, barhandle] = barplot_columns(dat, [other optional arguments])
%
% This function makes a barplot of columns of data, with standard error
% bars.  Optional arguments include removing continuous covariates before plotting,
% robust (IRLS) estimation of means and correlations with covariates, and
% within-subject error bars based on the subject x condition interaction
% (overall), which is not quite the standard error contrasts of interest,
% but is the standard error for a 1-way repeated measures ANOVA.
%
% plots circles around points at z >= 1.96
%
% plots individual points, unless you enter 4th argument
%
% if dat is a cell array, each entry becomes one "bar".  Useful if n
% observations is different for each column.
%
% :Examples: Just plot means and SE
% ::
%
%    h = barplot_columns(tmp,'Cluster 1',[],1);
%
% :Optional arguments:
%   1. Title for figure
%   2. covariates
%   3. String Arguments
%        - 'nofig' : do not make figure
%        - 'noind' : do not plot individual scores
%        - 'plotout': circle potential outliers at z>1.96 in red
%        - 'dorob' : do robust IRLS means and correlations
%        - 'dolines' : plot lines showing individual effects
%        - 'within' : within-subjects standard errors, followed by contrast
%                   matrix
%        - '95CI'   : error bars are 95% CI instead of SE
%        - 'line' : Make line plot instead of bar plot
%        - 'number' : plot case numbers instead of points
%        - 'x' : followed by x-axis values for bars
%        - 'color' : followed by color for bars (text: 'r' or [r g b]) OR
%               cell array with names of colors cell for each line/bar
%        - 'violin': add violin plot to each bar, with data points
%
%
% :Examples:
% ::
%
%    barplot_columns(ctmp,'RT effects by Switch Type',overall_sw,'nofig','dorob')
%
% Standard Errors ARE NOT Adjusted for covariate, right now.
%
% Example: within-subjects std. errors
% ::
%
%    barplot_columns(dat, 'Means', [], 'nofig', 'within', c);
%
% The example below uses color, width, and xposition arguments to make a grouped
% ::
%
%    barplot showing effects for two groups:
%    exp_dat = EXPT.error_rates(EXPT.group==1,:);
%    control_dat = EXPT.error_rates(EXPT.group==-1,:);
%    barplot_columns(exp_dat, 'Error rates', [], 'nofig', 'noind', 'color', 'r','width', .4);
%    barplot_columns(control_dat, 'Error rates', [], 'nofig', 'noind', 'color', 'b','width', .4, 'x', (1:9)+.5);
%    set(gca, 'XLim', [0 10], 'XTick', 1:9)
%
%    barplot_columns(nps_by_study, 'NPS by study', [], 'doind', 'colors', mycolors, 'nofig');
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
% See also: lineplot_columns, barplot_colored, line_plot_multisubject, violinplot

% ..
%    Defaults
% ..

dofig = 1;
doind = 1;
plotout = 0;
dorob = 0;
xdat = [];
dolines = 0;
dowithin = 0;
donumber = 0;
dojitter = 1; % jitter is for numbers only
mycolor = [.8 .8 .8];
barwidth = .8;
dolineplot = 0;
do95CI = 0;
nanwarningflag = 1;
doviolin = 1;
mytitle = [];
covs = [];

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
        if strcmp(varargin{i},'nofig'), dofig = 0;  end
        if strcmp(varargin{i},'noind'), doind = 0;  end
        if strcmp(varargin{i},'plotout'), plotout = 1;  end
        if strcmp(varargin{i},'dorob'), dorob = 1;  end
        if strcmp(varargin{i},'dolines'), dolines = 1;  end
        if strcmp(varargin{i},'number'), donumber = 1;  end
        if strcmp(varargin{i}, 'within'), dowithin = 1; end %cons = varargin{i + 1}; end
        if strcmp(varargin{i}, '95CI'), do95CI = 1; end
        if strcmp(varargin{i},'line'), dolineplot = 1;  end
        if strcmp(varargin{i}, 'x'), xvals = varargin{i + 1}; end
        if strcmp(varargin{i}, 'color') || strcmp(varargin{i}, 'colors')
            mycolor = varargin{i + 1};
            varargin{i + 1} = [];
        end
        if strcmp(varargin{i}, 'width'), barwidth = varargin{i + 1}; end
        if strcmp(varargin{i}, 'violin')
            doviolin = 1; 
            doind = 0;  % ind is already done in violin
        end
        
        if strcmp(varargin{i}, 'noviolin')
            doviolin = 0;
        end
        
        if strcmp(varargin{i}, 'covs')
            covs = varargin{i + 1};
            if ~isempty(covs), covs = scale(covs,1); end
        end
        
        if strcmp(varargin{i}, 'title'), mytitle = varargin{i + 1}; end
        
    end % for
end % varargin

% ----------------------------------------------------
% > Build design matrix X for controlling for covariates
% ----------------------------------------------------

dat = double(dat);

% delete nans casewise
%wh = find(any(isnan(dat),2));
%dat(wh,:) = [];
%if ~isempty(covs), covs(wh,:) = []; end

% replace nans with mean
for i = 1:size(dat, 2)
    if nanwarningflag && any(isnan(dat(:, i)))
        warning(sprintf('Some NaNs in Column %3.0f!', i))
    end
end

% find NaN columns
wh = find(all(isnan(dat),1));
dat(:, wh) = 0;

% get final design matrix, intercept is last column
[nn,ny] = size(dat);
k = size(covs,2);

% add intercept if not already in model
if ~isempty(covs)
    wh_oldintercept = find(all(diff(covs) < eps));
    covs(:, wh_oldintercept) = [];
end

X = [covs ones(nn,1)];
wh_intercept = k+1;

% ----------------------------------------------------
% > Get means and standard error of means
% With robust option, if specified, and removing
% covariates, if there are any.
% ----------------------------------------------------

stderr = [];

% key vars are :
% mymeans, stderr, mycor

wh_reg = 1; % regressor of interest

for i = 1:ny
    
    fprintf(1,'\nColumn %3.0f:\t',i);
    
    % ----------------------------------------------------
    % > Get [robust or non-robust] mean and standard error
    %   Return y, data from column, with nans in
    % ----------------------------------------------------
    
    % remove nans from this column
    tmpy = dat(:,i);
    tmpx = X;
    
    [wasnan, tmpx, tmpy] = nanremove(tmpx, tmpy);
    all_nan_index{i} = wasnan;
    
    n(i) = size(tmpy,1);
    
    % get mean and standard error of intercept (robust or OLS)
    % y is adjusted for all non-intercept covs
    % stats has weights, stats.w, which are all 1 for OLS
    [x, newy, r, p, stderr(i), mymeans(i), stats] = partialcor(tmpx, tmpy, wh_intercept, 1, dorob);
    
    %95% CI?
    if do95CI, stderr(i) = stderr(i) * 1.96; end
    
    y(:,i) = naninsert(wasnan, newy);
    
    %%%not needed y(:,i) = y(:,i) + mymeans(i);   % add mean
    myweights(:,i) = naninsert(wasnan, stats.w);
    
    % ----------------------------------------------------
    % > Use partialcor to remove covariates if requested
    %   Return y, adjusted y-values
    % ----------------------------------------------------
    
    if ~isempty(covs)
        % cov of interest here is fixed at 1 (see above)
        
        % if we have covs, leave in cov. of interest (cov1)
        % y is adjusted for all non-intercept covs
        [x,y(:,i),mycor(i),mycorrp(i)] = partialcor(tmpx,tmpy,wh_reg,1,dorob);
        
        %not needed %%% y(:,i) = y(:,i) + mymeans(i);   % add mean
        
    end
    
    %fprintf(1,'\n');
end

dat = y;  % adjusted data, for plot

if dowithin
    
    within_ste = barplot_get_within_ste(dat);
    
    stderr = repmat(within_ste, 1, size(dat, 2));
    
end

% ----------------------------------------------------
% > Make figure
% ----------------------------------------------------

if dofig
    f = figure('Color','w'); hout = gca; set(gca,'FontSize',18); %hold on; grid on;
else
    f = get(gcf); hout = gca; set(gca,'FontSize',18); hold on;
end


% ----------------------------------------------------
% > BARPLOT (or line plot)
% ----------------------------------------------------

if dolineplot
    h = plot(xvals, mymeans, 'o-', 'Color', mycolor, 'MarkerFaceColor', mycolor, 'MarkerSize', 8);
    h2 = errorbar(xvals, mymeans, stderr, stderr);
    set(h2, 'LineWidth', 2, 'Color', mycolor);
    
else
    
    h = bar(xvals, mymeans, barwidth);
    if iscell(mycolor)
        % each bar a different color
        for i = 1:length(xvals)
            bar(xvals(i), mymeans(i), 'FaceColor', mycolor{i});
            
            h2 = errorbar(xvals(i), mymeans(i), stderr(i), 'Color', mycolor{i} ./ 2, 'LineWidth', 3);
            
        end
        
    else %all bars the same color
        set(h,'FaceColor', mycolor); %,'LineWidth',2)
        
        for i = 1:length(xvals)
        h2 = errorbar(xvals(i), mymeans(i), stderr(i), 'Color', mycolor ./ 2, 'LineWidth', 3);
        end
        
    end
    
    %tor_bar_steplot(mymeans,stderr, {'k'}, xvals);

end

% add violin if entered
if doviolin
    
    if iscell(mycolor)
    mycolor = cat(1, mycolor{:});
    elseif size(mycolor, 1) < ny
        mycolor = repmat(mycolor, ny, 1);
    end
    
    Y = enforce_cell_array(y);
violinplot(Y, 'facecolor', mycolor, 'edgecolor', mycolor.*.75, 'mc', mycolor.*.5, 'x', xvals, 'medc', []);
legend off
    
    
end

set(gca,'XLim',[min(xvals)-.5 max(xvals) + .5], 'XTick', xvals, 'XTickLabel', xvals)
xlabel('Condition'), ylabel('Outcome value')
title(mytitle, 'FontSize', 24);


if doind
    
    % ----------------------------------------------------
    % > Plot individuals
    % ----------------------------------------------------
    
    % these are cells - get violin points
    xvalues = get_violin_points(xvals, y);
    
    for i = 1:ny
        xvalues{i} = naninsert(all_nan_index{i}, xvalues{i});
    end
    
    if dolines
        % this works because points are matched, equal numbers in each column...
        xvalues_for_lines = cat(2, xvalues{:});
    end
    
    % ----------------------------------------------------
    % > Sort individual scores by covariate if entered, for point plot
    % ----------------------------------------------------
    % sort individual scores by covariate, if covs
    if ~isempty(covs)
        
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
    for i = 1:ny % i is column
        
        % marker
        if mod(i,2)==0, mym='^'; myc=[.2 .2 .2]; else mym='o'; myc=[0 0 0]; end
        
        for j = 1:size(dat, 1) % j is data point
            
            % color by weight
            myc = [1 1 1] - ([1 1 1] .* sortedw(j,i));
            
            if i == 1
                %xdat(j,:) = x(1,j) + (1:size(dat,2)); % save x values for output (for line plotting)
                
                % plot parallel lines (if requested)
                if dolines
                    plot(xvalues_for_lines(j, :), dat(j, :),'k', 'LineWidth', .5, 'Color', [.7 .7 .7]);
                end
            end
            
            % plot marker
            %plot(x(:,j) + i,[dat(j,i) dat(j,i)]', 'w.'); % to set axis scale appropriately
            plot(xvalues{i}(j), dat(j, i), 'w.'); % to set axis scale appropriately
            
            if donumber && ~(any(isnan(x(:, j))) || isnan(dat(j, i)))
                %plot(x(:,j) + i,[dat(j,i) dat(j,i)]', 'w.'); % to set axis scale appropriately
                %text(x(1,j) + i, dat(j,i), num2str(j), 'Color', [0 0 0], 'FontSize', 14)
                text(xvalues{i}(j), dat(j,i), num2str(j), 'Color', [0 0 0], 'FontSize', 14)
                
            elseif ~(any(isnan(x(:, j))) || isnan(dat(j, i)))
                %plot(x(:,j) + i,[dat(j,i) dat(j,i)]',mym,'Color',[0 0 0],'LineWidth',1,'MarkerFaceColor',myc)
                plot(xvalues{i}(j), dat(j,i), mym, 'Color', [0 0 0], 'LineWidth', 1, 'MarkerFaceColor', myc);
                
            end
            
        end % j data points
        
        if plotout
            z = (dat(:,i) - mean(dat(:,i))) ./ std(dat(:,i));
            wh = find(abs(z) >= 1.96);
            
            if ~isempty(wh),plot(x(:,wh) + i,[dat(wh,i) dat(wh,i)]','ro','MarkerSize',14,'LineWidth',2),end
        end
        
    end
    
    %     if dolines
    %         for i = 1:size(dat,1)
    %             plot(xdat(i,:),dat(i,:),'k','LineWidth',.5,'Color',[.7 .7 .7]);
    %         end
    %     end
    
    % if donumber && dojitter
    %
    %     % Jitter positions
    %     han = findobj(gca, 'Type', 'text');
    %
    %     x = get(han, 'Position'); x = cat(1, x{:});
    %     n = size(x, 1);
    %     x(:, 1) = x(:, 1) + .2 * (rand(n, 1) - .5);
    %     for i = 1:size(x, 1)
    %         set(han(i), 'Position', x(i, :));
    %     end
    %
    %     % Bold
    %     set(han, 'FontWeight', 'b', 'Color', 'b');
    %
    % end % jitter
    
    
end % plot individuals

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
        
        if mod(length(my_xvals), 2)
            % odd number, use the midline point
            my_xvals = [myx; my_xvals(1:end-1)];
        end
        
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
