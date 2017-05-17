function [hout,dat,xdat, h] = barplot_columns3(dat,plottitle,varargin)
% :Usage:
% ::
%
%    [axishandle,adjusted data,x-data,barhandle] = barplot_columns(dat,title,[options])
%
% Makes a barplot of columns of data, with standard error bars.
% Optional arguments include removing continuous covariates before plotting,
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
% :Options:
%
%   **'cov':**
%        followed by matrix of covariates
%
%   **'labels':**
%        followed by cellstring of bar labels
%
%   **'nofig':**
%        do not make figure
%
%   **'ind':**
%        plot individual scores on top of bars
%
%   **'plotout':**
%        circle potential outliers at z>1.96 in red
%
%   **'robust':**
%        do robust IRLS means and correlations
%
%   **'indlines':**
%        plot lines showing individual effects
%
%   **'within':**
%        within-subjects standard errors, followed by contrast
%        matrix
%
%   **'line':**
%        Make line plot instead of bar plot
%
%   **'number':**
%        plot case numbers instead of points
%
%   **'x':**
%        followed by x-axis values for bars
%
%   **'color':**
%        followed by color for bars (text: 'r' or [r g b]) OR
%        cell array with names of colors cell for each line/bar
%
%   **'cons':**
%        followed by contrastXweights matrix
%
% To convert from long form to wide form
%
%   **'subcol':**
%        column numbers with subject numbers
%
%   **'propcols':**
%        vector of column numbers with properties
%
%   **'propnames':**
%        cell array of names of properties
%
% :Examples:
% ::
%
%    barplot_columns(ctmp,'RT effects by Switch Type',overall_sw,'nofig','robust')
%
% Standard Errors ARE NOT Adjusted for covariate, right now.
%
% Example: within-subjects std. errors
% ::
%
%    barplot_columns(dat, 'Means', 'nofig', 'within', 'ind');
%
% The example below uses color, width, and xposition arguments to make a grouped
% ::
%
%    barplot showing effects for two groups:
%    exp_dat = EXPT.error_rates(EXPT.group==1,:);
%    control_dat = EXPT.error_rates(EXPT.group==-1,:);
%    barplot_columns(exp_dat, 'Error rates', 'nofig', 'color', 'r', 'width', .4);
%    barplot_columns(control_dat, 'Error rates', 'nofig', 'color', 'b', 'width', .4, 'x', (1:9)+.5);
%    set(gca, 'XLim', [0 10], 'XTick', 1:9)
%
% ..
%    PROGRAMMERS' NOTES
%    4/2013 (Luka): changed input format (title and cov now varargin options, not sequential arguments)
%    4/2013 (Luka): changed input parsing
%    4/2013 (Luka): added labels option (adds diagonal labels)
%    4/2013 (Luka): moved labeling to after doind (axis is stable then)
%    4/2013 (Luka): added plabels and toplabels options
%    5/2013 (Luka): made no individuals default
%    need to sort out plabel option (including p values from glms in plot)
%    right now does not account for inclusion of cov
%    figure out how to add stars, legend?
% ..

% ..
%    Set up input arguments
% ..

disp('WARNING: DEPRECATED -- BARPLOT_COLUMNS.M IS RECOMMENDED.')

DO_fig = 1;
DO_ind = 0;
DO_plotout = 0;
DO_rob = 0;
DO_indlines = 0;
DO_within = 0;
DO_number = 0;
DO_lineplot = 0;
DO_plabels = 0;
DO_denan = 0;

mycolor = [.8 .8 .8];
barwidth = .8;
xdat = [];
covs = [];

%handle input of different lens -- passed in as cell array
if iscell(dat)
    maxlen=0;
    for i=1:length(dat)
        if length(dat{i}) > maxlen, maxlen = length(dat{i}); end
    end
    
    dat2 = repmat(NaN, maxlen, length(dat));
    for i=1:length(dat)
        dat2(1:length(dat{i}),i) = dat{i};
    end
    dat = dat2;
end

i=1;
while i <= numel(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            case 'subcol'
                subcol = varargin{i+1};
                i=i+1;
            case 'propcols'
                propcols = varargin{i+1};
                i=i+1;
            case 'propnames'
                propnames = varargin{i+1};
                i=i+1;
            case 'cons'
                cons = varargin{i+1};
                i=i+1;
            case 'connames'
                connames = varargin{i+1};
                i=i+1;
            case 'cov'                
                covs = scale(varargin{i+1},1);
                i=i+1;
            case 'labels'                
                xlabels = varargin{i+1};
                i=i+1;
            case 'ylabel'
                yl = varargin{i+1};
                i=i+1;
            case 'ylim',        ylimits = varargin{i+1}; i=i+1;
            case 'nofig',       DO_fig = 0;
            case 'ind',         DO_ind = 1;
            case 'plotout',     DO_plotout = 1;
            case 'robust',      DO_rob = 1;
            case 'denan',       DO_denan = 1;
            case 'indlines',    DO_indlines = 1;
            case 'number',      DO_number = 1;
            case 'within',      DO_within = 1;
            case 'line',        DO_lineplot = 1;
            case 'plabels',     DO_plabels = 1;
            case 'toplabels'
                toplabels = varargin{i+1};
                i=i+1;
            case 'x'
                xvals = varargin{i+1};                
                i=i+1;
            case 'color'
                mycolor = varargin{i+1};
                i=i+1;
            case 'width'                
                barwidth = varargin{i+1};
                i=i+1;
            otherwise
                error(['UNRECOGNIZED OPTION: ' varargin{i}])
        end
    else
        disp(varargin{i})
        error('Above option UNRECOGNIZED')
    end
    i=i+1;
end   

if exist('xlabels','var') && numel(xlabels) ~= size(dat,2)
    error('Number of labels given (%d) does not match the number of bars (%d)\n',numel(xlabels),size(dat,2))
end
if ~exist('plottitle','var')
    error('Second argument (title) is mandatory.')
end



if exist('subcol','var') && exist('propcols','var')    
    datacols = setdiff(1:size(dat,2),union(subcol,propcols));
    [propcombos ignore prowcode] = unique(dat(:,propcols),'rows');            
    [subjects ignore srowcode]= unique(dat(:,subcol));
    dat2 = [];
    for s = 1:numel(subjects)
        newrow = [];
        for c = 1:size(propcombos,1)
            newrow = [newrow nanmean(dat(srowcode==s & prowcode==c,datacols))];
        end
        dat2(s,:) = newrow;
    end
    
    if exist('xlabels','var')
        if numel(xlabels) ~= size(propcombos,1)
            error('number of labels must match number of property combos')
        end
    else
        for c = 1:size(propcombos,1)
            xlabels{c} = '';
            for p = 1:size(propcombos,2)
                xlabels{c} = [xlabels{c} '_' propnames{p} '(' num2str(propcombos(c,p)) ')'];
            end
        end
        xlabels = regexprep(xlabels,'^_','');
    end
end
clear s c p i j
dat=dat2;

if exist('cons','var')
    dat2 = [];
    xlabels2 = {};
    for i=1:size(cons,1)
        dat2(:,i) = sum(dat .* repmat(cons(i,:),size(dat,1),1),2);
        if exist('connames','var')
            xlabels2{i} = connames{i};
        else
            xlabels2{i} = sprintf('con%04d',num2str());
        end
    end
    
    dat = dat2;
    xlabels = xlabels2;
end

if ~exist('xvals','var'), xvals = 1:size(dat, 2); end

if DO_denan
    % delete nans casewise    
    wh = find(any(isnan(dat),2));
    if ~isempty(wh)
        fprintf('WARNING: The following rows are being excluded because they include NaNs:\n')
        fprintf('%d\n',wh)
        dat(wh,:) = [];
        if ~isempty(covs), covs(wh,:) = []; end
    end
end

% replace nans with mean
for i = 1:size(dat,2),
    if any(isnan(dat(:,i)))
        warning('Some NaNs!')
        %dat(find(isnan(dat(:,i))),i) = nanmean(dat(:,i));
    end
end

% find NaN columns
wh = find(all(isnan(dat),1));
dat(:,wh) = 0;

% get final design matrix, intercept is last column
[nn,ny] = size(dat);
k = size(covs,2);

% add intercept if not already in model
if ~isempty(covs)
    wh_oldintercept = find(all(diff(covs) < eps));
    covs(:,wh_oldintercept) = [];
end

X = [covs ones(nn,1)];
wh_intercept = k+1;

% ----------------------------------------------------
% > Get means and standard error of means
% With robust option, if specified, and removing
% covariates, if there are any.
% ----------------------------------------------------


% if DO_rob
%     % ROBUST FIT
stderr = [];

%[b,t,p,sig,f,fp,fsig,stat] = robust_reg_matrix(X,dat,1);

%for i = 1:k
%    [x,y,r,p,mycor(1,i),prob,serob] = partialcor(X,y,i)

% key vars are :
% mymeans, stderr, mycor

wh_reg = 1; % regressor of interest

for i = 1:ny
    fprintf(1,'\nColumn %3.0f:\n',i);
    
    % remove nans from this column
    tmpy = dat(:,i);
    tmpx = X;
    % %     wh = find(any(isnan(X),2) | any(isnan(tmpy),2));
    % %     tmpx(wh,:) = []; tmpy(wh) = [];
    [wasnan, tmpx, tmpy] = nanremove(tmpx, tmpy);
    
    n(i) = size(tmpy,1);
    
    % get mean and standard error of intercept (robust or OLS)
    % y is adjusted for all non-intercept covs
    % stats has weights, stats.w, which are all 1 for OLS
    [x,newy,r,p{i},stderr(i),mymeans(i),stats] = partialcor(tmpx,tmpy,wh_intercept,1,DO_rob);
    
    y(:,i) = naninsert(wasnan, newy);
    
    %%%not needed y(:,i) = y(:,i) + mymeans(i);   % add mean
    myweights(:,i) = naninsert(wasnan, stats.w);
    
    if ~isempty(covs)
        % cov of interest here is fixed at 1 (see above)
        
        % if we have covs, leave in cov. of interest (cov1)
        % y is adjusted for all non-intercept covs
        [x,y(:,i),mycor(i),mycorrp(i)] = partialcor(tmpx,tmpy,wh_reg,1,DO_rob);
        
        %not needed %%% y(:,i) = y(:,i) + mymeans(i);   % add mean
    end
    
    fprintf(1,'\n');
end
dat = y;  % adjusted data, for plot

if DO_within
    within_ste = barplot_get_within_ste(dat);    
    stderr = repmat(within_ste, 1, size(dat, 2));
end


% ----------------------------------------------------
% > Make figure
% ----------------------------------------------------

if DO_fig
    f = figure('Color','w'); hout = gca; set(gca,'FontSize',18); %hold on; grid on;
else
    f = get(gcf); hout = gca; set(gca,'FontSize',18); hold on;
end


% ----------------------------------------------------
% > BARPLOT (or line plot)
% ----------------------------------------------------

if DO_lineplot
    h = plot(xvals, mymeans, 'o-', 'Color', mycolor, 'MarkerFaceColor', mycolor, 'MarkerSize', 8);
    h2 = errorbar(xvals, mymeans, stderr, stderr);
    set(h2, 'LineWidth', 2, 'Color', mycolor);
else
    h = bar(xvals, mymeans, barwidth);
    if iscell(mycolor)
        % each bar a different color
        hold on
        for i = 1:length(xvals)      
            bar(xvals(i), mymeans(i), barwidth, 'FaceColor', mycolor{i}, 'LineWidth', 2);
        end
        hold off
    else %all bars the same color
        set(h,'FaceColor', mycolor); %,'LineWidth',2)
    end
    
    tor_bar_steplot(mymeans,stderr,{'k'}, xvals);
end

if exist('ylimits','var')
    ylim(ylimits);
end

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

%s = std(dat1);

if DO_ind
    % ----------------------------------------------------
    % > Plot individuals
    % ----------------------------------------------------

    x = [x; x]; 
    
    hold on
    for i = 1:size(dat,2)
        % marker
        if mod(i,2)==0, mym='^'; myc=[.2 .2 .2]; else mym='o'; myc=[0 0 0]; end

        for j = 1:size(dat,1)
            % color by weight
            myc = [1 1 1] - ([1 1 1] .* sortedw(j,i));
            
            if i == 1
                xdat(j,:) = x(1,j) + (1:size(dat,2)); % save x values for output (for line plotting)

                % plot lines (if requested)
                if DO_indlines
                    plot(xdat(j,:),dat(j,:),'k','LineWidth',.5,'Color',[.7 .7 .7]);
                end
            end
            
            % plot marker
            if DO_number && ~(any(isnan(x(:, j))) | isnan(dat(j, i)))
                plot(x(:,j) + i,[dat(j,i) dat(j,i)]', 'w.'); % to set axis scale appropriately
                text(x(:,j) + i,[dat(j,i) dat(j,i)]',num2str(j),'Color',[0 0 0],'FontSize', 14)
            elseif ~(any(isnan(x(:, j))) | isnan(dat(j, i)))
                plot(x(:,j) + i,[dat(j,i) dat(j,i)]',mym,'Color',[0 0 0],'LineWidth',1,'MarkerFaceColor',myc)
            end            
        end
        
        z = (dat(:,i) - mean(dat(:,i))) ./ std(dat(:,i)); wh = find(abs(z) >= 1.96);

        % print z-scores of potential outliers to text
        %if ~isempty(wh), fprintf(1,'\nCond %3.0f\t',i),for j = 1:length(wh), fprintf(1,'%3.2f\t',z(wh(j))); end,end

        if DO_plotout
            if ~isempty(wh), plot(x(:,wh) + i,[dat(wh,i) dat(wh,i)]','ro','MarkerSize',14,'LineWidth',2), end
        end
    end        
end % plot individuals


% now that axes are stabilized, add diagonal x labels
%set(gca,'XLim',[0 ny+1],'XTick',1:ny,'XTickLabel',1:ny)
if ~exist('xlabels','var')
    set(gca,'XLim',[min(xvals)-.5 max(xvals) + .5],'XTick',xvals,'XTickLabel',xvals)
    xlabel('Task Condition')
else
    set(gca,'XLim',[min(xvals)-.5 max(xvals) + .5])
    xaxis_labeling('Task Condition',xlabels);
end

if exist('yl','var')
    ylabel(yl);
end

title(plottitle,'FontSize',24)
%set(gcf,'Position',[464   283   930   827]), drawnow

if exist('toplabels','var') || (DO_plabels && size(dat,1)>1 && exist('p','var'))
    xd = get(get(h,'Children'),'XData');
    yl = get(get(h,'Parent'),'YLim');
    for i=1:size(xd,2)
        tx = (xd(2,i)+xd(3,i))/2;
        ty = yl(2);   
        if exist('toplabels','var')
            if iscell(toplabels)
                plabel = toplabels{i};
            else
                plabel = sprintf('%5.3f',toplabels(i));
            end
        else            
            plabel = sprintf('%5.3f',p{i});
        end
        text(tx,ty,plabel,'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',12);
    end
end

end


%% slightly edited code from mathworks.com (search "matlab angled tick labels")
% added by Luka 4/2013
function xaxis_labeling(xaxislabel,ticklabels)

% % reduce size of axis to fit labels
% pos = get(gca,'Position');
% set(gca,'Position',[pos(1), .2, pos(3) .65])

% set tick locations
Xt = 1:numel(ticklabels);
Xl = [0 numel(ticklabels)+1];
set(gca,'XTick',Xt,'Xlim',Xl);

ax = axis; % Current axis limits
axis(axis); % Set the axis limit modes to manual
Yl = ax(3:4); % Y-axis limits

% Place the text labels
t = text(Xt,Yl(1)*ones(1,numel(Xt)),ticklabels,'Interpreter','none');
set(t,'HorizontalAlignment','right','VerticalAlignment','top', ...
    'Rotation',45,'FontSize',15);

% Remove the default labels
set(gca,'XTickLabel','');

% get the Extent of each text object
for i = 1:length(ticklabels)
    ext(i,:) = get(t(i),'Extent'); %#ok
    set(t(i),'Units','pixels');
    extp(i,:) = get(t(i),'Extent'); %#ok
    set(t(i),'Units','data');
end

% Determine lower point for alignment of text
LowYPoint = min(ext(:,2));

% reduce size of axis to fit labels
pos = get(gca,'Position');
set(gca,'Position',[pos(1), .2, pos(3) .65])

% Place the axis label
XMidPoint = Xl(1) + abs(diff(Xl))/2;
t2 = text(XMidPoint,LowYPoint,xaxislabel,'VerticalAlignment','top','HorizontalAlignment','center');
set(t2,'FontSize',18);
t2ex = get(t2,'Extent');

% adjust size of window
pos = get(gcf,'Position');
apos = get(gca,'Position');
ss = get(0,'ScreenSize');
newpos4 = pos(4)-t2ex(2)+200;
set(gcf,'Position',[pos(1) ss(4)-newpos4-100 pos(3) newpos4]);
set(gca,'Position',[apos(1) apos(2) apos(3) 1-apos(2)-.1]);

end
