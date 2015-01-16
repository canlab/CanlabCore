function [hout,dat,xdat, h] = barplot_columns(dat,varargin)
% [axishandle,adjusted data,x-data, barhandle] = barplot_columns(dat,[title],[covs],[extra string args])
%
% This function makes a barplot of columns of data, with standard error
% bars.  Optional arguments include removing continuous covariates before plotting,
% robust (IRLS) estimation of means and correlations with covariates, and
% within-subject error bars based on the subject x condition interaction
% (overall), which is not quite the standard error contrasts of interest,
% but is the standard error for a 1-way repeated measures ANOVA.
%
% plots circles around points at z >= 1.96
% plots individual points, unless you enter 4th argument
%
% if dat is a cell array, each entry becomes one "bar".  Useful if n
% observations is different for each column.
%
% Examples:  Just plot means and SE
% h = barplot_columns(tmp,'Cluster 1',[],1);
%
% Optional arguments
% 1 - Title for figure
% 2 - covariates
% 3 - String Arguments
%       'nofig' : do not make figure
%       'noind' : do not plot individual scores
%       'plotout': circle potential outliers at z>1.96 in red
%       'dorob' : do robust IRLS means and correlations
%       'dolines' : plot lines showing individual effects
%       'within' : within-subjects standard errors, followed by contrast
%                   matrix
%       '95CI'   : error bars are 95% CI instead of SE
%       'line' : Make line plot instead of bar plot
%       'number' : plot case numbers instead of points
%       'x' : followed by x-axis values for bars
%       'color' : followed by color for bars (text: 'r' or [r g b]) OR
%               cell array with names of colors cell for each line/bar
%
%
% Examples:
% barplot_columns(ctmp,'RT effects by Switch Type',overall_sw,'nofig','dorob')
%
% Standard Errors ARE NOT Adjusted for covariate, right now.
%
% Example: within-subjects std. errors
% barplot_columns(dat, 'Means', [], 'nofig', 'within', c);
%
% The example below uses color, width, and xposition arguments to make a grouped
% barplot showing effects for two groups:
% exp_dat = EXPT.error_rates(EXPT.group==1,:);
% control_dat = EXPT.error_rates(EXPT.group==-1,:);
% barplot_columns(exp_dat, 'Error rates', [], 'nofig', 'noind', 'color', 'r','width', .4);
% barplot_columns(control_dat, 'Error rates', [], 'nofig', 'noind', 'color', 'b','width', .4, 'x', (1:9)+.5);
% set(gca, 'XLim', [0 10], 'XTick', 1:9)


% ----------------------------------------------------
% > Set up input arguments
% ----------------------------------------------------

dofig = 1; doind = 1; plotout = 0; dorob = 0; xdat = []; dolines = 0; 
dowithin = 0; donumber = 0; dojitter = 1; % jitter is for numbers only
mycolor = [.8 .8 .8];
barwidth = .8;
dolineplot = 0;
do95CI = 0;

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

xvals = 1:size(dat, 2);

if length(varargin) > 2
    for i = 3:length(varargin)
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
        if strcmp(varargin{i}, 'color'), mycolor = varargin{i + 1}; end
        if strcmp(varargin{i}, 'width'), barwidth = varargin{i + 1}; end
        
    end
end

dat = double(dat);

if length(varargin) > 1,
    covs = varargin{2};  if ~isempty(covs), covs = scale(covs,1); end
else covs = [];
end

% delete nans casewise
%wh = find(any(isnan(dat),2));
%dat(wh,:) = [];
%if ~isempty(covs), covs(wh,:) = []; end

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

%dat1 = dat;     % save original data for ind subj plot

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


% if dorob
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
    [x,newy,r,p,stderr(i),mymeans(i),stats] = partialcor(tmpx,tmpy,wh_intercept,1,dorob);
    
    %95% CI?
    if do95CI, stderr(i) = stderr(i) * 1.96; end
    
    y(:,i) = naninsert(wasnan, newy);
    
    %%%not needed y(:,i) = y(:,i) + mymeans(i);   % add mean
    myweights(:,i) = naninsert(wasnan, stats.w);
    
    if ~isempty(covs)
        % cov of interest here is fixed at 1 (see above)
        
        % if we have covs, leave in cov. of interest (cov1)
        % y is adjusted for all non-intercept covs
        [x,y(:,i),mycor(i),mycorrp(i)] = partialcor(tmpx,tmpy,wh_reg,1,dorob);
        
        %not needed %%% y(:,i) = y(:,i) + mymeans(i);   % add mean
        
    end
    
    fprintf(1,'\n');
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
        end
    else %all bars the same color
        set(h,'FaceColor', mycolor); %,'LineWidth',2)
    end
    tor_bar_steplot(mymeans,stderr,{'k'}, xvals);
end

%set(gca,'XLim',[0 ny+1],'XTick',1:ny,'XTickLabel',1:ny)
set(gca,'XLim',[min(xvals)-.5 max(xvals) + .5],'XTick',xvals,'XTickLabel',xvals)
xlabel('Task Condition'), ylabel('BOLD contrast')

if length(varargin) > 0, title(varargin{1},'FontSize',24),end
%set(gcf,'Position',[464   283   930   827]), drawnow

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

if ~doind
    % do nothing

else
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
                if dolines
                    plot(xdat(j,:),dat(j,:),'k','LineWidth',.5,'Color',[.7 .7 .7]);
                end
            end
            
            % plot marker
            if donumber && ~(any(isnan(x(:, j))) || isnan(dat(j, i)))
                plot(x(:,j) + i,[dat(j,i) dat(j,i)]', 'w.'); % to set axis scale appropriately
                text(x(1,j) + i, dat(j,i), num2str(j), 'Color', [0 0 0], 'FontSize', 14)
                
            elseif ~(any(isnan(x(:, j))) || isnan(dat(j, i)))
                plot(x(:,j) + i,[dat(j,i) dat(j,i)]',mym,'Color',[0 0 0],'LineWidth',1,'MarkerFaceColor',myc)
            end
            
        end
        
        z = (dat(:,i) - mean(dat(:,i))) ./ std(dat(:,i)); wh = find(abs(z) >= 1.96);

        % print z-scores of potential outliers to text
        %if ~isempty(wh), fprintf(1,'\nCond %3.0f\t',i),for j = 1:length(wh), fprintf(1,'%3.2f\t',z(wh(j))); end,end

        if plotout
            if ~isempty(wh),plot(x(:,wh) + i,[dat(wh,i) dat(wh,i)]','ro','MarkerSize',14,'LineWidth',2),end
        end

    end

%     if dolines
%         for i = 1:size(dat,1)
%             plot(xdat(i,:),dat(i,:),'k','LineWidth',.5,'Color',[.7 .7 .7]);
%         end
%     end

if donumber && dojitter
    
    % Jitter positions
    han = findobj(gca, 'Type', 'text');
    
    x = get(han, 'Position'); x = cat(1, x{:});
    n = size(x, 1);
    x(:, 1) = x(:, 1) + .2 * (rand(n, 1) - .5);
    for i = 1:size(x, 1)
        set(han(i), 'Position', x(i, :));
    end
    
    % Bold
    set(han, 'FontWeight', 'b', 'Color', 'b');
    
end % jitter

        
end % plot individuals