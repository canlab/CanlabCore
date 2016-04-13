function han = barplot_grouped(dat,X,xnames,seriesnames,varargin)
% :Usage:
% ::
%
%    han = barplot_grouped(dat,X,xnames,seriesnames, [optional args]);
%
% :Inputs:
%
%   **dat:**
%        n x 4 data matrix (2 x 2 grouping only for now, but easy to
%        expand later)
%
%   **X:**
%        covariates, no intercept; will be centered by this function
%
%   **xnames:**
%        x-axis labels
%
%   **seriesnames:**
%        series labels
%    
% First two bars are group, and last two bars are group
%    
% :Optional Inputs: (keywords)
%
%   **'within':**
%        within-error flag.  1 = errors based on subject x
%        condition interaction
%
%   **'stars':**
%        put stars for significance on graph (default)
%
%   **'nostars':**
%        do not plot stars
%
%   **'bars':**
%        followed by number of bars in group
%
%   **'pvals':**
%        followed by matrix of p-values (for stars; will calculate if missing)
%
%   **'inputmeans':**
%        input means and errors in first 2 inputs rather than data
%        and X
%
%   **'colors':**
%        followed by colors, e.g., mycol = {[1 0 0] [0 1 0] [1 0 1] [1 1 0] [0 0 1]}

    if isempty(dat), disp('Nothing to plot.'), return, end
    
    dowithin = 0;
    dostars = 1;
    nbars = 2;      % bars per group
    pvals = [];
    inputmeans = 0;

    for i = 1:length(varargin)
        if isstr(varargin{i})
            switch varargin{i}
                % reserved keywords
                case 'within', dowithin = 1;
                case 'stars', dostars = 1;
                case 'nostars', dostars = 0;
                case 'inputmeans', inputmeans = 1; m = dat; se = X;

                    % functional commands
                case 'bars', nbars = varargin{i+1};
                case 'pvals', pvals = varargin{i+1};

                case 'colors',mycolors = varargin{i+1};
                    
                otherwise, warning(['Unknown input string option:' varargin{i}]);
            end
        end
    end


    if inputmeans
        [ngroups,nbars] = size(m);
        
        n = 204; %input('Enter sample size: ');
        
        if isempty(pvals) && dostars
            disp('Must enter p-values for stars.'), dostars = 0; 
        else
            p = pvals'; p = p(:)'; pvals = p;
        end

        % rearrange into vector format for compatibility
        b = m'; b = b(:)';
        s = se'; s = s(:)'; se = s;

    else
        % we have data
        ncols = size(dat,2);
        if mod(ncols,nbars) ~= 0
            error('number of data columns must be divisible by number of bars per group.  try ''bars'' input option.');
        end

        ngroups = ncols ./ nbars;
        n = size(dat, 1);
        
        X = scale(X,1); X = [ones(n,1) X];    % design matrix
        b = pinv(X) * dat;                      % b(1,:) are means, removing centered covariates
        resid = dat - X*b; sigma = std(resid);

        Xerr = sqrt(diag(inv(X'*X)));
        sterr = Xerr * sigma;           % rows are predictors, cols are data columns

        m = b(1,:); se = sterr(1,:);    % for intercept; other rows are for covariates

        if isempty(pvals)
            pvals = 2 * (   1 - tcdf( abs(b(1,:)./se),size(X,1)-size(X,2) )   );
        end

        if dowithin
            se = barplot_get_within_ste(dat);
            se = repmat(se,1, ngroups .* nbars);
        end

        %tor_fig;
        m = reshape(m,nbars,ngroups)';    % rows are groups; plot bars in same order as original means
    end

    han = bar(m,'grouped');
    
    set(han(1),'FaceColor',[.2 .2 .2]);
    if length(han) > 1, set(han(2),'FaceColor',[.8 .8 .8]); end
    if exist('mycolors','var') && iscell(mycolors) && ~isempty(mycolors(1))
        for i =1:length(han), set(han(i),'FaceColor',mycolors{i}); end
    end
    
    if length(han) > 1
        set(gca,'XTick',1:ngroups,'XTickLabel',xnames);
        legend(han, seriesnames);
    else
        set(gca,'XTick',1:nbars,'XTickLabel',xnames);
    end
    
    xlocs = bar_centers(han);
    xlocs = sort(xlocs(:));
    
    tor_bar_steplot(b(1,:), se, {'k'}, xlocs)
    
    %tor_bar_steplot(b(1,:),se,{'k'},.35,.5,.2);

    if dostars
        star_plot(b(1,:),pvals,se,nbars,ngroups);
    end

    % set y-axis limits
    mx = max(b(1,:)+se); mi = min(b(1,:)-se);
    sc = .1 * (mx-mi);
    yl = [mi-sc mx+sc];
    set(gca,'YLim',yl);


    return


function star_plot(b,pvals,se,nbars,ngroups)

    xoffs = [-.15 .10];    % x offset, neg then pos
    starwid = .06;

    for mybar = 1:nbars
        y = b(1,mybar:nbars:end);
        p = pvals(mybar:nbars:end);
        s = se(mybar:nbars:end);

        yval = sign(y) .* ( abs(y) + abs(s) + abs(s) .* .1 );

        % adjust neg y-values for height of char
        yadj = zeros(size(y));
        yadj(yval < 0) = -.02;
        yval = yval + yadj;

        xval = (1:ngroups) + xoffs(mybar);

        % adjust for width of star
        xval = xval - .5*starwid;

        for i = 1:length(yval)
            if p(i) < .0015, mystr = '***';      xadj = -starwid*1.5;
            elseif p(i) < .015, mystr = '**';    xadj = -starwid;
            elseif p(i) < .055, mystr = '*';     xadj = 0;
            else mystr = ''; xadj = 0;
            end
            sh(i) = text(xval(i)+xadj,yval(i),mystr,'FontSize',24);
        end
    end

    return




