function [b,t,p,hh] = ancova(groups,x,y,varargin)
% :Usage:
% ::
%
%     b,t,p,pthandles] = ancova(groups,x,y,[plot],[covs of no interest])
%
% :Outputs:
%
%   **Elements of b, t, p:**
%        1st = intercept, 2nd = group effect, 3rd = slope, 4th = grp x slope
%        interaction
%
% recursive -- call ancova repeately if y is a matrix
% get pairwise standardized slopes (corrs)
% and group diffs and slope interaction for all 
% pairs of y vectors
%

hh = [];
doplot = 0; if length(varargin) > 0, doplot = varargin{1};,end

% covariates of no interest
if length(varargin) > 1, 
    covs = varargin{2};,

    [x,y,r,p,rrob,prob] = partialcor([x covs],y,1);
    % doesn't adjust groups now...
    
end





if size(y,2) > 1 & isempty(x)
    
    for i = 1:size(y,2)-1
        for j = i+1:size(y,2)
            
            [b,t,p] = ancova(groups,y(:,i),y(:,j));
            bb{1}(i,j) = b(2); tt{1}(i,j) = t(2); pp{1}(i,j) = p(2);  % group effect = cell 1
            bb{2}(i,j) = b(3); tt{2}(i,j) = t(3); pp{2}(i,j) = p(3);  % slope = cell 2
            bb{3}(i,j) = b(4); tt{3}(i,j) = t(4); pp{3}(i,j) = p(4);  % grp x slope = cell 3
        
        end
    end
    
    % clean up last row
    for i = 1:length(bb)
        bb{i}(end+1,:) = 0;  bb{i} = bb{i} + bb{i}';
        tt{i}(end+1,:) = 0;  tt{i} = tt{i} + tt{i}';
        pp{i}(end+1,:) = 0;  pp{i} = pp{i} + pp{i}';
    end
    
    b = bb;
    t = tt;
    p = pp;
    return
    
end


% -----------------------------------------
% basic ancova program
% -----------------------------------------

% center and format
x = [scale(groups,1) scale(x,1)];
%y = scale(y);

% add interaction in slopes (diff. slopes model)

x(:,3) = scale(x(:,1).*x(:,2));

[b,dev,stats]=glmfit(x,y);
bi = b(end); ti = stats.t(end); ppi = stats.p(end);

if ppi > .05
    
    % no sig. interaction, drop interaction and use parallel slopes
    
    [b,dev,stats]=glmfit(x(:,1:2),y);
    b(end+1) = bi; stats.t(end+1) = ti; stats.p(end+1) = ppi;
    
end

t = stats.t;
p = stats.p;


% -----------------------------------------
%  ancova plot
% -----------------------------------------


if doplot
    hh = [];
    colors = {'yo' 'm^'};
    [uni,b,grps] = unique(groups); 
    legstr = {['Group 1: ' num2str(uni(1))] ['Group 2: ' num2str(uni(2))]};
        uni = 1:max(grps);
    % median split, if continuous, high then low
    if length(uni) > 2,
        uni = [1 2]; grps = grps.*0;
        grps(groups>median(groups)) = 1; 
        grps(groups<median(groups)) = 2;
        legstr = {'High' 'Low'};
    end
    %tor_fig;
    for i = 1:length(uni)   % for each group, make a plot
        
        [tmp,tmp,tmp,hh(i)] = plot_correlation_samefig(x(find(grps==uni(i)),2),y(find(grps==uni(i))),[],colors{i},0,1);
    end
    xlabel('x'); ylabel('y'); 
    legend(hh,legstr)
end

    
return
