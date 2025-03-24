function [bestpval,bestmyclass,bestnames,bestX,where,clustnames,stats]=testclustnew(X,clust,varargin)
% This function tests for the significance of the cluster solution
%
% Usage
% ::
%
%     [bestpval,bestmyclass,bestnames,bestX,where,clustnames]=testclustnew(X,clust,[r],[nperm],[names],[remove],[linkagetype]);
%
% :Inputs:
%
%   **X:**
%        is the Group Space
%
%   **clust:**
%        is the cluster solutions (between how many and how many solutions is reasonable?)
%        the default for this option is 2:r/2 where r is the number of
%
%   **REGIONS:**
%        r is the number of DIMENSIONS in the solution (from choose_ndims)
%
%   **nperm:**
%        is the number of permutations for nonparametric testing
%         (default 1000)
%
%   **names:**
%        specifies the names of each region (object to be clustered)
%        if empty cell {}, will assign names
%
%   **remove:**
%        specifies what to do about elements which fit the cluster
%        solution poorly.  There are 3 possibilities:
%
%        'keep' - keeps all elements regardless of quality.  If you choose
%        this option, you are assuming that every region you enter
%        contributes something to your solution
%        'thresh' - this removes elements which fail to reach 95%
%        confidence as determined by random permutation testing
%
%        'iter'  - this option removes elements which fall in the bottom
%        5% of the permuted distribution, and remcomputes the solution
%        without these elements. This iterative recomputation means that
%        you are pruning your solution *until* you get a good one.  The
%        p-values which accompany each solution are thus difficult to
%        intepret if you use this option
%
%   **linkagetype:**
%        is the input to linkage ('single','average',etc) 
%
%        see help linkage_t
%
% :Outputs:
%
%   **bestpval:**
%        is the overall significance of the derived solution.
%        significance is calculated by permuting the group space in all
%        dimensions to derive a completely new configuration of elements in
%        n-dimensional space.  These elements are then assigned to clusters
%        exactly as for the point estimate solution, and the mean quality
%        (silhouette value) across the entire plot is used to form a
%        distribution against which to test the sihouette value of the true
%        solution
%
%   **bestmyclass:**
%        is the assignment of REGIONS to CLUSTERS for the best
%        soluton
%
%   **bestnames:**
%        contains the names of the REGIONS included in the
%        best solution (if remove==keep, this will be the same as names);
%
%   **bestX:**
%        is the best group space
%
%   **where:**
%        is a vector indexing which of the original regions (in
%        order) are included in the new group space.
%        clustnames is a cell containing the names of the elements in each
%        cluster
%
% :Calls: getmeanquality,clustquality,pdist1,linkage_t,cluster_t,makebinary
%
% Examples:
% -------------------------------------------------------------------------
% Generate some random data with 3 clusters in 10 variables:
% n_per_cluster = 33; n_vars = 10;
% 
% rng('default');  % For reproducibility
% X = [gallery('uniformdata',[n_per_cluster n_vars],12); ...
%     gallery('uniformdata',[n_per_cluster n_vars],13)+1.2; ...
%     gallery('uniformdata',[n_per_cluster n_vars],14)+2.5];
% y = [ones(n_per_cluster,1); 2*(ones(n_per_cluster,1)); 3*(ones(n_per_cluster,1))]; % Actual classes
% 
% % Test k = 2 - 7 clusters with 100 random permutations and make a plot of the results: 
% [bestpval,bestmyclass,bestnames,bestX,where,clustnames,stats]=testclustnew(X, [2:7], [], 100);



if length(varargin)>0,r=varargin{1}; else r = 2; end % inputs
if length(varargin)>1,nperm=varargin{2}; else nperm = 1000; end
if length(varargin)>2,names=varargin{3}; else names = cell(1,15); end
if length(varargin)>3,remove=varargin{4}; else remove='keep'; end
if length(varargin)>4,linkagetype=varargin{5}; else linkagetype = 'single'; end



XX=X;   %save a copy of the original group space in XX
Xnames=names;
num=size(X,1);  %thus is the number of elements at the start
bestpval = 1; beststdval = -Inf;
bestc = 0;

% initalize
bestnames = [];
bestmyclass = ones(num, 1);

fprintf(1,'testclustnew.m, linkage type is %s, Permutations: %d\n',linkagetype, nperm);

for c=clust  %loop through clusters
    fprintf(1,'permuting cluster %3.0f . ',c);
    X=XX;  % for each cluster solution, reset X
    names=Xnames; %and names;
    cont=0;     % if cont is 1, we stop iterating 
        
    while cont==0;
    [cq mcq myclass] = getmeanquality(X,c,linkagetype);   %point estimate values
    % cq: quality for each REGION
    % mcq: mean of cq
    % myclass: cluster assigments 
    % disp(['mean quality for cluster  ',num2str(c),' ......',num2str(mcq)]);    
    
    for p=1:nperm
        pX=reshape(shuffles(X(:)),[size(X)]);  % permute the group space values (shuffled across dimensions)
        [pcq pmcq(p)]=getmeanquality(pX,c,linkagetype);   %permuted values
        % pcq is quality  (discarded)
        % mean quality for each solution (pmcq) is entered into a
        % distribution
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% keep all REGIONS, including outliers %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if strcmp(char(remove),'keep')==1;
            cont=1;
            fprintf(1,' not removing elements . ');
            
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
    %%%% remove REGIONS in the bottom 5%, and recompute until none more removed %%%%%%%%%%%%%            
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif strcmp(char(remove),'iter')==1;
        try
            cutoff=prctile(pmcq,5);  %define 5% CI for permuted distribution
        catch
            cutoff=prctile_t(pmcq,5);  %define 5% CI for permuted distribution
        end
        
        if sum(cq>cutoff)==length(cq)
        disp(['found best solution for ',num2str(c),' clusters']);
        disp(['cutoff ',num2str(cutoff),' lowest quality ',num2str(min(cq))]);
        cont=1;
        else
        X=X(cq>=cutoff,:);
        fprintf(1,[' rejecting....',names{cq<cutoff}]);
        names=names(cq>=cutoff);
        myclass=myclass(cq>=cutoff,:);
        fprintf(1,[' keeping....',names{:}]);
        end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    %%%% 95% thresholding - remove REGIONS which don't reach this level %%%%%%%%%%%%%%%%%%%%    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif strcmp(char(remove),'thresh')==1;        
        disp('using 95% thresholding');
      try
          cutoff=prctile(pmcq,95);
      catch
          cutoff=prctile_t(pmcq,95);
      end
      cq=(cq>cutoff);
      fprintf(1,['rejecting....',names{cq<=cutoff}]);
      X=X(cq>cutoff,:);
      myclass=myclass(cq>cutoff,:);     
      names=names(cq>=cutoff);
      fprintf(1,['keeping....',names{:}]);
      cont=1;
    end     %ending if remove;
    
    %nonparametric significance for this solution
    pval = (nperm-(sum(mcq>pmcq)))/nperm;  %the p value for the point estimate
    fprintf(1,'pval = %3.4f\t',pval);
    stdval = (mcq-mean(pmcq))/std(pmcq);   % the number of st errors above the mean for the point estimate
    fprintf(1,'stdval = %3.2f\t',stdval);
    
    % if it's the best solution so far...
    if  pval < bestpval || pval==0     % if it's a better solution
        if stdval>beststdval       % if p==0, choose on the basis of SE above mean
            bestpval=pval;  %pval
            bestX=X;        %group space    
            
            % only replace classes if p < .05
            if pval < .05, bestmyclass = myclass; end     %cluster assignments
            
            bestnames = names; %names
            bestpmcq = pmcq;  % distribution of mean cluster quality for best
            bestcq = cq;      % quality * element
            bestmcq = mcq;    % mean quality
            bestc = c;        % which cluster
            beststdval = stdval;  %best stdval
         end   
    end
    
    end     %end while
    
    %save the best mean clustdr solution for each cluster for plotting
    cmcq(c)=mcq;  std_pcq(c)=std(pmcq);  cpval(c)=pval; cstdval(c)=stdval;
    indx = c - clust(1) + 1;
    try
        perm95(indx)=prctile(pmcq,95);
        perm05(indx)=prctile(pmcq,5);       
    catch
        perm95(indx)=prctile_t(pmcq,95);
        perm05(indx)=prctile_t(pmcq,5);
    end
    permuted_quality(indx)=mean(pmcq);
    
    fprintf(1,'\n');
    
end     %end loop through clusters

%%% use new and old versions of 'names' to find which ones have been
%%% rejected, if any.  This information is stored in 'where'.
where=zeros(1,length(Xnames));
for e1=1:length(Xnames);
    each=0;
    for e2=1:length(bestnames);
        if length(Xnames{e1})==length(bestnames{e2});
            if Xnames{e1}==bestnames{e2};
                each=1;
            end
        end
    end
    if each==1;
        where(e1)=1;
    end
end
where=find(where==1);
            

%generate clustnames;
if sum(cellfun(@isempty,Xnames)) == numel(Xnames)
    for i = 1:length(bestmyclass)
        Xnames{i} = ['v' num2str(i)];
    end
end

clustnames=cell(1,max(bestmyclass));
for n=1:max(bestmyclass)
    clustnames{n}=Xnames(bestmyclass==n);
end

% testclust_plot_plugin;

% -----------------------------------------------------------------------
% plot output
% -----------------------------------------------------------------------
%%%%%%%%%%bar graph
% scnsize = get(0,'ScreenSize');
% myposition = [50 50 scnsize(3)-100 scnsize(4)/2];
% myposition(3) = min(myposition(3), 1200);
% myposition(4) = min(myposition(4), 500);
% f1 = figure('position', myposition,'color','white'); set(gca,'FontSize',14),hold on
% subplot(1,3,1); set(gca,'FontSize',12);

f1 = create_figure('TestClust_Quality',1,3);


%%%------------------------- panel 1 ---------------------------------%%%
cmcq = cmcq(clust); %return to original format
cstdval = cstdval(clust);
bestc = bestc - clust(1)+1;
hold on;
plot(cmcq,'ko-','LineWidth',3,'MarkerFaceColor','k')
if exist('bestmcq', 'var')
    plot(bestc,bestmcq,'ks','LineWidth',3,'MarkerFaceColor',[.5 .5 .5],'MarkerSize',10)
end

plot(permuted_quality,'o-','Color',[.5 .5 .5],'MarkerFaceColor','k','LineWidth',2);
plot(perm95,'--','Color',[.5 .5 .5],'LineWidth',2)
plot(perm05,'--','Color',[.5 .5 .5],'LineWidth',2)

legend({'Real-data solutions' 'Chosen solution' 'Permuted solutions' '95% confidence'});
legend boxoff;
set(gca,'XTick',1:length(cmcq),'XTickLabel',clust,'XLim',[.5 length(clust)+.5])

%if ~isnan(coq)
%    set(gca,'Ylim',[0 max(coq)+max(cstd_poq)+.5]);
%end

title('Cluster quality')
xlabel('Number of clusters in solution')
ylabel('Mean silhouette value')
%%%------------------------- panel 1 ---------------------------------%%%


%%%------------------------- panel 2 ---------------------------------%%%
subplot(1,3,2); %set(gca,'FontSize',12);
hold on;
plot(cstdval,'ko-','LineWidth',3,'MarkerFaceColor','k')
plot(bestc,beststdval,'ks','LineWidth',3,'MarkerFaceColor',[.5 .5 .5],'MarkerSize',10)

%plot([bestc bestc],[min(cstdval) max(cstdval)-.01],'k:','Linewidth',2);

legend({'Real-data solutions' 'Chosen solution'});
legend boxoff;
set(gca,'XTick',1:length(cmcq),'XTickLabel',clust,'XLim',[.5 length(clust)+.5])

xlabel('Number of clusters in solution')
ylabel('Improvement over permuted data (s.d.)')
%%%------------------------- panel 2 ---------------------------------%%%


%%%------------------------- panel 3 ---------------------------------%%%
%%%%%%%%%%%%plot best histogram
subplot(1,3,3); %set(gca,'FontSize',12); hold on;
%for c=1:length(clust);

[h,x] = hist(bestpmcq,min(10,round(nperm ./ 40)));      %plot hist of best distribution
hh = bar(x,h);
set(hh,'FaceColor',[.3 .3 .3],'EdgeColor',[.3 .3 .3])

if exist('bestmcq', 'var')
    x = [bestmcq bestmcq];      %0:max(h);
    y = [0 max(h)];                  %ones(1,max(h)+1); %   *coq(best);

    plot(x,y,'color','k','LineWidth',3)
end

%set(gca,'XLim',min(h) - .1
title(['Permuted distribution for ',num2str(clust(bestc)),' clusters'])
ylabel('Mean silhouette value')
ylabel('Frequency')
%%%------------------------- panel 3 ---------------------------------%%%










stats.number_tested = clust;
stats.dims_of_space = size(X,2);
stats.linkagetype = linkagetype;
stats.silhouette = cmcq;
stats.permuted_sil = permuted_quality;     % average permuted quality
stats.perm95 = perm95;
stats.perm05 = perm05;
stats.permstd = std_pcq;
stats.pval = cpval;
stats.stdval = stdval;
stats.numclusters = clust(bestc);
stats.classes = bestmyclass;



end % function



% SUBFUNCTIONS

function Y = pdist(X,s,varargin)
%PDIST Pairwise distance between observations.
%   Y = PDIST(X) returns a vector Y containing the Euclidean distances
%   between each pair of observations in the M-by-N data matrix X.  Rows of
%   X correspond to observations, columns correspond to variables.  Y is an
%   (M*(M-1)/2)-by-1 vector, corresponding to the M*(M-1)/2 pairs of
%   observations in X.
%
%   Y = PDIST(X, DISTANCE) computes Y using DISTANCE.  Choices are:
%
%       'euclidean'   - Euclidean distance
%       'seuclidean'  - Standardized Euclidean distance, each coordinate
%                       in the sum of squares is inverse weighted by the
%                       sample variance of that coordinate
%       'cityblock'   - City Block distance
%       'mahalanobis' - Mahalanobis distance
%       'minkowski'   - Minkowski distance with exponent 2
%       'cosine'      - One minus the cosine of the included angle
%                       between observations (treated as vectors)
%       'correlation' - One minus the sample correlation between
%                       observatons (treated as sequences of values).
%       'hamming'     - Hamming distance, percentage of coordinates
%                       that differ
%       'jaccard'     - One minus the Jaccard coefficient, the
%                       percentage of nonzero coordinates that differ
%       function      - A distance function specified using @, for
%                       example @DISTFUN
%
%   A distance function must be of the form
%
%         function D = DISTFUN(XI, XJ, P1, P2, ...),
%
%   taking as arguments two L-by-N matrices XI and XJ each of which
%   contains rows of X, plus zero or more additional problem-dependent
%   arguments P1, P2, ..., and returning an L-by-1 vector of distances D,
%   whose Kth element is the distance between the observations XI(K,:)
%   and XJ(K,:).
%
%   Y = PDIST(X, DISTFUN, P1, P2, ...) passes the arguments P1, P2, ...
%   directly to the function DISTFUN.
%
%   Y = PDIST(X, 'minkowski', P) computes Minkowski distance using the
%   positive scalar exponent P.
%
%   The output Y is arranged in the order of ((1,2),(1,3),..., (1,M),
%   (2,3),...(2,M),.....(M-1,M)), i.e. the upper right triangle of the full
%   M-by-M distance matrix.  To get the distance between the Ith and Jth
%   observations (I < J), either use the formula Y((I-1)*(M-I/2)+J-I), or
%   use the helper function Z = SQUAREFORM(Y), which returns an M-by-M
%   square symmetric matrix, with the (I,J) entry equal to distance between
%   observation I and observation J.
%
%   Example:
%
%      X = randn(100, 5);                 % some random points
%      Y = pdist(X, 'euclidean');         % unweighted distance
%      Wgts = [.1 .3 .3 .2 .1];           % coordinate weights
%      Ywgt = pdist(X, @weucldist, Wgts); % weighted distance
%
%      function d = weucldist(XI, XJ, W) % weighted euclidean distance
%      d = sqrt((XI-XJ).^2 * W');
%
%   See also SQUAREFORM, LINKAGE, SILHOUETTE.

%   An example of distance for data with missing elements:
%
%      X = randn(100, 5);     % some random points
%      X(unidrnd(prod(size(X)),1,20)) = NaN; % scatter in some NaNs
%      D = pdist(X, @naneucdist);
%
%      function d = naneucdist(XI, XJ) % euclidean distance, ignoring NaNs
%      sqdx = (XI-XJ).^2;
%      nk = sum(~isnan(sqdx),2); nk(nk == 0) = NaN;
%      d = sqrt(nansum(sqdx')' .* size(XI,2) ./ nk); %correct for missing coords

%   Copyright 1993-2002 The MathWorks, Inc.
%   $Revision: 1.15 $

if nargin < 2
    s = 'euc';
    distfun = @distcalc;
    distargs = {s};
else
    if ischar(s)
        methods = strvcat('euclidean','seuclidean','cityblock','mahalanobis','minkowski','cosine','correlation','hamming','jaccard');
        i = strmatch(lower(s), methods);
        if length(i) > 1
            error(sprintf('Ambiguous ''DISTANCE'' argument:  %s.', s));
        elseif isempty(i)
            % error(sprintf('Unknown ''DISTANCE'' argument:  %s.', s));
            distfun = str2func(s);
            distargs = varargin;
            s = 'usr';
        else
            s = lower(methods(i,1:3));
            distfun = @distcalc;
            distargs = {s};
        end
    elseif isa(s, 'function_handle') |  isa(s, 'inline')
        distfun = s;
        distargs = varargin;
        s = 'usr';
    else
        error('The ''DISTANCE'' argument must be a string or a function.');
    end
end

[m, n] = size(X);
if any(imag(X(:))) & ~isequal(s,'usr')
   error('PDIST does not accept complex inputs.');
end

switch s
case 'seu' % Standardized Euclidean weights by coordinate variance
   distargs{end+1} = 1 ./ var(X)';
case 'mah' % Mahalanobis
   distargs{end+1} = inv(cov(X));
case 'min' % Minkowski distance needs a third argument
   if nargin < 3  % use default value for exponent
      distargs{end+1} = 2;
   elseif varargin{1} > 0
      distargs{end+1} = varargin{1}; % get exponent from input args
   else
      error('The exponent for the Minkowski metric must be positive.');
   end
case 'cos' % Cosine
   Xnorm = sqrt(sum(X.^2, 2));
   if min(Xnorm) <= eps * max(Xnorm)
       error(['Some points have small relative magnitudes, making them ', ...
              'effectively zero.\nEither remove those points, or choose a ', ...
              'distance other than cosine.'], []);
   end
   X = X ./ Xnorm(:,ones(1,n));
case 'cor' % Correlation
   X = X - repmat(mean(X,2),1,n);
   Xnorm = sqrt(sum(X.^2, 2));
   if min(Xnorm) <= eps * max(Xnorm)
       error(['Some points have small relative standard deviations, making ', ...
              'them effectively constant.\nEither remove those points, or ', ...
              'choose a distance other than correlation.'], []);
   end
   X = X ./ Xnorm(:,ones(1,n));
end

if m < 2
   % Degenerate case, just return an empty of the proper size.
   Y = zeros(1,0);
   return;
end

% Create (I,J) defining all pairs of points
p = (m-1):-1:2;
I = zeros(m*(m-1)/2,1);
I(cumsum([1 p])) = 1;
I = cumsum(I);
J = ones(m*(m-1)/2,1);
J(cumsum(p)+1) = 2-p;
J(1)=2;
J = cumsum(J);

% For large matrices, process blocks of rows as a group
n = length(I);
ncols = size(X,2);
blocksize = 1e4;                     % # of doubles to process as a group
M = max(1,ceil(blocksize/ncols));    % # of rows to process as a group
nrem = rem(n,M);
if nrem==0, nrem = min(M,n); end

Y = zeros(1,n);
ii = 1:nrem;
try
    Y(ii) = feval(distfun,X(I(ii),:),X(J(ii),:),distargs{:})';
catch
    if isa(distfun, 'inline')
        error(['The inline distance function generated the following ', ...
               'error:\n%s'], lasterr);
    elseif strfind(lasterr, ...
                   sprintf('Undefined function ''%s''', func2str(distfun)))
        error('The distance function ''%s'' was not found.', func2str(distfun));
    else
        error(['The distance function ''%s'' generated the following ', ...
               'error:\n%s'], func2str(distfun),lasterr);
    end
end;
for j=nrem+1:M:n
    ii = j:j+M-1;
    try
        Y(ii) = feval(distfun,X(I(ii),:),X(J(ii),:),distargs{:})';
    catch
        if isa(distfun, 'inline')
            error(['The inline distance function generated the following ', ...
                    'error:\n%s'], lasterr);
        else
            error(['The distance function ''%s'' generated the following', ...
                    'error:\n%s'], func2str(distfun),lasterr);
        end;
    end;
end

% ----------------------------------------------
function d = distcalc(XI,XJ,s,arg)
%DISTCALC Perform distance calculation for PDIST.
switch s
case 'euc',   d = sqrt(sum((XI-XJ).^2,2));            % Euclidean
case 'seu',   d = sqrt(((XI-XJ).^2) * arg);           % Standardized Euclidean
case 'cit',   d = sum(abs((XI-XJ)),2);                % City Block
case 'mah',   Y = XI - XJ;
              d = sqrt(sum((Y*arg).*Y,2));            % Mahalanobis
case 'min',   d = sum(abs((XI-XJ)).^arg,2).^(1/arg);  % Minkowski
case 'cos',   d = 1 - sum(XI.*XJ,2);                  % Cosine
case 'cor',   d = 1 - sum(XI.*XJ,2);                  % Correlation
case 'ham',   d = sum(XI ~= XJ,2) / size(XI,2);       % Hamming
case 'jac',   nz = XI ~= 0 | XJ ~= 0;
              ne = XI ~= XJ;
              d = sum(ne&nz,2) ./ sum(nz,2);          % Jaccard
end

end % distcalc

end % pdist1






function [cq mcq Xc]=getmeanquality(X,c,linkagetype)

    %point estimate cluster solution
    Y = pdist1(X);

    if exist('linkage.m', 'file') && exist('cluster.m', 'file')
        Z = linkage(Y,linkagetype);
        Xc = cluster(Z,c);
    else
        Z = linkage_t(Y,linkagetype);
        Xc = cluster_t(Z,c);
    end

    %convert X to X, binary input to doquality
    Xcx=makebinary(Xc)';

    %point estimate quality of each element;
    [cq center]=clustquality(Xcx,X);      %cq is cluster quality of each element

    %mean cq
    if all(cq == 1)
        mcq = 0;
    else
        mcq = mean( cq(cq~=1) );    %don't include single element clusters in your quality estimate
    end

end % getmeanquality







function [equality center] = clustquality(Xcx,X)
%[equality] = clustquality(Xcx,X)
%
% Xcx : binary indicator matrix of cluster assignments, 
% X   : stimulus coordinates in group space
% also: takes group spaces with zeros;


if size(Xcx>1);    
    for i = 1:size(Xcx,2)    % for each class
        tmp = mean(X(find(Xcx(:,i)),:),1);     % get center of this class
        center(i,:) = tmp;    
        % dist of all points in X from each center
        % number of cols of X is the number of dimensions; sums squared
        % vals across dims, takes sqrt to get Euclidean distance in N-d
        % space
        d(:,i) = sum((repmat(tmp,size(X,1),1) - X).^2,2).^0.5;
    end
    
    % rows of d are objects, columns are classes, values are dist from
    % class center
    for i=1:size(d,1)                   % i is the object (point)
        
        myclass = find(Xcx(i,:)==1);    % index of which class it is
        
        edist(i) = d(i,myclass);        %distance to center of own class
        
        otherclass = find(Xcx(i,:)==0);     % indices of columns for other classes
        
        otherdist(i) = min(d(i,otherclass));  %distance to center of cluster
        
    end

    % for each point, quality is distance to nearest neighbor - dist to own
    % cl / mak of those two
    equality = (otherdist - edist) ./ max([edist;otherdist]);
    
else
end

end % clustquality



function x=makebinary(y);
y=y(:)';
x=zeros(max(y),length(y));

for n=1:max(y);
    x(n,:)=y==n;
end
    
end

