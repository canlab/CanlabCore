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

testclust_plot_plugin;

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



