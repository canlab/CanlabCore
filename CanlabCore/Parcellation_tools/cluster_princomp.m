function [clusters,subclusters] = cluster_princomp(clusters,varargin)
% :Usage:
% ::
%
%     function [clusters,subclusters] = cluster_princomp(clusters,[behavioral score vector],[corr flag],[plotflag],[locflag])
%
% ALSO TRY: subcluster_montage(subclusters{1}) % to plot the output
%
% :Inputs:
%
%   **clusters:**
%        is structure of clusters from tor_extract_rois.m
%
%   behavioral vector is row vector of behavioral or other scores to correlate
%
%   corr flag:  *1 = work on correlations among voxels, 2 = work on covariance
%
%   plotflag:   *1 = yes, 0 = no.  plots.
%
%   locflag:    1 yes, *0 no; add XYZ voxel locations (scaled) to data submitted to clustering
%             pushes voxels closer in space to be classified in the same cluster
%
% try this to test the program on random data:
% ::
%
%    cl(1).all_data = randn(23,30);cl(1).numVox = 30;cl = cluster_princomp(cl,EXPT.behavior,1,1);
%    cl(1).all_data(:,1:10) = cl(1).all_data(:,1:10) + 10; cl = cluster_princomp(cl,EXPT.behavior,1,1);
%    cl(1).all_data(:,25:30) = cl(1).all_data(:,25:30) + repmat((EXPT.behavior .* 3)',1,6);
%    cl(1).all_data(:,21:24) = cl(1).all_data(:,21:24) + repmat((1:23)',1,4);
%    cl = cluster_princomp(cl,EXPT.behavior,1,1);
%
% mean-center everything now:
% ::
%
%    cl.PCA = [];
%    cl.all_data - cl.all_data - repmat(mean(cl.all_data),size(cl.all_data,1),1);
%    cl = cluster_princomp(cl,EXPT.behavior,1,1);
%
% add another correlated group:
% ::
%
%    cl.all_data(:,1:5) = cl.all_data(:,1:5) + repmat(rand(23,1)*5,1,5);
%    cl = cluster_princomp(cl,EXPT.behavior,1,1);
%
% if component scores are used and correlated with behavior, this means that the subjects
% tend to show the behavioral effect who also show the pattern associated with comp. x.  
% this may mean high on a number of voxels, or high on some and low on others.  
% the weights may be used to interpret what the components mean, and this can be done
% graphically.  
%
% t-tests on component scores have ambiguous interpretations, because a high t-score
% may indicate negative values or close-to-zero values on some voxels.
% a component could have the interpretation, "high on this component means high on V1
% and low on V2."  
%
% classifying voxels is done using cluster analysis (hierarchical, centroid linkage)
% on the voxels (observations) using the PCA weights (eigenvectors) as variables.
% This lets the clustering algorithm work in the reduced variable space with dimensionality
% equal to the number of components.  
% The max number of clusters is restricted based on the gradient of the eigenvalues in the PCA
% maxclusters = 1 + the number of eigenvalues with gradient at least 20% of the initial drop 
% from 1 to 2 eigenvalues.
% 
% Requires clustering library in Matlab.
% Robust option also uses the robust PCA algorithm RAPCA,
% created by:
% Hubert, M., Rousseeuw, P.J., Verboven, S. (2002),
%  "A fast method for robust principal components with applications to chemometrics", by Mia Hubert, Peter J. Rousseeuw, 
%  Chemometrics and Intelligent Laboratory Systems, 60, 101-111.
%
%

corrflag = 1; plotflag = 1; robustflag = 1; locflag = 0;
if length(varargin) > 1, corrflag = varargin{2};, end
if length(varargin) > 2, plotflag = varargin{3};, end
if length(varargin) > 3, locflag = varargin{4};, end

for i = 1:length(clusters)
    subclusters{i} = []; clusters(i).PCA = [];
    
    % check for NaNs and remove those voxels
    tst = any(isnan(clusters(i).all_data),1);
    clusters(i).all_data(:,tst) = [];
    clusters(i).XYZ(:,tst) = [];
    clusters(i).XYZmm(:,tst) = [];
    clusters(i).Z(:,tst) = [];
    if any(tst),disp(['Warning! Removed ' num2str(sum(tst)) ' voxels with NaN values.']),end
    
    a = clusters(i).all_data;
    
    % scale here, because robust pca doesn't use correlations, does it?
    % but robust PCA seems to be unaffected by scale changes on some variables
    if corrflag, a = scale(a);,end
        
    % if we choose to add the XYZ flag to add locations to clustering criteria
    if locflag,
        wfactor = round(size(a,1) ./ 3);                       % weight for loc; higher = more weight on location
        xyztmp = repmat(scale(clusters(i).XYZ')',wfactor,1);   % scale to make comparable to img values,
                                                                % but multiply by weighting factor
        a = [a; xyztmp];
    end
    
 
    if size(a,2) > 2    % must have 3 voxels to try clustering
    
    % -------------------------------------------------------------------------------
    % * All the real work is done here.  Compute pc's 
    % -------------------------------------------------------------------------------
    
    if ~robustflag
        [clusters(i).PCA.pcomps,clusters(i).PCA.weights,clusters(i).PCA.eigval,clusters(i).PCA.class] = pc(a,corrflag);
        
        % automatically pick number of clusters, based on gradient in eigenvalues
        g = abs(gradient(clusters(i).PCA.eigval));
        maxclusters = sum(g > g(1).*.2) + 1;    
    else
        
        % pick number of dimensions by hand
        disp(' ')
        disp([num2str(clusters(i).numVox) ' voxels in main cluster'])
        fprintf(1,'%3.2f Observations on %3.2f voxels\n',size(a,1),size(a,2))
        fprintf(1,'Save at least 2 eigenvectors to do clustering');
        % doing this on a' means voxels are observations, conditions/subj scores are variables
        % we'll use the scores, which has voxels as rows and components as columns, to classify
        out = rapca(a');
        clusters(i).PCA.pcomps = out.T;
        clusters(i).PCA.weights = out.P;
        clusters(i).PCA.eigval = out.L;
        maxclusters = length(out.L);
    end

    % -------------------------------------------------------------------------------
    % * All the real work is done here.  Classify
    % -------------------------------------------------------------------------------
    
    if size(clusters(1).all_data,1) > 12,
        disp('More than 12 dimensions (observations per voxel) in original data - using eigenvectors to classify')
        close all; try, pack, catch, end
        % we pick the number of CLASSES separately by hand, because # components not a good indicator
        % of how many classes there are
        clusters(i).PCA.class = docluster(out.T',[],plotflag);
    else
        clusters(i).PCA.class = docluster(out.T',[],plotflag);
    end
    
    % clean up if using locflag
    clusters(i).PCA.locflag = locflag;
    if locflag,
        clusters(i).PCA.pcomps = clusters(i).PCA.pcomps(1:size(clusters(i).all_data,1),:);
        %clusters(i).PCA.weights = clusters(i).PCA.weights(1:size(clusters(i).all_data,1),:);
        %clusters(i).PCA.avgs = clusters(i).PCA.avgs(1:size(clusters(i).all_data,1),:);
    end

    
    % -------------------------------------------------------------------------------
    % for each group, separate into contiguous clusters
    % -------------------------------------------------------------------------------
        
    grps = unique(clusters(i).PCA.class(clusters(i).PCA.class~=0)); % values are component of origin
    for j = 1:length(grps),
        wh = find(clusters(i).PCA.class==grps(j));
        XYZ = clusters(i).XYZ(:,wh);
        cl_index = spm_clusters(XYZ) ./ 100;
        clusters(i).PCA.class(wh) = clusters(i).PCA.class(wh) + cl_index;
    end
    ngrps = length(unique(clusters(i).PCA.class));
    fprintf(1,'%3.0f contiguous clusters separated by class',ngrps)
    
    % -------------------------------------------------------------------------------
    % average within classes / contiguous regions
    % -------------------------------------------------------------------------------
    
    grps = unique(clusters(i).PCA.class(clusters(i).PCA.class~=0)); % values are component of origin
    for j = 1:length(grps), 
        clusters(i).PCA.avgs(:,j) = mean(a(:,find(clusters(i).PCA.class==grps(j))),2);, 
        freq(j) = sum(clusters(i).PCA.class == grps(j));    
    end
    
    disp(['Cluster ' num2str(i) ', ' num2str(clusters(i).numVox) ' voxels: ' num2str(size(clusters(i).PCA.pcomps,2)) ' components'])
    fprintf(1,'\tMean\tEigval\tcorrel\t')
    
    % -------------------------------------------------------------------------------
    % * display each component and correlation with behavior
    % -------------------------------------------------------------------------------
    
    for j = 1:size(clusters(i).PCA.pcomps,2)
        
        %[H,P,CI,STATS] = TTEST(clusters(i).PCA.pcomps(:,j),0,.05,0);
        % skip the t-test.  t-tests on component scores don't make a lot of sense.
        fprintf(1,'\n\t%3.3f\t%3.3f\t',mean(clusters(i).PCA.pcomps(:,j)),clusters(i).PCA.eigval(j))
    
        if length(varargin) > 0
            if ~isempty(varargin{1})
                co = corrcoef(clusters(i).PCA.pcomps(:,j),varargin{1});
                co = co(1,2);
                fprintf(1,'%3.3f\t',co)
            end
        end
          
    end
    fprintf(1,'\n')
    
    % -------------------------------------------------------------------------------
    % * display classification info
    % -------------------------------------------------------------------------------
    disp(['Classified into ' num2str(ngrps) ' groups:'])
    fprintf(1,'\tClass\tVoxels\tMean\tt\tp\tcorrect. p\tcorrel\t')
    
    % for each component, test mean value and correlation with behavior
    for j = 1:length(grps)
        
        [H,P,CI,STATS] = TTEST(clusters(i).PCA.avgs(:,j),0,.05,0);
        fprintf(1,'\n\t%3.0f\t%3.0f\t%3.3f\t%3.3f\t%3.3f\t%3.3f\t',j,freq(j),mean(clusters(i).PCA.avgs(:,j)),STATS.tstat,P,P .* size(clusters(i).PCA.avgs,2))
    
        if length(varargin) > 0
            if ~isempty(varargin{1})
                co = corrcoef(clusters(i).PCA.avgs(:,j),varargin{1});
                co = co(1,2);
                fprintf(1,'%3.3f\t',co)
            end
        end
          
    end
    fprintf(1,'\n')
    
    % -------------------------------------------------------------------------------
    % * Plot, if requested
    % -------------------------------------------------------------------------------
    
    if plotflag,
        figure('Color','w'), subplot(1,3,1), imagesc(a), title(['Cl ' num2str(i) ': Data']), xlabel('Voxels'),ylabel('Subjects')
        subplot(1,3,2), imagesc(clusters(i).PCA.weights'), title(['Weights (eigenvectors)']), xlabel('Voxels'),ylabel('Eigenvectors')
        subplot(1,3,3), imagesc(clusters(i).PCA.pcomps), title(['Component scores (predictions)']), xlabel('Voxels'),ylabel('Subjects')
        
        a = [clusters(i).PCA.class' a']; a=sortrows(a,1); a = a(:,2:end)';
        figure;subplot 131; imagesc(a),title(['Cl ' num2str(i) ':Data sorted by class']), xlabel('Class'),ylabel('Subjects'),
        xlab = [sort(clusters(i).PCA.class(clusters(i).PCA.class~=0)) clusters(i).PCA.class(clusters(i).PCA.class==0)]; 
        set(gca,'XTick',1:length(clusters(i).PCA.class)); set(gca,'XTickLabel',xlab)
        subplot 132; imagesc(clusters(i).PCA.avgs),title('Class averages'), xlabel('Class'),ylabel('Subjects'),
        subplot 133; if length(varargin) > 0, if ~isempty(varargin{1}), imagesc(varargin{1}'), title('Behavior'),end,end
        
    end

    else    
        disp(['Cluster ' num2str(i) ' has less than 3 voxels.'])
        clusters(i).PCA.class = ones(1,clusters(i).numVox);
        clusters(i).PCA.avgs = clusters(i).timeseries;
        grps = 1;
        if ~isfield(clusters,'correl'), clusters(1).correl = [];, end
    end
    
    
    % -------------------------------------------------------------------------------
    % * separate into subclusters, based on class membership
    % -------------------------------------------------------------------------------

    if ~isfield(clusters,'correl'), clusters(i).correl = [];, end
    disp('Recomputing correlations and z-scores and saving in subclusters')
    if ~isfield(clusters,'XYZ'), clusters(i).XYZ = ones(3,size(clusters(i).all_data,2));, end
    if ~isfield(clusters,'Z'), clusters(i).Z = ones(1,size(clusters(i).XYZ,2));, end
    clear subc
    
    for j = 1:length(grps)
        
        try
            subc(j) = clusters(i);
        catch
            warning('clusters does not have all required fields.'); break
        end
        
        wh = find(clusters(i).PCA.class == grps(j));
        
        if length(varargin) > 0
            if ~isempty(varargin{1})
                co = corrcoef(clusters(i).PCA.avgs(:,j),varargin{1});
                subc(j).correl = co(1,2);
            end
        end
                 
        if isfield(subc(j),'title'), subc(j).title = [subc(j).title '_SUBCL_' num2str(j)];, end
        if isfield(subc(j),'name'), subc(j).name = [subc(j).name '_SUBCL_' num2str(j)];, end
        if isfield(subc(j),'numVox'), subc(j).numVox = length(wh);, end
        if isfield(subc(j),'Z'), subc(j).Z = subc(j).Z(wh);, end
        if isfield(subc(j),'XYZmm'), subc(j).XYZmm = subc(j).XYZmm(:,wh);, end
        if isfield(subc(j),'XYZ'), subc(j).XYZ = subc(j).XYZ(:,wh);, end
        
        if isfield(subc(j),'timeseries'), subc(j).timeseries = nanmean(subc(j).all_data(:,wh)')';, end
        %if isfield(subc(j),'snr'), subc(j).snr = subc(j).snr(wh);, end
        if isfield(subc(j),'center'), subc(j).center = center_of_mass(subc(j).XYZ,subc(j).Z);, end
        if isfield(subc(j),'mm_center'), subc(j).mm_center = center_of_mass(subc(j).XYZmm,subc(j).Z);, end
        
        if isfield(subc(j),'all_data'), 
            subc(j).all_data = subc(j).all_data(:,wh);, 
            for k = 1:size(subc(j).all_data,2)
                [H,P,CI,STATS] = ttest(subc(j).all_data(:,k),0,.05,0);
                subc(j).Z(k) = spm_t2z(STATS.tstat,STATS.df);
            end
        end
            
            
    end
    
    subclusters{i} = subc;
   
    
end


return



function [b,v,d,class] = pc(a,corrflag)
% a is original matrix, b is principal components, v is eigenvectors 
% (weights on columns, which = weights on voxels)
% class is classification of voxels into groups based on component loadings

if corrflag, [v,d]=eig(corrcoef(a));, else, [v,d]=eig(cov(a));,end
b = (pinv(v) * a')' ./ repmat((diag(d)').^.5,size(a,1),1);
% i made this up: think of rptating each subject's scores (in cols of a')
% by the rotation matrix pinv(v), and normalizing by the sqrt of the eigenvalues
% pinv(v) and v are rotation matrices because det = 1, no shearing or dilation
%
% this appears to work to give scores as well
% both methods (above,below) are scaled versions of the splus factor scores
% the problem is that doing it two different ways in splus flips the signs
% of some components and not others (gui vs cmd line).

%X = a; R = corrcoef(a); A = v * (d^.5); B = inv(R) * A;
%scores = X * B;
% X is data, A is factor loading matrix, B is factor score coeff matrix
% this method, from the text, and the one giving b above produce identical results

A = v * (d^.5);

b = fliplr(b); v = fliplr(v); A = fliplr(A); %scores = fliplr(scores);

num = min(10,sum(diag(d) >= 1));
b = b(:,1:num); v = v(:,1:num); A = A(:,1:num); 
origd = diag(d);
d = diag(d)'; d= fliplr(d); d = d(1:num);

if num == 0, warning('No eigenvalues above 1!');, origd, class = [];
    
else
    % classify each voxel into a group based on loading
    % use A, which re-introduces the comp variance, because we
    % want relationships with more variance to count more.
    % This just doesn't work so hot.  See docluster, below.
    
    %wh = A' == repmat(max(A'),size(A,2),1);
    %for i = 1:size(wh,2), tmp = find(wh(:,i));, class(i) = tmp(1); end
    %class(max(A') < .3) = 0;
    
end


%figure;plot(b,'r'),hold on;plot(a,'k'), hold on; plot(mean(a,2),'g--'),legend({'eig' 'orig' 'avg'})

return



function class = docluster(a,maxclusters,doplot)

    Y = pdist(a','euclid');     % transpose so the voxels are observations, eigenvectors the variables
    Z = linkage(Y,'complete');
    if maxclusters > 1
        class = cluster(Z,maxclusters)';
        
        if doplot, 
            dendrogram(Z,0); title('Dendrogram for clustering')
        end
    
    else
        dendrogram(Z,0); title('Dendrogram for clustering')
        set(gcf,'Position',[10   601   800   500])
        maxclusters = input('Pick number of classes to save: ');
        
        if maxclusters == 1,
            class = ones(1,size(a,2));
        else
            class = cluster(Z,maxclusters)';
        end
        
    end
    
    
    
return
    
