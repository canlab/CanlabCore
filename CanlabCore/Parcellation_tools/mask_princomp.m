function [clout] = mask_princomp(clusters,varargin)
% function [clusters] = mask_princomp(clusters,[behavioral score vector],[corr flag],[plotflag],[saveflag])
% 
% this function is just like cluster_princomp, except that it works on the SET
% of activations in all clusters, rather than within each cluster.
%
% clusters is structure of clusters from tor_extract_rois.m
% 
% behavioral vector is row vector of behavioral or other scores to correlate
% corr flag:  1 = work on correlations among voxels, 2 = work on covariance
% plotflag:   1 = yes, 0 = no.  plots.
%
%
% try this to test the program on random data:
% cl(1).all_data = randn(23,30);cl(1).numVox = 30;cl = cluster_princomp(cl,EXPT.behavior,1,1);
% cl(1).all_data(:,1:10) = cl(1).all_data(:,1:10) + 10; cl = cluster_princomp(cl,EXPT.behavior,1,1);
% cl(1).all_data(:,25:30) = cl(1).all_data(:,25:30) + repmat((EXPT.behavior .* 3)',1,6);
% cl(1).all_data(:,21:24) = cl(1).all_data(:,21:24) + repmat((1:23)',1,4);
% cl = cluster_princomp(cl,EXPT.behavior,1,1);
% mean-center everything now:
% cl.PCA = []; cl.all_data - cl.all_data - repmat(mean(cl.all_data),size(cl.all_data,1),1);
% cl = cluster_princomp(cl,EXPT.behavior,1,1);
% add another correlated group:
% cl.all_data(:,1:5) = cl.all_data(:,1:5) + repmat(rand(23,1)*5,1,5);
% cl = cluster_princomp(cl,EXPT.behavior,1,1);
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

corrflag = 1; plotflag = 1; robustflag = 1; saveflag = 0;
if length(varargin) > 1, corrflag = varargin{2};, end
if length(varargin) > 2, plotflag = varargin{3};, end
if length(varargin) > 3, saveflag = varargin{4};, end

% get all the voxels and data together

CLU = clusters2CLU(clusters);
CLU.all_data = cat(2,clusters.all_data);

clear clusters; 
clusters(1) = CLU;
i = 1;
subclusters{1} = []; clusters(i).PCA = [];
    
disp('Output saved in mask_princomp.out')
diary mask_princomp.out

a = CLU.all_data;
    
    if size(a,2) > 2    % must have 3 voxels to try clustering
    
    % -------------------------------------------------------------------------------
    % * All the real work is done here.  Compute pc's and clustering
    % -------------------------------------------------------------------------------
    
    if ~robustflag
        [clusters(i).PCA.pcomps,clusters(i).PCA.weights,clusters(i).PCA.eigval,clusters(i).PCA.class] = pc(a,corrflag);
        
        % automatically pick number of clusters, based on gradient in eigenvalues
        g = abs(gradient(clusters(i).PCA.eigval));
        maxclusters = sum(g > g(1).*.2) + 1;    
    else
        
        % pick number of clusters by hand
        out = rapca(a);
        clusters(i).PCA.pcomps = out.T;
        clusters(i).PCA.weights = out.P;
        clusters(i).PCA.eigval = out.L;
        maxclusters = length(out.L);
        
        if saveflag, 
            try mkdir mask_pca_results, catch, end
            saveas(gcf,'mask_pca_results/robust_pca_diagnostic','fig'), close
            saveas(gcf,'mask_pca_results/robust_eigenvalues','fig'), close
        end
        
    end

    clusters(i).PCA.class = docluster(a,maxclusters,plotflag);
    
    grps = unique(clusters(i).PCA.class(clusters(i).PCA.class~=0)); % values are component of origin
    for j = 1:length(grps), 
        clusters(i).PCA.avgs(:,j) = mean(a(:,find(clusters(i).PCA.class==grps(j))),2);, 
        freq(j) = sum(clusters(i).PCA.class == grps(j));    
    end
    
    disp([clusters(i).title ', ' num2str(clusters(i).numVox) ' voxels: ' num2str(size(clusters(i).PCA.pcomps,2)) ' components'])
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
    disp(['Classified into ' num2str(max(clusters(i).PCA.class)) ' groups:'])
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
        
        if saveflag, 
            try mkdir mask_pca_results, catch, end
            disp('Images saved in mask_pca_results folder')
            saveas(gcf,'mask_pca_results/PCA1','fig'),
        end
       
        a = [clusters(i).PCA.class' a']; a=sortrows(a,1); a = a(:,2:end)';
        figure;subplot 131; imagesc(a),title(['Cl ' num2str(i) ':Data sorted by class']), xlabel('Class'),ylabel('Subjects'),
        xlab = [sort(clusters(i).PCA.class(clusters(i).PCA.class~=0)) clusters(i).PCA.class(clusters(i).PCA.class==0)]; 
        set(gca,'XTick',1:length(clusters(i).PCA.class)); set(gca,'XTickLabel',xlab)
        subplot 132; imagesc(clusters(i).PCA.avgs),title('Class averages'), xlabel('Class'),ylabel('Subjects'),
        subplot 133; if length(varargin) > 0, if ~isempty(varargin{1}), imagesc(varargin{1}'), title('Behavior'),end,end
        
        if saveflag, saveas(gcf,'mask_pca_results/PCA_with_clustering','fig'),end
    end

    else    
        disp(['Cluster ' num2str(i) ' has less than 3 voxels.'])
        clusters(i).PCA.class = ones(1,clusters(i).numVox);
        clusters(i).PCA.avgs = clusters(i).timeseries;
        grps = 1;
        if ~isfield(clusters,'correl'), clusters(1).correl = [];, end
    end
    
    diary off
    
    % -------------------------------------------------------------------------------
    % * separate into subclusters, based on class membership
    % -------------------------------------------------------------------------------
    
    for j = 1:length(grps)
        
        wh = find(clusters(i).PCA.class == grps(j));
        
        if length(varargin) > 0
            if ~isempty(varargin{1})
                co = corrcoef(clusters(i).PCA.avgs(:,j),varargin{1});
                subc(j).correl = co(1,2);
            end
        end
        
        CLUcomp = clusters(1);
        CLUcomp.XYZmm = CLUcomp.XYZmm(:,wh);
        CLUcomp.XYZ = CLUcomp.XYZ(:,wh);
        CLUcomp.Z = CLUcomp.Z(:,wh);
        CLUcomp.all_data = CLUcomp.all_data(:,wh);
        CLUcomp.numVox = size(CLUcomp.all_data,2);
        
        clout{j} = tor_extract_rois([],CLUcomp,CLUcomp);
        if plotflag
            montage_clusters([],clout{j})
            if saveflag, saveas(gcf,['mask_pca_results/montage_subset' num2str(j)],'fig'),end
        end
        
    end
    
   
    
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


    if size(a,2) > 300,
        Y = pdist_by_parts(a');
    else
        Y = pdist(a','euclid');     % transpose so the voxels are observations, eigenvectors the variables
    end
    Z = linkage(Y,'centroid');
    class = cluster(Z,maxclusters)';
    
    if doplot, 
        dendrogram(Z,0); title('Dendrogram for clustering')
    end
    
return


function pd =  pdist_by_parts(x)

disp(['Computing distances between observations'])

indx = 1;
vindx = 2;

for i = 1:size(x,1)
    indx = 1;   % faster this way
    for j = vindx:size(x,1)
        pd{i}(indx) = (sum((x(i,:) - x(j,:)).^2)).^.5;
        indx = indx+1;
    end
    vindx = vindx+1;
    if mod(i,100)==0,fprintf(1,'.'),end
end

pd = cat(2,pd{:});
fprintf(1,'done\n')

return


    