function [K, stats] = cluster_confusion_matrix(pred_label,true_label,varargin)
% Given a set of inferred and true multi-class category labels (e.g., from classification)
% group categories into clusters such that each cluster is statistically distinguishable
% from each other cluster.
% 
% Performs hierarchical clustering of confusions with 2 to N clusters (groups),
% where N is the number of categories used in classification.
%
% Usage:
% ::
%
%   [K, stats] = cluster_confusion_matrix(pred_label,true_label,varargin)
%
% ..
%     Author and copyright information:
%
%     Copyright (C) 2018 Philip Kragel
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
% ..
% :Inputs:
%   **pred_label** - integer vector indicating model predictions
%   **true_label** - integer vector indicating ground truth labels
%
% :Optional inputs:
%   **'dofig'** - make a plot of the dendrogram for clustering and confusion matrix
%   **'labels'** - followed by a cell array of strings specifying labels for each category
%   **'method'** - followed by a string specifying method for linkage; see
%   documentation for the linkage function
%   **'bootstrap'** - compute bootstrap confidence interval for the number of
%   clusters using ward method for linkage with 5000 samples
%   **'perm'** - permute labels to estimate p-value (two-sided) that the
%   optimal number of clusters is less than the full number of clusters
%   **'pairwise'** - instead of optimizing one-vs-all accuracy for each 
%   cluster, optimize on most confusable pair
%
%
% :Outputs:
% **k** - optimal number of clusters
% **stats** - structure with information about performance
%	Fields that contain performance metrics for different clustering solutions:
%       multi_way_accuracy - overall accuracy for N-way classification
%       accuracy_mean - vector of average accuracy for each clustering
%       ste_mean - vectore of average standard error for each clustering
%        p_vals -  cell array of p-values for all clusters for each clustering solution
%       clustered_multi_way_accuracy - overall accuracy for each clustering solution
%       sig - cell array indicating rejection of null hypothesis at alpha = .05 for each clustering solution
%       num_sig - number of significant clusters in each solution
%       num_clusters - number of clusters in each solution
%	Fields that contain metrics for optimal clustering solution:
%       optimalK - 'optimal' number of clusters out of all solutions (maximizes ratio of num_sig/num_clusters)
%       optimalY - clustering solution for optimal number of clusters
%       optimalAccuracy - multi-way accuracy on for k clusters (this is positively biased estimate)
%       optimalCM - confusion matrix for k clusters; rows are true cluster
%       labels, columns are predicted labels
%   Fields for making inference:
%       K_phat - permutation based p-value (two-sided) that k is less than
%       the number of categories in the full model
%       K_bootci - 95% bootstrap confidence interval on k
%
% :Examples:
% ::
%   % Basic analysis on random data
%   % Generate true class labels (tl), which are integers
%
%   tl=randi(20,300,1); %create ground truth labels

%   % Generate predicted class labels (pl), which is a random permutation
%   % of the true labels.  This is effectively what is done repeatedly during
%   % permutation testing.
%   pl=tl; %create perfect match predictions
%   randvec=randperm(300); %create a random ordering of predictions
%   num_scrambled=200; %number of labels to 'corrupt'
%   pl(randvec(1:num_scrambled))=randi(20,1,num_scrambled); %assign random values
%
%   % Create category labels:
%   categories = {}; for i = 1:20, categories{i} = sprintf('Cat%d', i); end
%
%   % Run the clustering and confusion matrix, generate a plot
%   [k,stats] = cluster_confusion_matrix(pl,tl,'dofig','method','complete');
%
%   % Run the full analysis with ward linkage - bootstrap CI, permutation test for inference, and plot of results:
%   [k,stats] = cluster_confusion_matrix(pl,tl,'labels', categories, 'dofig', 'method', 'ward', 'perm', 'bootstrap');
%
%
% :See also:
%
% linkage, dendrogram, confusionmat
% ..
%    Programmers' notes:
%    TODO: need to update bootstrap CI to accept other linkage methods, not just default of ward
%      add bootstrap confidence intervals for clustering
%
%
%

%% Parse optional inputs
% make a figure if specified
if any(strcmp(varargin,'dofig'))
    dofig = true;
else
    dofig = false; %by default
end

% get labels from user input or use defaults
label_input = strcmp(varargin,'labels');
if any(label_input)
    labels = varargin{label_input+1};
else
    labels = num2str(unique(pred_label));
end

% get linkage method
method_input = strcmp(varargin,'method');
if any(method_input)
    method = varargin{find(method_input)+1};
else
    method = 'ward';
end


perm_input = strcmp(varargin,'perm');
if any(perm_input)
    doperm = true;
else
    doperm = false;
end

boot_input = strcmp(varargin,'bootstrap');
if any(boot_input)
    doboot = true;
else
    doboot = false;
end


pair_input = strcmp(varargin,'pairwise');
if any(pair_input)
    dopairs = true;
else
    dopairs = false;
end

%% Compute confusion matrix, distance matrix, and related stats

% generate confusion matrix for clustering
cm=confusionmat(true_label,pred_label);
stats.multi_way_accuracy = sum(diag(cm))/sum(sum(cm));

cm = bsxfun(@rdivide,cm,sum(cm,2)); %normalize confusion matrix

% Create hierarchical cluster tree.
D = 1 - (squareform(tril(cm,-1)) + squareform(triu(cm,1)'))/2; %distances
Z = linkage(D,method); % tree

if dofig
leafOrder = optimalleaforder(Z,D); %reorder for plotting if necessary
stats.leafOrder=leafOrder;
end
%% Search possible clustering solutions and compute stats

for nc=2:size(cm,2)
    c = cluster(Z,'MaxClust',nc,'Criterion','distance'); %cluster labels based on confusion matrix
    
    % get new labels and compute new confusion matrix
    clustered_pred_label = c(pred_label);
    clustered_true_label = c(true_label);
    clustered_cm = confusionmat(clustered_true_label,clustered_pred_label);
    
    % do a binomial test based on each row of confusion matrix
    clear RES;
    for k=1:size(clustered_cm,2)
        
        if ~dopairs
            %overall accuracy
            RES(k) = binotest([zeros(sum(clustered_cm(k,:))-clustered_cm(k,k),1); ones(clustered_cm(k,k),1)],sum(clustered_true_label==k)/length(clustered_true_label)); %#ok<AGROW>
            stats.accuracy_null{nc-1}(k)=sum(clustered_true_label==k)/length(clustered_true_label);

        else
            %accuracy relative to the most confused cluster
            tv=clustered_cm(k,:);
            tv(k)=nan;
            [maxval,ind]=max(tv);
            RES(k) = binotest([zeros(maxval,1); ones(clustered_cm(k,k),1)],sum(clustered_true_label==k)/sum(clustered_true_label==k|clustered_true_label==ind)); %#ok<AGROW>
            stats.accuracy_null{nc-1}(k)=sum(clustered_true_label==k)/sum(clustered_true_label==k|clustered_true_label==ind);

        
        end 
    end
    
    % compute mean, ste, and p-value for each cluster
    stats.accuracy_mean(nc-1) = nanmean([RES(:).prop]);
    stats.ste_mean(nc-1) = nanmean([RES(:).SE]);
    stats.p_vals{nc-1}(:) = [RES(:).p_val];
    
    % compute overall accuracy
    stats.clustered_multi_way_accuracy(nc-1) = sum(diag(clustered_cm))/sum(sum(clustered_cm));
    
    % identify the number of significant categories (.05 uncorrected) - also test direction 
    stats.sig{nc-1} = stats.p_vals{nc-1} < .05 & [RES(:).prop]-stats.accuracy_null{nc-1} >0;
    stats.num_sig(nc-1) = sum(stats.sig{nc-1});
    
end

%% Compute performance metrics and store in output structure
% store number of clusters in stats output
stats.num_clusters = 2:size(cm,2);

% find the cluster solution that maximizes the proportion of significant
% clusters (largest number of clusters will be selected for ties)

if ~dopairs
[~,ind] = max(fliplr([stats(1).num_sig]./stats(1).num_clusters));
stats.optimalK = 1+nc-ind;
else
% find the cluster solution where all clusters are significant
ind=max(find(stats.num_sig==stats.num_clusters )); %#ok<MXFND>
% store optimal K in output
stats.optimalK = 1+ind;
end
% create variable foor bootstrapping calls (first output of this function)
K=stats.optimalK;

% store multi-way accuracy for optimal clustering
stats.optimalAccuracy=stats.clustered_multi_way_accuracy(stats.optimalK-1);

% store predictions for optimal clustering solution
stats.optimalY = cluster(Z,'MaxClust',stats.optimalK,'Criterion','distance'); %cluster labels based on confusion matrix
clustered_pred_label = stats.optimalY(pred_label);
clustered_true_label = stats.optimalY(true_label);

% create normalized confusion matrix for optimal clustering solution
stats.optimalCM = confusionmat(clustered_true_label,clustered_pred_label);
stats.optimalCM = bsxfun(@rdivide,stats.optimalCM,sum(stats.optimalCM,2)); %normalize confusion matrix

%% find threshold for producing optimal number of clusters
th = Z(end-stats.optimalK+2,3)-eps;

%% Perform permutation test if specified

if doperm
    stats.K_phat = permutation_test_num_clusters(pred_label,true_label,stats.multi_way_accuracy,stats.optimalK);
end

%% Compute 95% bootstrap confidence interval on k if specified

if doboot
    
    p = gcp();  %check if we can do things in parallel
    if isempty(p)
        poolsize = 0;
    else
        poolsize = p.NumWorkers;
    end
    
    if poolsize>1 % If we have a parallel pool, use it!
        opts=statset('useparallel',true);
        stats.K_bootci=bootci(5000,{@cluster_confusion_matrix,pred_label,true_label},'type','norm','Options',opts);
        
    else
        
        stats.K_bootci=bootci(5000,{@cluster_confusion_matrix,pred_label,true_label},'type','norm');
    end
end



%% Plot results if requested

if dofig
    figure;
    subplot(1,2,1)
    [H,~,perm] = dendrogram(Z,'labels',labels,'reorder',fliplr(leafOrder),'colorThreshold',th,'orientation','left');
    
    
    hold on;
    plot([th th],[min(get(gca,'YLim')) max(get(gca,'YLim'))],'linewidth',1,'linestyle','--','color',[0 0 0])
    set(H,'Linewidth',1.5,'color',[.5 .5 .5])
    set(gca,'XTickLabelRotation',90)
    xlabel Distance
    title 'Optimal Clustering Solution'
    
    clrs=seaborn_colors(stats.optimalK);
    cluster_labels=stats.optimalY(perm);
    
    [b,m1] = unique(cluster_labels,'first');
    [~,d1] =sort(m1);
    b = b(d1);
    
    
    
    x_locations=1:length(cluster_labels);
    for i=1:stats.optimalK
        plot([min(get(gca,'XLim')) min(get(gca,'XLim'))],[min(x_locations(cluster_labels==i))-.5 max(x_locations(cluster_labels==i))+.5],'color',clrs{i},'linewidth',8)
    end
    
    set(gca,'Ylim',[.5 length(cluster_labels)+.5])
    axis square
    b=flipud(b);
    subplot(1,2,2)
    imagesc(stats.optimalCM(b,b));colorbar;
    colormap(flipud(gray))
    set(gca,'XTick',[]);
    set(gca,'YTick',[]);
    
    locations=1:max(cluster_labels);
    hold on;
    for i=1:max(locations)
        plot([locations(i)-1 locations(i)]+.5,.5+[max(get(gca,'YLim')) max(get(gca,'YLim'))],'color',clrs{b(i)},'linewidth',8)
        plot(-.5+[min(get(gca,'XLim')) min(get(gca,'XLim'))],[locations(i)-1 locations(i)]+.5,'color',clrs{b(i)},'linewidth',8)
    end
    
    set(gca, 'XLim',[-.25 max(locations)+1])
    set(gca,'YLim',.5+[-.25 max(locations)+.75])
    set(gca,'visible','off')
    h=title('Normalized Confusion Matrix');
    text(h.Position(1),h.Position(2)-.45,'Normalized Confusion Matrix','FontWeight','Bold','FontSize',11,'HorizontalAlignment','center');
    axis square
    
    set(gcf,'Color',[1 1 1])
    set(gcf,'Units','inches')
    set(gcf,'Position',[2 2 11 3])
end

end

%% subfunction for doing permutation test
function phat = permutation_test_num_clusters(pred_label,test_label,overall_accuracy,k)

nit=5000;

num_samples=length(pred_label);
num_scrambled=floor(num_samples*(1-overall_accuracy)); %number of labels to 'corrupt'

parfor it=1:nit
    tl=randi(max(test_label),num_samples,1); %create ground truth labels
    pl=tl; %create perfect match predictions
    randvec=randperm(num_samples); %create a random ordering of predictions
    toscramble=pl(randvec(1:num_scrambled));
    pl(randvec(1:num_scrambled))=toscramble(randperm(length(toscramble))); %permute to get random values
    [~,sim_stats(it)]=cluster_confusion_matrix(pl,tl);
end

phat = 1-(sum(k<[sim_stats(:).optimalK]))/(nit+1);
phat = phat*2; %two sided
end

