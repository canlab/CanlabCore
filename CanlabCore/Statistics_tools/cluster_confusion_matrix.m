function stats = cluster_confusion_matrix(pred_label,true_label,varargin)
%Find an 'optimal' clustering of categories in a multi-class classification
%problem. Given a set of predicted and ground truth labels, this script
%will perform hierarchical clustering of confusions with 2 to N clusters,
%where N is the number of categories used in classification.
%
% Inputs:
%   pred_label - integer vector indicating model predictions
%   true_label - integer vector indicating ground truth labels
%
% Outputs:
% stats - array with information about performance; fields listed below
%     multi_way_accuracy - overall accuracy for N-way classification
%     accuracy_mean - vector of average accuracy for each clustering
%     ste_mean - vectore of average standard error for each clustering
%     p_vals -  cell array of p-values for all clusters for each clustering solution
%     clustered_multi_way_accuracy - overall accuracy for each clustering solution
%     sig - cell array indicating rejection of null hypothesis at alpha = .05 for each clustering solution
%     num_sig - number of significant clusters in each solution
%     num_clusters - number of clusters in each solution
%     optimalK - 'optimal' number of clusters out of all solutions (maximizes ratio of num_sig/num_clusters)
%     optimalY - clustering solution for optimal number of clusters
%
% Example:
%   tl=randi(20,300,1); %create ground truth labels
%   pl=tl; %create perfect match predictions
%   randvec=randperm(300); %create a random ordering of predictions
%   num_scrambled=150; %number of labels to 'corrupt'
%   pl(randvec(1:num_scrambled))=randi(20,1,num_scrambled); %assign random values
%   stats=optimal_clustering_confusions(pl,tl);
method = 'ward';
dofig = true;

label_input=strcmp(varargin,'labels');
if ~isempty(label_input)
    labels=varargin{label_input+1};
else
    labels=unique(pred_label);
end


%generate confusion matrix for clustering
cm=confusionmat(true_label,pred_label);
stats.multi_way_accuracy = sum(diag(cm))/sum(sum(cm));

cm = bsxfun(@rdivide,cm,sum(cm,2)); %normalize confusion matrix

% Create hierarchical cluster tree.
D = 1 - (squareform(tril(cm,-1)) + squareform(triu(cm,1)'))/2;
Z = linkage(D,method);

leafOrder = optimalleaforder(Z,D);



%%
for nc=2:size(cm,2)
    c = cluster(Z,'MaxClust',nc,'Criterion','distance'); %cluster labels based on confusion matrix
    
    clustered_pred_label = c(pred_label);
    clustered_true_label = c(true_label);
    clustered_cm = confusionmat(clustered_true_label,clustered_pred_label);
    
    clear RES;
    for k=1:size(clustered_cm,2)
        RES(k) = binotest([zeros(sum(clustered_cm(k,:))-clustered_cm(k,k),1); ones(clustered_cm(k,k),1)],sum(clustered_true_label==k)/length(clustered_true_label)); %#ok<AGROW>
    end
    
    stats.accuracy_mean(nc-1) = nanmean([RES(:).prop]);
    stats.ste_mean(nc-1) = nanmean([RES(:).SE]);
    stats.p_vals{nc-1}(:) = [RES(:).p_val];
    
    stats.clustered_multi_way_accuracy(nc-1) = sum(diag(clustered_cm))/sum(sum(clustered_cm));
    
    
    stats.sig{nc-1} = stats.p_vals{nc-1} < .05;
    stats.num_sig(nc-1) = sum(stats.sig{nc-1});
    
end
stats.num_clusters = 2:size(cm,2);
[~,ind] = max(fliplr([stats(1).num_sig]./stats(1).num_clusters));
stats.optimalK = 1+nc-ind;

% simply maximize number of sig clusters
%  [~,ind]=max([stats(1).num_sig]);
%  stats.optimalK=1+ind;


stats.optimalY = cluster(Z,'MaxClust',stats.optimalK,'Criterion','distance'); %cluster labels based on confusion matrix
clustered_pred_label = stats.optimalY(pred_label);
clustered_true_label = stats.optimalY(true_label);
stats.optimalCM = confusionmat(clustered_true_label,clustered_pred_label);
stats.optimalCM = bsxfun(@rdivide,stats.optimalCM,sum(stats.optimalCM,2)); %normalize confusion matrix

%% find threshold for producing optimal number of clusters
th = Z(end-stats.optimalK+2,3)-eps;




%%

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
