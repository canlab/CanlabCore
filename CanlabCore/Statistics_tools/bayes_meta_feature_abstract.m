function [reducedY, cl, cluster_indx] = bayes_meta_feature_abstract(Y, thresh, shrink, bayes_model, volInfo, doplot)
% Purpose: Instead of original matrix of activations in studies x voxels,
% we may want to work with a data matrix of activations in contiguous
% regions x voxels
%
% This function calculates Y = 0,1 for Y = 1 in any voxel in each
% contiguous cluster
%
% :Examples:
% ::
%
%    reducedY = bayes_meta_feature_abstract(Y, .1, bayes_model, volInfo);
%    xval = classify_naive_bayes('xval', reducedY, Xi, 0, 0, bestk, bestg );
%    xval.prop_correct_by_class
%
% needs bayes_model.pa1_given_t, nclasses, params.k

fprintf(1, 'Getting likelihood. ');

nfeatures_orig = length(bayes_model.whkeep);

pa1_given_t = zeros(nfeatures_orig, bayes_model.nclasses);
pa1_given_t(bayes_model.whkeep,:) = bayes_model.pa1_given_t;

evidence_yes = (sum(Y(:, bayes_model.whkeep_from_notempty)) ./ bayes_model.nobs)';
evidence_yes = evidence_yes(:, ones(1, bayes_model.nclasses));
tmp = zeros(nfeatures_orig, bayes_model.nclasses);
tmp(bayes_model.whkeep,:) = evidence_yes;
evidence_yes = tmp;

lr = ( (pa1_given_t + shrink) ./ (evidence_yes + shrink) ) - 1;


% % pt = zeros(nfeatures_orig, bayes_model.nclasses);
% % pt(bayes_model.whkeep,:) = bayes_model.pt_given_act1;


lr(~bayes_model.whkeep, :) = 0;

if doplot, create_figure('Hist', 1, 2); hist(lr); plot_vertical_line(thresh); plot_vertical_line(-thresh); end
    
lr(lr < thresh & lr > -thresh) = 0;

if doplot, subplot(1, 2, 2); 
    tmp = lr(:);
    tmp(tmp == 0) = [];
    hist(tmp, 100); plot_vertical_line(thresh); plot_vertical_line(-thresh); 
    sum(lr > 0)
    sum(lr < 0)
end


% -------------------------------------------------------
fprintf(1, 'Getting clusters. ');

% get all clusters
% % % cl = iimg_indx2clusters(lr(:, 1), volInfo);
% % % 
% % % for i = 2:size(lr, 2)
% % %     cl = [cl iimg_indx2clusters(lr(:, i), volInfo)];    
% % % end
% % % 
% % % % reduce to remove overlap
% % % clu = clusters2CLU(cl);
% % % cl = tor_extract_rois([], clu, clu);
% this code does the same as the code above:
[cl, cluster_indx] = iimg_indx2clusters(any(lr')', volInfo);

fprintf(1, '%3.0f regions. ', length(cl))      
%
%cl = classify_naive_bayes('plot', bayes_model, 'lr', thresh, {[1 .7 0] [0 0 1]});
  
fprintf(1, 'Getting data in clusters. \n');

reducedY = Meta_cluster_tools('getdata',cl, Y',volInfo);
  
end
