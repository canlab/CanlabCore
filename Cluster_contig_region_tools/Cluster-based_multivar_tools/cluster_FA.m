function [clusters,m,mnames] = cluster_FA(clusters,OPT)
% [clusters,m,mnames]  = cluster_FA(clusters,OPT)
% Gets independent components for canonical timeseries values from clusters
%
% Requires OPT structure with the following fields (examples shown here):
%   O.HP = 100;
%   O.TR = 1.5;
%   O.doHP = 1;
%   O.doLP = 0;
%   O.scanadjust = 1;
%   O.percent = 0;
% 	O.filtertype = 'spm';
% 	O.nruns = 2;
%   O.adjustmatrix = custom adjustment matrix to regress out (e.g., movement params)
%   O.trimts = 3;
%
% Recommended now: average, or FA.  ICA can sometimes reverse sign, so can
% PCA potentially.  Not sure about sign thing.  ICA and PCA enforce
% orthogonal canonical variates.  
%
% Returns: clusters, with pcscore, ica, and varimax fields
%          average of voxels is stored in timeseries
%
% Maybe the best thing is to take voxels that load highly on ICs or PCs and
% return selective averages of those subsets - so gives subclusters.
%
% m is a matrix of rows = time (or observations), cols = timecourses within
% and across regions
%
% see also cluster_princomp.m

mnames = {};

for i = 1:length(clusters)
    
    m = clusters(i).all_data;
    m2 = roi_timeseries(m,OPT);
    
    close all
    
    clusters(i).pcscore = m2.pcscore;
    clusters(i).ica = m2.ica;
    clusters(i).varimax = m2.varimax;
    
    tmp = [clusters(i).title '_' num2str(i) '_' num2str(clusters(i).mm_center(1)) num2str(clusters(i).mm_center(2)) num2str(clusters(i).mm_center(3))];
    tmp(tmp==' ' | tmp=='.') = ['_'];
        
    for j = 1:size(clusters(i).varimax,2)
        mnames = [mnames {tmp}];
    end
end

m = cat(2,clusters.varimax);

return
