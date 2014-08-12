function [clpos,clneg] = mask_fisher(clsize,outname,k,Pimg,Timg)
% function [clpos,clneg] = mask_fisher(clsize,outname,k,Pimg,Timg)
% by Tor Wager
%
% This function performs combination of p-values across images using the 
% Fisher method (see refs below).  Voxels are thresholded with FDR, and can
% be positive and negative at once if both positive and negative sig
% effects exist in input images.  FDR correction is based on 2-tailed
% p-values - so this function assumes 2-tailed p-values!!  (robfit does
% this, glmfit and robustfit both return 2-tailed p-values as well.)
%
% Assumes input p-values are 1-tailed, and divides by two to get 2-tailed
% p-value inputs!!
%
% Inputs:
%   Pimg    string mtx of p-value image names
%   Timg    t-value img names, used to ensure signs are same across all tests
%   k       num voxels/num comparisons **NOT USED NOW - FDR USED INSTEAD
%           if empty, uses # of non-zero, non-NaN values in images
%
% empty i1 prompts for graphic selection of filenames
% extra arguments are more file names for 3 - n-way intersection
% empty outname prompts for entry of output img file name
%
% Described in:
% Lazar, N. A., Luna, B., Sweeney, J. A., & Eddy, W. F. (2002). 
% Combining brains: a survey of methods for statistical pooling 
% of information. Neuroimage, 16(2), 538-550.
%
% Stouffer, S. A., Suchman, E. A., DeVinney, L. C., Star, S. A., and
% Williams, R. M. 1949. The American Soldier: Vol. I. Adjustment
% During Army Life. Princeton University Press, Princeton.
% 
% Pimg=get_filename2('rob*\*_p_*0002.img'); Timg = get_filename2('rob*\*_tmap_0002.img')
% 

% get file names and calcstr to evaluate
% ------------------------------------------------

if isempty(Pimg),Pimg = spm_get(Inf,'*.img','Select p-value images',pwd,0);,end
if isempty(Timg),Timg = spm_get(Inf,'*.img','Select t-value images',pwd,0);,end

V =spm_vol(Pimg); v = spm_read_vols(V);
Vt =spm_vol(Timg); vt = spm_read_vols(Vt);

% this was used to make sure all values were of the same sign
% but if we use 2-tailed p-values, this isn't necessary.
% so we include an omnibus, and then these two as supplements.
pos = all(vt>0,4);
neg = all(vt<0,4);

% just include voxels in which one value at least is nonzero
%pos = any(vt>0,4); neg = any(vt<0,4);

% this was used before FDR correction was implemented - for bonferroni only
if isempty(k)
    mask = all(~isnan(v) & v ~= 0,4);
    k = sum(mask(:));
end

[z,p,sig,pt] = fisherp(v,.05);       %.05 ./ k);

disp('---------------------------------------')
disp('Output in mask_fisher_output.txt')
disp('---------------------------------------')
diary mask_fisher_output.txt

disp([inputname(2) ' Size threshold: ' inputname(1)])
disp(' ')
p = p(:); p(p==0 | isnan(p)) = []; 
fprintf(1,'\nMinimum combined p-value: %3.4f\n', min(p))

omnisig = z .* sig;
possig = z .* pos .* sig;
negsig = z .* neg .* sig;

tmp=omnisig(:); s = sum(~isnan(tmp) & tmp ~= 0); 
fprintf(1,'Omnibus: %3.0f voxels\n',s);

tmp=possig(:); s = sum(~isnan(tmp) & tmp ~= 0); 
tmp=negsig(:); sn = sum(~isnan(tmp) & tmp ~= 0); 
fprintf(1,'All betas Positive: %3.0f voxels, Negative: %3.0f voxels\n',s,sn);

fprintf(1,'Correction with Bonferroni, %3.0f voxels',k)
fprintf(1,'\nCorrection with FDR, p threshold = %3.4f, Z = %3.2f',pt,norminv(pt))
mynum = 2*size(Pimg,1);
tmp=z .* pos; tmp=tmp(:); tmp(tmp==0 | isnan(tmp)) = []; s = length(tmp); s2=max(tmp);
fprintf(1,'\nAll betas Positive: %3.0f eligible voxels with same signs, Max chisq(%2.0f) = %3.2f, p = %3.4f (corrected = %3.4f)\n',s,mynum,s2,1-chi2cdf(s2,mynum),(1-chi2cdf(s2,mynum)).*k);

tmp=z .* neg; tmp=tmp(:); tmp(tmp==0 | isnan(tmp)) = []; s = length(tmp); s2=max(tmp);
fprintf(1,'All betas Negative: %3.0f eligible voxels with same signs, Max chisq(%2.0f) = %3.2f, p = %3.4f (corrected = %3.4f)\n',s,mynum,s2,1-chi2cdf(s2,mynum),(1-chi2cdf(s2,mynum)).*k);

diary off

[d,f,e]=fileparts(outname);

tmp = fullfile(d,[f 'omni_fisher.img']); VV = V(1); VV.fname = tmp;
spm_write_vol(VV,omnisig);

warning off

if any(omnisig(:))
    
% make a structure like SPM.mat which can be used to extract clusters
% ------------------------------------------------
CLU = mask2struct(VV.fname,0,clsize);		% should be 0 or 1, so .5 will get all 1's
clusters = tor_extract_rois(Timg,CLU,CLU);
save fisher_omni_clusters clusters
try cluster_orthviews(clusters,{[1 0 0]}), catch,end

end


tmp = fullfile(d,[f 'pos_fisher.img']); VV = V(1); VV.fname = tmp;
spm_write_vol(VV,possig);

if any(possig(:))
% make a structure like SPM.mat which can be used to extract clusters
% ------------------------------------------------
CLU = mask2struct(VV.fname,0,clsize);		% should be 0 or 1, so .5 will get all 1's
clusters = tor_extract_rois(Timg,CLU,CLU);
clusters = chi2z(clusters,mynum);
[d,f] = fileparts(VV.fname);
if ~isempty(clusters)
    eval(['save ' f '_clusters clusters CLU'])
        montage_clusters([],clusters,{'r'});
        clpos = clusters;
end

end



tmp = fullfile(d,[f 'neg_fisher.img']); VV = V(1); VV.fname = tmp;
spm_write_vol(VV,negsig);


if any(negsig(:))
% make a structure like SPM.mat which can be used to extract clusters
% ------------------------------------------------
CLU = mask2struct(VV.fname,0,clsize);		% should be 0 or 1, so .5 will get all 1's
clusters = tor_extract_rois(Timg,CLU,CLU);
clusters = chi2z(clusters,mynum);
[d,f] = fileparts(VV.fname);
if ~isempty(clusters)
    eval(['save ' f '_clusters clusters CLU'])
    montage_clusters([],clusters,{'b'});
    clneg = clusters;
end
end

% apply cluster size threshold
% ------------------------------------------------
%[vol,numClusters,XYZ] = clusterSizeMask(clsize,vol);
%disp(['Found ' num2str(numClusters) ' clusters at k = ' num2str(clsize)])
diary on
if exist('clpos') == 1, fprintf(1,'\nPOSITIVE %s',clpos(1).title), cluster_table(clpos);,end
if exist('clneg') == 1, fprintf(1,'\nNEGATIVE %s',clneg(1).title), cluster_table(clneg);,end
diary off

warning on
return




function c = chi2z(c,df)

for i = 1:length(c)
    
    tmp = c(i).Z;
    tmp = 1 - chi2cdf(tmp,df);
    c(i).Z = norminv(1 - tmp);
    
end

return

