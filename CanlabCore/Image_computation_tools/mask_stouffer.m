function [vol,V,XYZ,clusters,Q] = mask_stouffer(clsize,outname,k,Pimg,Timg)
% This function performs combination of p-values across images using the 
% Stouffer method (see refs below).  
%
% :Usage:
% ::
%
%     function [vol,V,XYZ,clusters,Q] = mask_stouffer(clsize,outname,k,Pimg,Timg)
%
% :Inputs:
%
%   **Pimg:**
%        string mtx of p-value image names
%
%   **Timg:**
%        t-value img names, used to ensure signs are same across all tests
%
%   **k:**
%        num voxels/num comparisons 
%        if empty, uses # of non-zero, non-NaN values in images
%
% empty i1 prompts for graphic selection of filenames
% extra arguments are more file names for 3 - n-way intersection
% empty outname prompts for entry of output img file name
%
% :Described in:
% Lazar, N. A., Luna, B., Sweeney, J. A., & Eddy, W. F. (2002). 
% Combining brains: a survey of methods for statistical pooling 
% of information. Neuroimage, 16(2), 538-550.
%
% Stouffer, S. A., Suchman, E. A., DeVinney, L. C., Star, S. A., and
% Williams, R. M. 1949. The American Soldier: Vol. I. Adjustment
% During Army Life. Princeton University Press, Princeton.
%
% :Examples:
% ::
%
%    Pimg=get_filename2('rob*\*_p_*0002.img');
%    Timg = get_filename2('rob*\*_tmap_0002.img')
%
% ..
%    by Tor Wager
% ..

% ..
%    get file names and calcstr to evaluate
% ..

if isempty(Pimg),Pimg = spm_get(Inf,'*.img','Select p-value images',pwd,0);,end
if isempty(Timg),Timg = spm_get(Inf,'*.img','Select t-value images',pwd,0);,end

V =spm_vol(Pimg); v = spm_read_vols(V);
Vt =spm_vol(Timg); vt = spm_read_vols(Vt);

pos = all(vt>0,4);
neg = all(vt<0,4);

if isempty(k)
    mask = all(~isnan(v) & v ~= 0,4);
    k = sum(mask);
end

[z,p,sig] = stouffer(v,.05 ./ k);

p = p(:); p(p==0 | isnan(p)) = []; 
fprintf('\nMinimum combined p-value: %3.4f\n', min(p))

possig = z .* pos .* sig;
negsig = z .* neg .* sig;

tmp=possig(:); s = sum(possig(~isnan(possig)));
tmp=negsig(:); sn = sum(negsig(~isnan(negsig)));
fprintf(1,'Positive: %3.0f voxels, Negative: %3.0f voxels\n',s,sn);


tmp=z .* pos; tmp=tmp(:); tmp(tmp==0 | isnan(tmp)) = []; s = length(tmp); s2=max(tmp);
fprintf(1,'\nPositive: %3.0f eligible voxels with same signs, Max Z = %3.2f, p = %3.4f\n',s,s2,normcdf(1-s2));

tmp=z .* neg; tmp=tmp(:); tmp(tmp==0 | isnan(tmp)) = []; s = length(tmp); s2=max(tmp);
fprintf(1,'Negative: %3.0f eligible voxels with same signs, Max Z = %3.2f, p = %3.4f\n',s,s2,normcdf(1-s2));


s = sum(possig(~isnan(possig)));
tmp=negsig(:); sn = sum(negsig(~isnan(negsig)));
fprintf(1,'\nPositive: %3.0f voxels, Negative: %3.0f voxels\n',s,sn);

[d,f,e]=fileparts(outname);
tmp = fullfile(d,[f '_stoufZ.img']); VV = V(1); VV.fname = tmp;
spm_write_vol(VV,possig);

% make a structure like SPM.mat which can be used to extract clusters
% ------------------------------------------------
CLU = mask2struct(VV.fname,0,0);		% should be 0 or 1, so .5 will get all 1's
clusters = tor_extract_rois([],CLU,CLU);
[d,f] = fileparts(VV.fname);
if ~isempty(clusters)
    eval(['save ' f '_clusters clusters CLU'])
end


tmp = fullfile(d,[f '_stoufZ_pos_thresh.img']); VV = V(1); VV.fname = tmp;
spm_write_vol(VV,negsig);

% make a structure like SPM.mat which can be used to extract clusters
% ------------------------------------------------
CLU = mask2struct(VV.fname,0,0);		% should be 0 or 1, so .5 will get all 1's
clusters = tor_extract_rois([],CLU,CLU);
[d,f] = fileparts(VV.fname);
if ~isempty(clusters)
    eval(['save ' f '_clusters clusters CLU'])
end


% apply cluster size threshold
% ------------------------------------------------
%[vol,numClusters,XYZ] = clusterSizeMask(clsize,vol);
%disp(['Found ' num2str(numClusters) ' clusters at k = ' num2str(clsize)])




return

