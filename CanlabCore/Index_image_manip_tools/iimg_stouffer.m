function cl = iimg_stouffer(Pimg,maskname,thr,clsize,outname)
% This function performs combination of p-values across images using the 
% Stouffer method (see refs below).  
%
% :Usage:
% ::
%
%     function cl = iimg_stouffer(Pimg,maskname,thr,clsize,outname)
%
% :Inputs:
%
%   **Pimg:**
%        string mtx of p-value image names
%
%   **k:**
%        num voxels/num comparisons 
%        if empty, uses # of non-zero, non-NaN values in images
%
% empty Pimg prompts for graphic selection of filenames
% empty outname prompts for entry of output img file name
%
% :Outputs:
%
%   **cl:**
%        clusters, cl{1} is stouffer, cl{2} is image 1, cl{3} is image 2
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
% ..
%    by Tor Wager
% ..

% ..
%    get file names and calcstr to evaluate
% ...

if isempty(Pimg),Pimg = spm_get(Inf,'*.img','Select p-value images',pwd,0);,end

if ~isempty(maskname)
    [dat, volInfo] = iimg_mask(maskname,Pimg);
else
    [volInfo,dat] = iimg_read_img(Pimg,1);
    dat = dat(volInfo.wh_inmask,:);
end

k = volInfo.n_inmask;
fprintf('\nSearch space is %3.0f voxels', k)

[z,p,bonf_sig] = stouffer(dat,.05 ./ k);  % returns bonferroni sig.


% individual images
nimgs = size(dat,2);
cl = cell(1,nimgs+1);

for i = 1:nimgs
    fprintf('\n\nImage %3.0f\n---------------------\n\n',i);
    [bonfp, fdrp, uncp] = image_stats(dat(:,i), k, thr);
    
    % output clusters
    psig = dat(:,i) <= uncp;
    cl{i+1} = iimg_indx2clusters(psig,volInfo);
end

fprintf('\n\nCombined p-values: Stouffer\n---------------------\n\n');
[bonfp, fdrp, uncp] = image_stats(p, k, thr);

% output clusters
psig = p <= uncp;
cl{1} = iimg_indx2clusters(psig,volInfo);

% old way of getting clusters
%voldata = iimg_reconstruct_3dvol(psig,volInfo);
%cl = mask2clusters(voldata,volInfo.mat);




return



function [bonfp, fdrp, uncp] = image_stats(p, k,uncp)

fprintf('\nMinimum p-value: %3.6f\n', min(p))

bonfp = .05 ./ k; bonfz = norminv(1-bonfp);
nbonf = sum(p <= bonfp);

fdrp = FDR(p,.05); 
if isempty(fdrp), fdrp = 0; end
fdrz = norminv(1-fdrp);
nfdr = sum(p <= fdrp);

uncz = norminv(1-uncp);
nunc = sum(p <= uncp);

fprintf('\nSignificant voxels:\tMethod\tZ thresh\tp thresh\tNumber\n')

fprintf('\tBonferroni\t%3.2f\t%3.6f\t%3.0f\n', bonfz, bonfp, nbonf)

fprintf('\tFDR\t%3.2f\t%3.6f\t%3.0f\n', fdrz, fdrp, nfdr)

fprintf('\tp < %3.4f\t%3.2f\t%3.6f\t%3.0f\n', uncp, uncz, uncp, nunc)
 
fprintf('\n')

return
