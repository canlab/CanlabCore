function [clusters,vol,V,XYZ,Q] = mask_intersection2(clsize,outname,P,calcstr)
% empty P prompts for graphic selection of filenames
% extra arguments are more file names for 3 - n-way intersection
% empty outname prompts for entry of output img file name
%
% :Usage:
% ::
%
%     function [clusters,vol,V,XYZ,Q] = mask_intersection2(clsize,outname,P,calcstr)
%
% this version allows contrasts
% last argument, calcstr, is the string to evaluate in imcalc
%
% e.g., intersection of 4 images: 
% ::
%
%    'i1 & i2 & i3 & i4 & ~(isnan(i1) | isnan(i2) | isnan(i3) | isnan(i4))'
%
% e.g.,
% ::
%
%    'i1 & i2 & ~i3 & ~i4 & ~(isnan(i1) | isnan(i2)) = intersection of i1 and i2 and not i3 and not i4
%
% or
%
% ::
%
%    'i1 & i2 & isnan(i3) & isnan(i4) & ~(isnan(i1) | isnan(i2))'
%  ... for nan masked images
%
%     cl = mask_intersection2(5,'intersect001.img',P,'i1 & i2 & i3 & i4 & ~(isnan(i1) | isnan(i2) | isnan(i3) | isnan(i4))');
%
% ..
%    by Tor Wager
% ..

% ..
%    get file names and calcstr to evaluate
% ..

if isempty(P)
	P = spm_get(Inf,'*.img','Select binary mask images to intersect',pwd,0);

	if size(P,1) < 1, error('Must select at least 2 images'), end

	for i = 3:size(P,1)
		calcstr = [calcstr ' .* i' num2str(i)];
	end
end

disp(calcstr);
disp(P);

% make the intersection
% ------------------------------------------------
Q = spm_imcalc_ui(P,outname,calcstr);


% load the intersection mask file
% ------------------------------------------------
V = spm_vol(Q);
vol = spm_read_vols(V);

% apply cluster size threshold
% ------------------------------------------------
if isempty(clsize) || clsize == 0
    vol = double(vol);
else
    if clsize > 0
        [vol,numClusters,XYZ] = clusterSizeMask(clsize,vol);
        disp(['Found ' num2str(numClusters) ' clusters at k = ' num2str(clsize)])
    end
end


% write the size-filtered results back to the file
% ------------------------------------------------
V = spm_write_vol(V,vol);


% make a structure like SPM.mat which can be used to extract clusters
% ------------------------------------------------
CLU = mask2struct(V.fname,.5,clsize);		% should be 0 or 1, so .5 will get all 1's
clusters = tor_extract_rois([],CLU,CLU);
[d,f] = fileparts(V.fname);
if ~isempty(clusters)
    eval(['save ' f '_clusters clusters CLU'])
end

return

