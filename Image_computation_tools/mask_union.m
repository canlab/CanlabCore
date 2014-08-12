function [vol,V,XYZ,clusters,Q] = mask_union(clsize,outname,i1,i2,varargin)
% function [vol,V,XYZ,clusters,Q] = mask_union(clsize,outname,i1,i2,varargin)
% by Tor Wager
%
% empty i1 prompts for graphic selection of filenames
% extra arguments are more file names for 3 - n-way intersection
% empty outname prompts for entry of output img file name
%

if isempty(clsize), clsize = 0;,end
XYZ = []; vol = []; V = []; clusters = [];

% get file names and calcstr to evaluate
% ------------------------------------------------

calcstr = '(i1 > 0) | (i2 > 0)';
%calcstr2 = ' & (~isnan(i1) & ~isnan(i2)';

if isempty(i1),
	P = spm_get(Inf,'*.img','Select binary mask images to get union of',pwd,0);

	if size(P,1) < 1, error('Must select at least 2 images'), end

	for i = 3:size(P,1)
		calcstr = [calcstr ' | (i' num2str(i) ' > 0)'];
        %calcstr2 = [calcstr2 ' & ~isnan(i' num2str(i) ')'];
	end
else
	P = str2mat(i1,i2);
end


for i = 1:length(varargin)
	P = str2mat(P,varargin{i});
	calcstr = [calcstr ' | (i' num2str(i+2) ' > 0)'];
    %calcstr2 = [calcstr2 ' & ~isnan(i' num2str(i+2) ')'];
end

% make the intersection (union, in this case)
% ------------------------------------------------
%calcstr = [calcstr ')'];
%calcstr2 = [calcstr2 ')'];

Q = spm_imcalc_ui(P,outname,[calcstr]);
V = spm_vol(Q);

% apply cluster size threshold
% ------------------------------------------------
if clsize > 0
    % load the intersection mask file
    % ------------------------------------------------
    [vol,hdr] = readim2(Q(1:end-4));
    [vol,numClusters,XYZ] = clusterSizeMask(clsize,vol);
    disp(['Found ' num2str(numClusters) ' clusters at k = ' num2str(clsize)])
    
    % write the size-filtered results back to the file
    % ------------------------------------------------
    vol = double(vol);
    V = spm_write_vol(V,vol);
end





% make a structure like SPM.mat which can be used to extract clusters
% ------------------------------------------------
%CLU = mask2struct(V.fname,.5,clsize);		% should be 0 or 1, so .5 will get all 1's
%clusters = tor_extract_rois([],CLU,CLU);
clusters = mask2clusters(V.fname);
[d,f] = fileparts(V.fname);
if ~isempty(clusters)
    eval(['save ' f '_clusters clusters'])
end

return

