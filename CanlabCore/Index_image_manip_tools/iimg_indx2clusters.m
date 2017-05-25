% :Usage:
% ::
%
%     [cl, dat] = iimg_indx2clusters(dat, volInfo, [u], [k])
%
% cl is clusters
% volInfo is vol info, see iimg_read_img
%
% input dat is data from iimg_read_img, or any vectorized image
%
% output dat is vector of voxels in volInfo.xyzlist( in mask )
% whose elements are the cluster in cl they belong to (useful for adding
% fields later)
%
% :Optional Input:
%
%   **u:**
%        is height threshold
%
%        if single element, only dat values > u will be saved in clusters
%
%        if two elements, only dat values BETWEEN u(1) and u(2) will be saved
%
%        This convention is the same as that in iimg_threshold.
%
% ..
%    Tor Wager, updated 4/07 (documentation only)
%    Updated March 16, 2009. Fixed bug in cluster extent threhshold dat output
%    and bug in returning z-scores for matching voxels when extent threshold +
%    abs. is used
% ..

function [cl, dat] = iimg_indx2clusters(dat, volInfo, u, k)

    cl = [];

    %% check for fields
    if ~isfield(volInfo,'n_inmask')
        error('volInfo.n_inmask and other fields not found.  Use iimg_read_image.m with extended output flag to prepare volInfo.');
    end

    %% Prepare data from volinfo
    if nargin < 3
        u = 'unknown';
    end
        
    
    if size(dat,1) == volInfo.nvox
        z = dat(volInfo.wh_inmask);
    elseif size(dat,1) == volInfo.n_inmask
        z = dat;
    else
        error('data vector does not match size of image in volInfo!')
    end
    
    if nargin >= 3 && ~isempty(u)
        % impose threshold
        if length(u) == 2
            z(z < u(1) | z > u(2)) = 0; 
        else
            z(abs(z) < u) = 0; %8/5/15 Scott Schafer - Fixed sign issue that elimated negative activations
        end
    end
    
    xyz = double(volInfo.xyzlist)';  % enforce_variable_types compatibility: enforce double.  5/24/17
    XYZmm = voxel2mm(xyz, volInfo.mat);


    %% eliminate zero or NaN data values
    % And make sure we get spm cluster indices
    z(isnan(z)) = 0;
    wh = ~z;
    if any(wh)
        % we eliminate some voxels that are in-mask but don't have valid values
        z(wh) = [];
        xyz(:,wh) = [];
        XYZmm(:,wh) = [];
        clusters = get_cluster_index(xyz);
    elseif isfield(volInfo, 'clusters')
        clusters = volInfo.cluster;
    else
        clusters = spm_clusters(double(volInfo.xyzlist)')'; % enforce_variable_types compatibility: enforce double.  5/24/17
    end

    sig_regions = ~wh;
    
    %% save dat output with cluster indices in elements: which cl for each vox
    if nargout > 1
        dat = zeros(volInfo.n_inmask, 1); %double(dat);
        dat(sig_regions) = clusters;
    end

    %% Cluster size threshold, if specified
    if nargin > 3 && k > 1
        [nvox, dat, wh_omit_cl, wh_omit_vox] = iimg_cluster_sizes(clusters, dat, k, sig_regions);

        clusters(wh_omit_vox) = [];
        xyz(:,wh_omit_vox) = [];
        XYZmm(:,wh_omit_vox) = [];
        z(wh_omit_vox) = [];
        
        %% 2nd output stuff
%         if nargout > 1
%             for i = 1:length(wh_omit)
%                 dat(dat == wh_omit(i)) = 0;
%             end
%         end
    end

    %% Return if no voxels
    if isempty(clusters), return; end

    %% define cluster structure
    clnumbers = unique(clusters);   % some could be missing because they're too small
    nclust = length(clnumbers);
    cl = struct('title',[],'threshold',[],'M',[],'dim',[],'voxSize',[],'name',[],'Z',[],'XYZmm',[],'XYZ',[]);
    cl = repmat(cl,1,nclust);

    %% Fill in fields
    for i = 1:nclust
        % voxels in this cluster
        wh_incluster = find(clusters == clnumbers(i));
        nvox = length(wh_incluster);

        % re-number dat (cluster index) according to final output
        dat(dat == clnumbers(i)) = i;
        
        cl(i).title = sprintf('Cluster of %3.0f voxels from %s',nvox, volInfo.fname);

        % some included for backward compatibility
        cl(i).threshold = u;
        cl(i).M = volInfo.mat;
        cl(i).dim = volInfo.dim;
        cl(i).voxSize = diag(volInfo.mat(1:3,1:3))';
        cl(i).name = '';
        cl(i).numVox = nvox;
        cl(i).Z = z(wh_incluster)';


        cl(i).XYZmm = XYZmm(:,wh_incluster);
        cl(i).XYZ = xyz(:,wh_incluster);

        cl(i).mm_center = center_of_mass(cl(i).XYZmm,cl(i).Z);
    end
end




% Sub-functions

function clusters = get_cluster_index(xyz)
nvox = size(xyz,2);
if nvox == 0
    clusters = [];
    return;
    
    %  WORKS IN SPM5 -- CHANGING
    %     elseif nvox < 50000
else
    clusters = spm_clusters(xyz)';
    %     else
    %         clusters = ones(nvox,1);
end
end


% Count voxels in each contiguous cluster and threshold, if additional
% args are entered
function [nvox, dat, wh_omit_cl, wh_omit_vox] = iimg_cluster_sizes(clindx, dat, k, sig_regions)
    nclust = max(clindx);
    nvox = zeros(1,nclust);
    if nargin > 2
        wh_omit_cl = false(nclust,1);
        wh_omit_vox = false(length(clindx),1);
    end

    for i = 1:nclust
        wh = find(clindx == i);
        nvox(i) = length(wh);

        if nargin > 2
            if nvox(i) < k
                wh_omit_cl(i) = 1;
                wh_omit_vox(wh) = 1;
            end
        end
    end
    
    sr = find(sig_regions);
    dat(sr(wh_omit_vox)) = 0;
end

