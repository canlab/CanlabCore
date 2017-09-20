function parcel_images(image_names, extract_mask, nuisance_covs)
% :Usage:
% ::
%
%     parcel_images(image_names, extract_mask, nuisance_covs)
%
%
% This function performs the following steps, in this order
%
%   - map mask to functional space
%   - extract data from in-mask voxels
%   - remove nuisance covariates (before assessing connectivity)
%   - data reduction (pca)
%   - plot cases (detect outliers)
%   - separate data into a priori anatomical regions (LBPA40 hard-coded
%     right now; downloadable from web; see wiki)
%     (save label image mapped to functional space)
%   - cluster voxels in each region to get parcels
%   - save parcels and images
%   - NMDS on the parcels to group them into "networks" (default = use rank data)
%   - Visualization of the networks
%
% :Outputs:
%
% Creates and goes to a new directory: parcel_images_output
%
% Outputs saved to disk include
%   1. An image with unique numerical codes for each parcel 
%   2. A 'clusters' structure containing the parcels, with image data extracted and
%      averaged over voxels within each parcel
%
% :Inputs:
%
%   **image_names:**
%        names of images to extract data from, and to use for
%        functional parcellation. SEE ALTERNATE FORMAT BELOW FOR DIRECT DATA INPUT
%
%   **extract_mask:**
%        a mask image
%
%   **nuisance_covs:**
%        columns of a matrix to remove from the data, so that this
%        subspace is not used to determine connectivity
%
%        e.g., nuisance_covs = SPM.xX.X(:, 1:3); if these are nuisance
%        covariates...
%
% Alternate input formats for image_names:
%
% If you have data already extracted, image_names can be a structure with
% these fields:
%
% image_names.V, spm_vol-style volume info for the mask volume and image space
% image_names.data, extracted data from all valid in-mask (non-zero,
% non-nan) voxels, one column per voxel, in standard matlab (:) order.
% The mask and the data must match!
%
% ..
%    Principal components + anatomical parcellation of data
%    Tor Wager, v. June 2010
% ..

%%%
% Notes:
% 06/19/2017 - SG implemented a check whether lansvd.m (sparse 
% SVD) is on the path. Download lansvd from https://github.com/areslp/matlab/tree/master/PROPACK
% lansvd is memory efficient, but estimation is not very stable. can result
% in NAN-estimations and the code will fail. If lansvd is not on the path,
% MATLAB's svd.m is called.
% Fixed voxel data extraction. Line 240 had set the current cluster data to zero,
% but then all the extracted data in the saved cluster structure were zero,
% too. commented out that line, as memory shouldn't be a problem nowadays.
% SG moved NMDS clustering block back into the main function
% to be run after the initial parcelation (Section 'Multivariate networks: cluster
% regions into networks'). This is a two-step clustering as in Kober et al
% 2008.


spm_defaults

docorrmtx = 0;  % 0 to work on cov matrix for PCA (good for meta-analysis), or 1 to work on correlation matrix (maybe better for other cases?)

if isempty(which('lansvd.m'))
    % can be downloaded from https://github.com/areslp/matlab/tree/master/PROPACK
    warning('lansvd.m not found on path. will use svd.m instead');
end

% Go to directory: output images will be written here
% ============================================================
[dd, currdir] = fileparts(pwd);
if strcmp(currdir, 'parcel_images_output')
    disp('Already in parcel_images_output directory. Proceeding...');
else
    if ~exist(fullfile(pwd, 'parcel_images_output'), 'file')
        disp('Creating directory: parcel_images_output');
        mkdir('parcel_images_output');
    end
    cd('parcel_images_output');
end

% ============================================================
% Map mask image to image space
% ============================================================
if isstruct(image_names)
    scn_map_image(extract_mask, image_names.V(1).fname, 'write', 'extract_mask.img'); 
elseif ischar(image_names)
    scn_map_image(extract_mask, image_names(1,:), 'write', 'extract_mask.img');
else
    error('image_names must be either a string array of names or a structure with vol info and data.  see the help page.');
end
extract_mask = fullfile(pwd, 'extract_mask.img');

% ============================================================
% Map anatomical label image to functional space and write
% ============================================================

labelimg = which('atlas_labels_combined.img'); %anat_lbpa_thal.img'); %% 'lpba40.spm5.avg152T1.label.nii');
disp(labelimg)

if ~exist(labelimg, 'file')
    cd ..
    fprintf('Cannot find label image: %s\n', labelimg)
    error('Must be on path');
end

% note: must be nearest neighbor interp to work!!
scn_map_image(labelimg, extract_mask, 'write', 'label_image_resliced.img');
labelimg = fullfile(pwd, 'label_image_resliced.img');  

% ============================================================
% Extract data
% ============================================================

 [dat, maskInfo] = extract_image_data(extract_mask, image_names);
 
 % memory saver: don't need image_names.data anymore
 if isstruct(image_names), image_names.data = []; end
 
% ============================================================
% Remove nuisance covs
% ============================================================
if ~isempty(nuisance_covs)
    disp('Removing nuisance covariates')
    dat = dat - nuisance_covs * pinv(nuisance_covs) * dat;
else
    disp('No nuisance covariates specified.')
end

%%
% ============================================================
% Data reduction
% ============================================================
mysize = prod(size(dat));
isbigmatrix = mysize > 10^ 8;

if ~isbigmatrix && any(isnan(dat(:)))
    disp('warning - replacing nans with zeros');
    dat(isnan(dat))=0;
end

if isbigmatrix
    disp('Running PCA within each anatomical parcel instead of at the start.');
    %          [U, sigma, pc] = lansvd(dat, nrun, 'L');
    %          sigma = diag(sigma);
    %          score = U .* repmat(sigma',nrun,1); % == x0*coeff ***test to make sure output variables are same as svds and debug
else

    [pc, score, latent] = run_pca(dat, docorrmtx);
    % Eigenvalue plot
    figure('Color','w'), set(gca,'FontSize',18), bar(latent)
    xlabel('Components'),ylabel('Variance'); drawnow


    num_to_save = sum(latent>1);
    disp(num_to_save)

    % ============================================================
    % Plot cases : diagnostic
    % ============================================================

    % Component plot of images (are there groups of points?
    % this suggests clusters of images/subjects)
    create_figure('nmdsfig'); nmdsfig(score(:, 1:2)); drawnow

    % clustering of cases
    % unusual cases (outliers) will have different class from most subjects
    subj_classes = clusterdata(score(:,1:num_to_save),'maxclust',10,'linkage','average');
    create_figure('nmdsfig'); nmdsfig(score(:, 1:2), 'classes', subj_classes);
    title('Plot of cases: Unusual cases are possible outliers');

end

%%
% ============================================================
% Read all anatomical labels
% ============================================================

[labelInfo, labels] = iimg_read_img(labelimg, 2);
clear labelInfo
labels = round(labels); % small errors in resampling (?)
in_mask_labels = labels(maskInfo.wh_inmask);
all_labels = unique(in_mask_labels);

in_mask_labels = uint16(in_mask_labels);
clear labels

% treat zero as a non-label
all_labels(all_labels == 0) = [];

parcel_labels = zeros(maskInfo.n_inmask, 1);

%%
% ============================================================
% Cluster within anatomical regions
% ============================================================
% Loop through anatomical labels, and cluster data based on similarity
% within each region.  'region' refers to an anatomical region specified a
% priori, and 'parcel' refers to a functionally homogenous subset of voxels
% in that region.

fprintf('Parcellating %3.0f anatomical regions. 000', length(all_labels));

for i = 1:length(all_labels)
    fprintf('\b\b\b%3.0f', i);

    current_label = all_labels(i);
    wh_in_region = in_mask_labels == current_label;
    if sum(wh_in_region)>1

    %wh_in_region_image_space = maskInfo.wh_inmask(wh_in_region);

    % view this cluster
    current_cl = iimg_indx2clusters(double(wh_in_region), maskInfo);
    cluster_orthviews(current_cl, {[1 0 0]}, 'solid');

    % if huge matrix: PCA here
    if isbigmatrix
        clear pc score latent ans
        vox_dat = dat(:, wh_in_region);
        
        % memory saver: zero out dat for these voxels in sparse matrix;
        % re-load later - % removed, because parcel data are then all-zero,
        % too. SG. 06/19/2017
%         dat(:, wh_in_region) = 0;

        [pc, score, latent] = run_pca(vox_dat, docorrmtx);  %[pc, score, latent] = run_pca
        clear vox_dat score
        num_to_save = sum(latent>1);
        vox_pcs = double(pc(:, 1:num_to_save));

    else
        % select principal component weights for this region and cluster
        vox_pcs = double(pc(wh_in_region, 1:num_to_save));
    end

    max_parcels = max(2, round(sum(wh_in_region) ./ 50));
    vox_classes = clusterdata(vox_pcs,'maxclust',max_parcels,'linkage','average');

    % assign unique value for each parcel
    % This may NOT work will all labeling schemes!!
    current_unique_id = repmat(i * 100, sum(wh_in_region), 1) + vox_classes;
    parcel_labels(wh_in_region) = current_unique_id;


    % save clusters structure for each parcel
    wh_in_region_index = find(wh_in_region);

    for j = 1:max(vox_classes)
        wh_in_parcel = wh_in_region_index(vox_classes == j);
        in_parcel = false(length(wh_in_region), 1);
        in_parcel(wh_in_parcel) = 1;

        tmp_cl = iimg_indx2clusters(in_parcel, maskInfo);
        % if more than 1, just choose largest for cluster structure
        % this will not match label image exactly (!)
        howmany = cat(1, tmp_cl.numVox); 
        mymax = find(howmany == max(howmany));
        mymax = mymax(1);
        
        parcel_cl{i}(j) = tmp_cl(mymax);
    end
    % visualize the parcels for this region
    cluster_orthviews(parcel_cl{i}, 'unique', 'add');
    spm_orthviews('Reposition', current_cl.mm_center)
    else
        !echo skipping this region - all_labels
        num2str(i)
    end

end

% ============================================================
% Save output
% ============================================================

% Image of parcels
% =====================================
% write image with unique labels
labelimg_out_name = 'parcel_labels.img';
if size(parcel_labels, 1) == 1, parcel_labels = parcel_labels'; end

if isstruct(image_names)
    % if struct, we have changed mask, need to go back to original one
    [dummy, maskInfo_orig] = iimg_get_data(extract_mask, extract_mask, 'single', 'noexpand');
    clear dummy
    iimg_reconstruct_vols(parcel_labels, maskInfo_orig, 'outname', labelimg_out_name);
else
    iimg_reconstruct_vols(parcel_labels, maskInfo, 'outname', labelimg_out_name);
end
disp(['Written: ' labelimg_out_name]);


% Clusters of parcels with data
% =====================================
parcel_cl_flat = cat(2, parcel_cl{:});

% extent: 3 vox 
parcel_cl_flat(cat(1, parcel_cl_flat.numVox) < 3) = [];

save parcel_clusters_file parcel_cl 
clear parcel_cl parcel_labels wh_in_reg*

% Write parcel images
nvox = maskInfo.nvox; 
imgdat = zeros(nvox, 1);
for i = 1:length(parcel_cl_flat)
    pvec = iimg_clusters2indx(parcel_cl_flat(i), maskInfo);
    imgdat(pvec) = i;
end
iimg_reconstruct_vols(imgdat, maskInfo, 'outname', 'parcels.img', 'descrip', 'Clusters within combined anatomical mask boundaries');
disp('Written: parcels.img (4-D image with mask volume for each parcel');


% memory saver: re-load here, after clustering
%[dat, maskInfo] = extract_image_data(extract_mask, image_names);

% extract original image data again, and average over voxels within parcel
if isstruct(image_names)
    % Need to test data extraction and implement!! Use, for each parcel::
    for i = 1:length(parcel_cl_flat)
        [imgvec, maskvec] = iimg_clusters2indx(parcel_cl_flat(i), maskInfo);
        parcel_cl_flat(i).all_data = dat(:, maskvec);
        parcel_cl_flat(i).timeseries = nanmean(parcel_cl_flat(i).all_data')';
    end
else
    parcel_cl_flat = tor_extract_rois(image_names, parcel_cl_flat);
end
disp('Extracted and averaged image data for each parcel, attached to parcel_cl_flat variable');

save parcel_clusters_file -append parcel_cl_flat
disp('Saved parcel_cl_flat in parcel_clusters_file.mat');

% memory saver: view parcels at the end
% alternative method for large datasets

clear dat

%% view all the parcels in unique colors - commented out to save time, SG
% if isbigmatrix
%      % loop clusters and display each cluster separately
%     cluster_orthviews(parcel_cl_flat(1), {rand(1, 3)});
%     for i = 2:length(parcel_cl_flat)
%         cluster_orthviews(parcel_cl_flat(i), {rand(1, 3)});
%     end
% else
    % just show all the clusters at once in different colors.
    cluster_orthviews(parcel_cl_flat, 'unique');
% end

sizes = cat(1, parcel_cl_flat.numVox);
create_figure('Voxels per parcel'); hist(sizes, 100);


%% 
% ============================================================
% Multivariate networks: cluster regions into networks
% ============================================================
c = [];
doranks = 1;
c = nmdsfig_tools('get_data_matrix',c, parcel_cl_flat,'timeseries',1,[],doranks);
c = nmdsfig_tools('get_correlations',c);
[c.GroupSpace,c.obs,c.implied_dissim] = shepardplot(c.D,[]);

disp('saving key info in c variable in nmds_c_structure.mat');
save nmds_c_structure c
% at least n parcels/2, and if n parcels/2 > 5, then at least 5 but up to
% nparcels/10
max_cl_to_test = min(max(5, round(length(parcel_cl_flat) ./ 10)), round(length(parcel_cl_flat)./2));
max_cl_to_test = round(length(parcel_cl_flat)./2);
c = nmdsfig_tools('cluster_solution',c, c.GroupSpace, 2:max_cl_to_test, 5000, {});
scn_export_papersetup(550);
saveas(gcf,'clusterPermutationTesting','png');

%
c.colors = cluster_orthviews_classes(parcel_cl_flat,c.ClusterSolution.classes, [], [], 1);
set(gcf,'position',[1 1 1680 1000]);
saveas(gcf,'network_montage','fig');
scn_export_papersetup(550);
saveas(gcf,'network_montage','png');
%%
% apply and make network plot
c = nmdsfig_tools('apply_clusters',c);

saveas(gcf,'network_nmds_fig','fig');
scn_export_papersetup(450);
saveas(gcf,'network_nmds_fig','png');

%nmdsfig_tools('nmdsfig_plot',c, 0, 0, 'fill');

disp('saving updated key info in nmds_c_structure.mat');
save nmds_c_structure c


% name networks
c.APPLY_CLUSTER.names = cell(1, length(c.APPLY_CLUSTER.classes));
for i = 1:length(c.APPLY_CLUSTER.classes)
    wh = find(c.ClusterSolution.classes == i);
    cluster_orthviews(parcel_cl_flat(wh), c.colors(i), 'solid');
    c.APPLY_CLUSTER.names{i} = input(sprintf('Name network %3.0f: ', i), 's');
end

disp('saving updated key info in nmds_c_structure.mat');
save nmds_c_structure c

% re-make to leave open for display
cluster_orthviews_classes(parcel_cl_flat,c.ClusterSolution.classes, [], [], 1);
scn_export_papersetup(700); saveas(gcf, 'cluster_montage_networks.png');

cd('..');

%%



end % main function




function [dat, maskInfo] = extract_image_data(extract_mask, image_names)
    % Extract, using single-precision array and sequential image loading to save space
    % for large datasets.  'noexpand' is a bit faster and is compatible with SPM5 and above.
    if ischar(image_names)
        [dat, maskInfo] = iimg_get_data(extract_mask, image_names, 'verbose', 'single', 'noexpand');
    else
        % struct
        maskindex = logical(iimg_get_data(image_names.V.fname, extract_mask, 2, 'noexpand')); % 1/0 for in-mask/not
%         [~,maskInfo]  = iimg_get_data(extract_mask,extract_mask,'noverbose', 'single', 'noexpand'); % read maskInfo too. SG 6/19/2017
        if issparse(image_names.data)
            dat = image_names.data(:, maskindex); % single(full(image_names.data(:, maskindex))); %image_names.data(:, maskindex);
        else
            dat = single(image_names.data(:, maskindex));
        end
        image_names.data = [];
        maskInfo = image_names.V;
        maskInfo.wh_inmask = maskInfo.wh_inmask(maskindex);
        maskInfo.n_inmask = length(maskInfo.wh_inmask);
        maskInfo.xyzlist = maskInfo.xyzlist(maskindex, :);
    end

end


function  [pc, score, latent] = run_pca(dat, docorrmtx)
    
    [nobs, nvars] = size(dat);
    ncomponents = min( [min(size(dat)) 50 ] );

    if docorrmtx, dat = zscore(dat); end
    
%     if nvars > 4000
        warning off % lots of "slow matlab code" warnings
        if isempty(which('lansvd.m')) % check for lansvd.m in path (SG)
            [U, sigma, pc] = svd(dat,0);
        else
            [U, sigma, pc] = lansvd(dat, ncomponents, 'L');
        end
        sigma = diag(sigma);
        score = U .* repmat(sigma', nobs, 1); % == x0*coeff
        sigma = sigma ./ sqrt(nobs-1);
        latent = sigma.^2;
        warning on

%     elseif issparse(dat)
%         % note: this can be slower and return different values from princomp. actually different? i think so
%         % also: needs testing/debugging.
%         nrun = min(n, 50);
%         [U, sigma, pc] = svds(dat, nrun); % could be too big -- no memory!! [U, sigma, pc] = svds(dat, n);
%         score = U .* repmat(diag(sigma)',n,1); % == x0*coeff
%         sigma = sigma ./ sqrt(n-1);
%         latent = sigma.^2;

%     else
% 
%         [pc, score, latent] = princomp(single(full(dat)), 'econ');
% 
%     end

end % run_pca

%%
