function [mask_thresh, cl, inmaskvox, dat, outputname] = fmri_mask_thresh_canlab(fmri_file, outputname, implicit_masking_method, plotfigs)
% Implicit determination of which voxels are in-brain, based on the intensities of
% functional images.  Assumes much (most) of the image has near-zero
% background noise values, and the in-brain values are substantially
% higher.
%
% :Usage:
% ::
%
%     [mask_thresh, cl, inmaskvox, in_mask_logical_vector, maskfilename] = fmri_mask_thresh_canlab(fmri_file, outputname)
%
% :Inputs:
%
%   **fmri_file:**
%        is either a list of file names or an fmri_data object
%
%   **File names:**
%        a (preferably) 4-D file of imaging data, Analyze .img or .nii
%
%   **fmri_data object:**
%        With multiple images loaded with *no* mask
%
%   **outputname:**
%        is a mask file output name, e.g., 'mask.img', with .img
%        extension. Empty [] means do not write output image.
%
% :Implicit_masking_method:
%
%   **mean:**
%        take the top 95% of voxels above the mean value.  used by
%    default if no value is entered 
%
%   **dip:**
%        smooth the histogram and take the top 95% of values above the
%    first positive gradient
% 
%
%   **plotfigs
%     [1/0]: enable or suppress mask display and orthviews
%
% :Outputs:
%
%   **mask_thresh:**
%        signal-value above which voxels are considered in brain
%
%   **c1:**
%        clusters, from iimg_indx2clusters
%
%   **inmaskvox:**
%        number of inmask voxels
%
%   **dat:**
%        binary matrix of voxels that are in (1) or out (0) of mask
%
%
% Note: we want to be more inclusive than not at this stage.
%
% ..
%    Tor Wager
%    See code for revision notes
%
%    Replacement for Worsley's fmri_mask_thresh, which was not behaving well
% ..


% last edited Oct 2011 - add support for fmri_data/image_vector objects
% added figure suppression, SG 12/14/15

% defaults
if nargin < 4, plotfigs = 1; end
if nargin < 3, implicit_masking_method = 'mean'; end
if nargin < 2, outputname = []; end


if isa(fmri_file, 'image_vector')
    
    fmri_file = replace_empty(fmri_file);
    
    datin = mean(fmri_file.dat, 2);
    
    volInfo = fmri_file.volInfo;
else
    
    % Load up to 20 images spaced evenly across the (likely) 4-D volume.
    image_list = expand_4d_filenames(fmri_file);
    nvols = size(image_list, 1);
    wh_images = round(linspace(1, nvols, min(20, nvols)));
    
    V = spm_vol(image_list(wh_images, :));
    numframes = length(V);
    
    % read the 20 images and average data
    datin = zeros(prod(V(1).dim(1:3)), numframes);
    
    for i = 1:numframes
        dat = spm_read_vols(V(i));
        
        dattmp = dat(:); dattmp = dattmp(1:1:end, :); % every 1th voxel
        datin(:, i) = dattmp;
    end
    
    datin = mean(datin, 2);
    datin(isnan(datin)) = 0;
    
    V = spm_vol(fmri_file(1,:));
    %if (size(V,1) == 1)
        fname = [V(1).fname ',1']; %passed in 4D image
    %else
    %    fname = V{1}.fname % passed in 3D images
    %end
    volInfo = iimg_read_img(fname, 2, [], 1);
    
end


% display histogram, which is the voxels of the mean values for the first
% 20 images

nbin=100;
[freq, mask]=hist(datin, nbin);

if plotfigs
    create_figure('Implicit Mask');
    plot(mask, freq, 'k', 'LineWidth', 2); hold on;
end  

% get threshold

% this could be done using a mixture model, but the mixture models i tried
% did not coverge well on actual data. this seems more robust.
%
% This will likely work well SPECIFICALLY for images with a very
% low-intensity, homogeneous background and then some signal at much higher
% values.
%


%  OLD CODE THAT WASN'T DOING WHAT IT WAS SUPPOSED TO
%[dummy, minaccel] = min(gradient(freq));
%[dummy, whfreq] = min(freq(minaccel+1:end));
%whfreq = whfreq + minaccel;  % whfreq is now the index of the minimum derivative
% following the max deceleration
%partition_point = mask(whfreq);
%hh = plot_vertical_line(partition_point); set(hh, 'LineStyle', ':');
%lowermean = mean(datin(datin < partition_point));
%uppermean = mean(datin(datin > partition_point));
%plot_vertical_line(lowermean);
%plot_vertical_line(uppermean);

% Yoni's new in-brain threshold, 12/2011
% smooth the histogram, then take the first positive gradient
% -- not properly tested yet
sm = smooth(freq, 20); %big smooth of  the histogram
grad = gradient(sm); 
i = 15; %skip the first few bins
while grad(i) < 0
    i = i+ 1;
end
dip_partition_point = mask(i);

if plotfigs
    hh = plot_vertical_line(dip_partition_point); set(hh, 'LineStyle', '-');
end

% Basic way of calculating a basic, generous implicit mask
% take 95% of voxels whose value is above the mean value
mean_pp = prctile(datin(datin > mean(datin)), 5);
if plotfigs
    hh = plot_vertical_line(mean_pp);
    set(hh, 'LineStyle', '--');
end

if (exist('implicit_masking_method', 'var') && strcmp(implicit_masking_method, 'dip'))
    mask_thresh = dip_partition_point;
else
    mask_thresh = mean_pp;
end
% need fix if no in-mask voxels
inmaskvox = sum(datin > mask_thresh);

if inmaskvox < 1000
    warning('fmri_mask_thresh_canlab:uh-oh', 'Too few in-mask voxels!  Attempting to fix...check the data...');
    
    % simply assume 20% of image is in-brain.
    mask_thresh = prctile(datin, 80);
end

inmaskvox = sum(datin > mask_thresh);


if plotfigs
    xlabel('Mask value');
    ylabel('Frequency');
    
    t = title({'Implicit Mask, lines are in-brain threhold'; 'Dashed is just above mean value (actually used), Solid is "dip-based" (currently unused).'});
    axpos = get(gca,'pos');
    extent = get(t,'extent');
    set(gca,'pos',[axpos(1) axpos(2) axpos(3) .75])
end
 % Write mask and display
dat = datin > mask_thresh;

cl = iimg_indx2clusters(dat, volInfo);


if ~isa(fmri_file, 'image_vector') && plotfigs
    
    disp('EPI should be underlay, green should be mask');
    %montage_clusters([V(1).fname ',1'], cl, {[0 1 0]});
    % Tor checked: montage_clusters reverses x for anatomical vs. clusters
    % for radiological images, 
    % because clusters are not flipped in montage display, and so are
    % displayed in native radiological, but anatomical underlay is
    % displayed in neurological because flipping specified in the image is applied.
    % this is a "feature" (does not flip for display) but confusing.
    spm_check_registration([V(1).fname ',1']);
    cluster_orthviews(cl, {[0 1 0]}, 'trans', 'add');
    cluster_orthviews_montage(10, 'axial', [], 'onerow');
    fh = findobj('Tag', 'Graphics');
    set(fh, 'Visible', 'off');
else
    % could change image_vector.montage to handle this.
    
end

if ~isempty(outputname)
    fprintf('Writing mask file: %s\n', outputname);
    iimg_reconstruct_vols(dat, volInfo, 'outname', outputname);
end


% try to use a gaussian mixture model 10 times to segment
% isok = 0; cntr = 0;
% options = statset('Display','final');
% while ~isok && cntr < 20
%     try
%         obj = gmdistribution.fit(datin,2,'Options',options);
%         isok = 1;
%         ComponentMeans = obj.mu;
%     catch
%         cntr = cntr + 1;
%
%     end
% end
%
% if cntr >= 20, disp('Tried 20 times: Failed Gaussian mixture model'); end
%
% keyboard

% from Keith Worsley, but modified so alternative way is more stable.
% THis is not working well, need to replace with above.
% nbin=100;
% dbin=10;
% %dm=fmris_read_image(fmri_file(1,:),0,0);
%
% V = spm_vol(fmri_file(1,:));
% dat = spm_read_vols(V(1));
%
% numslices=size(dat, 3);
%
% [freq, mask]=hist(dat(:),nbin);
% fm=freq(1:nbin-2*dbin);
% f0=freq(1+dbin:nbin-dbin);
% fp=freq(1+2*dbin:nbin);
% h=(abs(f0-fm)+abs(f0-fp)).*(f0>fp).*(f0>fm);
% if any(h)
%     mh=min(find(h==max(h))+dbin);
%
%     mask_thresh = max(mask(find(freq==min(freq(1:mh)))));
%
% else
%
%     mask_thresh = prctile(dat(:), 20);
%
% end



end
