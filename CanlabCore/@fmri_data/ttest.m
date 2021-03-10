function statsimg = ttest(fmridat, pvalthreshold, thresh_type)
% T-test on fmri_data class object
%
% :Usage:
% ::
%
%    statsimg = ttest(fmridat, pvalthreshold, thresh_type)
%
% :Inputs:
%
%   **p-value threshold:**
%        p-value, e.g., .05 or .001 or [.001 .01 .05]
%
%   **thresh_type:**
%        'uncorrected', 'fwe', or 'fdr'
%
% :Examples:
% ::
%
% % Group one-sample t-test: 
% % A simple, complete example of a group analysis
% % For more, see walkthroughs on canlab.github.io 
% % --------------------------------------------------------------------
% % Load sample images, creating and fmri_data class object with 30 images     
%     imgs = load_image_set('emotionreg');
%
% % Display a slice from each image in a montage:
%     slices(imgs);
%
% % Display some useful summary plots of the dataset:
%     plot(imgs);
%
% % Perform a t-test on each voxel, returning a statistic_image object
% % containing t-stats and p-values:
%     t = ttest(imgs);
%
% % Display the unthresholded results in a quick-render montage of image
% % values only:
%     display_slices(t, 'axial'); colormap summer; colorbar;
% 
% % Display the unthresholded results over an anatomical underlay, 
% % on a combination of slices and surfaces, 
% % returning an fmridisplay class object with registered handles
%     o2 = canlab_results_fmridisplay(t, 'full');
%
% % Remove colors from slices and surfaces registered in the o2 object:
%     o2 = removeblobs(o2);
%
% % Threshold the t-statistic_image object at p < 0.005
%     t = threshold(t, .005, 'unc');
%
% % Re-display the thresholded images on slices/surfaces registered in o2:
%     o2 = addblobs(o2, region(t), 'nolegend');
%
% % Display the thresholded t-map with orthviews:
%     orthviews(t);
%
% % display on a slice montage:
%    create_figure('slices'); axis off; 
%    montage(t);
%
% % Re-threshold at q < 0.05 FDR, and re-display on orthviews:
%     t = threshold(t, .05, 'fdr');
%     orthviews(t);
%
% % Print a table of results, and return a region-class object r with labels
% % from a default atlas (an atlas-class object):
%     r = table(t);
%
% % Display a slice montage showing each activation blob, with labels:
%   montage(r, 'regioncenters', 'colormap');
%
% :Note: for two-sample T-test, use fmri_data.regress

fprintf('One-sample t-test\n')
fprintf('Calculating t-statistics and p-values\n');

err = ste(fmridat.dat')';

t = nanmean(fmridat.dat')' ./ err;

N = single(sum(~(isnan(fmridat.dat') | fmridat.dat' == 0) , 1));
N = N'; % make column

df = size(fmridat.dat, 2) - 1;

% two-tailed - now done in statistic_image constructor
% p = 2 * (1 - tcdf(abs(t), df));

% t = t';
% p = p';

% removed: done in constructor: 'p', p, 'p_type', 'two-tailed',
statsimg = statistic_image('type', 't', 'dat', t, 'ste', err, 'N', N, 'volInfo', fmridat.mask.volInfo, 'dfe', df);

% get rid of NaNs for empty vals
statsimg.p(isnan(t)) = 1;
statsimg.dat(isnan(t)) = 0;

statsimg.removed_voxels = fmridat.removed_voxels;

if nargin > 2
    
    statsimg = threshold(statsimg, pvalthreshold, thresh_type);
    
end



end
