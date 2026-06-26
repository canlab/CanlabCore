%% Regression Walkthrough

%%
% This walkthrough works through multiple aspects of a basic 2nd-level (group)
% analysis. We are interested in activity related to emotion regulation
% for a contrast comparing reappraisal of negative images vs. looking at
% negative images without reappraisal, i.e., a [Reappraise Negative - Look Negative]
% contrast. This walkthrough will regress brain
% contrast values in each voxel (y) onto individual differences in
% "Reappraisal Success", defined as the drop in negative emotion ratings
% for the [Reappraise Negative - Look Negative] contrast.
%
% We are interested in both the brain correlates of individual differences
% and the group-average [Reappraise Negative - Look Negative] effect in the
% brain, which is the intercept in our regression model.

%%
% This walkthrough includes multiple considerations you'll need when doing
% a practical data analysis. This goes far beyond just running a model!
% It's also important to understand the data and results. 
% Understanding the data allows you to make informed a priori choices about
% the analysis and avoid "P-Hacking". Understanding the results helps
% ensure they are valid, and helps you interpret them.
%
% A) Understanding the data
%
% # Loading the dataset into an fmri_data object
% # Examining the brain data and checking 
%
%   outliers
%   global signal
%   which voxels are analyzed (coverage) 
%   homogeneity across subjects (images)
%
% # Examining the behavioral (individual differences) data and leverages
% # Use the info gained to make data preparation and analysis decisions
%
%   Masking 
%   Transformation of the y variable
%   Rescaling/normalizing image data
%   Robust vs. OLS estimation
%
% B) Fitting the model
%
% # Running the regression and writing results images to disk
%
% C) Understanding the results
%
% # Explore variations in statistical threshold
% # Compare with other studies by comparing with a meta-analysis mask
% # Explore the impact of variations in data scaling/transformation
% # Localizing results, labeling regions, and making tables using a brain atlas
% # Extract and plot data from regions of interest
% # Explore the impact of plotting data in biased vs. unbiased regions of interest
%
% D) Multivariate prediction (overlaps with later walkthroughs)
% # Extra: Multivariate prediction within a set of unbiased ROIs

%% Part 1: Load sample data
% Load a sample dataset, including a standard brain mask
% And a meta-analysis mask that provides a priori regions of interest.
% ---------------------------------------------------------------

% Load sample data using load_image_set(), which produces an fmri_data
% object. Data loading exceeds the scope of this tutorial, but a more
% indepth demosntration may be provided by canlab_help_2_load_a_sample_dataset.m

% These are [Reappraise - Look Neg] contrast images, one image per person

[image_obj, networknames, imagenames] = load_image_set('emotionreg');

% Summarize and print a table with characteristics of the data:
desc = descriptives(image_obj);

% Load behavioral data
% This is "Reappraisal success", one score per person, in our example
% If you do not have the file on your path, you will get an error.

beh = importdata('Wager_2008_emotionreg_behavioral_data.txt')
success = beh.data(:, 2);           % Reappraisal success

% Load a mask that we would like to apply to analysis/results

mask = which('gray_matter_mask.img')

maskdat = fmri_data(mask, 'noverbose');

% Map of emotion-regulation related regions from neurosynth.org 
metaimg = which('emotion regulation_pAgF_z_FDR_0.01_8_14_2015.nii')

metamask = fmri_data(metaimg, 'noverbose');

%% Part 2: Examine the setup and coverage

% Visualize the mask
% ---------------------------------------------------------------
% BLOCK 2: Check mask

% This is an underlay brain with three separate montages so we can compare images:

o2 = canlab_results_fmridisplay([], 'noverbose', 'multirow', 3);
drawnow, snapnow;

% This is a basic gray-matter mask we will use for analysis:
% It can help reduce multiple comparisons relative to whole-image analysis
% but we should still look at what's happening in ventricles and
% out-of-brain space to check for artifacts.

o2 = addblobs(o2, region(maskdat), 'wh_montages', 1:2);

% Add a title to montage 2, the axial slice montage:
[o2, title_han] = title_montage(o2, 2, 'Mask image: Gray matter with border');
set(title_han, 'FontSize', 18)

% Visualize summary of brain coverage
% ---------------------------------------------------------------
% Check that we have valid data in most voxels for most subjects
% Always look at which voxels you are analyzing. The answer may surprise
% you!  Some software programs mask out voxels in individual analyses,
% which may then be excluded in the group analysis.

% Show summary of coverage - how many images have non-zero, non-NaN values in each voxel

o2 = montage(region(desc.coverage_obj), o2, 'trans', 'transvalue', .3, 'wh_montages', 3:4);

[o2, title_han] = title_montage(o2, 4, 'Coverage: Voxels included in images');
set(title_han, 'FontSize', 18)

% Visualize meta-analysis mask
% ---------------------------------------------------------------
o2 = montage(region(metamask), o2, 'colormap', 'wh_montages', 5:6);

[o2, title_han] = title_montage(o2, 6, 'Meta-analysis mask');
set(title_han, 'FontSize', 18)

% There are many other alternatives to the montage - e.g.,:
% orthviews(desc.coverage_obj, 'continuous');

%% Part 3: Examine the brain data 
% This is a really important step to understand what your images look like!
% Good data tend to produce good results.  

% The method plot() for an fmri_data object shows a number of plots that
% are helpful in diagnosing some problems.  

plot(image_obj);

% Check histograms of individual subjects for global shifts in contrast values

% The 'histogram' object method will allow you to examine a series of
% images in one panel each.  See help fmri_data.histogram for more options,
% including breakdown by tissue type.

hist_han = histogram(image_obj, 'byimage', 'color', 'b');

% This shows us that some of the images do not have the same mean as the
% others. This is fairly common, as individual subjects can often have
% global artifacts (e.g., task-correlated head motion or outliers) that
% influence the whole contrast image, even when baseline conditions are
% supposed to be "subtracted out".  
%
% It suggests that we may want to do an outlier analysis and/or standardize 
% the scale of the images. We'll return to this below.

% If we want to visualize this separated by tissue class, try this:

create_figure('histogram by tissue type')
hist_han = histogram(image_obj, 'byimage', 'by_tissue_type');

% Our previous peek at the data suggested that subjects 6 and 16 have
% whole-brain shifts towards acftivation and deactivation, respectively.

% This plot shows us that there is a substantial shift in gray matter values
% across the brain for some individuals.
% Also, the means across gray and white matter are correlated.  So are the standard deviations.
% We might consider a rescaling or nonparametric procedure (e.g., ranking)
% to standardize the images.  

% Some options are here:
help image_obj.rescale


%% Part 4: Examine the behavioral/task predictor(s) 
% Examine predictor distribution and leverages
% Leverage is a measure of how much each point influences the regression
% line. The more extreme the predictor value, the higher the leverage.
% Outliers will have very high leverage. High-leverage behavioral observations 
% can strongly influence, and sometimes invalidate, an analysis.

X = scale(success, 1); X(:, end+1) = 1;         % A simple design matrix, behavioral predictor + intercept
H = X * inv(X'* X) * X';                        % The "Hat Matrix", which produces fits. Diagonals are leverage

create_figure('levs', 2, 1); 
plot(success, 'o', 'MarkerFaceColor', [0 .3 .7], 'LineWidth', 3); 
set(gca, 'FontSize', 24); 
xlabel('Subject number'); 
ylabel('Reappraisal success');

subplot(2, 1, 2);
plot(diag(H), 'o', 'MarkerFaceColor', [0 .3 .7], 'LineWidth', 3); 
set(gca, 'FontSize', 24); 
xlabel('Subject number'); 
ylabel('Leverage');

drawnow, snapnow

% The distribution suggests that there are some high-leverage values to
% watch out for. They will influence the results disproportionately at every voxel in the brain!
% Ranking predictors and/or robust regression are good
% mitigation/prevention procedures in many cases.

%% Part 5: Make data prep and analysis decisions
% Before we've "peeked" at the results, let's make some sensible decisions
% about what to do with the analysis. We'll compare the results to the
% "standard analysis" if we hadn't examined the data at all at end.

% Choices:

% 1. Masking: Apply a gray-matter mask. This will limit the area subjected to multiple-comparisons 
%    correction to a meaningful set.  I considered {gray matter, no mask}

image_obj = apply_mask(image_obj, maskdat);         % Mask - analyze gray-matter

% 2. Image scaling: I considered arbitrary rescaling of images, e.g., by
% the l2norm or z-scoring images. Z-scoring can cause artifactual
% de-activations when we are looking at group effects. L2-norming is
% sensible if there are bad scaling problems... but I'd rather minimize ad
% hoc procedures applied to the data. Ranking each voxel is also sensible, and 
% if the predictor is ranked as well, this would implement Spearman's correlations.
% I considered {none, l2norm images, zscoreimages, rankvoxels}, and decided
% on ranking voxels.  For the group average (intercept), this is not  -- I'll use robust regression to help mitigate issues.

image_obj = rescale(image_obj, 'rankvoxels');     % rescale images within this mask (relative activation)

% 3. Predictor values: I considered dropping subject 16, which is a global
% outlier and also has high leverage, or windsorizing values to 3 sd. 
% I also considered ranking predictor values, a very sensible choice. 
% So I considered {do nothing, drop outlier, windsorize outliers, rank values}
% But robust regression (below) is partially redundant with this. I ultimately decided
% on ranking.

% This is done below:
% image_obj.X = scale(rankdata(image_obj.X), 1);

% 4. Regression method: I considered {standard OLS and robust IRLS
% regression}.  I chose robust regression because it helps mitigate the
% effects of outliers.  However, this is partly redundant with ranking, so
% it would have been sensible to choose ranking OR robust regression.

%% Part 6: Run the regression analysis and get results
% The regress method takes predictors that are attached in the object's X
% attribute (X stands for design matrix) and regresses each voxel's
% activity (y) on the set of regressors.  
% This is a group analysis, in this case correlating brain activity with
% reappraisal success at each voxel.

% .X must have the same number of observations, n, in an n x k matrix.
% n images is the number of COLUMNS in image_obj.dat

% mean-center success scores and attach them to image_obj in image_obj.X
image_obj.X = scale(success, 1);

image_obj.X = scale(rankdata(image_obj.X), 1);

% runs the regression at each voxel and returns statistic info and creates
% a visual image.  regress = multiple regression.
% % Track the time it takes (about 45 sec for robust regression on Tor's 2018 laptop)
tic

% out = regress(image_obj);             % this is what we'd do for standard OLS regression
out = regress(image_obj, 'robust');

toc

% out has statistic_image objects that have information about the betas
% (slopes) b, t-values and p-values (t), degrees of freedom (df), and sigma
% (error variance).  The critical one is out.t.
% out = 
% 
%         b: [1x1 statistic_image]
%         t: [1x1 statistic_image]
%        df: [1x1 fmri_data]
%     sigma: [1x1 fmri_data]

% This is a multiple regression, and there are two output t images, one for
% each regressor.  We've only entered one regressor, why two images?  The program always
% adds an intercept by default.  The intercept is always the last column of the design matrix

% Image   1  <--- brain contrast values correlated with "success"
% Positive effect:   0 voxels, min p-value: 0.00001192
% Negative effect:   0 voxels, min p-value: 0.00170529
% Image   2 FDR q < 0.050 threshold is 0.003193
% 
% Image   2 <--- intercept, because we have mean-centered, this is the
% average group effect (when success = average success).  "Reapp - Look"
% contrast in the whole group.
% Positive effect: 3133 voxels, min p-value: 0.00000000
% Negative effect:  51 voxels, min p-value: 0.00024068

% Now let's try thresholding the image at q < .05 fdr-corrected.
t = threshold(out.t, .05, 'fdr');
 
% Select the predictor image we care about: (the 2nd/last image is the intercept)
t = select_one_image(t, 1);  

% ...and display on a series of slices and surfaces
% There are many options for displaying results.
% montage and orthviews methods for statistic_image and region objects provide a good starting point.

o2 = montage(t, 'trans', 'full');
o2 = title_montage(o2, 2, 'Behavioral predictor: Reappraisal success  (q < .05 FDR)');
snapnow

% Or:
% orthviews(t)


%% Write the thesholded t-image to disk

% % t.fullpath = fullfile(pwd, 'Reapp_Success_005_FDR05.nii');
% % write(t)
% % 
% % disp(['Written to disk: ' t.fullpath])

%% Part 7. Explore: Display regions at a lower threshold, alongside meta-analysis regions
% Display the results on slices another way:
% multi_threshold lets us see the blobs with significant voxels at the
% highest (most stringent) threshold, and voxels that are touching
% (contiguous) down to the lowest threshold, in different colors.
% This shows results at p < .001 one-tailed (p < .002 two-tailed), with 10 contiguous voxels.
% Contiguous blobs are more liberal thresholds, down to .02 two-tailed =
% .01 one tailed, are shown in different colors. 

% Set up a display object with two separate brain montages:
o2 = canlab_results_fmridisplay([], 'multirow', 2);

% Use multi-threshold to threshold and display the t-map:
% Pass the fmridisplay object with montages (o2) in as an argument to use
% that for display.  Pass out the object with updated info.
% 'wh_montages', 1:2 tells it to plot only on the first two montages (one sagittal and one axial slice series) 

o2 = multi_threshold(t, o2, 'thresh', [.002 .01 .02], 'sizethresh', [10 1 1], 'wh_montages', 1:2);
o2 = title_montage(o2, 2, 'Behavioral predictor: Reappraisal success (p < .001 one-tailed and showing extent');

% Get regions from map from neurosynth.org
% Defined above: metaimg = which('emotion regulation_pAgF_z_FDR_0.01_8_14_2015.nii')

r = region(metaimg);

% Add these to the bottom montages
o2 = addblobs(o2, r, 'maxcolor', [1 0 0], 'mincolor', [.7 .2 .5], 'wh_montages', 3:4);

o2 = title_montage(o2, 4, 'Neurosynth mask: Emotion regulation');

% Here, the regions predictive of success are largely different from those
% that respond to reappraisal demand in past literature.

%% Part 8: Compare to a standard analysis with no scaling or variable transformation

% Start over and re-load the images:
[image_obj_untouched] = load_image_set('emotionreg', 'noverbose');

% attach success scores to image_obj in image_obj.X
% Just mean-center them so we can interpret the group result
image_obj_untouched.X = scale(success, 1);

% Do a standard OLS regression and view the results:
out = regress(image_obj_untouched, 'nodisplay');        

% Treshold at q < .05 FDR and display
t_untouched = threshold(out.t, .05, 'fdr');

o2 = montage(t_untouched);
o2 = title_montage(o2, 2, 'No scaling: Behavioral predictor: Reappraisal success (q < .05 FDR)');
o2 = title_montage(o2, 4, 'No scaling: Intercept: Group average activation');

% Since we have no results for the predictor, re-threshold at p < .005
% uncorrected and display so we can see what's there:

t_untouched = threshold(t_untouched, .005, 'unc');

o2 = montage(t_untouched);
o2 = title_montage(o2, 2, 'No scaling: Behavioral predictor: Reappraisal success (p < .005 uncorrected)');
o2 = title_montage(o2, 4, 'No scaling: Intercept: Group average activation');
snapnow

% select just the image for reappraisal success:
t_untouched = select_one_image(t_untouched, 1);

%% Make plots comparing the analysis with our choices and a naive analysis directly

o2 = canlab_results_fmridisplay([], 'multirow', 4);

t = threshold(t, .05, 'fdr');
o2 = montage(region(t), o2, 'wh_montages', 1:2, 'colormap');
o2 = title_montage(o2, 2, 'Reappraisal success with informed choices (q < .05 corrected)');

t_untouched = threshold(t_untouched, .05, 'fdr');
o2 = montage(region(t_untouched), o2, 'wh_montages', 3:4, 'colormap');
o2 = title_montage(o2, 4, 'Reappraisal success with no scaling (q < .05 corrected)');

% Here we'll use "multi-threshold" to display blobs at lower thresholds
% contiguous with those at higher thresholds. We'll use uncorrected p <
% .005 so we can see results in both maps:

o2 = multi_threshold(t, o2, 'thresh', [.002 .01 .02], 'k', [10 1 1], 'wh_montages', 5:6);
o2 = title_montage(o2, 6, 'Reappraisal success with informed choices (p < .001 one-tailed unc, showing extent to .01 one-tailed)');

o2 = multi_threshold(t_untouched, o2, 'thresh', [.002 .01 .02], 'k', [10 1 1], 'wh_montages', 7:8);
o2 = title_montage(o2, 8, 'Reappraisal success with no scaling (p < .001 unc, showing extent to .01 one-tailed)');

% One-tailed results are the default in SPM and FSL - so this is what you'd
% get with 0.001 in SPM.



%% Part 9: Localize results using a brain atlas
% The table() command is a simple way to label regions
% In addition, we can visualize each region in table to check the accuracy
% of the labels and the extent of the region.

% Select the Success regressor map
r = region(t);

% Autolabel regions and print a table
r = table(r);

% Make a montage showing each significant region
montage(r, 'colormap', 'regioncenters');

%% Part 10: Extract and plot data from regions of interest
% Here, we'll explore data extracted from two kinds of regions:

% 1 - biased regions based on the significant results. This is good for
% visualizing the distribution of the data in a region and checking
% assumptions, but examining the correlation in these regions is
% "circular", as the correlation was used to select voxels in the first
% place. The observed correlations in these regions will be upwardly
% biased as a result. In general, you cannot estimate the effect size
% (correlation or other metric) in a region or voxel after selecting it
% using criteria that are non-independent of the effect of interest.

% 2 - unbiased regions from the meta-analysis.  These are fair.



%% Extract and plot from (biased) regions of interest based on results
% Here, we'll compare the implications of picking BIASED ROIs (which is
% often done in the literature) vs. picking UNBIASED ROIs.
% 
% First, let's visualize the correlation scatterplots in the areas we've
% discovered as related to Success

% Extract data from all regions
% r(i).dat has the averages for each subject across voxels for region i
r = extract_data(r, image_obj);

% Select only regions with 3+ voxels
% wh = cat(1, r.numVox) < 3;
% r(wh) = [];

% Set up the scatterplots
nrows = floor(sqrt(length(r)));
ncols = ceil(length(r) ./ nrows);
create_figure('scatterplot_region', nrows, ncols);

% Make a loop and plot each region
for i = 1:length(r)
    
    subplot(nrows, ncols, i);

    % Use this line for non-robust correlations:
    %plot_correlation_samefig(r(i).dat, datno16.X);
    
    % Use this line for robust correlations:
    plot_correlation_samefig(r(i).dat, image_obj.X, [], 'k', 0, 1);
  
    set(gca, 'FontSize', 12);
    xlabel('Reappraise - Look Neg brain response');
    ylabel('Reappraisal success');
    
    % input('Press a key to continue');
    
end


%% Extract and plot data from unbiased regions of interest
% Let's visualize the correlation scatterplots in some meta-analysis
% derived ROIs

% Select the Neurosynth meta-analysis map
r = region(metaimg);

% Extract data from all regions
r = extract_data(r, image_obj);

% Select only regions with 20+ voxels
wh = cat(1, r.numVox) < 20;
r(wh) = [];

r = table(r);

% Make a montage showing each significant region
montage(r, 'colormap', 'regioncenters');

% Set up the scatterplots
nrows = floor(sqrt(length(r)));
ncols = ceil(length(r) ./ nrows);
create_figure('scatterplot_region', nrows, ncols);

% Make a loop and plot each region
for i = 1:length(r)
    
    subplot(nrows, ncols, i);

    % Use this line for non-robust correlations:
    %plot_correlation_samefig(r(i).dat, datno16.X);
    
    % Use this line for robust correlations:
    plot_correlation_samefig(r(i).dat, image_obj.X, [], 'k', 0, 1);
  
    set(gca, 'FontSize', 12);
    xlabel('Brain');
    ylabel('Success');
    title(' ');
   
end

%% Part 11: (bonus) Multivariate prediction from unbiased ROI averages
% Predict reappraisal success using brain images
% with 5-fold cross-validation

% Since we are predicting outcome data from brain images, we have to attach 
% data in the .Y field of the fmri_data object.  
% .Y is the outcome to be explained, in this case success scores
image_obj.Y = image_obj.X;  

% 5-fold cross validated prediction, stratified on outcome

[cverr, stats, optout] = predict(image_obj, 'algorithm_name', 'cv_lassopcrmatlab', 'nfolds', 5);

% Plot cross-validated predicted outcomes vs. actual

% Critical fields in stats output structure:
% stats.Y = actual outcomes
% stats.yfit = cross-val predicted outcomes
% pred_outcome_r: Correlation between .yfit and .Y
% weight_obj: Weight map used for prediction (across all subjects).
                 
% Continuous outcomes:

create_figure('scatterplot');
plot(stats.yfit, stats.Y, 'o')
axis tight
refline
xlabel('Predicted reappraisal success');
ylabel('Observed reappraisal success');

% Though many areas show some significant effects, these are not strong
% enough to obtain a meaningful out-of-sample prediction of Success

