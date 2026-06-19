%% Basic Image Visualization

%%
% Using Basic Plot Methods (plot, orthviews, and montage) to Examine a Dataset
%
% Objects can have methods with intuitive names, some of which overlap
% with names of functions in Matlab or other toolboxes. fmri_data.plot()
% is one of these. When you call plot() and pass in an fmri_data object, 
% you invoke the fmri_data object method, and a special plot for fmri_data
% objects is produced. 
%
% You can list methods for an object class (e.g., fmri_data) by typing:
%
% |methods(fmri_data)|
%
% You can get help for a method by typing
% |help <object class name>.<method name>| 
% e.g., |help fmri_data.plot|
%
% The plot() method takes an fMRI data object as a parameter and produces 
% an SPM Orthview presentation along with 6 plots of the data.
%
% The plot method expects an fMRI data object to be passed in. We can
% create an fMRI data object using the emotion regulation dataset
% via the following code:

[image_obj, networknames, imagenames] = load_image_set('emotionreg', 'noverbose');

%%
% Once created, we can pass this data object to the plot function to get
% the entire set of outputs, including Matlab console output regarding
% outliers and corresponding data visualizations, using the simple command:

plot(image_obj);

%% Explaining the output
% 
% The help file for fmri_data.plot did object has more information about the plots:
%
% e.g., |help fmri_data.plot|
%
% e.g., |help image_obj.plot|

help fmri_data.plot

%% Use Orthviews to visualize the mean image
%
% Some other commonly used methods to display images are the orthviews()
% and montage() methods. fmri_data objects have these methods, and other
% object classes do to, like the statistic_image, atlas, and region
% classes.
%
% First, let's create a mean image for the dataset:

m = mean(image_obj);
orthviews(m);
drawnow, snapnow

%% Threshold and display
% let's threshold it using the threshold() method.
% Here we'll exclude values between -1 and 1 and view the extreme values:

m = threshold(m, [-1 1], 'raw-outside');
orthviews(m);
drawnow, snapnow

%% Create and display a region object
% We can create a region class object, another type of object, from the
% thresholded image. This has additional info and options about each
% contiguous 'blob' in the suprathreshold map:

r = region(m);
orthviews(r);
drawnow, snapnow

%% Orthviews options
% Orthviews methods have a range of options.
% They use the cluster_orthviews function, which uses spm_orthviews
% See |help cluster_orthviews| for options.
%
% Let's try one that visualizes each contiguous blob in a different, solid
% color:

orthviews(r, 'unique', 'solid');
drawnow, snapnow

%% Use montage to visualize the thresholded mean image
% 
% Sometimes, we want to view map that shows a canonical range of slices.
% This is really useful for producing standard output for papers
% Arguably, one should *always* view and publish montage maps showing all
% slices, so as to show the "whole picture" and not omit any results.
%
% You can customize this a lot, as it uses the fmridisplay() object class,
% which allows you to add custom montages (in axial, saggital, and coronal 
% orientations, add blobs of various types, and remove them and re-plot.
% See |help fmridisplay| and |help fmridisplay.addblobs| for more details.
%
% For now, we'll just stick to a basic plot. We'll first create an empty
% figure,then plot the montage on it.

create_figure('montage'); axis off; 
montage(m);
drawnow, snapnow

% We've already thresholded it, so it'll use the previous threshold.
% however, we can re-threshold the image and redisplay it as well.

%% Use montage to visualize each blob in a thresholded map
%
% A really useful thing to do is to take a region object, often from a 
% thresholded map, and visualize each region.
% the montage() methods also have a number of options.
% Let's try one, |'regioncenters'|, that plots each blob (region) on a
% separate slice.
%
% Furthermore, we can use the |'colormap'| option to view the regions with
% colors mapped to their associated values, e.g., hot colors for positive
% values and cool colors for negative values.
%
% Finally, we might want to assign names to each region based on an atlas.
% We'll do that before plotting, so that the names appear on the plots.
% These names are saved in the r(i).shorttitle field, for each region i.
% The region.table() method automatically labels them as well.
% You can customize the atlas used; the default is the 'canlab2018_2mm'
% atlas (see |help load_atlas| for more info.)

r = autolabel_regions_using_atlas(r);

montage(r, 'regioncenters', 'colormap');
drawnow, snapnow

% Note that some large regions may span multiple areas.
% This can happen if the various regions are connected by contiguous 
% suprathreshold voxels.

%% Explore on your own
%
% 1. Try to re-threshold the image using some values you choose and re-plot.
% Look at the help for more options on thresholding. and pick one. What do you see?
% Try to plot only voxels with positive values in the mean image.
%
% 2. Try to bring up only *one* region from the region object (r) in orthviews.
% Can you visualize it in all three views?
%
% 3. Try using a couple of other display options in the montage() and orthviews() 
% methods. What do you see?

% That's it for this section!!