%% Perform a voxel-wise t-test
%
% In this example we perform a second level analysis on first level
% statistical parametric maps. Specifically we use a t-test to obtain an
% ordinary least squares estimate for the group level parameter. 
%
% The example uses the emotion regulation data provided with
% CANlab_Core_Tools. This dataset consists of a set of contrast images for 
% 30 subjects from a first level analysis. The contrast is [reappraise neg vs. look neg],
% reappraisal of negative images vs. looking at matched negative images. 
%
% These data were published in:
% Wager, T. D., Davidson, M. L., Hughes, B. L., Lindquist, M. A., 
% Ochsner, K. N.. (2008). Prefrontal-subcortical pathways mediating 
% successful emotion regulation. Neuron, 59, 1037-50.
%
% By the end of this example we will have regenerated the results of figure
% 2A of this paper.

%% Executive summary of the whole analysis
%
% As a summary, here is a complete set of commands to load the data and run
% the entire analysis:
% 
%     img_obj = load_image_set('emotionreg');         % Load a dataset
%     t = ttest(img_obj);                             % Do a group t-test
%     t = threshold(t, .05, 'fdr', 'k', 10);          % Threshold with FDR q < .05 and extent threshold of 10 contiguous voxels
%  
%     % Show regions and print a table with labeled regions:
%     montage(t);  drawnow, snapnow;                  % Show results on a slice display
%     r = region(t);                                  % Turn t-map into a region object with one element per contig region
%     table(r);                                       % Print a table of results using new region names
%  
% Here a graphical description of what it does:
%
% <<CANlab_ttest_flowchart.png>>
%
% Now, let's walk through it step by step, with a few minor additions.

%% Load sample data

% We can load this data using load_image_set(), which produces an fmri_data
% object. Data loading exceeds the scope of this tutorial, but a more
% indepth demonstration may be provided by load_a_sample_dataset.m

[image_obj, networknames, imagenames] = load_image_set('emotionreg');

%% Do a t-test

% Voxel-wise tests estimate parameters for each voxel independently. The
% fmri_data/ttest function performs this just like the native matlab
% ttest() function, and similarly returns a p-value and t-stat. Unlike the
% matlab ttest() function, fmri_data/ttest performs this for every voxel in
% the fmri_data object. All voxels are evaluated independently, but all are
% evaluated nonetheless. A statistic_image is returned.

t = ttest(image_obj);

%% Visualize the results
% There are many options. See methods(statistic_image) and methods(region)

orthviews(t)
drawnow, snapnow; 

% As can be seen, the ttest() result is unthresholded. The threshold
% function can be used to apply a desired alpha level using any of a number
% of methods. Here FDR is used to control for alpha=0.05. Note that no
% information is erased when performing thresholding on a statistic_image.

t = threshold(t, .05, 'fdr');
orthviews(t)
drawnow, snapnow; 

% Many neuroimaging packages (e.g., SPM and FSL) do one-tailed tests 
% (with one-tailed p-values) and only show you positive effects 
% (i.e., relative activations, not relative deactivations).  
% All the CANlab tools do two-sided tests, report two-tailed p-values. 
% By default, hot colors (orange/yellow) will show activations, and cool
% colors (blues) show deactivations.

% We can also apply an "extent threshold" of > 10 contiguous voxels:

t = threshold(t, .05, 'fdr', 'k', 10);
orthviews(t)
drawnow, snapnow; 

% montage is another visualization method. This function may require a
% relatively large amount of memory, depending on the resolution of the 
% anatomical underlay image you use. We recommend around 8GB of free memory. 
%
create_figure('montage'); axis off
montage(t)
drawnow, snapnow; 

%% Print a table of results

% First, we'll have to convert to another object type, a "region object".
% This object groups voxels together into "blobs" (often of contiguous
% voxels). It does many things that other object types do, and inter-operates
% with them closely.  See methods(region) for more info.

% Create a region object called "r", with contiguous blobs from the
% thresholded t-map:

r = region(t);

% Print the table:

table(r);

%% More on getting help

% Help info is available on https://canlab.github.io/
% and in the CANlab help examples repository: https://github.com/canlab/CANlab_help_examples
% This contains a series of walkthroughs, including this one.

% You can get help for the main object types by typing this in Matlab:
% doc fmri_data (or doc any other object types)

% Also, more detailed help is available for each function using the "help" command.
% You can get help and options for any object method, like "table". But
% because the method names are simple and often overlap with other Matlab functions 
% and toolboxes (this is OK for objects!), you will often want to specify
% the object type as well, as follows:

% help region.table

%% Write the t-map to disk

% Now we need to save our results. You can save the objects in your
% workspace or you can write your resulting thresholded map to an analyze
% file. The latter may be useful for generating surface projections using 
% Caret or FreeSurfer for instance.
%
% Thresholding did not actually eliminate nonsignificant voxels from our 
% statistic_image object (t). If we  simply write out that object, we will 
% get t-statistics for all voxels. 

t.fullpath = fullfile(pwd, 'example_t_image.nii');
write(t, 'overwrite')

% If we use the 'thresh' option, we'll write thresholded values:
% (write() does not overwrite an existing file unless asked, and we just
% wrote this path above, so pass 'overwrite' here too.)
write(t, 'thresh', 'overwrite')

t_reloaded = statistic_image(t.fullpath, 'type', 'generic');
orthviews(t_reloaded)


%% Some useful optional extensions

figure; surface(t);
figure; surface(t, 'foursurfaces');
figure; surface(t, 'coronal_slabs_5');

table_of_atlas_regions_covered(t);


%% Explore on your own
%
% 1. Click around the thresholded statistic image. What lobes of the brain
% are most of the results in? What can we say about areas that are not
% active, if anything?
%
% 2. Why might some of the results appear to be outside the brain? What
% does this mean for the validity of the analysis? Should we consider these
% areas in our interpretation, or should they be "masked out" (excluded)?
%
% That's it for this section!!

%% Explore More: CANlab Toolboxes
% Tutorials, overview, and help: <https://canlab.github.io>
%
% Toolboxes and image repositories on github: <https://github.com/canlab>
%
% <html>
% <table border=1><tr>
% <td>CANlab Core Tools</td>
% <td><a href="https://github.com/canlab/CanlabCore">https://github.com/canlab/CanlabCore</a></td></tr>
% <td>CANlab Neuroimaging_Pattern_Masks repository</td>
% <td><a href="https://github.com/canlab/Neuroimaging_Pattern_Masks">https://github.com/canlab/Neuroimaging_Pattern_Masks</a></td></tr>
% <td>CANlab_help_examples</td>
% <td><a href="https://github.com/canlab/CANlab_help_examples">https://github.com/canlab/CANlab_help_examples</a></td></tr>
% <td>M3 Multilevel mediation toolbox</td>
% <td><a href="https://github.com/canlab/MediationToolbox">https://github.com/canlab/MediationToolbox</a></td></tr>
% <td>M3 CANlab robust regression toolbox</td>
% <td><a href="https://github.com/canlab/RobustToolbox">https://github.com/canlab/RobustToolbox</a></td></tr>
% <td>M3 MKDA coordinate-based meta-analysis toolbox</td>
% <td><a href="https://github.com/canlab/Canlab_MKDA_MetaAnalysis">https://github.com/canlab/Canlab_MKDA_MetaAnalysis</a></td></tr>
% </table>
% </html>
% 
% Here are some other useful CANlab-associated resources:
%
% <html>
% <table border=1><tr>
% <td>Paradigms_Public - CANlab experimental paradigms</td>
% <td><a href="https://github.com/canlab/Paradigms_Public">https://github.com/canlab/Paradigms_Public</a></td></tr>
% <td>FMRI_simulations - brain movies, effect size/power</td>
% <td><a href="https://github.com/canlab/FMRI_simulations">https://github.com/canlab/FMRI_simulations</a></td></tr>
% <td>CANlab_data_public - Published datasets</td>
% <td><a href="https://github.com/canlab/CANlab_data_public">https://github.com/canlab/CANlab_data_public</a></td></tr>
% <td>M3 Neurosynth: Tal Yarkoni</td>
% <td><a href="https://github.com/neurosynth/neurosynth">https://github.com/neurosynth/neurosynth</a></td></tr>
% <td>M3 DCC - Martin Lindquist's dynamic correlation tbx</td>
% <td><a href="https://github.com/canlab/Lindquist_Dynamic_Correlation">https://github.com/canlab/Lindquist_Dynamic_Correlation</a></td></tr>
% <td>M3 CanlabScripts - in-lab Matlab/python/bash</td>
% <td><a href="https://github.com/canlab/CanlabScripts">https://github.com/canlab/CanlabScripts</a></td></tr>
% </table>
% </html>
%
% *Object-oriented, interactive approach*
% The core basis for interacting with CANlab tools is through object-oriented framework.
% A simple set of neuroimaging data-specific objects (or _classes_) allows you to perform
% *interactive analysis* using simple commands (called _methods_) that
% operate on objects. 
%
% Map of core object classes:
%
% <<CANlab_object_types_flowchart.png>>
%%
% % Bogdan Petre on 9/14/2016, updated by Tor Wager July 2018, Jan 2020
%