%% Using the atlas object for Region of Interest (ROI) analyses

%% About the atlas object
% An atlas-class object is a specialized subclass of fmri_data that stores
% information about a series of parcels, or pre-defined regions, and in
% some cases the probalistic maps that underlie the final parcellation.
% 
% Some common uses of atlas objects include:
% * Labeling regions by best-matching atlas regions, as in region.table()
% * Analysis within specific ROIs, or on ROI averages
% * Defining regions for connectivity and graph theoretic analyses
%
% A list of pre-defined atlases is contained in the function |load_atlas|.
% The default atlas for some CANlab functions is the "CANlab combined 2018"
% This was created by Tor Wager from other published atlases. It includes:
%
% * Glasser 2016 Nature 180-region cortical parcellation (in MNI space, not indivdidualized)
% * Pauli 2016 PNAS basal ganglia subregions
% * Amygdala/hippocampal and basal forebrain regions from SPM Anatomy Toolbox
% * Thalamus regions from the Morel thalamus atlas
% * Subthalamus/Basal forebrain regions from Pauli "reinforcement learning" HCP atlas
% * Diedrichsen cerebellar atlas regions
% * Multiple named brainstem nuclei localized based on individual papers
% * Additional Shen atlas parcels to cover areas (esp. brainstem) not otherwise named
%
% References for the corresponding papers are stored in the atlas object, and 
% printed in tables generated with the region.table() method.
%
% There are many methods for atlas, including all the fmri_data methods.
% You can see those with |methods(atlas)|. For example:
%
% Using the |select_atlas_subset()| method:
% You can select any subset of atlas regions by name or number, and return
% them in a new atlas object. You can also 'flatten' regions, combining
% them into a single new region. 
%
% % Using the |extract_data()| method:
% You can extract average data from every image in an atlas, for a series of target images.
%
% We will explore these here.

%% load an atlas
% 'atlas' objects are a class of objects specially designed for brain
% atlases. Here is more information on this class (also try |doc atlas|)

help atlas

% The function load_atlas in the CANlab toolbox loads a number of named
% atlases included with the toolbox.  Here is a list of named atlases:

help load_atlas

% Now load the "CANlab combined 2018" atlas:
atlas_obj = load_atlas('canlab2018_2mm');


%% visualize the atlas regions

orthviews(atlas_obj);

o2 = montage(atlas_obj);

drawnow, snapnow
%% Selecting regions of interest
% There are several ways of selecting regions of interest. You can do it:
% * By name or part of name (in the .labels, .labels_2, or other fields)
% * By number
% * Graphically, or by entering a coordinate and a search radius

%% 
% *Select regions graphically or by coordinate*

% First bring up the orthviews display:
orthviews(atlas_obj);

% Now click on a coordinate in the center of vmPFC
% We'll use the SPM code below to do this here:

spm_orthviews('Reposition', [0 38 -11])

% Running this will create a subset atlas containing only regions with centers 
% within 20 mm of that coordinate.

[obj_within_20mm,index_vals] = select_regions_near_crosshairs(atlas_obj, 'thresh', 20);

orthviews(obj_within_20mm);
drawnow, snapnow

% To reproduce this later in a script, you can enter the coordinates as
% input instead:

[obj_within_20mm,index_vals] = select_regions_near_crosshairs(atlas_obj, 'coords', [0 38 -11], 'thresh', 20);

% This function returns index values, so you can use these directly to
% specify the regions in select_atlas_subset as well.

find(index_vals)

test_obj = select_atlas_subset(atlas_obj, find(index_vals));

% % This should yield the same map of vmPFC regions: orthviews(test_obj)

%% 
% *Select regions by name*
% Let's select all the regions in the thalamus. All regions are labeled in
% the atlas object, so we can select them by name.

% Select all regions with "Thal" in the label:
thal = select_atlas_subset(atlas_obj, {'Thal'})

% Print the labels:
thal.labels

% Select a few thalamus/epithalamus regions of interest:
thal = select_atlas_subset(atlas_obj, {'Thal_Intra', 'Thal_VL', 'Thal_MD', 'Thal_LGN', 'Thal_MGN', 'Thal_Hb'});
thal.labels

% Select all the regions with "Thal" in the label, and collapse them into a single region:
whole_thal = select_atlas_subset(atlas_obj, {'Thal'}, 'flatten');

%% Extract and analyze data from atlas regions
% 

%%
% First, we'll load a dataset to extract data from for an ROI analysis
% The dataset contains data from 33 participants, with brain responses to six levels
% of heat (non-painful and painful).  
% 
% Aspects of this data appear in these papers:
% Wager, T.D., Atlas, L.T., Lindquist, M.A., Roy, M., Choong-Wan, W., Kross, E. (2013). 
% An fMRI-Based Neurologic Signature of Physical Pain. The New England Journal of Medicine. 368:1388-1397
% (Study 2)
%
% Woo, C. -W., Roy, M., Buhle, J. T. & Wager, T. D. (2015). Distinct brain systems 
% mediate the effects of nociceptive input and self-regulation on pain. PLOS Biology. 13(1): 
% e1002036. doi:10.1371/journal.pbio.1002036
%
% Lindquist, Martin A., Anjali Krishnan, Marina López-Solà, Marieke Jepma, Choong-Wan Woo, 
% Leonie Koban, Mathieu Roy, et al. 2015. ?Group-Regularized Individual Prediction: 
% Theory and Application to Pain.? NeuroImage. 
% http://www.sciencedirect.com/science/article/pii/S1053811915009982.
%
% This dataset is shared on figshare.com, under this link:
% https://figshare.com/s/ca23e5974a310c44ca93
%
% Here is a direct link to the dataset file with the fmri_data object:
% https://ndownloader.figshare.com/files/12708989
%
% The key variable is image_obj
% This is an fmri_data object from the CANlab Core Tools repository for neuroimaging data analysis.
% See https://canlab.github.io/
%
% image_obj.dat contains brain data for each image (average across trials)
% image_obj.Y contains pain ratings (one average rating per image)
%
% image_obj.additional_info.subject_id contains integers coding for which
% Load the data file, downloading from figshare if needed
%
% Alternative sample datasets:
% --------------------------------------------
% This dataset will take time to download from figshare if you haven't yet
% downloaded it. An alternative is to use the dataset included with the
% CANlab Core toolbox (though the results are less interesting for this
% example).
% image_obj = load_image_set('emotionreg');

fmri_data_file = which('bmrk3_6levels_pain_dataset.mat');

if isempty(fmri_data_file)
    
    % attempt to download
    disp('Did not find data locally...downloading data file from figshare.com')
    
    fmri_data_file = websave('bmrk3_6levels_pain_dataset.mat', 'https://ndownloader.figshare.com/files/12708989');
    
end

load(fmri_data_file);

descriptives(image_obj);


%% 
% *Now, we'll extract data from each atlas region*
% - "r" is a region-class object. see >> help region
% - r(i).dat contains averages over voxels within region i. for n images, it contains an n x 1 vector with average data.
% - r(i).data contains an images x voxels matrix of all data within region i.

r = extract_roi_averages(image_obj, thal);

%% 
% You could also extract data from the whole thalamus into a single
% variable using:
% |r = extract_roi_averages(data_obj, whole_thal);|

%% 
% An alternative is the extract_data method:
% |region_avgs = extract_data(thal, image_obj);|
% The values may differ slightly based on which image is being resampled to
% which. Resampling matches the voxels in two sets of images that may have different 
% voxel sizes and voxel boundaries by interpolating values from one to
% match the space of the other.

%% Do a group ROI analysis and make a "violin plot" (or barplot) for each region
%
% We'll use the |barplot_columns| function, which is a versatile function
% for plotting columns of matrices. It has many options for displaying or
% hiding individual data points, "violins" (distribution estimates),
% standard errors, and more. 
%
% The |barplot_columns| function also does a basic one-sample t-test
% analysis on each column, which corresponds to the average of each region
% here. This ROI analysis can be used for inference about a region.
%
% With the pain dataset, for a real analysis we might worry that multiple
% images come from the same subject. We can't do stats across all the
% iamges without accounting for the correlations among images belonging to
% the same subject. But first, let's plot the averages and see waht we get:

% Concatenate the region averages into a matrix:
roi_averages = cat(2, r.dat);

% Get the labels for each region from the atlas object, \
% and replace some characters to make them look a bit nicer:
thal_labels = format_strings_for_legend(thal.labels);

% Now make the figure:
create_figure('Thalamus regions', 1, 2);

barplot_columns(roi_averages, 'nofig', 'colors', scn_standard_colors(length(r)), 'names', thal_labels);
xlabel('Region')
ylabel('ROI activity')

subplot(1, 2, 2)

barplot_columns(roi_averages, 'nofig', 'colors', scn_standard_colors(length(r)), 'names', thal_labels, 'noind', 'noviolin');
xlabel('Region')
ylabel('ROI activity')

%% 
% Look at the results. Do they make sense for a pain dataset?
% The sensory/nociceptive regions of the thalamus (VL and intralaminar) respond to painful
% stimuli, but not the regions in visual and auditory pathways  (LGN, MGN). The
% habenula also responds.

%%
% To do a proper ROI analysis, let's just select a subset of the images
% responding to one stimulus intensity, with one image per person.
% the stimulus temperatures are stored in
% |image_obj.additional_info.temperatures|
% So we can use basic Matlab code to select one temperature.

wh = image_obj.additional_info.temperatures == 48;  % A logical vector for 48 degrees
indx_list = find(wh);                               % a list of which images

image_obj_48 = get_wh_image(image_obj, indx_list);  % Get the images we want

% You can try to extract data from the thalamus and repeat the t-test.

%% Selecting multiple regions using custom label fields
% In addition to the .labels field, atlas objects have fields for additional 
% labels (.labels_2, etc.). These can be used to group, and select, regions 
% at multiple levels of a spatial hierarchy.  For example, ?canlab2018_2mm? 
% atlas has macro-network and anatomical group labels in the ?labels_2? field. 
% You can use these to select all the regions associated with a particular 
% resting-state network or group like the brainstem or basal ganglia. 
% Here?s a code snippet to locate/pull out all the default mode network 
% regions using the select_atlas_subset() method. This method can also return 
% a vector that you can use to select elements of associated connectivity matrices.
%
% You can create and/or annotate your own atlas labels and add them to the
% .labels, labels_2, or other fields, and use select_atlas_subset to
% extract and manipulate them.

atlas_obj = load_atlas('canlab2018_2mm');

% List all the macro-level labels in .labels_2 field
disp(unique(atlas_obj.labels_2)')

% Pull out all regions with ?Def? (Default) as part of the label in labels_2 
atlas_obj = select_atlas_subset(atlas_obj, {'Def'}, 'labels_2');
create_figure('isosurface');
surface_handles = isosurface(atlas_obj);
% let's add a cortical surface for context
% We'll make the back surface (right) opaque, and the front (left) transparent
han = addbrain('hires right');
set(han, 'FaceAlpha', 1);
han2 = addbrain('hires left');
set(han2, 'FaceAlpha', 0.1);
axis image vis3d off
material dull
view(-88, 31);

lightRestoreSingle
drawnow, snapnow

%% Explore on your own
%
% 1. One of the ideas about emotion regulation is that positive appraisal involves
% activity increases in the nucleus accumbens (NaC) in the basal ganglia.
% Load the emotion regulation dataset from previous walkthroughs.
% Try to locate and extract data from this region, and do an ROI analysis.
% What do you see?
% Hint: The NAc is part of the ventral striatum. 
% To extract the ventral striatum from the canonical CANlab atlas, try: 
% NAc = select_atlas_subset(atlas_obj, {'V_Str'});
%
% 2. try the t-test on 48 degree only images above. What do you see? Which
% regions show sigificant activation?
%
% 3. Try to identify which prefrontal regions show the greatest emotion
% regulation effect, and which show the greatest response during pain.
% Are they the same?

% That's it for this section!!




