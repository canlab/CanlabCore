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
