%% Datasets used in CANlab tutorials
%
% Note: this report was generated from |canlab_help_2c_loading_datasets.m|
% in the repository <https://github.com/canlab/CANlab_help_examples>
%

%% The Neuroimaging_pattern_masks Github repository and website
%
% For these walkthroughs, we'll load several image datasets.
% Some are stored in Github, if the files are small enough.
% The Wager 2008 emotion regulation dataset, for example, is
% in the CANlab Core repository. Others are on figshare or Neurovault.
%
% Many tutorials apply pre-trained patterns and masks. 
% These are stored in this Github repository:
% 
% <https://github.com/canlab/Neuroimaging_Pattern_Masks>
% 
% In addition, you can explore these and find more information here:
% <https://sites.google.com/dartmouth.edu/canlab-brainpatterns/home>
%
% To run this tutorial and many others in this series, you'll need to
% download the Neuroimaging_Pattern_Masks Github repository and add it to
% your Matlab path with subfolders.
% The Github site has three main types of datasets, shown here:
% 
% <<Neuroimaging_pattern_masks_site.png>>
%
%% Registries for easy loading
% There are several functions that contain sets of images that you can load
% by name. For example:
%
% |load_image_set()| includes a registry of datasets that you can access by
% name. Try |help load_image|set| for a list of images you can load by
% keyword.  
%
% |load_atlas()| also has a named set of atlases/brain parcellations you can load.
%
% |canlab_load_ROI| has a registry of many named regions derived from previous
% studies. This is particularly useful for loading a subcortical ROI and
% visualizing it or applying it to data.

%% Emotion regulation dataset
% "Wager_et_al_2008_Neuron_EmotionReg"
% The dataset is a series of contrast images from N = 30 participants.
% Each image is a contrast image for [reappraise neg vs. look neg]
% 
% These data were published in:
% Wager, T. D., Davidson, M. L., Hughes, B. L., Lindquist, M. A., 
% Ochsner, K. N.. (2008). Prefrontal-subcortical pathways mediating 
% successful emotion regulation. Neuron, 59, 1037-50.
%
% To load, try the following:

[data_obj, subject_names, image_names] = load_image_set('emotionreg', 'noverbose');

montage(mean(data_obj), 'trans', 'compact2');
drawnow, snapnow

% data_obj is an fmri_data object containing all 30 images.
% subject_names is a list of short names for each image.
% image_names is a list of the full image names with their path.

%% Buckner lab resting-state connectivity maps
%
% A popular series of masks is a set of "networks" developed from resting-state fMRI
% connectivity in 1,000 people. For our purposes, the "network" maps
% consist of a set of voxels that load most highly on each of 7 ICA
% components from the 1,000 Functional Connectomes project. 
% These were published by Randy Buckner's lab in 2011 in three papers:
% Choi et al., Buckner et al., and Yeo et al.  We'll use the cortical
% networks from Yeo et al. here.
%

[bucknermaps, mapnames] = load_image_set('bucknerlab', 'noverbose'); % loads 7 masks from Yeo et al.
disp(char(mapnames))

% Create a montage plot
o2 = montage(get_wh_image(bucknermaps, 1), 'trans', 'compact2', 'color', rand(1, 3));
for i = 2:7, o2 = addblobs(o2, region(get_wh_image(bucknermaps, i)), 'trans', 'color', rand(1, 3)); end
drawnow, snapnow

%% Pain, Cognition, Emotion balanced N = 270 dataset
%
% This dataset is an fmri_data object class object created using CANlab tools (canlab.github.io). 
% It contains 270 single-participant fMRI contrast maps across 18 studies (with 15 participants per study). 
% Studies are grouped into three domains: Pain, Cognitive Control, and Negative Emotion, with 9 studies each. 
% Each domain includes 3 subdomains, with 3 studies each.
%
% The dataset is from Kragel et al. 2018, Nature Neuroscience. 
% It is too large for Github, and it's stored on Neurovault.org
% as collection #504. You could get it using CANlab tools like this:
%
% |[files_on_disk, url_on_neurovault, mycollection, myimages] = retrieve_neurovault_collection(504);|
%
% It is also available on Figshare, with DOI 10.6084/m9.figshare.24033402
% https://figshare.com/ndownloader/files/42143352
%
% If you download from Neurovault, you'd
% have to add metadata for the study category labels to use it.
% Therefore, we suggest you use the CANLab load_image_set() function, as below. 
% It saves a file on your hard drive the first time you run it:
% 
% |kragel_2018_nat_neurosci_270_subjects_test_images.mat|
%
% This file also has metadata that is not necessarily on Neurovault.
% If you save this file somewhere on your Matlab path, you'll be able to
% reload and reuse the dataset easily.

[test_images, names] = load_image_set('kragel18_alldata', 'noverbose');

% This field contains a table object with metadata for each image:
metadata = test_images.metadata_table;

disp('Metadata available in test_images.metadata_table:')
metadata(1:5, :)

% Make a plot of the spatial correlation of the average image for each study

imgs = cellstr(test_images.image_names);
m = mean(test_images);

for i = 1:length(names)
    
    % Create a mean image for each study and store in "m" object.
    
    wh = metadata.Studynumber == i;
    studymean = mean(get_wh_image(test_images, find(wh)));
    m.dat(:, i) = studymean.dat;
    
end

disp('Map of spatial correlations across the mean images for each study');

plot_correlation_matrix(m.dat, 'names', names, 'partitions', [ones(1, 6) 2*ones(1, 6) 3*ones(1, 6)], 'partitionlabels', {'Pain', 'Cognition', 'Emotion'});
drawnow, snapnow

%% Pain across six intensity levels per person (BMRK3)
%
% The dataset contains data from 33 participants, with brain responses to six levels
% of heat (non-painful and painful). Each image is the average over several
% (4-8) trials of heat delivered at a single stimulus intensity, ranging
% from 44.3 - 49.3 degrees C in one-degree increments. Each image is also
% paired with an average reported pain value for that set of trials, rated
% immmediately after heat experience. 
%
% This dataset is interesting for mixed-effects and predictive analyses, as
% it has both within-person and between-person sources of variance.
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
% The key variable is image_obj, an fmri_data class object.
% See https://canlab.github.io/
%
% image_obj.dat contains brain data for each image (average across trials)
% image_obj.Y contains pain ratings (one average rating per image)
%
% image_obj.additional_info.subject_id contains integers coding for which
% Load the data file, downloading from figshare if needed

fmri_data_file = which('bmrk3_6levels_pain_dataset.mat');

if isempty(fmri_data_file)
    
    % attempt to download
    disp('Did not find data locally...downloading data file from figshare.com')
    
    fmri_data_file = websave('bmrk3_6levels_pain_dataset.mat', 'https://ndownloader.figshare.com/files/12708989');
    
end

load(fmri_data_file);

descriptives(image_obj);

subject_id = image_obj.additional_info.subject_id;
ratings = reshape(image_obj.Y, 6, 33)';
temperatures = image_obj.additional_info.temperatures;

% Plot the ratings

create_figure('ratings');
title('Pain ratings by stimulus intensity in BMRK3 dataset')
hold on; plot(ratings', '-', 'Color', [.7 .7 .7], 'LineWidth', .5);
lineplot_columns(ratings, 'color', [.7 .3 .3], 'markerfacecolor', [1 .5 0]);
xlabel('Temperature');
ylabel('Rating');
set(gca, 'XTickLabel', [44.3:49.3], 'FontSize', 18);
drawnow, snapnow

%% CANlab 2018 Combined Atlas object
% An atlas-class object is a specialized subclass of fmri_data that stores
% information about a series of parcels, or pre-defined regions, and in
% some cases the probalistic maps that underlie the final parcellation.
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

atlas_obj = load_atlas('canlab2018');

% Create a figure showing cortical parcels
create_figure('iso'); isosurface(atlas_obj, 'alpha', .5);
drawnow, snapnow

% Create a figure of the thalamus regions only
thal = select_atlas_subset(atlas_obj, {'Thal'});
create_figure('iso'); isosurface(thal);
hh = addbrain('hires');
set(hh, 'FaceAlpha', 0.05);
set(gca, 'XLim', [-40 40], 'YLim', [-40 20]);
drawnow, snapnow


%% Kragel 2015 emotion category-predictive patterns
% 
% This loads a series of 7 multivariate patterns from Kragel et al. 2015.
% These patterns were developed to predict emotion categories.
% For our purposes here, we'll just treat them as data images.

[test_dat, names] = load_image_set('kragelemotion', 'noverbose'); % loads 7 masks from Kragel et al.
disp(char(names))

% Plot spatial correlations
plot_correlation_matrix(test_dat.dat, 'names', names);
drawnow, snapnow

%% Subcortical and brainstem regions
% |canlab_load_ROI| has a registry of many named regions derived from previous
% studies. This is particularly useful for loading a subcortical ROI and
% visualizing it or applying it to data.

help canlab_load_ROI

vmpfc = canlab_load_ROI('vmpfc');
nac =  canlab_load_ROI('nac');
pag =  canlab_load_ROI('pag');

create_figure('brain');
isosurface(vmpfc, 'color', {[1 0 1]});
isosurface(nac, 'color', {[1 .5 0]});
isosurface(pag, 'color', {[1 0 0]});
bhan = addbrain('hires right');
view(-119, 5);
drawnow, snapnow

