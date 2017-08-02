% Note: You need the CANlab Core Tools to run this.
% https://github.com/canlab

[cortical_images, linenames, imgnames] = load_image_set('pauli_cortex');

%% POLAR PLOT OF REALATION WITH BUCKNERLAB MAPS

% !gunzip *gz
% f = filenames('*nii')
% 
% cortical_images = fmri_data(char(f{2:end})); orthviews(cortical_images);

%colors = {[.3 .3 1] [.3 1 .3] [.5 1 .7] [1 .3 .3] [.7 .2 .7]};
colors = {[.2 1 .2] [.65 1 .85] [.3 .3 1] [.8 .3 .8] [1 .3 .3] };

create_figure;
image_similarity_plot(cortical_images, 'bucknerlab', 'colors', colors, 'cosine_similarity');

linenames = {'Ca' 'Cp' 'Pa' 'Pp' 'VS'};
makelegend(linenames, colors);

%% PLOT ON LATERAL SURFACES

cortical_images = threshold(cortical_images, [90 Inf], 'raw-between');

create_figure('surfaces', 1, 2)
hh = addbrain('hires');
view(-85, 10); lightRestoreSingle;
      
wh = 1:5;

n = length(wh);

for i = 1:n
    
    r = region(get_wh_image(cortical_images, wh(i)));
    
    surface_handles{i} = surface_cutaway('cl', r, 'surface_handles', hh, 'existingfig', 'pos_colormap', colors{wh(i)});
    
end

subplot(1, 2, 2)

hh = addbrain('hires');
view(85, 10); lightRestoreSingle;
   
for i = 1:n
    
    r = region(get_wh_image(cortical_images, wh(i)));
    
    surface_handles{i} = surface_cutaway('cl', r, 'surface_handles', hh, 'existingfig', 'pos_colormap', colors{wh(i)});
    
end

%% PLOT ON MONTAGES

% orthviews(cortical_images)

o2 = canlab_results_fmridisplay('montagetype', 'compact2');

for i = 1:5
    my_img = get_wh_image(cortical_images, i);
    my_regions = region(my_img);
    o2 = addblobs(o2, my_regions, 'color', colors{i}, 'trans');
end

%% Or: Differentiate the greens a bit more


o2 = removeblobs(o2);

for i = 1:5
    my_img = get_wh_image(cortical_images, i);
    my_regions = region(my_img);
    o2 = addblobs(o2, my_regions, 'color', colors{i}, 'trans');
end

%or 'striatum' directory