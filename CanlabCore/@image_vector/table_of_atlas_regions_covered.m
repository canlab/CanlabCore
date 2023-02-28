function [results_table_pos, results_table_neg, r, excluded_region_table, r_excluded, region_list] = table_of_atlas_regions_covered(obj, varargin)
% Make a table of which atlas parcels are covered by a set of regions or map identified in a study 
%
% [results_table_pos, results_table_neg, r, excluded_region_table, r_excluded, region_list] = table_of_atlas_regions_covered(obj, [atlas_obj])
%
% fMRI activation maps often include large suprathreshold areas that span
% multiple brain regions.  Traditional tables divide areas by contiguous
% regions ("blobs"), but this type of division is not useful if blobs cannot 
% be easily summarized by a single anatomical label. A complementary approach
% is to identify which areas defined in a standard brain atlas ("parcels")
% are covered by the activation map. Parcels can be defined based on anatomy 
% (e.g., cytoarchitecture) or function (e.g., resting-state fMRI).  
% For example, a large blob may span the insula, claustrum, and putamen. 
% This can be described in a table that lists all of these regions, along with 
% the percentage of each parcel covered by the activation map. 
% 
% This function uses the object method subdivide_by_atlas() to create such
% a table, which is dividided into two tables with positive and negative average
% statistic values (results_table_pos, results_table_neg). These are Matlab
% table-class objects.  It also returns a region-class object for reference 
% and rendering the subdivided blobs on brain slices or surfaces.
%
% You can enter any atlas-class object as a 2nd argument. If you don't,
% a default atlas object will be used.
%
% Examples:
% -----------------------------------------------------------------
% % Load a thresholded map:
% % (Note: You must have Neuroimaging_Pattern_Masks repo with subfolders on your path)
% % This is a predictive map for drug craving (Koban et al. 2022, Nat Neurosci)
%
% ncsthr = fmri_data(which('NCS_multithr_001k5_005_05_pruned.nii'), 'noverbose');
% 
% % Convert to region object:
% r = region(ncsthr);
%
% [results_table_pos, results_table_neg, r, excluded_region_table, table_legend_text] = table_of_atlas_regions_covered(r);
% montage(r, 'regioncenters', 'colormap');
%
% % Also see large blobs that didn't sufficiently cover any one atlas parcel:
% montage(r_excluded(cat(1, r_excluded.numVox) > 20), 'regioncenters', 'colormap');
%
% % You can use the image_vector method to apply this to a single-image
% % thresholded map:
% [results_table_pos, results_table_neg, r, excluded_region_table, r_excluded, region_labels] = table_of_atlas_regions_covered(ncsthr);
%
% Note: For a key to labels when using the Glasser 2016 Nature atlas, or the 
% CANlab combined atlas (which uses Glasser), see GlasserTableS1.pdf in the
% original paper or CANlab github, or http://braininfo.rprc.washington.edu/
%
% see also region.table, for autolabeling of regions with an atlas and more
% detailed table output.
%
% Tor Wager, Feb 2023


if (size(obj.dat, 2) > 1)
    error('Use this function with fmri_data or image_vector objects containing only a single image. Use get_wh_image() to select one.')
end

r = region(obj);

[results_table_pos, results_table_neg, r, excluded_region_table, r_excluded, region_list] = table_of_atlas_regions_covered(r, varargin{:});

end


