function [results_table_pos, results_table_neg, r, excluded_region_table, r_excluded, region_list] = table_of_atlas_regions_covered(r, varargin)
% Make a table of which atlas parcels are covered by a set of regions or map identified in a study 
%
% [results_table_pos, results_table_neg, r, excluded_region_table, r_excluded, region_list] = table_of_atlas_regions_covered(r, [atlas_obj])
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

n_cols = 140;                       % 140 good for HTML reports
sep_str = repmat('_', 1, n_cols);   % see textwrap

dosep = true;
dolegend = true;
table_legend_text = '';

[results_table_pos, results_table_neg] = deal(table());

if isempty(r) || length(r) == 0
    disp('No regions to display');
    return
end

if ~isempty(varargin)

    atlas_obj = varargin{1};
    if ~isa(atlas_obj, 'atlas'), error('2nd input must be an atlas-class object'), end
else

    atlas_obj = load_atlas('canlab2018_2mm');

end

[r, percent_coverage] = subdivide_by_atlas(r, atlas_obj);

% attach this so we can deal with reordering of regions
for i = 1:length(r)
    r(i).custom_info2 = percent_coverage(i);
    r(i).custom_info2_descrip = 'atlas parcel coverage';
end

% Separate subregions with positive and negative values if requested
% -------------------------------------------------------------------------
if dosep
    % separate pos and neg
    [poscl, negcl] = posneg_separate(r);

    r = [poscl negcl];
    ispos = [true(1, length(poscl)) false(1, length(negcl))]; % logical for splitting combined cl later

    clear poscl negcl

    fprintf('\n%s\nPositive Effects\n', sep_str)

else

    fprintf('\n%s\nTable of all regions\n', sep_str)
end


% Volume in mm^3
% count per cubic mm - voxel count * voxel volume
get_vox_volume_in_mm = @(r) size(r.XYZ, 2) .* prod(abs(diag(r.M(1:3, 1:3))));

Volume = zeros(length(r), 1);

for i = 1:length(r)

    Volume(i, 1) = get_vox_volume_in_mm(r(i));

end

results_table = table({r.shorttitle}', 'VariableNames', {'Region'});

results_table.Volume = Volume;

results_table.XYZ = round(cat(1, r.mm_center));

results_table = [results_table get_signed_max(r, 'Z', 'maxZ')];  % use function because may be empty, handle if so

results_table.Properties.VariableNames{4} = 'maxVal';

atlas_region_coverage = cat(1, r.custom_info2);

results_table = addvars(results_table, round(atlas_region_coverage), 'NewVariableNames', 'Coverage');

% Exclude
wh_exclude = cat(1, r.numVox) < 2 | atlas_region_coverage < 25;

excluded_region_table = results_table(wh_exclude, :);

results_table(wh_exclude, :) = [];

ispos(wh_exclude) = [];
r_excluded = r(wh_exclude);

r(wh_exclude) = [];

if dosep

    results_table_pos = results_table(ispos, :);
    results_table_neg = results_table(~ispos, :);

end

% Legend info

voxvol = prod(abs(diag(r(1).M(1:3, 1:3))));

table_legend_text = {'Note: XYZ values are Montreal Neurologic Institute coordinates in mm.'};
table_legend_text{2} = sprintf('Voxel volume is %3.2fmm^3', voxvol);
table_legend_text{3} = 'maxVal is maximum positive or negatively-valued value in the .Z field. Its meaning varies across analyses. For a t-test, it usually reflects the maximum t-value.';
table_legend_text{4} = 'Coverage: Percentage of voxels in the labeled atlas parcel covered by the identified activation region.';
table_legend_text{5} = 'Regions shown are those with >1 voxel and >=25% coverage of the atlas parcel. See excluded_region_table for additional regions excluded in this way.';

% Now split into positive and neg sub-tables and display
if dosep
    % Tables for pos and neg results separately

    if any(ispos)
        disp(results_table_pos)
    else
        disp('No regions to display');
    end

    fprintf('\nNegative Effects\n')
    if any(~ispos)
        disp(results_table_neg)
    else
        disp('No regions to display');
    end

else
    % No separation

    if size(results_table, 1) > 0
        disp(results_table)
    else
        disp('No regions to display');
    end

end % dosep

% Get region list
% -------------------------------------------------------------
region_list = cell(4, 1);
for i = 1:4, region_list{i} = ''; end
region_list{1} = 'Positive effects:';
region_list{3} = 'Negative effects:';

for i = 1:size(results_table_pos, 1)
    region_list{2} = [region_list{2} sprintf('%s (%2.0f%%), ', results_table_pos.Region{i}, results_table_pos.Coverage(i))];
end

for i = 1:size(results_table_neg, 1)
    region_list{4} = [region_list{4} sprintf('%s (%2.0f%%), ', results_table_neg.Region{i}, results_table_neg.Coverage(i))];
end


if dolegend == false || isempty(table_legend_text)

    return

else
    % Print legend text

    canlab_print_legend_text(table_legend_text{:});
end

disp('')
disp('List of regions covered:')
for i = 1:4
    disp(region_list{i});
end


% clean up text for legend
% if length(table_legend_text) > 2
%     table_legend_text(2:3) = [];
% end




end % main function






function val_table = get_signed_max(cl, myfield, tablevarname)
% Returns a table with var "MaxZ" nregions x 1, or empty if cl.Z is empty

%maxZ = @(i) max(double(cl(i).Z));

smax = @(i) signedmax(cl(i).(myfield));

for i = 1:length(cl)

    if ~isempty(cl(i).(myfield))
        myZ(i, 1) = smax(i);
    else
        myZ(i, 1) = NaN;
    end

    %if isinf(smax(i)), keyboard, end

end

% Fix infinite vals - only for .Z . so this is not generalizable beyond
% this function without modifications:
maxZ = norminv(1 - 1E-12);
myZ(myZ > maxZ) = maxZ;

if all(isnan(myZ))
    val_table = [];
else
    val_table = table(myZ,  'VariableNames', {tablevarname});
end

end % function



function val = signedmax(vals)

vals = double(vals);
[maxabs, wh] = max(abs(vals));

val = sign(vals(wh)) .* maxabs;

end

