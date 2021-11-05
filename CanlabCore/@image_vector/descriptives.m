function desc = descriptives(dat, varargin)
% Get descriptives for an fmri_data or other image_vector object
% - Returns a structure with useful numbers: min/max, percentiles
% - Vectors for nonempty and complete voxels and images
% - Returns summary fmri_data objects showing coverage: numbers of images with valid data in each voxel 
%
% Image_vector (and subclass fmri_data) objects are 4-d datasets in which
% 3-D images are vectorized into columns in a 2-D matrix.
% - For data object dat, dat.dat contains the data
% - By convention, zero indicates missing (empty) data and is not a valid value
%
% :Usage:
% ::
%
%     desc = descriptives(dat, ['noverbose'])
%
% For objects: Type methods(object_name) for a list of special commands
%              Type help object_name.method_name for help on specific
%              methods.
%
% ..
%     Author and copyright information:
%
%     Copyright (C) 2018 Tor Wager
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
% ..
%
% :Optional Inputs:
%
%   **'noverbose':**
%        Suppress printing of output summary
%
% :Outputs:
% **desc.coverage_obj**
% a summary fmri_data object showing coverage: numbers of images with valid data in each voxel   
%
% **desc.coverage_obj_binned**
% fmri_data object showing coverage binned by percentiles of images:
% all images have data = 100
% 80% - 99.9% (all but one) = 80
% 50% - 99.9% (all but one) = 50
% 1 image - 50% = 50
% Rationale: Sometimes it's hard to see if there is only one or a few
% images missing voxels/areas out of a large set.
%
% **desc.coverage_obj_complete** 
% fmri_data object showing voxels with complete coverage
%
% Example:
% % Load a standard dataset and create a color map of the coverage
% obj = load_image_set('emotionreg');
% desc = descriptives(obj);
% o2 = montage(desc.coverage_obj_binned, 'trans', 'maxcolor', [.5 1 .5], 'mincolor', [1 0 0], 'cmaprange', [1 100], 'transvalue', 0.8);
%
% Show areas with valid data for all voxels
% o2 = montage(desc.coverage_obj_complete, 'trans', 'maxcolor', [.5 1 .5], 'mincolor', [0 0 0], 'cmaprange', [1 100], 'transvalue', 0.8);
%
% Show histogram of missing voxels per image
% figure; hist(desc.percent_missing_per_image, 30)
% xlabel('Percentage of voxels missing'); ylabel('Number of images');

% :See also:
%   - methods(image_vector) and methods(fmri_data)
%

% ..
%    Programmers' notes:
%    List dates and changes here, and author of changes
% ..

% ..
%    DEFAULTS AND INPUTS
% ..

doverbose = true;
% initalize optional variables to default values here.


% optional inputs with default values
for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            
            case 'noverbose', doverbose = false;
                
            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end

desc.n_images = size(dat.dat, 2);
desc.n_vox = size(dat.dat, 1);

desc.wh_zero = dat.dat == 0;
desc.wh_nan = isnan(dat.dat);

% Object with summary of coverage - how many images have non-zero, non-NaN
% values in each voxel
m = mean(dat);
m.dat = sum(~isnan(dat.dat) & dat.dat ~= 0, 2);
desc.coverage_obj = m;

% Object with coverage summary binned by percentage categories for easy viewing.
cutoffs = [floor(desc.n_images .* .5) floor(desc.n_images .* .8)]; 
obj2 = m;
obj2.dat(m.dat == desc.n_images) = 100;
obj2.dat(m.dat > 0 & m.dat <= cutoffs(1)) = 1;
obj2.dat(m.dat > cutoffs(1) & m.dat <= cutoffs(2)) = 50;
obj2.dat(m.dat > cutoffs(2) & m.dat < desc.n_images) = 80;

desc.coverage_obj_complete = m;
desc.coverage_obj_complete.dat(m.dat == desc.n_images) = 100;
desc.coverage_obj_complete.dat(m.dat ~= desc.n_images) = 0;

% orthviews(desc.coverage_obj_binned, 'continuous')

% By convention, zero indicates missing (empty) data and is not a valid value.

desc.nonempty_vox_descrip = '.nonempty_voxels: Voxels with non-zero, non-NaN data values for at least one image';
desc.nonempty_voxels = ~all(desc.wh_zero | desc.wh_nan, 2);
desc.n_nonempty_vox = sum(desc.nonempty_voxels);
desc.n_in_mask = dat.volInfo.n_inmask;

% percentage of missing voxels per image
desc.percent_missing_per_image_descrip = sprintf('.percent_missing_per_image: percentage of coverage area in desc.nonempty_voxels (%3.0f voxels) with valid data for each image.', desc.n_nonempty_vox);  
desc.percent_missing_per_image = 100 .* sum((desc.wh_zero | desc.wh_nan) & desc.nonempty_voxels) ./ sum(desc.nonempty_voxels);
desc.images_missing_over_50percent = desc.percent_missing_per_image >= 49.5;
desc.images_missing_over_10percent = desc.percent_missing_per_image >= 9.5;

desc.complete_voxels_descrip = '.complete_voxels: Voxels with non-zero, non-NaN data values for all images';
desc.complete_voxels = ~any(desc.wh_zero | desc.wh_nan, 2);
desc.n_complete_vox = sum(desc.complete_voxels);

desc.nonempty_image_descrip = 'Images with non-zero, non-NaN data values for at least one nonempty voxel';
desc.nonempty_images = ~all(desc.wh_zero(desc.nonempty_voxels, :) | desc.wh_nan(desc.nonempty_voxels, :), 1);
desc.n_nonempty_images = sum(desc.nonempty_images);

desc.complete_image_descrip = 'Images with non-zero, non-NaN data values for all nonempty voxels';
desc.complete_images = desc.percent_missing_per_image == 0; % ~any(desc.wh_zero(desc.nonempty_voxels, :) & desc.wh_nan(desc.nonempty_voxels, :), 1);
desc.n_complete_images = sum(desc.complete_images);

datavec = dat.dat(desc.nonempty_voxels, desc.nonempty_images);
datacat = double(datavec(:));
datacat(datacat == 0 | isnan(datacat)) = [];  % still need to remove invalid voxels

desc.min = min(datacat);
desc.max = max(datacat);
%desc.quartiles_25_50_75 = prctile(datacat, [25 50 75]);

desc.unique_vals = unique(datacat);
desc.num_unique_vals = length(desc.unique_vals);

desc.prctiles = [.1 .5 1 5 25 50 75 95 99 99.5 99.9];
desc.prctile_vals = prctile(datacat, desc.prctiles);

Percentiles = desc.prctiles';
Values = desc.prctile_vals';
desc.prctile_table = table(Percentiles, Values);

desc.mean = mean(datacat);
desc.std = std(datacat);

if doverbose
    
    try
        % try...catch here because these fields may not contain expected
        % data types.
        
        disp(' ')
        
        if ~isempty(dat.source_notes)
            
            fprintf('Source: %s\n', dat.source_notes);
            
        end
        
        if ~isempty(dat.dat_descrip)
            
            fprintf('Data: %s\n', dat.dat_descrip);
            
        end
        
        disp(' ')
        % References (ad hoc, not obligatory field in fmri_data
        
        if isstruct(dat.additional_info) && isfield(dat.additional_info, 'references')
            myrefs = dat.additional_info.references;
            if isempty(myrefs)
                % do nothing
            else
                if iscell(myrefs)
                    canlab_print_legend_text(dat.additional_info.references{:})
                else
                    canlab_print_legend_text(dat.additional_info.references);
                end
            end
        end
        
    catch
        disp('Source and dat_info fields do not have expected char arrays; skipping printout')
        
    end
    
    disp(' ')
    disp('Summary of dataset')
    disp('______________________________________________________')
    
    fprintf('Images: %3.0f\tNonempty: %3.0f\tComplete: %3.0f\n', desc.n_images, desc.n_nonempty_images, desc.n_complete_images);
    
    fprintf('  Images missing >50%% of voxels: %3.0f\t', sum(desc.images_missing_over_50percent));
    if any(desc.images_missing_over_50percent)
        fprintf('%3.0f ', find(desc.images_missing_over_50percent));
    end
    fprintf('\n');
    
    fprintf('  Images missing >10%% of voxels: %3.0f\t', sum(desc.images_missing_over_10percent));
    if any(desc.images_missing_over_10percent)
        fprintf('%3.0f ', find(desc.images_missing_over_10percent));
    end
    fprintf('\n');
    
    fprintf('Voxels: %3.0f\tNonempty (1+ images have valid data): %3.0f\tComplete  (all images have data): %3.0f\n', desc.n_vox, desc.n_nonempty_vox, desc.n_complete_vox);
    
    fprintf('Unique data values: %3.0f\n', desc.num_unique_vals);
    
    disp(' ')
    
    fprintf('Min: %3.3f\tMax: %3.3f\tMean: %3.3f\tStd: %3.3f\n', desc.min, desc.max, desc.mean, desc.std);
    
    disp(' ');
    
    disp(desc.prctile_table);
    
    disp(' ');
    
    disp('Saved desc.coverage_obj_binned with maps of number of valid images,')
    disp('and desc.coverage_obj_binned with values of 100=100% valid images, 80=80-99.9% valid, 50=50-80% valid, and 1=>50% valid images');
    
    % warnings
    if desc.n_nonempty_vox > desc.n_complete_vox
        disp('Warning: Some voxels have data for 1+ images but not for all images')
    end
    
    if any(desc.images_missing_over_50percent)
            disp('Warning: Some images are missing over 50% of voxels')
    end
    
end


end % function

