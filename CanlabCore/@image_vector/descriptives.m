function desc = descriptives(dat, varargin)
% Get descriptives for an fmri_data or other image_vector object
% - Returns a structure with useful numbers: min/max, percentiles
% - Vectors for nonempty and complete voxels and images
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

% By convention, zero indicates missing (empty) data and is not a valid value.

desc.nonempty_vox_descrip = 'Voxels with non-zero, non-NaN data values for at least one image';
desc.nonempty_voxels = ~all(desc.wh_zero | desc.wh_nan, 2);
desc.n_nonempty_vox = sum(desc.nonempty_voxels);
desc.n_in_mask = dat.volInfo.n_inmask;

desc.complete_voxels = ~any(desc.wh_zero | desc.wh_nan, 2);
desc.n_complete_vox = sum(desc.complete_voxels);

desc.nonempty_image_descrip = 'Images with non-zero, non-NaN data values for at least one nonempty voxel';
desc.nonempty_images = ~all(desc.wh_zero(desc.nonempty_voxels, :) | desc.wh_nan(desc.nonempty_voxels, :), 1);
desc.n_nonempty_images = sum(desc.nonempty_images);

desc.complete_image_descrip = 'Images with non-zero, non-NaN data values for all nonempty voxels';
desc.complete_images = ~any(desc.wh_zero(desc.nonempty_voxels, :) & desc.wh_nan(desc.nonempty_voxels, :), 1);
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
    
    fprintf('Voxels: %3.0f\tNonempty: %3.0f\tComplete: %3.0f\n', desc.n_vox, desc.n_nonempty_vox, desc.n_complete_vox);
    
    fprintf('Unique data values: %3.0f\n', desc.num_unique_vals);
    
    disp(' ')
    
    fprintf('Min: %3.3f\tMax: %3.3f\tMean: %3.3f\tStd: %3.3f\n', desc.min, desc.max, desc.mean, desc.std);
    
    disp(' ');
    
    disp(desc.prctile_table);
    
    disp(' ');
    
end


end % function

