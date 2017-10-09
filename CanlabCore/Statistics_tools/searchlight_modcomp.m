function [r_corr, dat, sl_size] = searchlight_modcomp(mask1, rdm, varargin)
% Calculate the local pattern similarity between two pattern maps using
% the searchlight approach. 
%
% :Usage:
% ::
%
%     [r_corr, dat, sl_size] = searchlight_correlation(mask1, mask2, [additional_inputs])
% 
%
% :Inputs:
%
%   **mask1:**
%        pattern or activation maps 1
%
%   **mask2:**
%        pattern or activation mediaps 2
%
% :Optional inputs: 
%
%   **'r':**
%        searchlight sphere radius (in voxel) (default: r = 3 voxels)
%
%   **'type':**
%        This calls corr.m, and can take 'type' option.
%
%        'Pearson' (default), 'Kendall', 'Spearman'.
%
% :Outputs:
%
%   **r_corr:**
%        Correlation between weights of two pattern masks
%
%   **dat:**
%        This contains a statistic_image object that contain 
%        correlation values between weights of two pattern masks 
%        (=r_corr; in .dat) and p values for the correlation values 
%        (in .p). 
%
%   **sl_size:**
%        The number of voxels within each searchlight. Based on this 
%        number, you can detect searchlights on the edge (searchlights 
%        with low sl_size should be on the edge of the brain.
%
% ..
%     Author and copyright information:
%
%     Copyright (C) 2017  Phil Kragel
%     Copyright (C) 2014  Wani Woo
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
% :Examples:
% ::
%
%    mask1 = which('weights_NSF_grouppred_cvpcr.img');
%    mask2 = which('nonnoc_v6_109subjmap_mean.nii');
%
%    [r, dat] = searchlight_correlation(mask1, mask2, 'r', 5);

 
r = 4; % default radius (in voxel)
corr_type = 'Pearson';

% parsing varargin

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            % functional commands
            case 'type'
                corr_type = varargin{i+1};
            case 'r' % radius
                r = varargin{i+1};
        end
    end
end

dat1 = fmri_data(mask1); dat1 = remove_empty(dat1);


dat1 = remove_empty(dat1);
n = size(dat1.dat,1);

r_corr = NaN(n,1);
p_corr = NaN(n,1);
sl_size = zeros(n,1);

fprintf('\n Calculating corrleation for voxel                 ');
xyzlist=dat1.volInfo.xyzlist(~dat1.removed_voxels,:);
dat=dat1.dat;
for i = 1:n %(1):vox_to_run(10)
    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%07d/%07d', i, n);
    searchlight_indx = searchlight_sphere_prep(xyzlist, i, r);
    [r_corr(i), p_corr(i)] = corr(pdist(dat(searchlight_indx,:)','correlation')',rdm', 'type', corr_type);
    sl_size(i) = sum(searchlight_indx);
end

dat = statistic_image;
eval(['dat.volInfo = dat' num2str(1) '.volInfo;']);
eval(['dat.removed_voxels = dat' num2str(1) '.removed_voxels;']);
dat.dat = r_corr;
dat.p = p_corr;

end

% ========== SUBFUNCTION ===========

function indx = searchlight_sphere_prep(xyz, i, r)
seed = xyz(i,:);
indx = sum([xyz(:,1)-seed(1) xyz(:,2)-seed(2) xyz(:,3)-seed(3)].^2, 2) <= r.^2;
end


