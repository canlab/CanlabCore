function dat_comb = searchlight_applymask_collate(dat2, file_list)
% Estimate local pattern weight on train data using SVR and searchlight and
% apply weights to the test dataset
%
% :Usage:
% ::
%
%     dat_comb = searchlight_applymask(dat2, file_list)
%
% :Inputs:
%
%   **dat2:**
%        fmri_data test object
%
%   **file_list:**
%        Cellarray of list of file distributed in parallel.  Make sure
%        it is sorted correctly (e.g., sort_nat())
%
% :Output:
%
%   **dat_comb:**
%        An fmri_data object with searchlight weights applied to
%        test dataset.
%
% ..
%     Author and copyright information:
%
%     Copyright (C) 2015  Luke Chang & Wani Woo
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
% :Example:
% ::
%
%    dat_comb = searchlight_applymask_collate(dat2, file_list);


if ~isa(dat2,'fmri_data') % Check inputs
    dat1 = fmri_data(dat2); dat2 = remove_empty(dat2);
end

if ~iscell(file_list)
    error('Make sure searchlight file list is a cellarray and that files are on the path')
end

% Collate searchlight .mat files
vox_num = size(dat2.dat,1);
dist_n = length(file_list);
dat_comb = dat2;
for i = 1:dist_n
    load(file_list{i});
    dist_indx = select_voxels(i, dist_n, vox_num);
    dat_comb.dat(dist_indx,:) = dat.dat;
end
end

% ========== SUBFUNCTION ===========


function dist_indx = select_voxels(dist_id, dist_n, vox_num)
% preparation of distribution indx

dist_indx = false(vox_num,1);

unit_num = ceil(vox_num/dist_n);
unit = true(unit_num,1);

start_point = (unit_num.*(dist_id-1)+1);
end_point = min(start_point+unit_num-1, vox_num);
dist_indx(start_point:end_point) = unit(1:end_point-start_point+1);

dist_indx = logical(dist_indx);
end

