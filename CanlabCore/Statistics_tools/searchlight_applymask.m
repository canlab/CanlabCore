function [dat, sl_size] = searchlight_applymask(dat1, dat2, varargin)
% Estimate local pattern weight on train data using SVR and searchlight and 
% apply weights to the test dataset 
%
% :Usage:
% ::
%
%     [dat, sl_size] = searchlight_applymask(train, test, varargin)
%
%
% :Inputs:
%
%   **dat1:**
%        fmri_data object with train.Y == size(train.dat,2)
%
%   **dat2:**
%        Data to apply local weight map.
%
% :Optional inputs:
%
%   **'r':**
%        searchlight sphere radius (in voxel) (default: r = 3 voxels)
%
%   **'parallel':**
%        run subset of voxels to distribute on a cluster.  flag must
%               be followed by array specifing id and total number of jobs
%               (e.g., 'parallel',[1,10]);
%
% :Outputs:
%
%   **dat:**
%        This contains an fmri_data object that contain
%        correlation pattern expression values
%
%   **sl_size:**
%        The number of voxels within each searchlight. Based on this
%        number, you can detect searchlights on the edge (searchlights
%        with low sl_size should be on the edge of the brain.
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
% :Examples:
% ::
%
%    [r, dat] = searchlight_applymask(train, test, 'r', 5);
%
%    [r, dat] = searchlight_applymask(train, test, 'r', 5,'parallel',[1,10]);


r = 4; % default radius (in voxel)
doParallel = 0;

% parsing varargin
for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            % functional commands
            case 'r' % radius
                r = varargin{i+1};
                varargin{i} = []; varargin{i+1} = [];
            case 'parallel'
                doParallel = 1;
                tmp = varargin{i+1};
                dist_id = tmp(1);
                dist_n = tmp(2);
                varargin{i} = []; varargin{i+1} = [];
        end
    end
end

if ~isa(dat1,'fmri_data')
    dat1 = fmri_data(dat1); dat1 = remove_empty(dat1);
end
if ~isa(dat2,'fmri_data')
    dat2 = fmri_data(dat2); dat2 = remove_empty(dat2);
end

isdiff = compare_space(dat1, dat2);

if isdiff == 1 % diff space, not just diff voxels
    % resample to the image with bigger voxel size
    [~, idx] = max([abs(dat1.volInfo.mat(1,1)), abs(dat2.volInfo.mat(1,1))]);
    
    if idx == 1
        dat2 = resample_space(dat2, dat1);
    else
        dat1 = resample_space(dat1, dat2);
    end
    
end

dat1 = remove_empty(dat1);
dat2 = remove_empty(dat2);

% if any(sum(sum(dat1.volInfo.xyzlist ~= dat2.volInfo.xyzlist))), error('dat and dat2 should be in the same space.'); end

if doParallel
    fprintf(['Running in Parallel: id' num2str(dist_id) ' of ' num2str(dist_n)])
    vox_num = size(dat1.dat,1);
    dist_indx = select_voxels(dist_id, dist_n, vox_num);
    dat1.dat = dat1.dat(dist_indx, :);
    dat1.removed_voxels(~dist_indx) = true;
    dat2.dat = dat2.dat(dist_indx, :);
    dat2.removed_voxels(~dist_indx) = true;
end

[~, idx] = min([size(dat1.dat,1) size(dat2.dat,1)]);

eval(['n = size(dat' num2str(idx) '.dat,1);']);

r_corr = NaN(n,size(dat2.dat,2));
sl_size = zeros(n,1);

dd1 = dat1;
dd2 = dat2;

fprintf('\n Calculating correlation for voxel                   ');
for i = 1:n %(1):vox_to_run(10)
    clc
    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%07d/%07d', i, n);
    eval(['searchlight_indx = searchlight_sphere_prep(dat' num2str(idx) ...
        '.volInfo.xyzlist(~dat' num2str(idx) '.removed_voxels,:), i, r);']);
    dd1.dat = dat1.dat(searchlight_indx,:);
    [yfit, wt] = alg_svr(dd1);
    
    dd2.dat = dat2.dat(searchlight_indx,:);
    for j = 1:size(dd2.dat, 2)
        r_corr(i,j) = corr(dd2.dat(:,j),wt.dat);
    end
    sl_size(i) = sum(searchlight_indx);
end

dat = fmri_data;
eval(['dat.volInfo = dat' num2str(idx) '.volInfo;']);
eval(['dat.removed_voxels = dat' num2str(idx) '.removed_voxels;']);
dat.dat = r_corr;

end

% ========== SUBFUNCTION ===========

function [yfit, wt, int] = alg_svr(dat, varargin);
% run SVR on data object

dat.dat = double(dat.dat);
dat.Y = double(dat.Y);

dataobj = data('spider data', dat.dat', dat.Y);

% slack parameter
slackstr = 'C=1';
wh = find(strcmpi('C=', varargin));
if ~isempty(wh), slackstr = varargin{wh(1)}; end

svrobj = svr({slackstr, 'optimizer="andre"'});
[res, svrobj] = train(svrobj, dataobj);
yfit = res.X;
wt = dat;
wt.dat = get_w(svrobj)';
end

function indx = searchlight_sphere_prep(xyz, i, r)
seed = xyz(i,:);
indx = sum([xyz(:,1)-seed(1) xyz(:,2)-seed(2) xyz(:,3)-seed(3)].^2, 2) <= r.^2;
end

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

