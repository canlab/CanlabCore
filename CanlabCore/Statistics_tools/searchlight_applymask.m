function [dat, sl_size] = searchlight_applymask(train, test, varargin)

% Calculate the local pattern similarity between two pattern maps using
% the searchlight approach. 
%
% Usage:
% -------------------------------------------------------------------------
% [r_corr, dat, sl_size] = searchlight_correlation(mask1, mask2, [additional_inputs])
% 
%
% Inputs:
% -------------------------------------------------------------------------
% mask1         pattern or activation maps 1
% mask2         pattern or activation mediaps 2
% 
% Optional inputs: 
% -------------------------------------------------------------------------
% 'r'           searchlight sphere radius (in voxel) (default: r = 3 voxels)
% 'type'        This calls corr.m, and can take 'type' option.
%               'Pearson' (default), 'Kendall', 'Spearman'.
% Outputs:
% -------------------------------------------------------------------------
% r_corr        Correlation between weights of two pattern masks
% dat           This contains a statistic_image object that contain 
%               correlation values between weights of two pattern masks 
%               (=r_corr; in .dat) and p values for the correlation values 
%               (in .p). 
% sl_size       The number of voxels within each searchlight. Based on this 
%               number, you can detect searchlights on the edge (searchlights 
%               with low sl_size should be on the edge of the brain.
%
% Author and copyright information:
% -------------------------------------------------------------------------
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
%
% Example 
% -------------------------------------------------------------------------
% mask1 = which('weights_NSF_grouppred_cvpcr.img');
% mask2 = which('nonnoc_v6_109subjmap_mean.nii');
%
% [r, dat] = searchlight_correlation(mask1, mask2, 'r', 5);
%

% Programmers' notes:
% 

%% set-up variables
 
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
dat2 = fmri_data(mask2); dat2 = remove_empty(dat2);

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

if any(sum(sum(dat1.volInfo.xyzlist ~= dat2.volInfo.xyzlist))), error('dat and dat2 should be in the same space.'); end

[~, idx] = min([size(dat1.dat,1) size(dat2.dat,1)]);

eval(['n = size(dat' num2str(idx) '.dat,1);']);

r_corr = NaN(n,1);
p_corr = NaN(n,1);
sl_size = zeros(n,1);

fprintf('\n Calculating corrleation for voxel                 ');
for i = 1:n %(1):vox_to_run(10)
    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%07d/%07d', i, n);
    eval(['searchlight_indx = searchlight_sphere_prep(dat' num2str(idx) ...
        '.volInfo.xyzlist(~dat' num2str(idx) '.removed_voxels,:), i, r);']);
    [r_corr(i), p_corr(i)] = corr(dat1.dat(searchlight_indx), dat2.dat(searchlight_indx), 'type', corr_type);
    sl_size(i) = sum(searchlight_indx);
end

dat = statistic_image;
eval(['dat.volInfo = dat' num2str(idx) '.volInfo;']);
eval(['dat.removed_voxels = dat' num2str(idx) '.removed_voxels;']);
dat.dat = r_corr;
dat.p = p_corr;

end

% ========== SUBFUNCTION ===========

function [yfit, w, int] = alg_svr(dat, varargin);


dataobj = data('spider data', dat.dat, dat.Y);

% slack parameter
slackstr = 'C=1';
wh = find(strcmpi('C=', varargin));
if ~isempty(wh), slackstr = varargin{wh(1)}; end

svrobj = svr({slackstr, 'optimizer="andre"'});
[res, svrobj] = train(svrobj, dataobj);
w = get_w(svrobj);

end

function [yfit, w, dummy] = cv_svr(xtrain, ytrain, xtest, cv_assignment, varargin)

dataobj = data('spider data', xtrain, ytrain);

% Define algorithm
%svrobj = svr({'C=1', 'optimizer="andre"', kernel({'rbf', 1})})

% slack parameter
slackstr = 'C=1';
wh = find(strcmpi('C=', varargin));
if ~isempty(wh), slackstr = varargin{wh(1)}; end

svrobj = svr({slackstr, 'optimizer="andre"'});

% Training
fprintf('Training...')
t1 = tic;
[res, svrobj] = train(svrobj, dataobj);
w = get_w(svrobj);
t2 = toc(t1);
fprintf('Done in %3.2f sec\n', t2)

% Testing
fprintf('Testing...')
res2 = test(svrobj, data('test data', xtest, []));
yfit = res2.X; % this is proportional to res.X*w (perfectly correlated), but scale is different
fprintf('\n');
% vec = [res.X res.Y dataobj.X * w' res2.X]
% corrcoef(vec)

w = w';
dummy = 0;

end % function

function indx = searchlight_sphere_prep(xyz, i, r)
seed = xyz(i,:);
indx = sum([xyz(:,1)-seed(1) xyz(:,2)-seed(2) xyz(:,3)-seed(3)].^2, 2) <= r.^2;
end


