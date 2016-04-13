function out = searchlight_disti(dat, mask, dist_i, additional_inputs)
% Run the actual searchlight analysis on each brain chunck 
% searchlight_dream.m will generate codes to run this funtion.
%
% :Usage:
% ::
%
%     out = searchlight_disti(dat, mask, dist_i, [additional_inputs])
%
% ..
%     Author and copyright information:
%
%     Copyright (C) 2014  Wani Woo, Tor Wager
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
% :See Also:
% searchlight_dream.m, xval_cross_classify.m, fmri_data.predict.m, 


do_cross = false;
do_saveweight = false;
r = 3; % default radius (in voxel)

% parsing varargin (additional_inputs)

for i = 1:length(additional_inputs)
    if ischar(additional_inputs{i})
        switch additional_inputs{i}
            % functional commands
            case 'dat2'
                do_cross = true;
                dat2 = additional_inputs{i+1};
                additional_inputs{i} = []; additional_inputs{i+1} = [];
            case 'algorithm_name'
                algorithm_name = additional_inputs{i+1};
                additional_inputs{i} = []; additional_inputs{i+1} = [];
            case 'r' % radius
                r = additional_inputs{i+1};
                additional_inputs{i} = []; additional_inputs{i+1} = [];
            case 'cv_assign'
                cv_assign = additional_inputs{i+1};
                additional_inputs{i} = []; additional_inputs{i+1} = [];
            case 'save_weights'
                do_saveweight = true;
                additional_inputs{i} = [];
        end
    end
end

% prep dat2 if do_cross

if do_cross
    [dat2, mask] = apply_mask(dat2, mask);
    dat2 = trim_mask(dat2);
end

%% searchlight prep

% data check
if do_cross, if any(sum(sum(dat.volInfo.xyzlist ~= dat2.volInfo.xyzlist))), error('dat and dat2 should be in the same space.'); end, end

% get voxel index to run
vox_to_run = find(dist_i)';

index_within_this_run = 0;
for i = vox_to_run %(1):vox_to_run(10)
    searchlight_indx = searchlight_sphere_prep(dat.volInfo.xyzlist, i, r);
    
    % data prep
    
    if do_cross
        results = xval_cross_classfy_wrapper(algorithm_name, dat, dat2, cv_assign, searchlight_indx, do_saveweight, additional_inputs);
    else
        results = predict_wrapper(algorithm_name, dat, cv_assign, searchlight_indx, do_saveweight, additional_inputs);
    end
    
    % output prep
    if ~exist('out', 'var')
        for jj = 1:numel(results{1})
            f = fields(results{1}{jj});
            for ii = 1:numel(f)
                eval(['out.test_results{jj}.' f{ii} ' = NaN(dat.volInfo.n_inmask, size(results{1}{jj}.' f{ii} ',2));']);
            end
        end
    end
    
    index_within_this_run = index_within_this_run+1;
    out = parse_results(out, results, i, index_within_this_run);
    
end

end


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%
% Sub-functions
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% searchlight_sphere_prep

function indx = searchlight_sphere_prep(xyz, i, r)
seed = xyz(i,:);
indx = sum([xyz(:,1)-seed(1) xyz(:,2)-seed(2) xyz(:,3)-seed(3)].^2, 2) <= r.^2;
end

% -------------------------------------------------------------------------
% xval_cross_classification 

function results = xval_cross_classfy_wrapper(algorithm_name, dat1, dat2, cv_assign, searchlight_indx, do_saveweight, additional_inputs)

whempty = cellfun(@isempty, additional_inputs);
additional_inputs(whempty)=[];

for i = 1:2
    eval(['dat' num2str(i) '.dat = dat' num2str(i) '.dat(searchlight_indx,:);']);
    eval(['dat' num2str(i) '.removed_voxels(~searchlight_indx) = true;']);
end

if do_saveweight
    [results{1}, results{2}, results{3}] = xval_cross_classfy(algorithm_name, dat1, dat2, cv_assign, additional_inputs{:});
else
    results{1} = xval_cross_classfy(algorithm_name, dat1, dat2, cv_assign, additional_inputs{:});
end

end

% -------------------------------------------------------------------------
% predict wrapper

function results = predict_wrapper(algorithm_name, dat, cv_assign, searchlight_indx, do_saveweight, additional_inputs)

whempty = cellfun(@isempty, additional_inputs);
additional_inputs(whempty)=[];

[out_method, test_Y, dat] = setup_testvar(dat, additional_inputs{:});

dat.dat = dat.dat(searchlight_indx, :);
dat.removed_voxels(~searchlight_indx) = true;

[dummy, stats] = predict(dat, 'algorithm_name', algorithm_name, 'nfolds', cv_assign, 'verbose', 0);

if do_saveweight
    stats = rmfield(stats, {'weight_obj', 'teIdx', 'trIdx'});
    results{1} = get_test_results(stats, test_Y, out_method, algorithm_name);
    results{2} = stats;
else
    stats = rmfield(stats, {'weight_obj', 'teIdx', 'trIdx'});
    results{1} = get_test_results(stats, test_Y, out_method, algorithm_name);
end

end

% -------------------------------------------------------------------------
% parse test_results

function out = parse_results(out, results, i, index_within_this_run)

for j = 1:numel(results{1})
    
    f = fields(results{1}{j});
    
    for ii = 1:numel(f)
        eval(['out.test_results{j}.' f{ii} '(i,:) = results{1}{j}.' f{ii} ';']);
    end
end

if numel(results) == 2
    out.stats{index_within_this_run} = results{2};
    out.stats{index_within_this_run}.vox_indx = i;
elseif numel(results) == 3
    for ii = 1:2
        eval(['out.stats' num2str(ii) '{index_within_this_run} = results{' num2str(ii+1) '};']);
        eval(['out.stats' num2str(ii) '{index_within_this_run}.vox_indx = i;']);
    end
end

end


% -------------------------------------------------------------------------
% set up test variables

function [out_method, test_Y, dat] = setup_testvar(dat, varargin)

% default for binary
out_method_input = 'twochoice';

% get outcome method input
for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            case {'outcome_method'}
                out_method_input = varargin{i+1};
        end
    end
end

% if dat.Y is the test variable
if size(dat.Y,2) == 1
    
    if numel(unique(dat.Y(dat.Y~=0)))==2 
        out_method{1} = out_method_input;
        test_Y{1} = dat.Y; 
        
    elseif numel(unique(dat.Y(dat.Y~=0))) > 2 
        out_method{1} = 'correlation';
        test_Y{1} = dat.Y; 
    end
    
% if data.Y has more than one column
elseif size(dat.Y,2) > 1
    
    for i = 2:size(dat.Y,2)

        if numel(unique(dat.Y(dat.Y(:,i)~=0,i)))==2 
            out_method{i-1} = out_method_input;
            test_Y{i-1} = dat.Y(:,i); 
        elseif numel(unique(dat.Y(dat.Y(:,i)~=0,i))) > 2 
            out_method{i-1} = 'correlation';
            test_Y{i-1} = dat.Y(:,i); 
        end
        
    end
end
    
% delete test variables from data
dat.Y = dat.Y(:,1);

end

function test_results = get_test_results(stats, test_Y, out_method,algorithm_name)

for i = 1:numel(out_method)
    
    whempty = test_Y{i} == 0;
    
    switch out_method{i}
        case 'correlation'
            [test_results{i}.r, test_results{i}.p] = corr(stats.yfit(~whempty), test_Y{i}(~whempty));
        case 'singleinterval'
            try
                if ~isempty(strfind(algorithm_name, 'svm')) 
                    ROC = roc_plot(stats.dist_from_hyperplane_xval, test_Y{i}==1, 'include', ~whempty, 'threshold', 0, 'balanced');
                elseif ~isempty(strfind(algorithm_name, 'lassopcr'))
                    roc_temp = roc_plot(stats.Y, test_Y{i}==1, 'include', ~whempty, 'balanced');
                    thr = roc_temp.class_threshold;
                    ROC = roc_plot(stats.yfit, test_Y{i}==1, 'include', ~whempty, 'threshold', thr, 'balanced');
                end
                test_results{i}.acc = ROC.accuracy;
                test_results{i}.se = ROC.accuracy_se;
                test_results{i}.p = ROC.accuracy_p;
                test_restuls{i}.thr = ROC.class_threshold;
            catch dummy
                test_results{i}.acc = NaN;
                test_results{i}.se = NaN;
                test_results{i}.p = NaN;
                test_restuls{i}.thr = NaN;
            end
            
            test_results{i}.p(test_results{i}.p > 1) = 1;
            
        case 'twochoice'
            try
                if ~isempty(strfind(algorithm_name, 'svm')) 
                    ROC = roc_plot(stats.dist_from_hyperplane_xval, test_Y{i}==1, 'twochoice', 'include', ~whempty);
                elseif ~isempty(strfind(algorithm_name, 'lassopcr'))
                    ROC = roc_plot(stats.yfit, test_Y{i}==1, 'twochoice', 'include', ~whempty);
                end
                test_results{i}.acc = ROC.accuracy;
                test_results{i}.se = ROC.accuracy_se;
                test_results{i}.p = ROC.accuracy_p;
            catch dummy
                test_results{i}.acc = NaN;
                test_results{i}.se = NaN;
                test_results{i}.p = NaN;
            end
            
            test_results{i}.p(test_results{i}.p > 1) = 1;
    end
end

close all;

end
