function [results_obj, stats, indx] = searchlight(dat, varargin)
% Run searchlight multivariate prediction/classification on an image_vector
% or fmri_data object OR two objects, for cross-prediction.
%
% :Usage:
% ::
%
%    [list outputs here] = function_name(list inputs here, [optional inputs])
%    [results_obj, stats, indx] = searchlight(dat, [optional inputs])
%
%
% :Features:
%  - Runs searchlight with standard, pre-defined algorithms
%  - Custom-entry definition of holdout sets
%  - Can re-use searchlight spheres after initial definition
%  - Custom-entry definition of any spheres/regions of interest
%  - Uses Matlab's parallel processing toolbox (parfor)
%
% Type help image_vector.searchlight to display this help information
%
% ..
%     Author and copyright information:
%
%     Copyright (C) 2014  Tor Wager and...
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
%
% :Inputs:
%
%   **dat:**
%        image_vector or fmri_data object with data
%
%   **dat.Y:**
%        required: true outcomes for each observation (image) in dat
%
% :Optional Inputs:* Keyword followed by input variable:
%
%   **r:**
%        searchlight radius, voxels
%   **dat2:**
%        second dataset, for cross-prediction
%   **indx:**
%        sparse logical matrix. each COLUMN is index of inclusion sets for each region/sphere in searchlight
%        This takes a long time to calculate, but can be saved and
%        re-used for a given mask
%   **holdout_set:**
%        Followed by integer vector of which observations belong to which
%        holdout set, for cross-validation. This is passed into fmri_data.predict.m.  Default is
%        empty.
%
% :Outputs:
%
%   **results_obj:**
%        fmri_data object with results maps
%
%   **stats:**
%        selected statistics for each sphere in searchlight
%
%   **indx:**
%        sparse logical matrix. each COLUMN is index of inclusion sets for each region/sphere in searchlight
%        * this can be re-used for all data with the same mask/structure. *
%
%
% :Examples:
% ::
%
%    % Define a sensible gray-matter mask:
%    dat = fmri_data(which('scalped_avg152T1_graymatter.img'));
%    dat = threshold(dat, [.8 Inf], 'raw-between');
%    dat = trim_mask(dat);
%
%    % Create fake data and holdout indicator index vector
%    dat.dat = randn(dat.volInfo.n_inmask, 30);
%    dat.Y = dat.dat(111111, :)' + .3 * randn(30, 1);
%    holdout_set = ones(6, 1); for i = 2:5, holdout_set = [holdout_set; i*ones(6, 1)]; end
%
%    % Run, and run again with existing indx
%    pool = parpool(12);  % initialize parallel processing (12 cores)
%    [results_obj, stats, indx] = searchlight(dat, 'holdout_set', holdout_set);
%    results_obj = searchlight(dat, 'holdout_set', holdout_set, 'indx', indx);
%
% :See also:
% region.m, fmri_data.predict.m
%
% ..
%    Programmers' notes:
%    List dates and changes here, and author of changes
% ..

% ..
% DEFAULTS AND INPUTS
% ..

r = 5; % For defining regions
indx = [];      % cell; with logical index of inclusion sets for each region/sphere in searchlight

% For prediction/classification
algorithm_name = 'cv_lassopcr';
holdout_set = [];
results_obj = [];

% process variable arguments
for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            % functional commands
            case {'r', 'indx', 'algorithm_name', 'dat2', 'holdout_set'}
                str = [varargin{i} ' = varargin{i + 1};'];
                eval(str)
                varargin{i} = [];
                varargin{i + 1} = [];
                
            case {'cross_predict'}
                algorithm_name = varargin{i};
                
                
            otherwise
                % use all other inputs as optional arguments to specific
                % functions, to be interpreted by them
                predfun_inputs{end + 1} = varargin{i};
                
                if (i+1) <= length(varargin) && ~ischar(varargin{i + 1})
                    predfun_inputs{end + 1} = varargin{i + 1};
                end
                
        end
    end
end

n = dat.volInfo.n_inmask;

%% Set up indices for spherical searchlight, if not entered previously

if isempty(indx)
    
    indx = searchlight_sphere_prep(dat, r);
    
else
    
    fprintf('Using input indx to define regions/spheres...\n');
    
end


%% Check for data and return if empty

if ~isa(dat, 'fmri_data') || isempty(dat.Y)
    
    fprintf('Returning indx only: No data in dat.Y to predict or data is not an fmri_data object');
    
    return
    
end


%% Run cross-predict (or other) function

% Get rough time estimate first
predict_time_estimate(dat, indx, algorithm_name, holdout_set);

t = tic;
fprintf('Running prediction in each region...');

% parfor i = 1:n
%
%     indxcell{i} = indx(:, i);
%
% end

output_val = cell(1, n);

parfor i = 1:n
    
    output_val{i} = predict_wrapper(dat, indx(:, i), algorithm_name, holdout_set);
    
end

e = toc(t);

[hour, minute, second] = sec2hms(e);
fprintf(1,'Done in %3.0f hours %3.0f min %2.0f sec\n',hour, minute, second);


%% Save results map(s) in object

results_obj = mean(dat);
results_obj.dat = cat(1, output_val.cverr{:});


end % function



% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%
% Sub-functions
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


% SPHERE PREP
% -------------------------------------------------------------------------

function indx = searchlight_sphere_prep(dat, r)

n = dat.volInfo.n_inmask;
indx = cell(1, n);

t = tic;
fprintf('Setting up seeds...');

parfor i = 1:n
    seed{i} = dat.volInfo.xyzlist(i, :);
end

e = toc(t);
fprintf('Done in %3.2f sec\n', e);

% Set up indices for spherical searchlight
% -------------------------------------------------------------------------
% These could be indices for ROIs, user input, previously saved indices...

% First, a rough time estimate:
% -------------------------------------------------------------------------
fprintf('Searchlight sphere construction can take 20 mins or more! (est: 20 mins with 8 processors/gray matter mask)\n');
fprintf('It can be re-used once created for multiple analyses with the same region definitions\n');
fprintf('Getting a rough time estimate for how long this will take...\n');

n_to_run = min(500, n);
t = tic;
parfor i = 1:n_to_run
    
    mydist = sum([dat.volInfo.xyzlist(:, 1) - seed{i}(1) dat.volInfo.xyzlist(:, 2) - seed{i}(2) dat.volInfo.xyzlist(:, 3) - seed{i}(3)] .^ 2, 2);
    indx{i} = mydist <= r.^2;
    
end
e = toc(t);
estim = e * n / n_to_run;

[hour, minute, second] = sec2hms(estim);
fprintf(1,'\nEstimate for whole brain = %3.0f hours %3.0f min %2.0f sec\n',hour, minute, second);

% Second, do it for all voxels/spheres:
% -------------------------------------------------------------------------

t = tic;

fprintf('Constructing spheres for each seed...');

%infdist = inf * ones(size(dat.volInfo.xyzlist(:, 1), 1), 1);

parfor i = 1:n
    
    mydist = sum([dat.volInfo.xyzlist(:, 1) - seed{i}(1) dat.volInfo.xyzlist(:, 2) - seed{i}(2) dat.volInfo.xyzlist(:, 3) - seed{i}(3)] .^ 2, 2);
    indx{i} = mydist <= r.^2;
    
end

e = toc(t);
fprintf('Done in %3.2f sec\n', e);

% Make sparse matrix
t = tic;
fprintf('Sparsifying region indices...');

indxcat = cat(2, indx{:});
indxcat = sparse(indxcat);

e = toc(t);
fprintf('Done in %3.2f sec\n', e);


% NOTE: IdenTICAL but twice as slow
% indx1 = cell(1, n);
% tic
% parfor i = 1:1000
%     indx1{i} = logical(iimg_xyz2spheres(seed{i}, dat.volInfo.xyzlist, r));
% end
% toc


end % searchlight_sphere_prep




% PREDICTION
% -------------------------------------------------------------------------

function output_val = predict_wrapper(dat, indx, algorithm_name, holdout_set)
% This function wraps predict.m to use its functionality.
% It could be extended to handle optional inputs, etc.
% Right now it is basic.
% Note: pass in indx{i}, which becomes indx here, for ONE
% region/searchlight sphere

dat.dat = dat.dat(indx, :);
dat.removed_voxels(~indx) = true;

[cverr, stats] = predict(dat, 'algorithm_name', algorithm_name, 'nfolds', holdout_set, 'useparallel', 0, 'verbose', 0);

% Parse output
% Only storing important info and full data weight map - we will probably
% want xval weight maps too at some point. We will also probably want a
% p-value and standard errors for accuracy too.

output_val = struct;
switch algorithm_name
    case 'cv_svm'
        output_val.cverr = stats.cverr;
        output_val.dist_from_hyperplane_xval = stats.dist_from_hyperplane_xval;
        output_val.weight_obj = stats.weight_obj;
        output_val.intercept = stats.other_output{3};
        
    case 'cv_lassopcr'
        output_val.yfit = stats.yfit;
        output_val.pred_outcome_r = stats.pred_outcome_r;
        output_val.weight_obj = stats.weight_obj;
        output_val.intercept = stats.other_output{3};
end

end % subfunction



% PREDICTION time estimate
% -------------------------------------------------------------------------

function predict_time_estimate(dat, indx, algorithm_name, holdout_set)

t = tic;
fprintf('Getting rough time estimate: Running prediction for up to 500 voxels...');

n = dat.volInfo.n_inmask;
output_val = cell(1, n);

n_to_run = min(500, n);

parfor i = 1:n_to_run
    
    output_val{i} = predict_wrapper(dat, indx(:, i), algorithm_name, holdout_set);
    
end

e = toc(t);
fprintf('Done in %3.2f sec\n', e);

estim = e * n / n_to_run;

[hour, minute, second] = sec2hms(estim);
fprintf(1,'Estimate for whole brain = %3.0f hours %3.0f min %2.0f sec\n',hour, minute, second);

end




% Convert time
% -------------------------------------------------------------------------


function [hour, minute, second] = sec2hms(sec)
%SEC2HMS  Convert seconds to hours, minutes and seconds.
%
%   [HOUR, MINUTE, SECOND] = SEC2HMS(SEC) converts the number of seconds in
%   SEC into hours, minutes and seconds.

hour   = fix(sec/3600);      % get number of hours
sec    = sec - 3600*hour;    % remove the hours
minute = fix(sec/60);        % get number of minutes
sec    = sec - 60*minute;    % remove the minutes
second = sec;
end

% % WANI'S CROSS-CLASSIFICATION FUNCTION - to edit and integrate
% % -------------------------------------------------------------------------
% 
% function [acc, p, se, stats1 stats2] = cv_svm_cross_subf(data1, data2, cv_assign, dobalanced, balanced_ridge, doscale)
% 
% % run cv_svm on train data and test data
% 
% if doscale
%     data1 = rescale(data1, 'zscoreimages');
%     data2 = rescale(data2, 'zscoreimages');
% end
% 
% % preallocate some variables
% predicted11 = NaN(size(data1.Y));
% predicted12 = NaN(size(data2.Y));
% predicted21 = NaN(size(data1.Y));
% predicted22 = NaN(size(data2.Y));
% 
% dist_from_hyper11 = NaN(size(data1.Y));
% dist_from_hyper12 = NaN(size(data2.Y));
% dist_from_hyper21 = NaN(size(data1.Y));
% dist_from_hyper22 = NaN(size(data2.Y));
% 
% % get cv assignment
% u = unique(cv_assign);
% nfolds = length(u);
% [trIdx, teIdx] = deal(cell(1, nfolds));
% 
% for i = 1:length(u)
%     teIdx{i} = cv_assign == u(i);
%     trIdx{i} = ~teIdx{i};
% end
% 
% stats1.Y = data1.Y;
% stats2.Y = data2.Y;
% 
% % Use all
% trainobj1 = data(double(data1.dat'),data1.Y);
% trainobj2 = data(double(data2.dat'),data2.Y);
% 
% svmobj1 = svm({'optimizer="andre"','C=1','child=kernel'});
% if dobalanced, svmobj1.balanced_ridge = balanced_ridge; end
% [~, svmobj1] = train(svmobj1, trainobj1);
% 
% svmobj2 = svm({'optimizer="andre"','C=1','child=kernel'});
% if dobalanced, svmobj2.balanced_ridge = balanced_ridge; end
% [~, svmobj2] = train(svmobj2, trainobj2);
% 
% % 1) weights
% stats1.all{1} = get_w(svmobj1)';
% stats2.all{1} = get_w(svmobj2)';
% % 2) intercepts
% stats1.all{2} = svmobj1.b0;
% stats2.all{2} = svmobj2.b0;
% 
% stats1.all_descrip = {'1:weights 2:intercept 3:distance from hyperplane for test data1 4: dist for data2'};
% stats2.all_descrip = {'1:weights 2:intercept 3:distance from hyperplane for test data2 4: dist for data1'};
% 
% % Cross-val loop starts here
% 
% for i = 1:numel(teIdx)
%     
%     % Prepare the hot-warm and Rej-friend data in Spider format
%     trainobj1 = data(double(data1.dat'),data1.Y);
%     trainobj2 = data(double(data2.dat'),data2.Y);
%     
%     % Select training and test data - trainobj, testobj
%     testobj1 = trainobj1;
%     testobj2 = trainobj2;
%     
%     trainobj1.X = trainobj1.X(trIdx{i}, :);
%     trainobj1.Y = trainobj1.Y(trIdx{i}, :);
%     
%     trainobj2.X = trainobj2.X(trIdx{i}, :);
%     trainobj2.Y = trainobj2.Y(trIdx{i}, :);
%     
%     testobj1.X = testobj1.X(teIdx{i}, :);
%     testobj1.Y = testobj1.Y(teIdx{i}, :);
%     
%     testobj2.X = testobj2.X(teIdx{i}, :);
%     testobj2.Y = testobj2.Y(teIdx{i}, :);
%     
%     
%     %-------------------------------------------------------------------
%     % Train on data1 and test on data1 and data2
%     
%     % Set up an SVM object
%     svmobj1 = svm({'optimizer="andre"','C=1','child=kernel'});
%     
%     if dobalanced
%         svmobj1.balanced_ridge = balanced_ridge;
%     end
%     
%     % train on hw data set
%     [~, svmobj1] = train(svmobj1, trainobj1);
%     
%     % Test on both HW and RF test set
%     
%     testobj11 = test(svmobj1, testobj1);
%     testobj12 = test(svmobj1, testobj2);
%     
%     w = get_w(svmobj1)';
%     b0 = svmobj1.b0;
%     
%     dist11 = testobj1.X * w + b0;
%     dist12 = testobj2.X * w + b0;
%     
%     if ~isequal((dist11>0)*2-1,testobj11.X)
%         error('something is wrong');
%     elseif ~isequal((dist12>0)*2-1,testobj12.X)
%         error('something is wrong');
%     end
%     
%     % Save predictions (cross-val)
%     predicted11(teIdx{i}) = testobj11.X;
%     predicted12(teIdx{i}) = testobj12.X;
%     
%     
%     %-------------------------------------------------------------------
%     % Train on RF and test on HW and RF
%     
%     % Set up an SVM object
%     svmobj2 = svm({'optimizer="andre"','C=1','child=kernel'});
%     
%     if dobalanced
%         svmobj2.balanced_ridge = balanced_ridge;
%     end
%     % train on rf data set
%     [~, svmobj2] = train(svmobj2, trainobj2);
%     
%     % Test on both HW and RF test set
%     
%     testobj21 = test(svmobj2, testobj1);
%     testobj22 = test(svmobj2, testobj2);
%     
%     w = get_w(svmobj2)';
%     b0 = svmobj2.b0;
%     
%     dist21 = testobj1.X * w + b0;
%     dist22 = testobj2.X * w + b0;
%     
%     if ~isequal((dist21>0)*2-1,testobj21.X)
%         error('something is wrong');
%     elseif ~isequal((dist22>0)*2-1,testobj22.X)
%         error('something is wrong');
%     end
%     
%     % Save predictions (cross-val)
%     predicted21(teIdx{i}) = testobj21.X;
%     predicted22(teIdx{i}) = testobj22.X;
%     
%     % get stats from svmobj
%     % 1) weights
%     stats1.cvoutput{i,1} = get_w(svmobj1)';
%     stats2.cvoutput{i,1} = get_w(svmobj2)';
%     % 2) intercepts
%     stats1.cvoutput{i,2} = svmobj1.b0;
%     stats2.cvoutput{i,2} = svmobj2.b0;
%     % 3) distance from hyperplane for test data 1 and 2
%     stats1.cvoutput{i,3} = testobj1.X * stats1.cvoutput{i,1} + stats1.cvoutput{i,2};
%     stats1.cvoutput{i,4} = testobj2.X * stats1.cvoutput{i,1} + stats1.cvoutput{i,2};
%     
%     stats2.cvoutput{i,3} = testobj1.X * stats2.cvoutput{i,1} + stats2.cvoutput{i,2};
%     stats2.cvoutput{i,4} = testobj2.X * stats2.cvoutput{i,1} + stats2.cvoutput{i,2};
%     
%     stats1.cvoutput_descrip = {'1:weights 2:intercept 3:distance from hyperplane for test data1 4: dist for data2'};
%     stats2.cvoutput_descrip = {'1:weights 2:intercept 3:distance from hyperplane for test data2 4: dist for data1'};
%     
%     dist_from_hyper11(teIdx{i}) = stats1.cvoutput{i,3};
%     dist_from_hyper12(teIdx{i}) = stats1.cvoutput{i,4};
%     dist_from_hyper21(teIdx{i}) = stats2.cvoutput{i,3};
%     dist_from_hyper22(teIdx{i}) = stats2.cvoutput{i,4};
% end
% 
% stats1.yfit = predicted11;
% stats2.yfit = predicted22;
% stats1.testfit = predicted12;
% stats2.testfit = predicted21;
% 
% stats1.all{3} = dist_from_hyper11;
% stats1.all{4} = dist_from_hyper12;
% stats2.all{3} = dist_from_hyper22;
% stats2.all{4} = dist_from_hyper21;
% 
% % RO on H vs. W
% % res_temp = binotest(predicted11(1:(length(predicted11)/2)) ==
% %data1.Y(1:(length(predicted11)/2)), 0.5);
% % test_results.accuracy(1,1) = res_temp.prop;
% % test_results.accuracy_se(1,1) = res_temp.SE;
% % test_results.accuracy_p(1,1) = res_temp.p_val;
% %
% % res_temp = binotest(predicted12(1:(length(predicted12)/2)) ==
% %data2.Y(1:(length(predicted12)/2)), 0.5);
% % test_results.accuracy(1,2) = res_temp.prop;
% % test_results.accuracy_se(1,2) = res_temp.SE;
% % test_results.accuracy_p(1,2) = res_temp.p_val;
% %
% % res_temp = binotest(predicted21(1:(length(predicted21)/2)) ==
% %data1.Y(1:(length(predicted21)/2)), 0.5);
% % test_results.accuracy(1,3) = res_temp.prop;
% % test_results.accuracy_se(1,3) = res_temp.SE;
% % test_results.accuracy_p(1,3) = res_temp.p_val;
% %
% % res_temp = binotest(predicted22(1:(length(predicted22)/2)) ==
% %data2.Y(1:(length(predicted22)/2)), 0.5);
% % test_results.accuracy(1,4) = res_temp.prop;
% % test_results.accuracy_se(1,4) = res_temp.SE;
% % test_results.accuracy_p(1,4) = res_temp.p_val;
% 
% subjn = numel(teIdx);
% outcome = [true(subjn,1); false(subjn,1)];
% 
% hvsw_dist11 = stats1.all{3}(1:(subjn*2));
% rvsf_dist12 = stats1.all{4}(1:(subjn*2));
% hvsw_dist21 = stats2.all{4}(1:(subjn*2));
% rvsf_dist22 = stats2.all{3}(1:(subjn*2));
% 
% 
% % ROC = roc_plot(hvsw_dist11, outcome, 'threshold', 0);
% % test_results.accuracy_thresh0(1) = ROC.accuracy;
% % test_results.accuracy_thresh0_se(1) = ROC.accuracy_se;
% % test_results.accuracy_thresh0_p(1) = ROC.accuracy_p;
% 
% try
%     ROC = roc_plot(hvsw_dist11, outcome, 'twochoice');
%     acc(1) = ROC.accuracy;
%     se(1) = ROC.accuracy_se;
%     p(1) = ROC.accuracy_p;
% catch err
%     % ROC = roc_plot(hvsw_dist11, outcome, 'twochoice');
%     acc(1) = NaN;
%     se(1) = NaN;
%     p(1) = NaN;
% end
% 
% %
% % ROC = roc_plot_wani(rvsf_dist12, outcome, 'threshold', 0);
% % test_results.accuracy_thresh0(2) = ROC.accuracy;
% % test_results.accuracy_thresh0_se(2) = ROC.accuracy_se;
% % test_results.accuracy_thresh0_p(2) = ROC.accuracy_p;
% %
% try
%     ROC = roc_plot(rvsf_dist12, outcome, 'twochoice');
%     acc(2) = ROC.accuracy;
%     se(2) = ROC.accuracy_se;
%     p(2) = ROC.accuracy_p;
% catch err
%     acc(2) = NaN;
%     se(2) = NaN;
%     p(2) = NaN;
% end
% %
% %
% % ROC = roc_plot(hvsw_dist21, outcome, 'threshold', 0);
% % test_results.accuracy_thresh0(3) = ROC.accuracy;
% % test_results.accuracy_thresh0_se(3) = ROC.accuracy_se;
% % test_results.accuracy_thresh0_p(3) = ROC.accuracy_p;
% %
% try
%     ROC = roc_plot(hvsw_dist21, outcome, 'twochoice');
%     acc(3) = ROC.accuracy;
%     se(3) = ROC.accuracy_se;
%     p(3) = ROC.accuracy_p;
% catch err
%     acc(3) = NaN;
%     se(3) = NaN;
%     p(3) = NaN;
% end
% 
% % ROC = roc_plot(rvsf_dist22, outcome, 'threshold', 0);
% % test_results.accuracy_thresh0(4) = ROC.accuracy;
% % test_results.accuracy_thresh0_se(4) = ROC.accuracy_se;
% % test_results.accuracy_thresh0_p(4) = ROC.accuracy_p;
% try
%     ROC = roc_plot(rvsf_dist22, outcome, 'twochoice');
%     acc(4) = ROC.accuracy;
%     se(4) = ROC.accuracy_se;
%     p(4) = ROC.accuracy_p;
% catch err
%     acc(4) = NaN;
%     se(4) = NaN;
%     p(4) = NaN;
% end
% 
% close;

% end
