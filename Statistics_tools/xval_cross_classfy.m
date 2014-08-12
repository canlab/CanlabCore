function [test_results, stats1, stats2] = xval_cross_classfy(algorithm_name, data1, data2, cv_assign, varargin)

% Run SVM or LASSOPCR on two image_vector or fmri_data objects for cross-classification.
%
% Usage:
% -------------------------------------------------------------------------
% [test_results, stats1, stats2] = xval_cross_classfy(algorithm_name, dat1, dat2, cv_assign, [optional inputs])
%
% Features:
% - ...
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
% Inputs:
% -------------------------------------------------------------------------
% algorithm_name: 'cv_svm' (linear svm) or 'cv_lassopcr' 
% -------------------------------------------------------------------------
% Data:
%
% dat1          image_vector or fmri_data object with data
% dat1.Y(:,1)   for cv_svm: true(1) or false(-1) for each observation (image) in Y(:,1)
%               for cv_lassopcr: continuous value for Y(:,1)
%
% dat2          image_vector or fmri_data object with data
% dat2.Y        for cv_svm: true(1) or false(-1) for each observation (image) in Y(:,1)
%               for cv_lassopcr: continuous value for Y(:,1),
%
% dat1.Y(:,[2:n]), dat2.Y(:,[2:n]) 
%               Test sets: could be binary: and true(1), false(-1),
%               ignore(0) or continuous values
%
% cv_assign     vector of integers for membership in custom holdout set of each fold
%               If two datasets have different order of cross-validation
%               folds, you can put two different cv_assign vector in cell
%               arrays. e.g.) cv_assign{1} = whfolds1, cv_assign{2} = whfolds2
% 
% Optional inputs: 
% -------------------------------------------------------------------------
% 'scale'       z-scored input data in image_vector or fmri_data object
% 'balanced'    use the balanced ridge option - balanced ridge value should
%               be followed.
% 'outcome_method'  followed by the following options
%       'correlation' - "default" for for continuous measures
%       'twochoice'- "default" for binary outcome
%       'singleinterval'  - for binary outcome
% 
%
% Outputs:
% -------------------------------------------------------------------------
% test_results:     for binary classification - four values of accuracy, p, se      
%                   for continuous prediction - four values of r(pearson's 
%                    correlation), p-value. 
%                   The order of the test results are [dat1-on-dat1, dat1-on-dat2,
%                    dat2-on-dat1, dat2-on-dat2]. All results are cross-validated 
%                    results.
%
% stats1, stats2    stats1 and stats2 are similar to the outputs of predict
%                   function. help fmri_data.predict
%
% Examples:
% % Data preparation
% % -------------------------------------
% dat1 = fmri_data(which('scalped_avg152T1_graymatter.img'));
% dat1 = threshold(dat, [.8 Inf], 'raw-between');
% dat1 = trim_mask(dat);
% dat2 = dat1;
% 
% % Create fake data and holdout indicator index vector
% % -------------------------------------
% dat1.dat = randn(dat.volInfo.n_inmask, 30);
% dat2.dat = randn(dat.volInfo.n_inmask, 30);
% 
% % for svm
% % dat1.Y = ((dat.dat(111111, :)' + .3 * randn(30, 1))>0).*2-1;
% % dat2.Y = ((dat.dat(111, :)' + .3 * randn(30, 1))>0).*2-1;
% % cv_assign = ones(6, 1); for i = 2:5, cv_assign = [cv_assign; i*ones(6, 1)]; end
% 
% % for cv_lassopcr
% dat1.Y = dat.dat(111111, :)' + .3 * randn(30, 1);
% dat2.Y = dat.dat(111, :)' + .3 * randn(30, 1);
% dat1.Y(:,2) = [ones(10,1); zeros(10,1); -ones(10,1)];
% dat2.Y(:,2) = [-ones(10,1); zeros(10,1); ones(10,1)];
% cv_assign = ones(6, 1); for i = 2:5, cv_assign = [cv_assign; i*ones(6, 1)]; end
% 
% % Run, and run again with existing indx
% % -------------------------------------
% % [test_results, stats1, stats2] = xval_cross_classfy('cv_svm', dat1, dat2, cv_assign, 'outcome_method', 'singleinterval')
% [test_results, stats1, stats2] = xval_cross_classfy('cv_lassopcr', dat1, dat2, cv_assign, 'outcome_method', 'singleinterval')
%
% See also:
% fmri_data.predict.m

% Programmers' notes:
% 


% ----- defaults and input ----- 

doscale = false;
dobalanced = false; balanced_ridge = [];

% ----- organize test variables ----- 
[out_method, test_Y1, test_Y2, data1, data2] = setup_testvar(data1, data2, varargin);

% ----- optional inputs ----- 
for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            % functional commands
            case {'scale'}
                doscale = true;
            case {'balanced'}
                dobalanced = true;
                balanced_ridge = varargin{i+1};
        end
    end
end

% ----- z-score data if doscale ----- 
if doscale
    data1 = rescale(data1, 'zscoreimages');
    data2 = rescale(data2, 'zscoreimages');
end

% ----- run algorithms ----- 

switch algorithm_name    
    case 'cv_svm'
        [test_results, stats1, stats2] = cv_svm_cross_classfy(data1, data2, cv_assign, test_Y1, test_Y2, out_method, dobalanced, balanced_ridge);
    case 'cv_lassopcr'        
        [test_results, stats1, stats2] = cv_lassopcr_cross_classfy(data1, data2, cv_assign, test_Y1, test_Y2, out_method);
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
% CV-ASSIGNMENT

function [teIdx, trIdx] = do_cv_assign(cv_assign)

% assign cross-validation
if ~iscell(cv_assign) 
    for j = 1:size(cv_assign,2)
        temp_cv{j} = cv_assign(:,j);
    end 
    cv_assign = temp_cv;
end

if numel(cv_assign) == 1, cv_assign{2} = cv_assign{1}; end

for j = 1:numel(cv_assign)
    u = unique(cv_assign{j});
    nfolds = length(u);
    [trIdx{j}, teIdx{j}] = deal(cell(1, nfolds));

    for i = 1:length(u)
        teIdx{j}{i} = cv_assign{j} == u(i);
        trIdx{j}{i} = ~teIdx{j}{i};
    end
end

end

% -------------------------------------------------------------------------
% CHECK SVM ERRORS, IF ANY

function test_error(testobj1, testobj2, testobj3, testobj4, w, b0)

% briefly check everything is right
test_dist1 = testobj1.X * w + b0;
test_dist2 = testobj2.X * w + b0;

if ~isequal((test_dist1>0)*2-1,testobj3.X)
    error('Something is wrong.');
elseif ~isequal((test_dist2>0)*2-1,testobj4.X)
    error('Something is wrong.');
end

end

% -------------------------------------------------------------------------
% set up test variables

function [out_method, test_Y1, test_Y2, data1, data2] = setup_testvar(data1, data2, varargin)

% default for binary
out_method_input = 'twochoice';

% get outcome method input
for i = 1:length(varargin{1})
    if ischar(varargin{1}{i})
        switch varargin{1}{i}
            case {'outcome_method'}
                out_method_input = varargin{1}{i+1};
        end
    end
end

% first, check the data
if size(data1.Y,2) ~= size(data2.Y,2)
    error('The column size of Y should be same between dat1 and dat2.')
    
% if data.Y is the test variable
elseif size(data1.Y,2) == 1
    
    if numel(unique(data1.Y(data1.Y~=0)))==2 && numel(unique(data2.Y(data2.Y~=0)))==2
        out_method{1} = out_method_input;
        test_Y1{1} = data1.Y; test_Y2{1} = data2.Y;
        
    elseif numel(unique(data1.Y(data1.Y~=0))) > 2 && numel(unique(data2.Y(data2.Y~=0))) > 2
        out_method{1} = 'correlation';
        test_Y1{1} = data1.Y; test_Y2{1} = data2.Y;
        
    else
        error('dat1.Y and dat2.Y are not consistent.');
    end
    
% if data.Y has more than one column
elseif size(data1.Y,2) > 1
    
    for i = 2:size(data1.Y,2)

        if numel(unique(data1.Y(data1.Y(:,i)~=0,i)))==2 && numel(unique(data2.Y(data2.Y(:,i)~=0,i)))==2
            out_method{i-1} = out_method_input;
            test_Y1{i-1} = data1.Y(:,i); test_Y2{i-1} = data2.Y(:,i);
        elseif numel(unique(data1.Y(data1.Y(:,i)~=0,i))) > 2 && numel(unique(data2.Y(data2.Y(:,i)~=0,i))) > 2
            out_method{i-1} = 'correlation';
            test_Y1{i-1} = data1.Y(:,i); test_Y2{i-1} = data2.Y(:,i);
        else
            error('dat1.Y and dat2.Y are not consistent.');
        end
        
    end
end
    
% delete test variables from data
data1.Y = data1.Y(:,1);
data2.Y = data2.Y(:,1);

end

% -------------------------------------------------------------------------
% get test results

function test_results = get_test_results(stats1, stats2, test_Y1, test_Y2, out_method)

out_data = {'stats1.all{3}', 'stats1.all{4}', 'stats2.all{4}', 'stats2.all{3}'};
out_test = {'test_Y1', 'test_Y2', 'test_Y1', 'test_Y2'};
out_whempty = {'whempty1', 'whempty2', 'whempty1', 'whempty2'};

for i = 1:numel(out_method)
    
    whempty1 = test_Y1{i} == 0;
    whempty2 = test_Y2{i} == 0;
    
    switch out_method{i}
        case 'correlation'
            
            for j = 1:4
                eval(['[test_results{i}.r(' num2str(j) '), test_results{i}.p(' num2str(j) ')] = corr(' out_data{j} '(~' out_whempty{j} '), ' out_test{j} '{i}(~' out_whempty{j} '));']);
            end
            
        case 'singleinterval'
            
            roc_temp1 = roc_plot(stats1.Y, test_Y1{i} == 1, 'include', ~whempty1, 'balanced');
            thr{1} = roc_temp1.class_threshold; thr{2} = thr{1};
            roc_temp2 = roc_plot(stats2.Y, test_Y2{i} == 1, 'include', ~whempty2, 'balanced');
            thr{3} = roc_temp2.class_threshold; thr{4} = thr{3};
            
            for j = 1:4
                try
                    eval(['ROC = roc_plot(' out_data{j} ', ' out_test{j} '{i}==1, ''include'', ~' out_whempty{j} ', ''balanced'', ''threshold'', thr{' num2str(j) '});']);
                    eval(['test_results{i}.acc(' num2str(j) ') = ROC.accuracy;']);
                    eval(['test_results{i}.se(' num2str(j) ') = ROC.accuracy_se;']);
                    eval(['test_results{i}.p(' num2str(j) ') = ROC.accuracy_p;']);
                    eval(['test_results{i}.thr(' num2str(j) ') = ROC.class_threshold;']);
                    
                catch dummy
                    eval(['test_results{i}.acc(' num2str(j) ') = NaN;']);
                    eval(['test_results{i}.se(' num2str(j) ') = NaN;']);
                    eval(['test_results{i}.p(' num2str(j) ') = NaN;']);
                    eval(['test_results{i}.thr(' num2str(j) ') = NaN;']);
                end
                
            end
            
            test_results{i}.p(test_results{i}.p > 1) = 1;
            
        case 'twochoice'
            
            for j = 1:4
                try
                    eval(['ROC = roc_plot(' out_data{j} ', ' out_test{j} '{i}==1, ''twochoice'', ''include'', ~' out_whempty{j} ');']);
                    eval(['test_results{i}.acc(' num2str(j) ') = ROC.accuracy;']);
                    eval(['test_results{i}.se(' num2str(j) ') = ROC.accuracy_se;']);
                    eval(['test_results{i}.p(' num2str(j) ') = ROC.accuracy_p;']);
                catch dummy
                    eval(['test_results{i}.acc(' num2str(j) ') = NaN;']);
                    eval(['test_results{i}.se(' num2str(j) ') = NaN;']);
                    eval(['test_results{i}.p(' num2str(j) ') = NaN;']);
                end
            end
            
            test_results{i}.p(test_results{i}.p > 1) = 1;
    end
end

close all;

end

% -------------------------------------------------------------------------
% CROSS-CLASSIFY USING CV-SVM

function [test_results, stats1, stats2] = cv_svm_cross_classfy(data1, data2, cv_assign, test_Y1, test_Y2, out_method, dobalanced, balanced_ridge)

% ----- activate the spider toolbox. ----- 
% %% this is slow. Please activate spider toolbox before running this function
% if isempty(which('use_spider.m'))
%     error('Error: You need to have the Spider toolbox in your path. Check your path.');
% else
%     use_spider;
% end

% ----- preallocate variables -----
for j = 1:2
    eval(['predicted1' num2str(j) ' = NaN(size(data' num2str(j) '.Y));']);
    eval(['predicted2' num2str(j) ' = NaN(size(data' num2str(j) '.Y));']);
    eval(['dist_from_hyper1' num2str(j) ' = NaN(size(data' num2str(j) '.Y));']);
    eval(['dist_from_hyper2' num2str(j) ' = NaN(size(data' num2str(j) '.Y));']);
end

% ----- get cv assignment ----- 
[teIdx, trIdx] = do_cv_assign(cv_assign);

%  ----- setup the variables for SVM ----- 
for j = 1:2, eval(['trainobj' num2str(j) ' = data(double(data' num2str(j) '.dat''),data' num2str(j) '.Y);']); end

% ----- run it with all data (not cv)----- 
for j = 1:2
    eval(['svmobj' num2str(j) ' = svm({''optimizer="andre"'',''C=1'',''child=kernel''});']);
    eval(['if dobalanced, svmobj' num2str(j) '.balanced_ridge = balanced_ridge; end']);
    eval(['[~, svmobj' num2str(j) '] = train(svmobj' num2str(j) ', trainobj' num2str(j) ');']);
end

% ----- outputs I: using all data 1 & 2 -----
for j = 1:2
    eval(['stats' num2str(j) '.Y = data' num2str(j) '.Y;']);
    % output: weights
    eval(['stats' num2str(j) '.all{1} = get_w(svmobj' num2str(j) ')'';']);
    % output: intercepts
    eval(['stats' num2str(j) '.all{2} = svmobj' num2str(j) '.b0;']);
    % description of out outputs: 3 and 4 have not been collected yet, though.
    eval(['stats' num2str(j) '.all_descrip = {''1:weights 2:intercept 3:distance from hyperplane for testing data' num2str(j) '(xval) 4: dist for the other data(xval)''};']);
end

% ----- Cross-validation loop starts here -----

if numel(teIdx)>1, if ~isequal(numel(teIdx{1}), numel(teIdx{2})), error('The number of folds should be same between dat1 and dat2.'); end, end 

for i = 1:numel(teIdx{1})
    
    % Prepare the data in Spider format
    for j = 1:2
        % prepare train data 1 & 2
        eval(['trainobj' num2str(j) ' = data(double(data' num2str(j) '.dat''),data' num2str(j) '.Y);']);
        eval(['testobj' num2str(j) ' = trainobj' num2str(j) ';']);
        
        eval(['trainobj' num2str(j) '.X = trainobj' num2str(j) '.X(trIdx{j}{i}, :);']);
        eval(['trainobj' num2str(j) '.Y = trainobj' num2str(j) '.Y(trIdx{j}{i}, :);']);
        % prepare test data 1 & 2
        eval(['testobj' num2str(j) '.X = testobj' num2str(j) '.X(teIdx{j}{i}, :);']);
        % eval(['testobj' num2str(j) '.Y = testobj' num2str(j) '.Y(teIdx{i}, :);']);
    end
    
    
    % ----- train on data1 & 2 and test on data1 and data2 -----
    
    for j = 1:2
        % -- set-up an SVM object
        eval(['svmobj' num2str(j) ' = svm({''optimizer="andre"'',''C=1'',''child=kernel''});']);
        eval(['if dobalanced, svmobj' num2str(j) '.balanced_ridge = balanced_ridge; end']);
    
        % -- train on data1
        eval(['[~, svmobj' num2str(j) '] = train(svmobj' num2str(j) ', trainobj' num2str(j) ');'])
    
        % -- test on data1 and data2
        eval(['testobj' num2str(j) '1 = test(svmobj' num2str(j) ', testobj1);']);
        eval(['testobj' num2str(j) '2 = test(svmobj' num2str(j) ', testobj2);']);
    
        eval(['w = get_w(svmobj' num2str(j) ')'';']);
        eval(['b0 = svmobj' num2str(j) '.b0;']);
    
        eval(['test_error(testobj1, testobj2, testobj' num2str(j) '1, testobj' num2str(j) '2, w, b0);']);
        
        % -- save predictions (cross-val)
        eval(['predicted' num2str(j) '1(teIdx{1}{i}) = testobj' num2str(j) '1.X;']);
        eval(['predicted' num2str(j) '2(teIdx{2}{i}) = testobj' num2str(j) '2.X;']);
    end
    
    % ----- outputs II: cross-val data -----
    
    for j = 1:2
        % output: cv-weights
        eval(['stats' num2str(j) '.cvoutput{i,1} = get_w(svmobj' num2str(j) ')'';']);
        % output: cv-intercepts
        eval(['stats' num2str(j) '.cvoutput{i,2} = svmobj' num2str(j) '.b0;']);
        % output: cv-distance-from-hyperplane for test data1 and data2
        if j == 1
            stats1.cvoutput{i,3} = testobj1.X * stats1.cvoutput{i,1} + stats1.cvoutput{i,2};
            stats1.cvoutput{i,4} = testobj2.X * stats1.cvoutput{i,1} + stats1.cvoutput{i,2};
            stats1.cvoutput_descrip = {'1:weights 2:intercept 3:distance from hyperplane for testing data2(xval) 4: dist for data1(xval)'};
            dist_from_hyper11(teIdx{j}{i}) = stats1.cvoutput{i,3};
            dist_from_hyper12(teIdx{j}{i}) = stats1.cvoutput{i,4};
        else
            stats2.cvoutput{i,3} = testobj2.X * stats2.cvoutput{i,1} + stats2.cvoutput{i,2};
            stats2.cvoutput{i,4} = testobj1.X * stats2.cvoutput{i,1} + stats2.cvoutput{i,2};
            stats2.cvoutput_descrip = {'1:weights 2:intercept 3:distance from hyperplane for testing data2(xval) 4: dist for data1(xval)'};
            dist_from_hyper21(teIdx{j}{i}) = stats2.cvoutput{i,4};
            dist_from_hyper22(teIdx{j}{i}) = stats2.cvoutput{i,3};
        end
    end
    
end

% ----- outputs I again: cross-val data -----

stats1.yfit = predicted11;
stats2.yfit = predicted22;
stats1.crosstestfit = predicted12;
stats2.crosstestfit = predicted21;

stats1.all{3} = dist_from_hyper11;
stats1.all{4} = dist_from_hyper12;
stats2.all{3} = dist_from_hyper22;
stats2.all{4} = dist_from_hyper21;

test_results = get_test_results(stats1, stats2, test_Y1, test_Y2, out_method);

end

% -------------------------------------------------------------------------
% CROSS-CLASSIFY USING CV-LASSOPCR

function [test_results, stats1, stats2] = cv_lassopcr_cross_classfy(data1, data2, cv_assign, test_Y1, test_Y2, out_method)

% ----- check lassopcr is in the path ----- 
if isempty(which('lasso_rocha.m'))
    error('Error: You need to have the lasso (rocha) toolbox in your path. Check your path.');
end

% ----- preallocate variables -----
for j = 1:2
    eval(['predicted1' num2str(j) ' = NaN(size(data' num2str(j) '.Y));']);
    eval(['predicted2' num2str(j) ' = NaN(size(data' num2str(j) '.Y));']);
end

% ----- get cv assignment ----- 
[teIdx, trIdx] = do_cv_assign(cv_assign);

% ----- run it with all data (not cv)----- 
for j = 1:2, eval(['[dummy, res' num2str(j) '] = predict(data' num2str(j) ', ''algorithm_name'', ''cv_lassopcr'', ''nfolds'', 1);']); end

% ----- outputs I: using all data -----
for j = 1:2
    eval(['stats' num2str(j) '.Y = res' num2str(j) '.Y;']);
    % output: weights and intercepts
    eval(['for i = 1:2, stats' num2str(j) '.all{i} = res' num2str(j) '.other_output{i}; end']);
    % description of out outputs: 3 and 4 have not been collected yet, though.
    eval(['stats' num2str(j) '.all_descrip = {''1:weights 2:intercept 3:pexp from testing data' num2str(j) '(xval) 4: pexp for the other data(xval)''};']);
end

% ----- Cross-validation loop starts here -----

if numel(teIdx)>1, if ~isequal(numel(teIdx{1}), numel(teIdx{2})), error('The number of folds should be same between dat1 and dat2.'); end, end 

for i = 1:numel(teIdx{1})
    
    for j = 1:2
        % prepare training data 1 & 2
        eval(['trainobj' num2str(j) ' = data' num2str(j) ';']); 
        eval(['trainobj' num2str(j) '.dat = data' num2str(j) '.dat(:,trIdx{j}{i});']);
        eval(['trainobj' num2str(j) '.Y = data' num2str(j) '.Y(trIdx{j}{i},1);']);
        % prepare test data 1 & 2
        eval(['testobj' num2str(j) ' = data' num2str(j) ';']);
        eval(['testobj' num2str(j) '.dat = data' num2str(j) '.dat(:, teIdx{j}{i});']);
    end
    
    % --- train on data1 and test on data1 and data2 ---
    
    for j = 1:2
        % -- train on data1
        eval(['[dummy, res' num2str(j) '] = predict(trainobj' num2str(j) ', ''algorithm_name'', ''cv_lassopcr'', ''nfolds'', 1);']);
        
        % -- save predictions (cross-val)
        eval(['pexp' num2str(j) '1 = testobj1.dat'' * res' num2str(j) '.other_output{1} + res' num2str(j) '.other_output{2};']);
        eval(['pexp' num2str(j) '2 = testobj2.dat'' * res' num2str(j) '.other_output{1} + res' num2str(j) '.other_output{2};']);
        eval(['predicted' num2str(j) '1(teIdx{1}{i}) = pexp' num2str(j) '1;']);
        eval(['predicted' num2str(j) '2(teIdx{2}{i}) = pexp' num2str(j) '2;']);
    end
    
    
    % ----- outputs II: cross-val data -----
    
    for j = 1:2
        % output: cv-weights
        eval(['stats' num2str(j) '.cvoutput{i,1} = res' num2str(j) '.other_output{1};']);
        % output: cv-intercepts
        eval(['stats' num2str(j) '.cvoutput{i,2} = res' num2str(j) '.other_output{2};']);
        % output: cv-pexp for test data1 and data2
        if j == 1
            stats1.cvoutput{i,3} = pexp11;
            stats1.cvoutput{i,4} = pexp12;
            stats1.cvoutput_descrip = {'1:weights 2:intercept 3:pexp for testing data1(xval) 4: dist for data2(xval)'};
        else
            stats2.cvoutput{i,3} = pexp22;
            stats2.cvoutput{i,4} = pexp21;
            stats2.cvoutput_descrip = {'1:weights 2:intercept 3:pexp for testing data2(xval) 4: dist for data1(xval)'};
        end
    end
end

stats1.yfit = predicted11;
stats2.yfit = predicted22;
stats1.crosstestfit = predicted12;
stats2.crosstestfit = predicted21;
stats1.test_Y = test_Y2;
stats2.test_Y = test_Y1;

stats1.all{3} = predicted11;
stats1.all{4} = predicted12;
stats2.all{3} = predicted22;
stats2.all{4} = predicted21;

test_results = get_test_results(stats1, stats2, test_Y1, test_Y2, out_method);

end