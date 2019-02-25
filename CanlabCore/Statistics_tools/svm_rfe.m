function out = svm_rfe(dat, varargin)

% This function performs Recursive Feature Elimination (RFE) using
% 'predict' function for model training and cross-validation. It returns 
% details of each elimination step and a model trained on the number of 
% features that shows the best accuracy.
%
%
% :Features:
%       - repeatedly trains a desired model (SVM recommended) using the 'predict'
%       function and eliminates least significant weights until a desired
%       number of features is reached
%       - elimination step and target number of features can be defined by
%       user
%       - evaluates accuracy in each iteration based on CV error
%
%
% :Usage:
%       out = svm_rfe(dat, varargin)
%
%
% :Input variables:
%       **n_removal**
%           The number of features to be removed in each iteration. Default
%           is 10000.
%           e.g., 'n_removal', 15000
%
%       **n_finalfeat**
%           The final number of features after elimination. Default is
%           50000.
%           e.g., 'n_finalfeat', 30000
%
%       **n_initialfeat**
%           The number of features, from which the algorithm starts the elimination
%           by the defined step. The excessive features are eliminated in one step.
%           Default is the whole feature set.
%           e.g., 'n_initialfeat', 100000
%
%       **for other inputs see fmri_data.predict
%
% :Example:
%
%       out = svm_rfe(dat, 'algorithm_name', 'cv_svm', 'nfolds', 5, 'error_type', 'mcr');
%       out = svm_rfe(dat, 'n_finalfeat', 30000, 'algorithm_name', 'cv_svm', 'nfolds', 5, 'error_type', 'mcr);
%       out = svm_rfe(dat, 'n_removal', 12000, 'n_finalfeat', 30000, 'algorithm_name', 'cv_svm', 'nfolds', 5, 'error_type', 'mcr');
%       out = svm_rfe(dat, 'n_removal', 12000, 'n_finalfeat', 30000, 'n_initialfeat', 150000, 'algorithm_name', 'cv_svm', 'nfolds', 5, 'error_type', 'mcr');
%
%
% Created by Lada Kohoutova (ladakohoutova@gmail.com)
%            Cocoan lab (Wani Woo's lab), 2/19/2019
%

% default

n_removal = 10000;
n_finalfeat = 50000;
init_feat = 0;

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            % functional commands
            case 'n_removal'
                n_removal = varargin{i+1};
                varargin{i+1} = [];
            case 'n_finalfeat'
                n_finalfeat = varargin{i+1};
                varargin{i+1} = [];
            case 'n_initialfeat'
                n_initialfeat = varargin{i+1};
                init_feat = 1;
                varargin{i+1} = [];
        end
    end
end

orig_indx = 1:size(dat.dat,1);
data_dim = size(dat.dat,1);
out.descending_indx = [];
out.whkeep_orginal_idx{1} = orig_indx;
out.removed_index{1} = [];
dat_loop = dat;
i = 1;

if init_feat == 1 % the first iteration
    
    n_initremoval = data_dim - n_initialfeat;    % calculate how many features to remove
    
    str = sprintf('Initial training with all features...');
    fprintf('%s\n', str)
    
    [~, stats, optout] = predict(dat_loop, varargin{:});
    
    out.cv_accuracy(i,1) = 1-stats.cverr;            % collecting outputs (weights, accuracy) from training
    out.n_features(i,1) = size(dat_loop.dat,1);
    w = optout{1};
    
    w = abs(w);
    [~, out.ascending_idx{i}] = sort(w);
    ordered_orig_indx = orig_indx(out.ascending_idx{i});
    out.descending_indx = [ordered_orig_indx(1:n_initremoval) out.descending_indx];  % keeping ordered indices
    out.removed_index{i+1} = ordered_orig_indx(1:n_initremoval);
    
    % removal
    remove_idx = out.ascending_idx{i}(1:n_initremoval);
    dat_loop.dat(remove_idx, :) = [];
    
    orig_indx(out.ascending_idx{i}(1:n_initremoval)) = []; % clear removed original indices
    out.whkeep_orginal_idx{i+1} = orig_indx;    % save kept original indices
    
    data_dim = data_dim - n_initremoval;  % reduce the dimension
    
    fprintf('\n')
    str = sprintf('Eliminated %d features.', n_initremoval);
    fprintf('%s\n', str)
    
    i = i+1;
    
end

while 1
    
    str = sprintf('Training with %d features...', data_dim);
    fprintf('%s\n', str)
    
    [~, stats, optout] = predict(dat_loop, varargin{:});
    
    out.cv_accuracy(i,1) = 1-stats.cverr;            % collecting outputs (weights, accuracy) from training
    out.n_features(i,1) = size(dat_loop.dat,1);
    w = optout{1};
    
    if out.n_features(i) <= n_finalfeat, break, end
    
    w = abs(w);
    [~, out.ascending_idx{i}] = sort(w);               % sorting weights in ascending order
    
    if data_dim-n_removal < n_finalfeat, n_removal = data_dim-n_finalfeat; end   % condition for the final removal
    
    ordered_orig_indx = orig_indx(out.ascending_idx{i});
    out.descending_indx = [ordered_orig_indx(1:n_removal) out.descending_indx];  % keeping ordered indices
    out.removed_index{i+1} = ordered_orig_indx(1:n_removal);
    
    % removal
    remove_idx = out.ascending_idx{i}(1:n_removal);
    dat_loop.dat(remove_idx, :) = [];
    
    orig_indx(out.ascending_idx{i}(1:n_removal)) = []; % clear removed original indices
    out.whkeep_orginal_idx{i+1} = orig_indx;    % save kept original indices
    
    data_dim = data_dim - n_removal;  % reduce the dimension
    
    fprintf('\n')
    str = sprintf('Eliminated %d features.', n_removal);
    fprintf('%s\n', str)
    
    % rerun predict
    
    i = i + 1;
    
end

out.descending_indx = [orig_indx out.descending_indx];

out.smallestnfeat_stats = stats;  % save stats of the smallest model
indxToRemove = [];
for i=1:length(out.cv_accuracy)
    indxToRemove = [out.removed_index{i} indxToRemove];
end

% [indx,~] = find(~out.smallestnfeat_stats.weight_obj.removed_voxels);
indx = (1:size(out.smallestnfeat_stats.weight_obj.volInfo.wh_inmask,1))';
out.smallestnfeat_stats.weight_obj.removed_voxels(indx(indxToRemove)) = true;

[max_acc, max_idx] = max(out.cv_accuracy);

out.best_n_features = out.n_features(max_idx);
out.best_accuracy = max_acc;

dat.dat = dat.dat(out.whkeep_orginal_idx{max_idx},:);

fprintf('\n')
str = sprintf('Training the final model...');
fprintf('%s\n', str)

[cverr, out.bestaccuracy_stats, optout] = predict(dat, varargin{:}); % out.bestaccuracy_stats

indxToRemove = [];
for i=1:max_idx
    indxToRemove = [out.removed_index{i} indxToRemove];
end

%[indx,~] = find(~out.bestaccuracy_stats.weight_obj.removed_voxels);
indx = (1:size(out.bestaccuracy_stats.weight_obj.volInfo.wh_inmask,1))';
out.bestaccuracy_stats.weight_obj.removed_voxels(indx(indxToRemove)) = true;   % update removed voxels in the weight_obj


end
