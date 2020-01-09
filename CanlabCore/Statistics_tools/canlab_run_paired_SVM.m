function [paired_d, stats, optout, cat_obj] = canlab_run_paired_SVM(varargin)
% Run 5-fold, leave-whole-subject out cross-validated SVM classification on a pair of fmri_data image objects.
%
% - Requires two fmri_data objects with paired images in each (one image per pair in each object)
% - Images to be paired must be in same order in each object
%
% :Usage:
% ::
%
% [paired_d, stats, optout, cat_obj] = canlab_run_paired_SVM(pos_image_obj, neg_image_obj)
%
% ..
%     Author and copyright information:
%     Copyright (C) Tor Wager, Jan 2020
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
% :Inputs:
%
%   **pos_image_obj:**
%        An fMRI-data object, often with one image per subject
%
%   **neg_image_obj:**
%         An fMRI-data object, often with one image per subject
%         Must have same number of images as pos_image_obj!
%         This function assumes the images are paired (e.g., same subjects)
%
% :Optional Inputs:
%   **nfolds:**
%        Followed by number of folds desired (default = 5)
%       
%   **dobootstrap:**
%        Followed by logical 1/0
%
%   **boot_n:**
%        Followed by number of boostrap samples desired (default = 100)
%
%   **parallelstr:**
%        Followed by string 'parallel' to use parallel processing in fmri_data.predict.m
%
%   **dol2norm:**
%        Normalize each image by the L2-norm
%
%   **dozscoreimages:**
%        Normalize each image by replacing values with z-scores across voxels
%
% :Outputs:
%
%   **paired_d:**
%        Simple Cohen's d effect size for forced-choice cross-validated classification
%        This is the main metric of interest for evaluating classification performance.
%
%   **stats:**
%        Stats output structure for SVM.
%        This function adds specific fields:
%        stats.paired_hyperplane_dist_scores: The distance from the hyperplane for the paired pos - neg images for each unit (e.g., subject). 
%                                             Positive indicates correct classification.
%        stats.paired_d:                      Simple Cohen's d effect size for forced-choice cross-validated classification
%                                             This is the main metric of interest for evaluating classification performance.
%        stats.paired_accuracy:               Accuracy of forced-choice (two choice) classification, cross-validated
%        stats.ROC:                           More output from roc_plot,
%        including sensitivity/specificity/PPV (will be identical for two-choice paired classification)
%
% :Examples:
% ::
%
%    [paired_d, stats] = canlab_run_paired_SVM(pos_images, neg_images);
%    paired_d
%
%    % Plot Receiver Operating Characteristic (ROC):
%    ROC = roc_plot(stats.dist_from_hyperplane_xval, stats.Y > 0, 'twochoice');
%
%    % Plot individual participants
%    create_figure('subjects'); 
%    plot(stats.paired_hyperplane_dist_scores, 'o');
%    plot_horizontal_line(0); 
%    xlabel('Participant'); ylabel('Classifier score');
%
%    % Re-run with z-scored images:
%    [paired_d, stats] = canlab_run_paired_SVM(pos_images, neg_images, 'dozscoreimages', true);
%
% :References:
%   Vapnik 1995 (SVM)
%
% :See also:
%   - canlab batch scripts on canlab.github.io
%   - fmri_data.predict()
%

% ..
%    Programmers' notes:
%    List dates and changes here, and author of changes
%    Created by Tor Wager, Jan 2020
% ..

% ----------------------------------------------------------------------
% Parse inputs
% ----------------------------------------------------------------------
ARGS = parse_inputs(varargin{:});

% Stack images 
% --------------------------------------------------------------------

n1 = size(ARGS.pos_image_obj.dat, 2);
n2 = size(ARGS.neg_image_obj.dat, 2);

if n1 ~= n2, error('Two image sets compared must be equal n for paired analysis!'); end

% Stack images 
% --------------------------------------------------------------------

cat_obj = cat(ARGS.pos_image_obj, ARGS.neg_image_obj);

% Define and add outcome
% --------------------------------------------------------------------

cat_obj.Y = [ones(n1, 1); -ones(n2, 1)];

% Define holdout set - 5-fold leave-whole-subject out
% --------------------------------------------------------------------

% Assign test set for each fold to all images in paired conditions
% If participants are crossed with conditions, this leaves one whole participant out (all images across conditions)

cvpart = cvpartition(n1,'k',ARGS.nfolds);
holdout_set = zeros(n1, 1); % container to store integers

for i = 1:ARGS.nfolds

    mytest = cvpart.test(i);
    holdout_set(mytest) = i;
    
end

holdout_set = [holdout_set; holdout_set]; % Replicate to leave out both conditions for each subject

% Optional scaling, masking, transformations
% --------------------------------------------------------------------
if ARGS.dol2norm
    
    cat_obj = rescale(cat_obj, 'l2norm_images');
    
end

if ARGS.dozscoreimages
    
    cat_obj = rescale(cat_obj, 'zscoreimages');
    
end


% Run SVM prediction model
% --------------------------------------------------------------------
if ARGS.dobootstrap
    [cverr, stats, optout] = predict(cat_obj, 'algorithm_name', 'cv_svm', 'nfolds', holdout_set, 'bootsamples', ARGS.boot_n, 'error_type', 'mcr', ARGS.parallelstr);
    % Threshold, if possible - can re-threshold later with threshold() method
    stats.weight_obj = threshold(stats.weight_obj, .05, 'unc');
    
else
    [cverr, stats, optout] = predict(cat_obj, 'algorithm_name', 'cv_svm', 'nfolds', holdout_set, 'error_type', 'mcr',  ARGS.parallelstr);
end

% Paired stats and effect size
% --------------------------------------------------------------------

% The distance from the hyperplane for the paired pos - neg images for each unit (e.g., subject). 
% Positive indicates correct classification.
stats.paired_hyperplane_dist_scores = stats.dist_from_hyperplane_xval(1:n1) - stats.dist_from_hyperplane_xval(n1+1:end);

% Effect size
stats.paired_d = mean(stats.paired_hyperplane_dist_scores) ./ std(stats.paired_hyperplane_dist_scores);

paired_d = stats.paired_d;

% Two-choice classification accuracy
stats.paired_accuracy = sum(stats.paired_hyperplane_dist_scores > 0) ./ length(stats.paired_hyperplane_dist_scores);

stats.ROC = roc_plot(stats.dist_from_hyperplane_xval, stats.Y > 0, 'twochoice', 'noplot');

end % function




function ARGS = parse_inputs(varargin)

p = inputParser;

% Validation functions - customized for each type of input
% ----------------------------------------------------------------------

valfcn_fmridata = @(x) validateattributes(x, {'fmri_data'}, {'nonempty'});

valfcn_scalar = @(x) validateattributes(x, {'numeric'}, {'nonempty', 'scalar'});
valfcn_logical = @(x) validateattributes(x, {'logical' 'numeric'}, {'nonempty', 'scalar'}); % could enter numeric 0,1 or logical

% Required inputs
% ----------------------------------------------------------------------
p.addRequired('pos_image_obj', valfcn_fmridata);
p.addRequired('neg_image_obj', valfcn_fmridata);

% Optional inputs
% ----------------------------------------------------------------------
% Pattern: keyword, value, validation function handle

p.addParameter('nfolds', 5, valfcn_scalar);
p.addParameter('dobootstrap', false, valfcn_logical);
p.addParameter('boot_n', 100, valfcn_scalar);
% p.addParameter('parallelstr', 'parallel', @(x) validateattributes(x, {'char'}));
p.addParameter('parallelstr', 'parallel');
p.addParameter('dol2norm', false, valfcn_logical);
p.addParameter('dozscoreimages', false, valfcn_logical);

% Parse inputs and distribute out to variable names in workspace
% ----------------------------------------------------------------------
p.parse(varargin{:});

ARGS = p.Results;

end

    