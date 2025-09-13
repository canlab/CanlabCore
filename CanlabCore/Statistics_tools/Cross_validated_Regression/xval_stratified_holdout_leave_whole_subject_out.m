% Create stratified k-fold holdout sets, balancing on outcome class and keeping observations from a grouping variable (e.g., subject) together
%
% - Works for classification (binary classes) and regression (continuous outcomes)
% - Balances (stratifies) folds on outcome class membership (for categorical outcomes Y) or
%   quartiles of outcome (for continuous outcomes Y)
%
% :Usage:
% ::
%
%     [list outputs here] = xval_stratified_holdout_leave_whole_subject_out(Y, id, [optional inputs])
%
% For objects: Type methods(object_name) for a list of special commands
%              Type help object_name.method_name for help on specific
%              methods.
%
% ..
%     Author and copyright information:
%
%     Copyright (C) 2020 Tor Wager
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
%   **Y:**
%        Class variable to stratify by, numeric vector, defines classes (not continuous
%        at the moment). This would be the outcome in SVM
%
%   **id:**
%        Grouping id variable to keep together (e.g., subject ID)
%
% :Optional Inputs:
%   **nfolds:**
%        Number of folds in k-fold; default = 10
%
%   **doverbose:**
%        Verbose output
%
%   **doplot:**
%        Plot output/holdout sets
%
%   **downsample_to_balanced_classes:**
%        True or false [default]. Remove training cases (undersample) to
%        keep training and test sets balanced on outcome, for SVMs and
%        other algorithms that work best with balanced training sets.
%
% :Outputs:
%
%   **out1:**
%        description of out1
%
%   **out2:**
%        description of out2
%
% :Examples:
% ::
%
%    % give examples of code here
%    param1 = abc();
%    param2 = xyz();
%    [out1,out2] = func_call(param1, param2)
%
% :See also:
%   - cvpartition, other CANlab holdout set tools (could be merged!)
%
% ..
%    Programmers' notes:
%
% - Matlab has many great tools, e.g. cvpartition, but they don't group obs from same subject together
% - It is desirable to stratify (balance) across outcome classes, so
% train/test proportions are similar across folds, and also hold out all
% observations from a group (subject) together to ensure independence
% across training folds.
% - Issue described here: https://www.mathworks.com/matlabcentral/answers/203155-how-to-manually-construct-or-modify-a-cross-validation-object-in-matlab
% ..

% ..
%    DEFAULTS AND INPUTS
% ..

function [trIdx, teIdx, holdout_integer_indic] = xval_stratified_holdout_leave_whole_subject_out(varargin)


% ----------------------------------------------------------------------
% Parse inputs
% ----------------------------------------------------------------------
% Uses the inputParser object. Older schemes are below.

ARGS = parse_inputs(varargin{:});

fn = fieldnames(ARGS);

for i = 1:length(fn)
    str = sprintf('%s = ARGS.(''%s'');', fn{i}, fn{i});
    eval(str)
end

% Logical flags
% ----------------------------------------------------------------------
if any(strcmp(varargin, 'noverbose')), doverbose = false; end
if any(strcmp(varargin, 'noplots')), doplot = false; end

% ----------------------------------------------------------------------
% Select holdout sets
% -------------------------------------------------------------------------

if doverbose
    
    fprintf('Holdout set:  %d-fold cross-validation, leaving whole subject out, stratifying on treatment group\n', nfolds);
    
end

% Keep subject together (all images in same train/test)
[group_together, wh] = unique(id,'stable');

 % Y is categorical (vs. continuous)...deal with stratification and plots differently if so
class_labels = unique(Y);  % Unique outcome labels or values
nclasses = length(class_labels);

is_cat = length(nclasses) < 3;

% Gid: group, subject-wise. Assumes all obs for a subject have same Y value
% for stratification purposes. If they do not (e.g., Y=1 and Y=-1 for each
% of a pair of obs per subject), the CV partition could be non-random
% depending on the order you have entered cases, so check carefully. Future
% versions may want to check for this and turn off stratification in
% cvpartition in this case.

% prep vectors of which obs belong to which group
class_by_group = count_class_by_group(id, Y, group_together);

Yid = Y(wh);


% Partition subject-wise, and then re-assign observations to folds
% This keeps all subjects together.
% nfolds defined in input parser, default = 10

if is_cat
    % by itself, this will keep cases in each class approximately
    % proportional to the base rate of the two classes. We need more below
    % to under/oversample cases if we want balanced classes in each
    % training fold.
    cvpart = cvpartition(Yid, 'k', nfolds, 'Stratify', true);
    
else
    % Group into 4 categories based on mean for each id, and stratify by this.
    % 
    Ybins = scores_to_bins(Yid);
    cvpart = cvpartition(Ybins, 'k', nfolds, 'Stratify', true);

end

% Reassign train/test in observation-wise list
[trIdx, teIdx] = deal(cell(1, nfolds));

for i = 1:nfolds
    
    trIdx{i} = ismember(id, group_together(cvpart.training(i)));
    
    teIdx{i} = ismember(id, group_together(cvpart.test(i)));
    
end

% Down-sample (or up, in future) largest class to make them equal, if
% requested
if downsample_to_balanced_classes
    trIdx = undersample_to_smallest_class(trIdx, nfolds, id, group_together, Y, doverbose);
end

testmat = cat(2, teIdx{:});
holdout_integer_indic = indic2condf(testmat);

% ----------------------------------------------------------------------
% Collect and plot diagnostic output
% -------------------------------------------------------------------------

% Tabulate proportions of each class in each fold
for i = 1:nfolds
    
    if is_cat
        % tabulate(Y(trIdx{1})) How many obs of each outcome class in fold?
        prop_table = tabulate(Y(trIdx{i}));
        class_proportions(i, :) = prop_table(:, 3)';
        
    else
        % use mean outcome for each fold
        outcome_means(i, 1) = mean(Y(trIdx{i}));
        outcome_std(i, 1) = std(Y(trIdx{i}));
        outcome_vals{i} = Y(trIdx{i});
    end
    
end

if doverbose
    fprintf('Groups (e.g., subjects): %d, Images: %d, Classes: %d\n', length(Yid), length(id), length(unique(Y)))
end

% Check that all ids appear in exactly 1 fold

    fold_sim = condf2indic(id)' * cat(2, teIdx{:});
    num_folds_for_each_grouping_id = sum(fold_sim' > 0);
    
    if all(num_folds_for_each_grouping_id) == 1
        if doverbose, fprintf('Observations from each id appear in exactly 1 fold\n'); end
    else
        disp('Observations from each id appear in > 1 fold! Check/debug this function.');
        keyboard
    end
    
% Plot holdout sets
if doplot
    
    % scaledid = 1 + (rankdata(id) ./ max(id));
    
    to_plot = [teIdx{:} scale(Y)+2];
    
    create_figure('Test sets', 2, 2);

    imagesc(to_plot); axis tight; set(gca, 'YDir', 'reverse');
    title('holdout sets, with outcome (Y)');
    xlabel('Fold');
    ylabel('Observation index');
    
    subplot(2, 2, 2);
    [~, wh] = sort(Y);
    imagesc(to_plot(wh, :)); axis tight; set(gca, 'YDir', 'reverse');
    title('holdout sets, sorted by outcome');
    ylabel('Observations (resorted)');
    
    subplot(2, 2, 3);
    if is_cat
        bar(class_proportions,'stacked');
        ylabel('Class proportions')
    else
        %         barplot_columns(outcome_vals, 'noviolin', 'nofig'); % slow
        for i = 1:length(outcome_vals)
            plot(i * ones(size(outcome_vals{i})), outcome_vals{i}, 'o', 'Color', [.6 .6 .6]);
        end
        bar(outcome_means);
        h = errorbar(outcome_means, outcome_std);
        set(h, 'Color', [.3 .3 .3]);
        ylabel('Outcome value');
    end
    
    colormap(cool)
    
    xlabel('Fold');
    axis tight;
    
    subplot(2, 2, 4);
	imagesc(fold_sim)
    set(gca, 'YDir', 'reverse');
    xlabel('Fold');
    ylabel('Unique grouping id');
    axis tight; 
    title('Images in each fold, by grouping id')
    colorbar

end

end % main function




% ======== This goes at the end of the file / after the main function ========
function ARGS = parse_inputs(varargin)

p = inputParser;

% Validation functions - customized for each type of input
% ----------------------------------------------------------------------

valfcn_scalar = @(x) validateattributes(x, {'numeric'}, {'nonempty', 'scalar'});

valfcn_number = @(x) validateattributes(x, {'numeric'}, {'nonempty'}); % scalar or vector

% valfcn_logical = @(x) validateattributes(x, {'numeric' 'logical'}, {'nonempty', 'scalar', '>=', 0, '<=', 1}); % could enter numeric 0,1 or logical


% Required inputs
% ----------------------------------------------------------------------
p.addRequired('Y', valfcn_number);
p.addRequired('id', valfcn_number);

% Optional inputs
% ----------------------------------------------------------------------
% Pattern: keyword, value, validation function handle

p.addParameter('doverbose', true);
p.addParameter('doplot', true);
p.addParameter('nfolds', 10, valfcn_scalar);
p.addParameter('downsample_to_balanced_classes', false);

% Parse inputs and distribute out to variable names in workspace
% ----------------------------------------------------------------------
% e.g., p.parse([30 1 0], [-40 0 10], 'bendpercent', .1);
p.parse(varargin{:});

ARGS = p.Results;

end % parse_inputs


function x_binned = scores_to_bins(x)

n_bins = 4;
% adapted from line_plot_multisubject

bins = prctile(x, linspace(0, 100, n_bins + 1));
bins(end) = Inf;

for j=1:n_bins %make the bins
    
    wh = x >= bins(j) & x < bins(j+1);
    
    x_binned(wh) = j; %nanmean(t(wh, :));
    
end

end % scores to bins



    function class_by_group = count_class_by_group(id, Y, group_together)

        class_labels = unique(Y);  % Unique outcome labels or values
        nclasses = length(class_labels);
        ngroups = length(group_together);

        group_indic = condf2indic(id, 'integers', ngroups);

        class_by_group = zeros(ngroups, nclasses);

        for j = 1:nclasses
            for i = 1:length(group_together)
                is_class = Y == class_labels(j);

                class_by_group(i, j) = sum(group_indic(:, i) & is_class);

            end
        end

    end % count_class_by_group


% ************
% function trIdx = undersample_to_smallest_class(trIdx, nfolds, id, group_together, Yid, Y)
% 
% for i = 1:nfolds
% 
%     ytrain = Y(trIdx{i});
%     classes = unique(ytrain);
%     idtrain = id(trIdx{i});
% 
%     % Count number obs in each class
%     for j = 1:length(classes), count(j) = sum(ytrain == classes(j)); end
% 
%     [maxCount, class_to_downsample] = max(count);
%     [minCount, smallest_class] = min(count);
% 
%     if maxCount == minCount, continue; end  % we are already balanced
% 
%     class_by_group = count_class_by_group(idtrain, ytrain, group_together);
% 
%     class_diffs_by_group = diff(class_by_group')';
% 
%     if class_to_downsample == 1
%         eligible_groups = group_together(class_diffs_by_group < 0); % indices of which groups have more obs for the larger class 
% 
%     else % 2 is greater
% 
%         eligible_groups = group_together(class_diffs_by_group > 0); 
%     end
% 
%     rp = randperm(length(eligible_groups)); % random permutation - choose among eligible groups randomly
%     eligible_groups = eligible_groups(rp);
% 
%     % cumulative sum of differences in cases by group, for eligible groups
%     diffcount = class_diffs_by_group(eligible_groups);
%     diffcount = cumsum(diffcount(rp));
% 
%     wh_groups_to_remove = abs(diffcount) <= abs(diff(count));
%     groups_to_remove = eligible_groups(wh_groups_to_remove);
% 
%     indx = ismember(id, groups_to_remove);
% 
%     trIdx{i}(indx) = 0;  % remove them
% 
%    % REPORT COUNT
%        ytrain = Y(trIdx{i});
%     classes = unique(ytrain);
%     idtrain = id(trIdx{i});
% 
%     % Count number obs in each class
%     for j = 1:length(classes), count(j) = sum(ytrain == classes(j)); end
%     [ count class_to_downsample]
% 
% 
% end % folds
% 
% end % function


function trIdx = undersample_to_smallest_class(trIdx, nfolds, id, group_together, Y, doverbose)
% undersample_to_smallest_class undersamples the largest class within each fold.
%
% :Usage:
% ::
%     trIdx = undersample_to_smallest_class(trIdx, nfolds, id, group_together, Y, doverbose)
%
% :Inputs:
%
%   **trIdx:**
%        Cell array of logical vectors indicating which observations are in-fold.
%
%   **nfolds:**
%        Number of folds.
%
%   **id:**
%        Numeric vector of group identifiers for each observation.
%
%   **group_together:**
%        A vector of group numbers that should be considered together.
%
%   **Yid:**
%        (Legacy input; not used in current implementation.)
%
%   **Y:**
%        Vector of class labels for each observation.
%
% :Outputs:
%
%   **trIdx:**
%        Updated cell array of logical vectors where, for each fold, groups
%        from the largest class have been removed until the classes are approximately balanced.
%
% :Examples:
% ::
%     trIdx = undersample_to_smallest_class(trIdx, nfolds, id, group_together, Yid, Y);
%
% :References:
%   See also count_class_by_group.
%
% :See also:
%   count_class_by_group
%
% -------------------------------------------------------------------------

if doverbose
    disp('Undersampling largest class to keep numbers of observations balanced in each training fold')
end


for i = 1:nfolds
    % Get indices for observations in the current fold.
    fold_idx = find(trIdx{i});
    ytrain = Y(fold_idx);
    classes = unique(ytrain);
    idtrain = id(fold_idx);

    % Reinitialize count for current fold.
    count = zeros(length(classes), 1);
    for j = 1:length(classes)
        count(j) = sum(ytrain == classes(j));
    end

    [maxCount, class_to_downsample] = max(count);
    [minCount, smallest_class] = min(count);

    % If the classes are already balanced, do nothing.
    if maxCount == minCount
        continue;
    end

    % Count observations by group for each class.
    % The function count_class_by_group should return a matrix with rows corresponding
    % to each group in group_together and columns corresponding to classes.
    class_by_group = count_class_by_group(idtrain, ytrain, group_together);

    % Ensure we have exactly two classes.
    if size(class_by_group, 2) ~= 2
        error('undersample_to_smallest_class currently supports only two classes.');
    end

    % Compute the difference between counts for class2 and class1 for each group.
    % A positive value means more observations in class2; negative means more in class1.
    class_diffs_by_group = class_by_group(:,2) - class_by_group(:,1);

    % Determine which groups are eligible for removal.
    % If the largest class is the first class, we want groups with negative differences.
    if class_to_downsample == 1
        eligible_idx = find(class_diffs_by_group < 0);
    else % largest class is the second class
        eligible_idx = find(class_diffs_by_group > 0);
    end

    % Map eligible indices to group numbers.
    eligible_groups = group_together(eligible_idx);

    % Randomize the order of eligible groups.
    rp = randperm(length(eligible_groups));
    eligible_groups = eligible_groups(rp);

    % Get the differences for these eligible groups (in randomized order)
    diffcount = class_diffs_by_group(ismember(group_together, eligible_groups));
    % Compute the cumulative sum of differences.
    diffcount = cumsum(diffcount);

    % Determine the threshold for removal based on the overall class imbalance.
    overall_diff = abs(diff(count));

    % Identify groups to remove: remove groups until the cumulative difference
    % is less than or equal to the overall difference.
    wh_groups_to_remove = abs(diffcount) <= overall_diff;
    groups_to_remove = eligible_groups(wh_groups_to_remove);

    % Remove observations in the current fold that belong to the groups to remove.
    remove_idx = fold_idx(ismember(idtrain, groups_to_remove));
    trIdx{i}(remove_idx) = false;


    if doverbose

        % (Optional) Report counts after removal.
        ytrain_after = Y(trIdx{i});
        new_count = zeros(length(classes), 1);
        for j = 1:length(classes)
            new_count(j) = sum(ytrain_after == classes(j));
        end
        fprintf('Fold %d: before removal: %s, after removal: %s\n', i, mat2str(count'), mat2str(new_count'));

    end

end

end % function undersample





