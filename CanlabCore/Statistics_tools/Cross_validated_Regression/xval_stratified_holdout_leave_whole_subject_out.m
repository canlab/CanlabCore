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


% Gid: group, subject-wise. Assumes all obs for a subject have same Y value.
Yid = Y(wh);

 % Y is categorical (vs. continuous)...deal with stratification and plots differently if so
is_cat = length(unique(Y)) < 3;

% Partition subject-wise, and then re-assign observations to folds
% This keeps all subjects together.
% nfolds defined in input parser, default = 10

if is_cat
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

end

