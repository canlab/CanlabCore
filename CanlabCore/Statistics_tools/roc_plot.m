function ROC = roc_plot(input_values, binary_outcome, varargin)
% This function makes a specific kind of ROC curve plot, based on input
% values along a continuous distribution and a binary outcome variable
% (logical)
%
% :Usage:
% ::
%
%     ROC = roc_plot(input_values, binary_outcome, ['include', include])
%
% Include is an optional logical variable of cases to include
%
% :Optional Inputs:
%
%   **'include':**
%        followed by logical vector of cases to include
%
%   **'threshold':**
%        followed by a priori threshold cutoff for determining misclassification
%
%   **'threshold_type':**
%        followed by thresh type: choices below:
%          - 'Optimal balanced error rate'
%          - 'Optimal overall accuracy' [default]
%          - 'Minimum SDT bias'
%          - [Enter threshold OR threshold_type]
%
%   **'color':**
%        followed by color, e.g., 'r' or [1 .5 0]
%
%   **'plotmethod':**
%        followed by 'deciles' [default] or 'observed'
%
%   **'nonormfit':**
%        suppress normal curve fitting to ROC
%
%   **'plothistograms':**
%        plot histograms of the signal present/absent distributions
%
%   **'writerscoreplus':**
%        Write text file for input into RScorePlus by Lew Harvey
%
%   **'boot':**
%        [default] Bootstrap 95% confidence intervals for sens, spec, PPV at threshold
%
%   **'noboot':**
%        Skip bootstrap
%
%   **'balanced':**
%        Balanced accuracy for single interval classification
%        THIS IS NOT COMPLETELY IMPLEMENTED BECAUSE IT AFFECTS ACCURACY
%        ESTIMATES, BUT NOT P-VALUES OR THRESHOLD AT WHICH TO EVALUATE SENS/SPEC
%
%   **'dependent':**
%        followed by vector of subject IDs, e.g., ('dependent',[1,1,2,2,3,3].
%
%        This will perform multilevel version of binomial test for single interval classification.
%
%   **'noplot':**
%        Skip generating plots
%
%   **'nooutput':**
%        Suppress text output
%
% :Outputs:
%
%   A structure containing the true and false pos rates (tpr, fpr) along the curve
%   and the criterion threshold values of the input variable (thr) corresponding to these rates.
%
%   Uses the function roc_calc.m
%
%   Also returns some information about misclassified observations
%   and line handle for ROC line plot and other statistics:
%     - area under ROC curve
%     - accuracy statistics based on binomial test
%     - PPV
%
% :Examples:
% ::
%
%    ROC = roc_plot(pattern_exp_values, ishot);
%    ROC = roc_plot(pattern_exp_values, ishot, 'threshold', 2.5);
%    ROC = roc_plot(pattern_exp_values, ishot, 'color', 'r', 'twochoice');
%    ROC = roc_plot(pattern_exp_values, ishot, 'color', 'r', 'twochoice', 'nonormfit');
%    ROC = roc_plot(pexp, logical(outcome), 'color', 'g', 'plothistograms', 'threshold', 0.3188);
%    ROC = roc_plot(pexp, logical(outcome), 'twochoice', 'color', 'b', 'plothistograms');
%    ROC = roc_plot(pexp, logical(outcome), 'writerscoreplus');
%    ROC = roc_plot(pexp, logical(outcome), 'color', 'r', 'plotmethod', 'observed', 'plothistograms');
%    ROC = roc_plot(pexp, logical(outcome), 'color', 'm', 'plotmethod', 'observed', 'plothistograms', 'Optimal overall accuracy');
%
% For a whole image with p-values, this may be helpful.
%
% Pre-specifies p-values you want to evaluate at.
% ::
%    rocout = roc_plot(1-t.p, truevals, 'plotmethod', 'observed', 'valuestoevaluate', 1 - [.5:-.1:.1 .05 .01 .005 .001 .0001 .00001 .000001 .0000001 .00000001]);
%
% ..
%    Tor Wager, Feb 2012
%    See Notes in text for more programming details.
%

%    Notes:
%    Tor Edited 3/17/2012 to add standard Gaussian signal detection fit curves,
%    effect size estimates based on Gaussian equal variance model
%
%    Tor Edited 3/20/2012 to fix AUC estimate and add/change output.
%
%    Tor Tested 3/20/12 against RScorePlus.  There are some consequences of
%    binning for input to RScorePlus, inclding that the "response" criteria are inferred
%    in RScorePlus, but they are given if we are selecting arbitrary criteria based on
%    continuous measures (e.g., brain activity). This influences the actual
%    estimates of the mean and std. signal distributions, as well as the
%    sens/spec estimates within response bins. This function does not use
%    arbitrary bins whenever possible, and uses the full ROC across thresholds
%    corresponding to each unique increment in specificity for AUC
%    calculation.
%
%    Edited 3/11/2014: Luke Chang to add balanced accuracy option for single
%    interval classification when classes are unbalanced
%
%    Edited 6/24/2014: Luke Chang to add multilevel binomial test for single
%    interval classification.  Assumes each subject has equal number of
%    trials.  Requires at least more than 20 subjects to ensure distribution
%    is reasonably approximated by normal distribution.  Uses one sample-test
%    across subjects.
%
%    Edited 1/13/2015: Luke Chang - added option to suppress plots for running
%    on cluster.
%
%    Edited 3/25/2015: Luke Chang - added option to suppress text output for
%    speeding up computations on cluster
%
%    Edited 8/2015: Tor Wager - reduced length of thr to speed computation
%    with large numbers of values (e.g., images with 50K+voxels)
%
%    12/8/2018: Tor Wager - made 2 substantive changes:
%     For 'Optimal balanced error rate', we don't want the threshold
%     that optimizes overall sensitivity+specificity; we want the threshold that
%     minimizes the discrepancy between sensitivity and specificity.
%     This is because with unbalanced classes, one can always choose
%     "yes" or "no" and the threshold will be at the bounds, yielding a
%     sens/spec of 0% and 100% or vice versa.
% 
%     Threshold optimization correction for binomial P-values
%     Chance values can be tricky, and optimizing a threshold can be
%     circular. For example, if only 33% of observations are true
%     positives, and the classifier always guesses "no", it will be right
%     66% of the time, which is better than the 50% expected under random guessing.
%     We want to build in a correction for threshold selection.
%     The expectation under chance is that we'll adopt a strategy that
%     gives us at least the base-rate of true positive or true negative
%     results. For 33% true-pos, this would be 66%.
%     We take the max p-value of the difference from 0.5 and the baserate
%     or 1 - baserate.
%     Not implemented for 'dependent' binomial test.
%
% ..

include = true(size(binary_outcome));
threshold_type = 'Optimal overall accuracy';
class_thr = [];
color = [.2 .2 .2];
plotmethod = 'deciles'; %'npoints'; % 'observed';
donormfit = 1;
istwochoice = 0;
reportstats90 = 0;
plothistograms = 0;
writerscoreplus = 0;
doboot = 1;
dobalanced = 0;
doDependent = 0;
doplot = 1;
doOutput = 1;
valuestoevaluate = 'auto';

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            
            case 'include', include = logical(varargin{i+1});
                
            case 'threshold'
                class_thr = varargin{i + 1};
                threshold_type = 'A priori threshold';
                
            case {'threshold_type'}
                threshold_type = varargin{i + 1}; varargin{i + 1} = [];
                
            case {'Optimal overall accuracy', 'Optimal balanced error rate', 'Minimum SDT bias'}
                threshold_type = varargin{i};
                
            case 'color'
                color = varargin{i + 1}; varargin{i + 1} = [];
                
            case 'plotmethod'
                plotmethod = varargin{i + 1}; varargin{i + 1} = [];
                
            case 'valuestoevaluate'
                valuestoevaluate = varargin{i + 1}; varargin{i + 1} = [];
                
            case 'nonormfit'
                donormfit = 0;
                
            case {'twochoice', 'forcedchoice', 'pairedobservations'}
                istwochoice = 1;
                
            case 'reportstats90', reportstats90 = 1;
                
            case 'plothistograms', plothistograms = 1;
                
            case 'boot', doboot = 1;
            case 'noboot', doboot = 0;
                
            case 'balanced'
                dobalanced = 1;
                error('THIS OPTION IS NOT IMPLEMENTED CORRECTLY; SEE CODE. CONSIDER USING ''Optimal balanced error rate'' OPTION');
                
            case 'noplot'
                doplot = 0;
                
            case 'nooutput'
                doOutput = 0;
                
            case 'dependent'
                doDependent = 1;
                subject_id = varargin{i + 1};
                if(~ismatrix(subject_id) || ~isnumeric(subject_id) || length(input_values)~=length(subject_id))
                    error('Make sure ''dependent'' flag is followed by valid subject_id vector')
                end
                
                disp('ROC for single interval classification of paired observations.')
                
            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end

if length(include) ~= length(binary_outcome) || ~any(include)
    error('Problem with include variable')
end

input_values = input_values(include);
binary_outcome = logical(binary_outcome(include));

% Deal with paired observations
% -------------------------------------------------------------------------

if istwochoice
    % Adjust input scores to reflect differences from the mean of each pair
    % This allows us to do forced-choice classification for pairs based on
    % which is higher. The higher one will always be above the mean.
    % The threshold used here should be zero.
    
    if doOutput
        disp('ROC for two-choice classification of paired observations.')
        disp('Assumes pos and null outcome observations have the same subject order.')
        disp('Using a priori threshold of 0 for pairwise differences.')
    end
    
    meanscores = (input_values(binary_outcome) + input_values(~binary_outcome)) ./ 2;
    
    input_values(binary_outcome) = input_values(binary_outcome) - meanscores;
    input_values(~binary_outcome) = input_values(~binary_outcome) - meanscores;
    
    threshold_type = 'A priori threshold';
    class_thr = 0;
end



% Get ROC values and main results
% -------------------------------------------------------------------------
if ischar(valuestoevaluate) && strcmp(valuestoevaluate, 'auto') % default
    
len = length(binary_outcome);
newlen = min(50*len, 500); % max 500 values to evaluate.

thr = linspace(min(input_values), max(input_values), newlen); %min(input_values):.01:max(input_values);

else
    thr = valuestoevaluate;
end

% AUC will be replaced with theoretical value in 2AFC case!
[dummy, tpr, fpr, auc, c_bias] = roc_calc(input_values, binary_outcome, thr);

% count signal present and signal absent
n1 = sum(~isnan(input_values(binary_outcome)) & ~isinf(input_values(binary_outcome)));
n0 = sum(~isnan(input_values(~binary_outcome)) & ~isinf(input_values(~binary_outcome)));

% Get criterion threshold for final stats output
% -------------------------------------------------------------------------

switch threshold_type
    case 'Optimal balanced error rate'
        
        % For 'Optimal balanced error rate', we don't want the threshold
        % that optimizes overall sensitivity+specificity; we want the threshold that
        % minimizes the discrepancy between sensitivity and specificity.
        % This is because with unbalanced classes, one can always choose
        % "yes" or "no" and the threshold will be at the bounds, yielding a
        % sens/spec of 0% and 100% or vice versa.
        
        % old:
%         avg = mean([tpr; 1-fpr]);
%         [dummy, wh] = max(avg);
%         class_thr = thr(wh);
        
        % new:
        tfdiff = diff([tpr; (1-fpr)]);
        [dummy, wh] = min(abs(tfdiff));
        class_thr = thr(wh);
        
        dobalanced = 1;
        
    case 'A priori threshold'
        
        % Report a priori threshold
        %         wh = find(thr <= class_thr);
        %         wh = wh(end);
        
    case 'Optimal overall accuracy'
        
        ncorrt = tpr .* n1;
        ncorrf = (1 - fpr) .* n0;
        
        mysum = sum([ncorrt; ncorrf]);
        [dummy, wh] = max(mysum);
        class_thr = thr(wh);
        
    case 'Minimum SDT bias'
        [dummy, wh] = min(abs(c_bias));
        class_thr = thr(wh);
        
    otherwise
        error('Unknown threshold type')
end


% Save stuff
% -------------------------------------------------------------------------
ROC.baserate = sum(binary_outcome) ./ length(binary_outcome);

ROC.all_vals.thr = thr;

ROC.class_threshold = class_thr;
ROC.sensitivity = sum(input_values(binary_outcome) >= class_thr) ./ n1; % go back to original thresh for precision when using specific input thresh
ROC.specificity = 1 - ( sum(input_values(~binary_outcome) >= class_thr) ./ n0 );
ROC.AUC = auc;
ROC.AUC_descrip = 'Numerically integrated, nonparametric area under curve';

% vectors of true/false positives and accuracy
if istwochoice
    falseneg = (input_values <= class_thr & binary_outcome); % Wani added this line to fix an error when two values are same (in the two-choice test)
else
    falseneg = (input_values < class_thr & binary_outcome);
end
falsepos = (input_values >= class_thr & ~binary_outcome);
misclass =  falseneg | falsepos;
truepos = binary_outcome & ~misclass;
trueneg = ~binary_outcome & ~misclass;

ROC.threshold_type_for_misclass = threshold_type;


if istwochoice
    % Reshape to reflect pairs for stats/output
    sz = [length(misclass) ./ 2 2];
    
    if any(sz ~= round(sz))
        disp('Two-choice classification assumes you enter paired observations in order:')
        disp('signal present for obs (1:n) followed by signal absent for (1:n) in the same order or vice versa.')
        error('The input has the wrong size.')
    end
    
    truepos = truepos(binary_outcome);
    trueneg = trueneg(~binary_outcome);
    falseneg = falseneg(binary_outcome);
    falsepos = falsepos(~binary_outcome);
    misclass = falsepos | falseneg;
    
    if any(truepos & falsepos) || any(trueneg & falseneg)
        disp('Two-choice classification assumes you enter paired observations in order:')
        disp('signal present for obs (1:n) followed by signal absent for (1:n) in the same order or vice versa.')
        error('Inconsistent output: The observations are likely not in order.')
    end
    
end

% Stuff for figuring out which points are misclassified

ROC.observations.truepos = truepos;
ROC.observations.trueneg = trueneg;
ROC.observations.falseneg = falseneg;
ROC.observations.falsepos = falsepos;
ROC.observations.misclass = misclass;

ROC.PPV = sum(ROC.observations.truepos) ./ (sum(ROC.observations.truepos) + sum(ROC.observations.falsepos));

% Accuracy stats
if ~dobalanced
    accuracy = 1 - (sum(misclass) ./ length(misclass));
else
    accuracy = (sum(truepos)/(sum(truepos) + sum(falseneg)) + sum(trueneg)/(sum(trueneg) + sum(falsepos)))/2;
end

ROC.accuracy = accuracy;

if ~doDependent
    
    % Threshold optimization correction for binomial P-values
    % Chance values can be tricky, and optimizing a threshold can be
    % circular. For example, if only 33% of observations are true
    % positives, and the classifier always guesses "no", it will be right
    % 66% of the time, which is better than the 50% expected under random guessing.
    % We want to build in a correction for threshold selection.
    % The expectation under chance is that we'll adopt a strategy that
    % gives us at least the base-rate of true positive or true negative
    % results. For 33% true-pos, this would be 66%. 
    % We take the max p-value of the difference from 0.5 and the baserate
    % or 1 - baserate.
    % Not implemented for 'dependent' binomial test.

    RES = binotest(double(~misclass), .5);
    RES2 = binotest(double(~misclass), ROC.baserate);
    RES3 = binotest(double(~misclass), 1 - ROC.baserate);
    RES.p_val = max([RES.p_val RES2.p_val RES3.p_val]);
    
else %Run hierarchical version of binomial test
    
    %Create new subject matrix
    subID = unique(subject_id);
    for i = 1:length(subID)
        sdat(i,:) = ~misclass(subject_id == subID(i));
    end
    [RES1, RES2, RES3, RES4] = binotest_dependent(sdat,.5);
    RES = RES4;
    
end
ROC.N = RES.n;
ROC.accuracy_p = RES.p_val;
ROC.accuracy_se = RES.SE;

% Stuff for re-plotting full ROC
ROC.all_vals.tpr = tpr;
ROC.all_vals.fpr = fpr;

% -------------------------------------------------------------------------
% Plot stuff
%
% Get ROC values for plot - default is to plot 10 points
% -------------------------------------------------------------------------

switch plotmethod
    case 'deciles'
        plotthr = prctile(input_values, [10 20 30 40 50 60 70 80 90]);
        [dummy, plottpr, plotfpr] = roc_calc(input_values, binary_outcome, plotthr);
        plotsymbol = 'o';
        linewid = 2;
        
    case 'observed'
        plotthr = thr;
        plottpr = tpr;
        plotfpr = fpr;
        plotsymbol = '-';
        linewid = 3;
end


if doplot
    han = plot(plotfpr, plottpr, plotsymbol, 'Color', color, 'LineWidth', linewid);
    
    set(gca, 'FontSize', 32, 'XLim', [-.02 1.02], 'XTick', 0:.2:1, 'YTick', 0:.2:1);
    xlabel('(1 - Specificity)');
    ylabel('Sensitivity')
    
    
    inline_plothistograms();
end

% effect size
% -------------------------------------------------------------------------

if istwochoice
    % Two-alternative forced choice
    % -------------------------------------------------------------
    
    diffscores = input_values(binary_outcome) - input_values(~binary_outcome);
    
    meandiff = mean(diffscores);
    d = meandiff ./ std(diffscores);
    
    % forced-choice is variance of the difference, which is 2 * pooled variance
    % std of difference therefore = pooledsd * sqrt(2), and pooledsd = std(diffs)/sqrt(2)
    pooledsd = std(diffscores) ./ sqrt(2);
    
    d_a_model = meandiff ./ pooledsd; % estimate of what d_a would be for single-interval
    
    % From Lew Harvey's notes - this should be the "observed" d_a based on
    % empirical accuracy. But it will be inaccurate as accuracy approaches
    % 1. Use "model" because it's closer to the data, no norminv inaccuracy
    d_a_obs = sqrt(2) * norminv(ROC.accuracy);
    
    if donormfit
        % standard equal-variance signal detection
        x = [-3:.1:3];
        tprn = 1 - normcdf(x, d, 1);
        fprn = 1 - normcdf(x, -d, 1);
        hold on;
        if doplot
            han = [han plot(fprn, tprn, '-', 'Color', color, 'LineWidth', linewid)];
        end
    end
    
    aucn = calc_auc(fprn, tprn);
    expected_acc = 1 - normcdf(0, d, 1); % expected to be = to the AUC!
    
    ROC.Gaussian_model.type = 'Two-alternative forced choice';
    ROC.Gaussian_model.names = {'diff_scores'};
    ROC.Gaussian_model.means = meandiff;
    ROC.Gaussian_model.n = length(diffscores);
    ROC.Gaussian_model.sd = std(diffscores);
    ROC.Gaussian_model.d_a = d_a_model;
    ROC.Gaussian_model.d_a_descrip = '''Observed'' d_a based on empirical accuracy.';
    ROC.Gaussian_model.sensitivity = expected_acc;
    ROC.Gaussian_model.specificity = expected_acc;
    
    % PPV is just the accuracy, too.
    ROC.Gaussian_model.PPV = ROC.Gaussian_model.sensitivity ./ (ROC.Gaussian_model.sensitivity + 1 - ROC.Gaussian_model.specificity);
    
    ROC.Gaussian_model.AUC_numerical_integration = aucn;
    ROC.Gaussian_model.AUC = normcdf(d_a_model/sqrt(2)); % should just invert the accuracy equation. ...and yes, it does. same as accuracy.
    
else
    % single-interval
    % ------------------------------------------------------------
    %     N1 = length(input_values(binary_outcome)); done above
    %     N2 = length(input_values(~binary_outcome));
    
    meanpres = mean(input_values(binary_outcome));
    meanabs = mean(input_values(~binary_outcome));
    
    v1 = var(input_values(binary_outcome));
    v0 = var(input_values(~binary_outcome));
    
    pooledsd = sqrt((v1.*(n1-1) + v0.*(n0-1)) ./ (n1 + n0 - 2));
    
    d = (meanpres - meanabs) ./ pooledsd;
    
    zpres = meanpres ./ pooledsd;
    zabs = meanabs ./ pooledsd;
    
    if donormfit
        % integrate across PDFs for each of signal absent, signal present distributions
        x = [zabs-3:.1:zpres+3]; % relative to null dist
        tprn = 1 - normcdf(x, zpres, 1);
        fprn = 1 - normcdf(x, zabs, 1);
        hold on;
        if doplot
            han = [han plot(fprn, tprn, '-', 'Color', color, 'LineWidth', linewid)];
        end
    end
    
    aucn = calc_auc(fprn, tprn);
    
    ROC.Gaussian_model.type = 'Single-interval';
    ROC.Gaussian_model.names = {'sig. abs.' 'sig. pres.'};
    ROC.Gaussian_model.means = [meanabs meanpres];
    ROC.Gaussian_model.n = [n0 n1];
    ROC.Gaussian_model.sd = [sqrt(v0) sqrt(v1)];
    ROC.Gaussian_model.mean_diff = meanpres - meanabs;
    ROC.Gaussian_model.pooledsd = pooledsd;
    ROC.Gaussian_model.d_a = d;
    ROC.Gaussian_model.sensitivity = 1 - normcdf(class_thr, meanpres, sqrt(v1));
    ROC.Gaussian_model.specificity = normcdf(class_thr, meanabs, sqrt(v0));
    ROC.Gaussian_model.PPV = ROC.Gaussian_model.sensitivity ./ (ROC.Gaussian_model.sensitivity + 1 - ROC.Gaussian_model.specificity);
    
    ROC.Gaussian_model.AUC_numerical_integration = aucn;
    ROC.Gaussian_model.AUC = normcdf(d/sqrt(2));
    
    
end

if doplot
    ROC.line_handle = han;
end

% Boostrap, if asked for
% -------------------------------------------------------------------------
if doboot
    
    [ci, names] = roc_boot(input_values, binary_outcome, ROC.class_threshold, 0);
    
    ROC.sensitivity_ci = ci{1};
    ROC.specificity_ci = ci{2};
    ROC.PPV_ci = ci{3};
    
end




% fprintf('\nROC_PLOT Output: %s, %s\n', ROC.Gaussian_model.type, threshold_type);
%
% fprintf('  Nonparametric AUC:\t%3.2f\tParametric d_a:\t%3.2f\n', ROC.AUC, ROC.Gaussian_model.d_a);
%
% fprintf('  Threshold:\t%3.2f\tSens:\t%3.0f%%\tSpec:\t%3.0f%%\tPPV:\t%3.0f%%\n', ...
%     ROC.class_threshold, 100*ROC.sensitivity, 100*ROC.specificity, 100*ROC.PPV);
%
% fprintf('  Accuracy:\t%3.0f%% +- %3.1f%% (SE), P = %3.6f\n', ...
%     100*ROC.accuracy, 100*ROC.accuracy_se, ROC.accuracy_p);

if doOutput
    % Single line format
    fprintf('\nROC_PLOT Output: %s, %s\n', ROC.Gaussian_model.type, threshold_type);
    
    if doboot
        fprintf('Threshold:\t%3.2f\tSens:\t%3.0f%% CI(%.0f%%-%.0f%%)\tSpec:\t%3.0f%% CI(%.0f%%-%.0f%%)\tPPV:\t%3.0f%% CI(%.0f%%-%.0f%%)\t', ...
            ROC.class_threshold, 100*ROC.sensitivity, 100*ROC.sensitivity_ci, 100*ROC.specificity, 100*ROC.specificity_ci, 100*ROC.PPV, 100*ROC.PPV_ci);
    else
        fprintf('Threshold:\t%3.2f\tSens:\t%3.0f%%\tSpec:\t%3.0f%%\tPPV:\t%3.0f%%\t', ...
            ROC.class_threshold, 100*ROC.sensitivity, 100*ROC.specificity, 100*ROC.PPV);
    end
    
    fprintf('Nonparametric AUC:\t%3.2f\tParametric d_a:\t%3.2f\t', ROC.AUC, ROC.Gaussian_model.d_a);
    
    fprintf('  Accuracy:\t%3.0f%% +- %3.1f%% (SE), P = %3.6f\n', ...
        100*ROC.accuracy, 100*ROC.accuracy_se, ROC.accuracy_p);
end

% Report stats (max sens) at 90% specificity
if reportstats90
    
    ROC = report_at_threshold_spec();
    
end


% Calculate statistics by criterion threshold value
% Parallels output of RSCOREplus by Lew Harvey
% -------------------------------------------------------------------------

bins = [-Inf plotthr Inf];   % bin edges: EDGES(k) <= X(i) < EDGES(k+1)
npoints = length(bins);
s0 = histc(input_values(~binary_outcome), bins)';
s1 = histc(input_values(binary_outcome), bins)';

for i = 1:length(bins)
    %bin_names{i} = sprintf('R_%3.1f', bins(i));
    bin_names{i} = sprintf('R_%d', i);
end

ROC.Binned_output.s0 = s0;
ROC.Binned_output.s1 = s1;
ROC.Binned_output.bin_edges = bins;
ROC.Binned_output.bin_names = bin_names;

ROC.Binned_output.spec_bins = cumsum(s0) ./ sum(s0);
ROC.Binned_output.sens_bins = 1 - (cumsum(s1) ./ sum(s1));
ROC.Binned_output.PPV_bins = (ROC.Binned_output.sens_bins .* sum(s1)) ./ (ROC.Binned_output.sens_bins .* sum(s1) + (1 - ROC.Binned_output.spec_bins) .* sum(s0));




if writerscoreplus
    
    write_rscoreplus_input()
    
end

% INLINE FUNCTIONS
% -------------------------------------------------------------------------



    function inline_plothistograms()
        
        if plothistograms
            
            if ~isempty(findobj(get(0, 'Children'), 'type', 'figure'))
                returncurrfig = 1;
            end
            
            if returncurrfig
                figh = gcf;
            end
            
            % plot histograms...
            
            create_figure('distributions');
            h = histfit(input_values(binary_outcome), 20);
            hbar = get(h(1), 'Children');
            set(hbar, 'FaceAlpha', .3, 'FaceColor', [0 0 1], 'EdgeColor', 'none');
            % changing in diff versions of matlab...
            if isempty(hbar)
                set(h(1), 'FaceAlpha', .3, 'FaceColor', [0 0 1], 'EdgeColor', 'none');
            end
            
            set(h(2), 'Color', [0 0 1]);
            
            hold on;
            h = histfit(input_values(~binary_outcome), 20);
            hbar = get(h(1), 'Children');
            set(hbar, 'FaceAlpha', .3, 'FaceColor', [0 1 0], 'EdgeColor', 'none');
            
            % changing in diff versions of matlab...
            if isempty(hbar)
                set(h(1), 'FaceAlpha', .3, 'FaceColor', [0 1 0], 'EdgeColor', 'none');
            end
            
            set(h(2), 'Color', [0 .4 0]);
            
            h = plot_vertical_line(class_thr);
            set(h, 'Color', 'k', 'LineStyle', '-', 'LineWidth', 2)
            
            if returncurrfig
                figure(figh)
            end
            
        end
        
        
    end % inline



    function write_rscoreplus_input()
        
        fname = 'roc_rscoreplus_input.txt';
        disp('Writing RScorePlus input .txt file: %s\n', fname);
        
        fid = fopen(fullfile(pwd, fname), 'w+');
        if fid == -1, disp('Cannot open rscoreplus_input.txt file. Skipping.'); return, end
        
        fprintf(fid, 'Heading\nroc_plot.m_output\n');
        
        % parameters => see rscore plus help on lew harvey's website
        fprintf(fid, '%01d\t%01d\t%01d\t%01d\t%01d\t%01d\t', length(bins), 2, 1, 0, 0, 1);
        
        if istwochoice, paradigmstr = 'MAFC'; else paradigmstr = 'SINT'; end
        fprintf(fid, '%s\n', paradigmstr);
        
        fprintf(fid, 'labels\t');
        fprintf(fid, '%s\t', bin_names{:});
        fprintf(fid, '\n');
        
        fprintf(fid, 's0\t');
        fprintf(fid, '%d\t', s0);
        fprintf(fid, '\n');
        
        fprintf(fid, 's1\t');
        fprintf(fid, '%d\t', s1);
        fprintf(fid, '\n');
        
        fprintf(fid, '%3.1f\t%3.1f\t%3.1f\t%3.1f\n', 0, 1, .5, 1); % params of reference dist?
        fprintf(fid, '%01d\t%01d\t%01d\t%01d\n', 0, 0, 1, 1); % params of reference dist?
        
        fprintf(fid, 'end of file\n-1\n');
        fclose(fid);
        
    end % inline



    function report_at_threshold_spec
        wh = find(fpr <= .1);
        [mymax, wh2] = max(tpr(wh));
        wh = wh(wh2);
        if ~isempty(wh), wh = wh(1); else wh = NaN; end
        
        npos = sum(tpr(wh) .* binary_outcome);
        nfp = sum(fpr(wh) .* ~binary_outcome);
        ppv = npos ./ (npos + nfp);
        
        if isnan(wh)
            fprintf('At 90+%% Spec: Thresh = %3.2f, Sensitivity = %3.0f%%, Specificity = %3.0f%%, PPV = %3.0f%%\n', NaN, NaN, NaN, NaN);
            
            [ROC.thresh_90_percent_spec, ROC.sens_90_percent_spec] = deal(NaN);
        else
            
            fprintf('At 90+%% Spec: Thresh = %3.2f, Sensitivity = %3.0f%%, Specificity = %3.0f%%, PPV = %3.0f%%\n', thr(wh), tpr(wh)*100, (1-fpr(wh))*100, ppv * 100);
            
            ROC.thresh_90_percent_spec = thr(wh);
            ROC.sens_90_percent_spec = tpr(wh);
        end
        
        wh = find(fpr <= .05);
        [mymax, wh2] = max(tpr(wh));
        wh = wh(wh2);
        if ~isempty(wh), wh = wh(1); else wh = NaN; end
        
        if isnan(wh)
            fprintf('At 95+%% Spec: Thresh = %3.2f, Sensitivity = %3.0f%%, Specificity = %3.0f%%, PPV = %3.0f%%\n', NaN, NaN, NaN, NaN);
            
            [ROC.thresh_95_percent_spec, ROC.sens_95_percent_spec] = deal(NaN);
            
        else
            fprintf('At 95+%% Spec: Thresh = %3.2f, Sensitivity = %3.0f%%, Specificity = %3.0f%%, PPV = %3.0f%%\n', thr(wh), tpr(wh)*100, (1-fpr(wh))*100, ppv * 100);
            
            ROC.thresh_95_percent_spec = thr(wh);
            ROC.sens_95_percent_spec = tpr(wh);
        end
        
    end % inline


end % main function

%
% function [xvals, tpr, fpr] = roc_calc(input_vals, binary_outcome, xvals)
% % Calculate Receiver Operating Characteristic plot (ROC) given P-values
% %
% %  function [xvals, tpr, fpr] = roc_calc(input_vals or input values, input_vals, [treshold vals to assess])
% %
% % Modified from roc_calc.m in scnlab tools
% % Tor Wager, 2012
%
%
% [tpr, fpr] = deal(zeros(size(xvals)));
%
% indx = 1;
% for x = xvals
%     wh = input_vals >= x;
%
%     tpr(indx) = sum(wh(binary_outcome)) ./ sum(binary_outcome);
%     fpr(indx) = sum(wh(~binary_outcome)) ./ sum(~binary_outcome);
%
%     indx = indx + 1;
% end
%
% end % function




function auc = calc_auc(fpr, tpr)

[u, wh] = unique(fpr);
u2 = tpr(wh);

% fix for AUC = 1 if no overlap; triangle method not perfectly accurate
% here.
if any(u == 0 & u2 == 1), auc = 1; return, end

for i = 2:length(u)
    
    xdiff = u(i) - u(i - 1);
    ydiff = u2(i) - u2(i - 1);
    a(i) = xdiff * u2(i - 1) + xdiff * ydiff / 2;  % area of rect + area of triangle
    
end


auc = sum(a);

end


