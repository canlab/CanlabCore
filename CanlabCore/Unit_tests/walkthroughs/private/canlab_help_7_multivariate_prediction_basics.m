%% Multivariate prediction of a continuous outcome

%%
% This walkthrough uses the |predict()| method for fmri_data objects to
% predict a continuous outcome using cross-validated principal component
% regression (PCR / LASSO-PCR). 
%
% The |predict()| method can run multiple algorithms with a range of options.
% The main ones used in the CANlab are SVM for two-choice classification and PCR for
% regression.

%% About the pain dataset
% -------------------------------------------------------------------------
%
% The dataset contains data from 33 participants, with brain responses to six levels
% of heat (non-painful and painful). Each image is the average over several
% (4-8) trials of heat delivered at a single stimulus intensity, ranging
% from 44.3 - 49.3 degrees C in one-degree increments. Each image is also
% paired with an average reported pain value for that set of trials, rated
% immmediately after heat experience. 
%
% This dataset is interesting for mixed-effects and predictive analyses, as
% it has both within-person and between-person sources of variance.  
% 
% Aspects of this data appear in these papers:
% Wager, T.D., Atlas, L.T., Lindquist, M.A., Roy, M., Choong-Wan, W., Kross, E. (2013). 
% An fMRI-Based Neurologic Signature of Physical Pain. The New England Journal of Medicine. 368:1388-1397
% (Study 2)
%
% Woo, C. -W., Roy, M., Buhle, J. T. & Wager, T. D. (2015). Distinct brain systems 
% mediate the effects of nociceptive input and self-regulation on pain. PLOS Biology. 13(1): 
% e1002036. doi:10.1371/journal.pbio.1002036
%
% Lindquist, Martin A., Anjali Krishnan, Marina López-Solà, Marieke Jepma, Choong-Wan Woo, 
% Leonie Koban, Mathieu Roy, et al. 2015. ?Group-Regularized Individual Prediction: 
% Theory and Application to Pain.? NeuroImage. 
% http://www.sciencedirect.com/science/article/pii/S1053811915009982.
%
%% Load dataset with images and print descriptives
% This dataset is shared on figshare.com, under this link:
% https://figshare.com/s/ca23e5974a310c44ca93
%
% Here is a direct link to the dataset file with the fmri_data object:
% https://ndownloader.figshare.com/files/12708989
%
% The key variable is image_obj
% This is an fmri_data object from the CANlab Core Tools repository for neuroimaging data analysis.
% See https://canlab.github.io/
%
% image_obj.dat contains brain data for each image (average across trials)
% image_obj.Y contains pain ratings (one average rating per image)
%
% image_obj.additional_info.subject_id contains integers coding for which
% Load the data file, downloading from figshare if needed

fmri_data_file = which('bmrk3_6levels_pain_dataset.mat');

if isempty(fmri_data_file)
    
    % attempt to download
    disp('Did not find data locally...downloading data file from figshare.com')
    
    fmri_data_file = websave('bmrk3_6levels_pain_dataset.mat', 'https://ndownloader.figshare.com/files/12708989');
    
end

load(fmri_data_file);

descriptives(image_obj);

%% Collect some variables we need

% subject_id is useful for cross-validation
%
% this used as a custom holdout set in fmri_data.predict() below will
% implement leave-one-subject-out cross-validation

subject_id = image_obj.additional_info.subject_id;

% ratings: reconstruct a subjects x temperatures matrix
% so we can plot it
%
% The command below does this only because we have exactly 6 conditions
% nested within 33 subjects and no missing data. 

ratings = reshape(image_obj.Y, 6, 33)';

% temperatures are 44 - 49 (actually 44.3 - 49.3) in order for each person.

temperatures = image_obj.additional_info.temperatures;

%% Plot the ratings
% -------------------------------------------------------------------------  

create_figure('ratings');
hold on; plot(ratings', '-', 'Color', [.7 .7 .7], 'LineWidth', .5);
lineplot_columns(ratings, 'color', [.7 .3 .3], 'markerfacecolor', [1 .5 0]);
xlabel('Temperature');
ylabel('Rating');
set(gca, 'XTickLabel', [44.3:49.3], 'FontSize', 18);

%% Prediction options and choices
% There are many choices for algorithm and process/parameter tuning 
% Here, we describe a constrained set of some of the most common choices
%
% Base model: Linear, whole brain prediction
% 
% Outcome distribution:
% Continuous: Regression (LASSO-PCR) Categorical/2 choice: SVM, logistic regression
%
% Model Options:
% - Feature selection (part of the brain)
% - Feature integration (new features)
%
% Testing options:
% - Input data (what level)
% - Cross-validation
% - What testing metric (classification accuracy? RMSE?) and at what level
% (single-trial, condition averages, one-map-per-subject)
%
% Inference options:
% - Component maps, voxel-wise maps
% - Thresholding (Bootstrapping, Permutation)
% - Global null hypothesis testing / sanity checks (Permutation)

%% Run the Base model
%
% *relevant functions:*
% predict method (fmri_data)                            
% predict_test_suite method (fmri_data) 
%

%% 
% *Define the holdout set for cross-validation*
% We want to define custom holdout set.  If we use subject_id, which is a vector of
% integers with a unique integer per subject, then we are doing
% leave-one-subject-out cross-validation. 
%
% Let's build five-fold cross-validation set that leaves out ALL the images
% from a subject together. That way, we are always predicting out-of-sample
% (new individuals). If not, dependence across images from the same
% subjects may invalidate the estimate of predictive accuracy.

holdout_set = zeros(size(subject_id));      % holdout set membership for each image
n_subjects = length(unique(subject_id));
C = cvpartition(n_subjects, 'KFold', 5);
for i = 1:5
    teidx = test(C, i);                     % which subjects to leave out
    imgidx = ismember(subject_id, find(teidx)); % all images for these subjects
    holdout_set(imgidx) = i;
end

%% 
% *Run the prediction*

algoname = 'cv_lassopcr'; % cross-validated penalized regression. Predict pain ratings

[cverr, stats, optout] = predict(image_obj, 'algorithm_name', algoname, 'nfolds', holdout_set);

%% Plot cross-validated predicted vs. actual outcomes
%
% Critical fields in stats output structure:
% stats.Y = actual outcomes
% stats.yfit = cross-val predicted outcomes
% pred_outcome_r: Correlation between .yfit and .Y
% weight_obj: Weight map used for prediction (across all subjects).
                 
% Continuous outcomes:

create_figure('scatterplots', 1, 2);
plot(stats.yfit, stats.Y, 'o')
axis tight
refline
xlabel('Predicted pain');
ylabel('Observed pain');

% Consider that each subject has their own series of multiple observations.
% Separate subjects into multiple cells, which will be plotted as separate
% individual lines:

n = max(subject_id);

for i = 1:n
    YY{i} = stats.Y(subject_id == i);
    Yfit{i} = stats.yfit(subject_id == i);
end

subplot(1, 2, 2);
line_plot_multisubject(Yfit, YY, 'nofigure');
xlabel('Predicted pain');
ylabel('Observed pain');
axis tight

% Yfit is a cell array, one cell per subject, with fitted values for that subject

%% Summarize within-person classification accuracy
% 
% Classification can turn a continuous outcome, such as the magnitude of
% activity in an integrated model/pattern, into a binary choice (Type A or
% B, Yes/No, Pain/No pain). Classification accuracy is an easy-to-interpret
% metric, which makes it appealing to report. 
% 
% But though it seems easy to understand, it is perhaps tricker to
% interpret than we might think.
%
% Classification accuracy is not one metric, but a family of them. The accuracy
% you obtain depends on what type of classification you're doing, and what the
% units of analysis are. For example, is your model based on a single person 
% (individualized for that person), or based on other participants' data? 
% Are you classifying images corresponding to a single trial, or an average across a group of trials? 
% Do you have knowledge about which conditions/images belong to the same person, or not?
% Here, we will:
% (a) use a model based on other people's data (cross-validated across participants)
% (b) classify images summarizing a group of trials, and
% (c) assume we have 2 images belonging to the same person, and we classify
% which is which
%
% On this last point, there are two basic senses of classification, which assume different 
% levels of background knowledge:
% - In _single-interval classification_, we do not know which images belong to which participants,
%   and we are classifying a single image as Type A or B based on whether
%   the response in the model is above or below a threshold. 
% - In _forced-choice_ classification, we have two images that we know come
%   from the same person, and we know one is Type A and one is Type B...or,
%   equivalently, that one is "more A" than the other (e.g., more painful).
%   The classifier tries to pick out which is the best answer.
%
% Forced-choice classification is an inherently easier problem, with higher
% resulting accuracy, because (a) each participant serves as their own control
% (many sources of noise are canceled out), and (b) you assume more
% background knowledge, conferring an advantage. 
%
% In this case, we want to summarize each subject with two images: One high
% pain, and one low pain. We'll use stimulus intensity as the variable to
% classify here, and ask if we can predict which intensity is highest given
% pairs of images one degree apart for each participant.

% If the stim intensity increases by 1 degree, how often do we correctly
% predict an increase?
diffs = cellfun(@diff, Yfit, 'UniformOutput', false); % successive increases in predicted values with increases in stim intensity

% Subjects (rows) x comparisons (cols), each comparison 1 degree apart
diffs = cat(2, diffs{:})';

acc_per_subject = sum(diffs > 0, 2) ./ size(diffs, 2);
acc_per_comparison = (sum(diffs > 0) ./ size(diffs, 1))';
cohens_d_per_comparison = (mean(diffs) ./ std(diffs))';

fprintf('Mean acc is : %3.2f%% across all successive 1 degree increments\n', 100*mean(acc_per_subject));

% Create and display a table:

comparison = {'45-44' '46-45' '47-46' '48-47' '49-48'}';
results_table = table(comparison, acc_per_comparison, cohens_d_per_comparison);
disp(results_table)

% 45-44 degrees, 46-45, 47-46, etc.

%% Visualize the classifier weight map
%
% The weight map is stored in |stats.weight_obj|, as a statistic_image
% class object. If we have requested bootstrapping, the image will be
% thresholded (0.05 uncorrected by default). Otherwise, it will be unthresholded,
% and we'll see weights everywhere within the analysis mask.

w = stats.weight_obj;  % This is an image object that we can view, etc.

orthviews(w)

create_figure('montage');
o2 = montage(w);
o2 = title_montage(o2, 5, 'Predictive model weights');


%% Bootstrap weights: Get most reliable weights and p-values for voxels
%
% Here is an exercise: 
% Re-run the predictive model, adding a 'bootstrap' flag. 
% For a final analysis, at least 5K bootstrap samples is a good idea.
% For now, try it with just 1,000 bootstrap samples (and be prepared to
% wait a bit). What does the resulting weight map look like?
%
% (This is not run in the example; see |help fmri_data.predict| for how to
% add the bootstrap option.)

%% Try normalizing/scaling data
%
% One of the biggest sources of noise can be whole-image (or
% large-spatial-scale shifts) in image values across participants. These can
% be apparent even in contrast or "subtraction" images where various
% sources of noise and artifacts are supposed to have been "subtracted
% out". This may result from different chance correlations between task
% regressors and physiological noise or artifacts/outliers across
% individuals. 
%
% Rescaling the data can remove whole-brain signal. This should be used
% with caution, as you're also removing signal from the images, and the
% resulting weight maps should be interpreted as effects *relative* to the mean
% signal across the image. i.e., a task may be associated with whole-brain
% increases, and *relative* decreases in a local area that appear only
% after the global increase has been accounted for.
%
% Another reason to rescale images is to ask whether the *pattern* of
% activity across the brain (or within a local area) is predictive of task
% state, after removing the overall intensity of activation. This relates
% to pattern information.
%
% Here, we'll rescale the images and ask if this improves our
% cross-validated predictions or not. If it does, the global signal we've
% removed is likely more noise than signal. If not, the global signal may
% contain information about pain.

% z-score every image, removing the mean and dividing by the std.
% Remove some scaling noise, but also some signal...
 image_obj = rescale(image_obj, 'zscoreimages');

% re-run prediction and plot

[cverr, stats, optout] = predict(image_obj, 'algorithm_name', algoname, 'nfolds', holdout_set);

create_figure('scatterplots');

for i = 1:n
    YY{i} = stats.Y(subject_id == i);
    Yfit{i} = stats.yfit(subject_id == i);
end

line_plot_multisubject(Yfit, YY, 'nofigure');
xlabel('Predicted pain');
ylabel('Observed pain');
axis tight

fprintf('Overall prediction-outcome correlation is %3.2f\n', stats.pred_outcome_r)

figure; o2 = montage(stats.weight_obj);
o2 = title_montage(o2, 5, 'Predictive model weights');

%% Try selecting an a priori 'network of interest'

% re-load
load('bmrk3_6levels_pain_dataset.mat', 'image_obj')

% load an a priori pain-related mask and apply it

mask = fmri_data(which('pain_2s_z_val_FDR_05.img.gz')); % in Neuroimaging_Pattern_Masks repository

image_obj = apply_mask(image_obj, mask);

% show mean data with mask
m = mean(image_obj);
figure; o2 = montage(m);
o2 = title_montage(o2, 5, 'Mean image data, masked');


%% re-run prediction and plot

[cverr, stats, optout] = predict(image_obj, 'algorithm_name', algoname, 'nfolds', holdout_set);

create_figure('scatterplots');

for i = 1:n
    YY{i} = stats.Y(subject_id == i);
    Yfit{i} = stats.yfit(subject_id == i);
end

line_plot_multisubject(Yfit, YY, 'nofigure');
xlabel('Predicted pain');
ylabel('Observed pain');
axis tight

fprintf('Overall prediction-outcome correlation is %3.2f\n', stats.pred_outcome_r)

figure; o2 = montage(stats.weight_obj);
o2 = title_montage(o2, 5, 'Predictive model weights');

%% Other ideas
% There are many ways to build from this basic analysis.

% Data
% ------------
% Examine/clean outliers
% Transformations of data: Scaling/centering to remove global mean,
% ventricle mean
% Identify unaccounted sources of variation/subgroups (experimenter sex,
% etc.)

% Features
% ------------
% A priori regions / meta-analysis
% Feature selection in cross-val loop: Univariate, recursive feature elim

