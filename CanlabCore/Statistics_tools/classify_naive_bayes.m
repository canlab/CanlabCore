function varargout = classify_naive_bayes(meth, varargin)
% Naive Bayes classifier
%
% :Usage:
% ::
%
%     % set up model structure
%     bayes_model = classify_naive_bayes('setup', Y, Xi, [activation_cutoff, selectivity_cutoff]);
%
% :Inputs:
%
%   **Y:**
%        is full (not sparse); empty features will be eliminated
%
%   **Xi:**
%        is obs x classes, an indicator matrix of 1's and 0's
%
% :TEST: test classifier; make a prediction about classes from data
% ::
%
%    [class_est log_joint best_log map p_obs_act_given_class] = classify_naive_bayes('test', Y, bayes_model);
%
% :EVAL: evaluate classification accuracy
% ::
%
%    [prop_correct, confusion_mtx, misclass, prop_correct_by_class, chance, chance_95_ci] = classify_naive_bayes('eval', true_class, class_est, wh_obs);
%
% :APPARENT: apparent classification; with full dataset
% ::
%
%     bayes_model = classify_naive_bayes('apparent', Y, bayes_model);
%
% :XVAL: crossvalidate
% ::
%
%    xval = classify_naive_bayes('xval', Y, Xi);
%
%    bayes_model = classify_naive_bayes('write', bayes_model, Y, volInfo, conditionnames)
%    bayes_model = classify_naive_bayes('write', bayes_model, Y, MC_Setup.volInfo, MC_Setup.Xinms);
%
% To add feature abstraction step within xval:
% ::
%
%    xval = classify_naive_bayes('xval', Y, Xi, 0, .9, .05, 1, volInfo );
%
% :PLOT:
% ::
%
%    classify_naive_bayes('plot', bayes_model, ['pa|t', 'lr', 'lr surface',  'map plot', or 'class plot']);
%
% :Optional Inputs: (any order)
%
% Threshold (abs. value), and colors in cell array
% ::
%
%    classify_naive_bayes('plot', bayes_model, 'lr', .10, {[1 .7 0] [0 0 1]});
%
%
% :Examples:
% ::
%
%    [bayes_model, Y] = classify_naive_bayes('setup', Y, Xi);
%
%    % Test obs. 2
%    tic, [class_est log_joint best_log map] = classify_naive_bayes('test',Y(2,:), bayes_model); toc
%
%    % Get apparent classification rate and look at confusion matrix
%    tic, bayes_model = classify_naive_bayes('apparent', Y, bayes_model); toc
%    bayes_model.apparent.confusion_mtx
%
%    % Cross-validate
%    bayes_model.xval = classify_naive_bayes('xval', Y, Xi);
%    bayes_model.xval.confusion_mtx
%
% :Example 2:
% ::
%
%     Select features, and do apparent and cross-validated classification
%    Y = MC_Setup.unweighted_study_data';
%    wh = sum(Y) > 5;
%    Y = Y(:, wh);
%    whos Y
%    [bayes_model, Y] = classify_naive_bayes('setup', Y, Xi);
%    bayes_model = classify_naive_bayes('apparent', Y, bayes_model);
%    bayes_model.apparent.prop_correct, bayes_model.apparent.confusion_mtx
%    bayes_model.xval = classify_naive_bayes('xval', Y, Xi);
%    bayes_model.xval.prop_correct, bayes_model.xval.confusion_mtx
%
% :Example 3:
%    create_figure('hist'); hist(bayes_model.pa1_given_t, 100);
%    xval = classify_naive_bayes('xval', Y, Xi, .05, .3);
%
% Get results from key regions and run classifier only on those regions:
% ::
%
%     cl = classify_naive_bayes('plot', bayes_model, 'lr', .10, {[1 .7 0] [0 0 1]});
%    [studybyroi,studybyset] = Meta_cluster_tools('getdata',cl{1},dat',volInfo);
%    bayes_model_regions = classify_naive_bayes('setup', studybyroi, Xi);
%    bayes_model_regions = classify_naive_bayes('apparent', studybyroi,bayes_model_regions);
%    disp('Apparent confusion')
%    disp(bayes_model_regions.apparent.confusion_mtx)
%
%    bayes_model_regions.xval = classify_naive_bayes('xval', studybyroi, Xi);
%    disp('Cross-validated confusion')
%    disp(bayes_model_regions.xval.confusion_mtx)
%
%    fprintf('Proportion correct: Apparent: %3.0f%%  Xval: %3.0f%%\n', 100*bayes_model_regions.apparent.prop_correct, 100*bayes_model_regions.xval.prop_correct);
%    fprintf('Proportion correct by class: \t'); fprintf('%3.0f%%\t', 100*bayes_model_regions.xval.prop_correct_by_class);
%    fprintf('\n');

%    sz = cat(1,cl{1}(:).numVox);
%    cl{1}(sz < 10) = [];
%
%    subcl = subclusters_from_local_max(cl{1}, 10);
%
% ..
%    Tor Wager, Sept 2007
% ..

switch meth


    case 'setup'
        Y = full(varargin{1});
        Xi = varargin{2};

        activation_cutoff = [];
        selectivity_cutoff = [];
        return_full_ps = 1;         % for faster processing if not all fields are needed (i.e., in xval)

        if length(varargin) > 2

            activation_cutoff = varargin{3};
            selectivity_cutoff = varargin{4};

        end

        if length(varargin) > 4
            return_full_ps = varargin{5};
        end

        k = 1; % default k regularization param
        if length(varargin) > 5
            k = varargin{6};
        end
        

        % setup
        [Y, whkeep] = eliminate_empty_Y(Y);

        [bayes_model, Y] = classify_naive_bayes_setup(Y, Xi, activation_cutoff, selectivity_cutoff, whkeep, return_full_ps, k);

        gam = 1; % default gam prior weighting param
        if length(varargin) > 6
            gam = varargin{7};
        end
        bayes_model.params.gam = gam;

        varargout{1} = bayes_model;
        varargout{2} = Y;


    case 'apparent'

        Y = full(varargin{1});
        bayes_model = varargin{2};

        bayes_model.true_class = indic2condf(bayes_model.Xi);

        % test apparent classification rate
        [bayes_model.apparent.class_est bayes_model.apparent.log_joint bayes_model.apparent.best_log bayes_model.apparent.map bayes_model.apparent.pact_given_task] = ...
            classify_naive_bayes_test(Y, bayes_model);

        [bayes_model.apparent.prop_correct, bayes_model.apparent.confusion_mtx, bayes_model.apparent.misclass, ...
            bayes_model.apparent.prop_correct_by_class, bayes_model.apparent.chance, bayes_model.apparent.chance_95_ci] = ...
            classify_naive_bayes_eval(bayes_model.true_class, bayes_model.apparent.class_est);


        varargout{1} = bayes_model;



    case 'test'
        Y = full(varargin{1});
        bayes_model = varargin{2};

        [class_est log_joint best_log map p_obs_act_given_class] = classify_naive_bayes_test(Y, bayes_model);

        varargout{1} = class_est;
        varargout{2} = log_joint;
        varargout{3} = best_log;
        varargout{4} = map;
        varargout{5} = p_obs_act_given_class;

    case 'eval'
        bayes_model = varargin{1};
        class_est = varargin{2};
        wh_obs = [];
        if length(varargin) > 2, wh_obs = varargin{3}; end

        [prop_correct, confusion_mtx, misclass, prop_correct_by_class, chance, chance_95_ci] = classify_naive_bayes_eval(bayes_model, class_est, wh_obs);

        varargout{1} = prop_correct;
        varargout{2} = confusion_mtx;
        varargout{3} = misclass;
        varargout{4} = prop_correct_by_class;
        varargout{5} = chance;
        varargout{6} = chance_95_ci;

        
        
        case 'abstract'
          %bayes_model = classify_naive_bayes('abstract', Y, thresh, bayes_model, volInfo);  
            
        Y = full(varargin{1});
        thresh = varargin{2};

        bayes_model = varargin{3};
        volInfo = varargin{4};
        
        reducedY = bayes_meta_feature_abstract(Y, thresh, bayes_model, volInfo);
        
        [bayes_model, Y] = classify_naive_bayes('setup', reducedY, bayes_model.Xi, bayes_model.params.activation_cutoff, bayes_model.params.selectivity_cutoff, 1, bayes_model.params.k);
        
        varargout{1} = bayes_model;
        varargout{2} = Y;
        
    case 'xval'
        Y = full(varargin{1});
        Xi = varargin{2};

        activation_cutoff = [];
        selectivity_cutoff = [];
        if length(varargin) > 2

            activation_cutoff = varargin{3};
            selectivity_cutoff = varargin{4};

        end

        k = .5;
        gam = 1;
        if length(varargin) > 4

            k = varargin{5};
            gam = varargin{6};

        end
        
        if length(varargin) > 6
            volInfo = varargin{7};
            do_feature_abstract = 1;
        end

    
        % initialize
        if do_feature_abstract, Yorig = Y; end
        
        [Y, whkeep] = eliminate_empty_Y(Y);

        Y = sparse(Y);

        nobs = size(Y, 1);
        fprintf('\nCross-validating %3.0f observations: 00000', nobs);

        log_joint = zeros(nobs, size(Xi, 2));
        map = zeros(nobs, size(Xi, 2));
        class_est = zeros(nobs, 1);
        best_log = zeros(nobs, 1);

        true_class = indic2condf(Xi);

        % Update volinfo for feature abstraction
% %         volInfo.wh_inmask = double(whkeep');
% %         volInfo.n_inmask = sum(whkeep);
% %         volInfo.image_indx(volInfo.image_indx) = whkeep';
% %         volInfo.xyzlist(~whkeep',:) = [];

        % run
        for i = 1:nobs

            if mod(i,10) == 0, fprintf(1,'\b\b\b\b\b%05d',i); end

            include_in_training = true(nobs, 1);
            include_in_training(i) = 0;

            % setup: does feature-selection as well, if requested; no
            % need to return Y, because ...
            bayes_model = classify_naive_bayes_setup(Y(include_in_training, :), Xi(include_in_training, :), activation_cutoff, selectivity_cutoff, whkeep, 0, k);

            bayes_model.params.k = k;
            bayes_model.params.gam = gam;

            testdata = Y(i, bayes_model.whkeep_from_notempty);
            
            % feature abstraction step
            if do_feature_abstract
                [reducedY, cl, cluster_indx] = bayes_meta_feature_abstract(Yorig(include_in_training, :), .1, .1, bayes_model, volInfo, 0);
                cluster_indx = cluster_indx(bayes_model.whkeep);
                
                bayes_model = classify_naive_bayes('setup', reducedY, Xi(include_in_training, :), 0, 0, gam, k);
                
                wh_cl = cluster_indx(testdata);
                wh_cl = unique(wh_cl); wh_cl(wh_cl == 0) = [];
                testdata = false(1, bayes_model.nfeatures);
                testdata(wh_cl) = 1;

            end
            
            [class_est(i) log_joint(i,:) best_log(i) map(i,:)] = classify_naive_bayes_test(testdata, bayes_model);
        end

        fprintf(1, '\n');

        [prop_correct, confusion_mtx, misclass, prop_correct_by_class, chance, chance_95_ci] = classify_naive_bayes_eval(true_class, class_est);

        varargout{1} = struct('class_est', class_est, 'log_joint', log_joint, 'best_log', best_log, 'map', map, ...
            'prop_correct', prop_correct, 'prop_correct_by_class', prop_correct_by_class, 'confusion_mtx', confusion_mtx, ...
            'misclass', misclass, 'chance', chance, 'chance_95_ci', chance_95_ci, 'params', bayes_model.params);


    case {'write', 'images', 'write_images'}
        bayes_model = varargin{1};
        Y = varargin{2};
        volInfo = varargin{3};
        conditionnames = varargin{4};

        bayes_model = write_images(bayes_model, Y, volInfo, conditionnames);
        varargout{1} = bayes_model;


    case {'plot', 'output'}
        bayes_model = varargin{1};
        disptype = varargin{2};

        cl = display_output(bayes_model, disptype, varargin{3:end});

        varargout{1} = cl;

    otherwise
        error('Unknown method (invalid first argument).');
end






end   % end main function










function [bayes_model, Y] = classify_naive_bayes_setup(Y, Xi, activation_cutoff, selectivity_cutoff, whkeep, return_full_ps, k)

% Y is data matrix, nobs observations x nfeatures variables
% in a meta-analysis, e.g., this is studies x voxels
%
% Each row of Y is an observation; if Y has one row, we're classifying an
% observation (e.g., one brain map)
%
% Xi are indicators for tasks (classes), ntasks columns of nobs elements
% (i.e., an nobs x ntasks matrix)

% eliminate features with no 'on' values (activations)


% k is regularization param, biases towards 0.5
if nargin < 7, k = 1; end

whkeep_from_notempty = whkeep;  % indicator of features after feature-selection, in index of not-empties

% feature selection, if requested
if ~isempty(selectivity_cutoff) && ~isempty(activation_cutoff)

    [Y, whkeep, whkeep_from_notempty] = select_features(Y, Xi, activation_cutoff, selectivity_cutoff, whkeep);

end

[nobs, nfeatures] = size(Y);
nclasses = size(Xi, 2);

if return_full_ps
    [priors, pa1_given_t, pa0_given_t, pt_given_act1, pt_given_act0, pa1_given_not_t] = ...
        bayes_get_probabilities(Y, Xi, k);

else
    % key results (not all) only; saves computation time in
    % crossvalidation.
    [priors, pa1_given_t, pa0_given_t] = bayes_get_probabilities(Y, Xi, k);

end

true_class = indic2condf(Xi);


bayes_model = struct('nobs', nobs, 'nfeatures', nfeatures, 'nclasses', nclasses, ...
    'whkeep', whkeep, 'whkeep_from_notempty', whkeep_from_notempty, ...
    'params', [], ...
    'priors', priors, 'pa1_given_t', pa1_given_t, 'pa0_given_t', pa0_given_t, ...
    'Xi', Xi, 'true_class', true_class);

bayes_model.params.gam = 1;  % gamma, weighting on prior vs. joint evidence
bayes_model.params.k = k;      % k, regularization of pa|t, bias towards 0.5.
bayes_model.params.activation_cutoff = activation_cutoff;
bayes_model.params.selectivity_cutoff = selectivity_cutoff;

if return_full_ps
    bayes_model.pt_given_act1 = pt_given_act1;
    bayes_model.pt_given_act0 = pt_given_act0;
    bayes_model.pa1_given_not_t = pa1_given_not_t;

end


end




function [class_est log_joint best_log map p_obs_act_given_class] = classify_naive_bayes_test(Y, bayes_model)

% Estimate classes, log joint p, and MAP estimates for a test data set

% P(task | act)
% The var. below is all that is needed to choose the best class.
% argmax log_joint
% do this for each observation

[nobs, nfeatures] = size(Y);

if nfeatures ~= bayes_model.nfeatures
    %  try eliminating zero voxels
    Y = Y(:, bayes_model.whkeep);
    [nobs, nfeatures] = size(Y);
end

if nfeatures ~= bayes_model.nfeatures
    error('Number of features (i.e., voxels/regions) in Y must equal that in bayes_model structure or that of original training data.');
end

log_joint = zeros(nobs, bayes_model.nclasses);
map = zeros(nobs, bayes_model.nclasses);

class_est = zeros(nobs, 1);
best_log = zeros(nobs, 1);

gam = bayes_model.params.gam;    % prior weighting factor; from 0 to nfeatures.

logpriors = log(bayes_model.priors);
%logpriors = logpriors(ones(nfeatures, bayes_model.nclasses));

for i = 1:nobs

    wh_active = Y(i,:)' > 0;

    % sum of log prob. of activation at each active voxel, plus sum of
    % log prob. of NOT activation at each inactive voxel (1 - p(activation))
    % Same as ln of product of all p(act state = observed state |
    % class) for all features

    % prob. that activity state is the observed one in the test vector,
    % given each class.
    % p(Fi = fi | C = c)
    % p(ACTi = acti | T = t)

    % % % %     p_obs_act_given_class(wh_active, :) = bayes_model.pt_given_act1(wh_active, :);         % p(active | t, obs active)
    % % % %     p_obs_act_given_class(~wh_active, :) = bayes_model.pt_given_act0(~wh_active, :);

    p_obs_act_given_class = zeros(nfeatures, bayes_model.nclasses);
    p_obs_act_given_class(wh_active, :) = bayes_model.pa1_given_t(wh_active, :);         % p(active | t, obs active)
    p_obs_act_given_class(~wh_active, :) = bayes_model.pa0_given_t(~wh_active, :);   % p(not active | t, obs not active)

    loglikelihood = sum(log(p_obs_act_given_class));
    log_joint(i,:) = gam * log(bayes_model.priors) + loglikelihood;      % log of joint prob. p(T, Y), p(task, act1, act2, act3, etc.)

    %log_joint(i,:) = sum(log(p_obs_act_given_class)); % proportional to the posterior.

    [best_log(i) class_est(i)] = max(log_joint(i,:));

    % Divide joint by evidence
    % get scaled estimate of p(task | activation)
    % These are Max. a posteriori estimates, and values may be > 1 because
    % we're assuming independence when that assumption is not likely to be true
    % * something still wrong with this
    % % % %
    % % % %         % evidence
    % % % %         log_evidence = sum(bayes_model.log_evidence_act(wh_active)) + sum(bayes_model.log_evidence_notact(~wh_active));
    % % % %         % because something is wrong, exp() gives Inf values a lot, so isn't
    % very sensible.
    % % % %         map(i,:) = log_joint(i,:) - log_evidence;  %exp(log_joint(i,:) - log_evidence);


    map(i,:) = exp(log_joint(i,:) ./ (gam + nfeatures) );  %exp(log_joint(i, :)); % will be very small
    %sum(log_joint) + nfeatures * logpriors;

    %map = map ./ sum(map);


end


end




function [prop_correct, m, misclass, prop_correct_by_class, chance, chance_95_ci] = classify_naive_bayes_eval(true_class, class_est, wh_obs)

if nargin < 3 || isempty(wh_obs)
    % assume we have all observations.
    wh_obs = 1:length(true_class);
end

[m,dp,corr,far,misclass] = confusion_matrix(true_class(wh_obs),class_est);

totaln = sum(m(:));
prop_correct = trace(m) ./ totaln;

nclasses = length(unique(true_class(true_class ~= 0)));

prop_correct_by_class = diag(m) ./ sum(m,2);


chance = 1 ./ nclasses;
chance_95_ci = [chance - sqrt(chance * (1 - chance) ./ totaln)   chance + sqrt(chance * (1 - chance) ./ totaln)];

end



function [Y, whkeep] = eliminate_empty_Y(Y)

fprintf('Eliminating any empty features. \n');

% using whkeep helps avoid out of memory errors.
whzero = (sum(Y) < eps);
whkeep = ~whzero;
Y = Y(:, whkeep);

end



function [Y, whkeep, whomsave] = select_features(Y, Xi, activation_cutoff, selectivity_cutoff, whkeep)

% select features
% -----------------------------------------------------------------

[nobs, nfeatures] = size(Y);
ntasks = size(Xi,2);

sumtasks = sum(Y); sumtasks(sumtasks==0) = 1;  % avoid warning for empty tasks

% p(task | activation) estimates
ptask = (Xi' * Y) ./ repmat(sumtasks,ntasks,1);

maxp = max(ptask);  % max prob of task across vars

pt = (Xi' * Y) ./ repmat(sum(Xi)',1,nfeatures);  % proportion of each task that activated in an area
ptmax = max(pt);
whomsave = maxp > (selectivity_cutoff) & ptmax > activation_cutoff;

if sum(whomsave) == 0
    disp('Warning: No variables meet feature selection criteria.');
    whomsave = (maxp == max(maxp));
end

whkeep(whkeep) = whomsave;  % indices of saved voxels in whole-brain (original) index vector
%whkeep = whkeep(whomsave);  % overall indices of saved voxels
Y = Y(:,whomsave);

%fprintf(1,'Selected: %3.0f ',nvars);
% get weights on voxels based on n activations.... ******

end




function bayes_model = write_images(bayes_model, Y, volInfo, conditionnames)
% bayes_model = write_images(bayes_model, MC_Setup.volInfo, MC_Setup.Xinms)

% Xi = bayes_model.Xi;

%[priors, pa1_given_t, pa0_given_t, varargout] = bayes_get_probabilities(Y, Xi, bayes_model.params.k);

nfeatures_orig = length(bayes_model.whkeep);

pa1_given_t = NaN * zeros(nfeatures_orig, bayes_model.nclasses);
pa1_given_not_t = NaN * zeros(nfeatures_orig, bayes_model.nclasses);

pt_given_act1 = NaN * zeros(nfeatures_orig, bayes_model.nclasses);
pt_given_act0 = NaN * zeros(nfeatures_orig, bayes_model.nclasses);

pa1_given_t(bayes_model.whkeep,:) = bayes_model.pa1_given_t;
pa1_given_not_t(bayes_model.whkeep,:) = bayes_model.pa1_given_not_t;

pt_given_act1(bayes_model.whkeep,:) = bayes_model.pt_given_act1;
pt_given_act0(bayes_model.whkeep,:) = bayes_model.pt_given_act0;


% now done in bayes_get_probabilities
% % pa1_given_not_t = (double(~(Xi)') * Y)' ./ repmat(sum(~Xi), size(Y, 2), 1);
% % pa1_given_not_t(pa1_given_not_t < eps) = 1 ./ bayes_model.nobs;
% % pa1_given_not_t(pa1_given_not_t == 1) = 1 - (1 ./ bayes_model.nobs);

%%%%pa1_given_not_t(~bayes_model.whkeep, :) = NaN;


disp('Writing images to disk in the current directory.');

for i = 1:bayes_model.nclasses
    pa1_given_tname{i} = ['pa|t_' conditionnames{i} '.img'];
    iimg_reconstruct_3dvol(pa1_given_t(:,i), volInfo, 'outname', pa1_given_tname{i});

    disp(pa1_given_tname{i})

    pa_given_nottname{i} = ['pa|nott_' conditionnames{i} '.img'];
    iimg_reconstruct_3dvol(pa1_given_not_t(:,i), volInfo, 'outname', pa_given_nottname{i});

    disp(pa_given_nottname{i})

    pt_given_aname{i} = ['pt|a_' conditionnames{i} '.img'];
    iimg_reconstruct_3dvol(pt_given_act1(:,i), volInfo, 'outname', pt_given_aname{i});

    disp(pt_given_aname{i})

    pt_given_notaname{i} = ['pt|nota_' conditionnames{i} '.img'];
    iimg_reconstruct_3dvol(pt_given_act0(:,i), volInfo, 'outname', pt_given_notaname{i});

    disp(pt_given_notaname{i})
end


% likelihood ratio at each voxel
% should penalize this in some way for low activations...
%pa1_given_t(pa1_given_t == 0) = 1;
%pa1_given_not_t(pa1_given_not_t == 0) = 1;

% higher numbers mean more association between activation and class.
% lr for 2 classes can both be positive IF not all obs. fall within class 1
% or class 2!!!
evidence_yes = (sum(Y) ./ bayes_model.nobs)';
evidence_yes = evidence_yes(:, ones(1, bayes_model.nclasses));

% % priors = NaN * zeros(nfeatures_orig, bayes_model.nclasses);
% % priors(bayes_model.whkeep,:) = bayes_model.priors(ones(bayes_model.nfeatures, bayes_model.nclasses));

% add .1 to regularize...this could be k param
lr = ( (pa1_given_t + .1) ./ (evidence_yes + .1) ) - 1;
%pt_given_act1 ./ pt_given_act0;
lr(lr < .05 & lr > -.05) = NaN;
%lr = log(lr);

for i = 1:bayes_model.nclasses
    lrname{i} = ['likeratio_' conditionnames{i} '.img'];
    iimg_reconstruct_3dvol(lr(:,i), volInfo, 'outname', lrname{i});

    disp(lrname{i})
end

bayes_model.image_output.volInfo = volInfo;
bayes_model.image_output.pa1_given_tname = pa1_given_tname;
bayes_model.image_output.pa_given_nottname = pa_given_nottname;

bayes_model.image_output.pt_given_aname = pt_given_aname;
bayes_model.image_output.pt_given_notaname = pt_given_notaname;

bayes_model.image_output.likerationame = lrname;

bayes_model.image_output.pa1_given_t = pa1_given_t;
bayes_model.image_output.pa1_given_not_t = pa1_given_not_t;
bayes_model.image_output.likeratio = lr;

end




function cl = display_output(bayes_model, disptype, varargin)
% display_output(bayes_model, disptype, varargin)

switch spm('Ver')
    case 'SPM2'
        % spm_defaults is a script
        disp('WARNING: spm defaults not set for spm2. Make sure your defaults are set correctly');

    case 'SPM5'
        % spm_defaults is a function
        spm_defaults()
     
    case 'SPM8'
        % spm_defaults is a function
        spm_defaults()
        
    otherwise
        % unknown SPM
        disp('Unknown version of SPM!');
        spm_defaults()
end

cl = [];
colors1 = [1 .7 0];
colors2 = [0 0 1];

for i = 1:length(varargin)
    if iscell(varargin{i})
        colors1 = varargin{i}{1};
        colors2 = varargin{i}{2};
    elseif length(varargin{i}) == 1
        thresh = varargin{i};
    end
end

switch disptype

    case {'pa|t', 'lr', 'pt|a'}

        overlay = which('spm2_single_subj_T1_scalped.img');
        spm_check_registration(repmat(overlay, bayes_model.nclasses, 1));

        for i = 1:bayes_model.nclasses

            switch disptype
                case 'pa|t'
                    cl{i} = mask2clusters(bayes_model.image_output.pa1_given_tname{i});
                    fprintf('Displayed: p(a | t) for each task.\n');

                    spm_orthviews_name_axis(bayes_model.image_output.pa1_given_tname{i}, i);


                case 'lr'
                    cl{i} = mask2clusters(bayes_model.image_output.likerationame{i});
                    fprintf('Displayed: regions with positive likelihood ratio for each task.\n');

                    spm_orthviews_name_axis(bayes_model.image_output.likerationame{i}, i);

                case 'pt|a'
                    cl{i} = mask2clusters(bayes_model.image_output.pt_given_aname{i});
                    fprintf('Displayed: p(t | a) for each task.\n');

                    spm_orthviews_name_axis(bayes_model.image_output.pt_given_aname{i}, i);

            end

            % threshold
            if exist('thresh', 'var')
                for j = 1:length(cl{i})
                    wh_omit = abs( cat(2, cl{i}(j).Z) ) < thresh;
                    cl{i}(j).XYZ(:, wh_omit) = [];
                    cl{i}(j).XYZmm(:, wh_omit) = [];
                    cl{i}(j).Z(:, wh_omit) = [];
                end

                clu = clusters2CLU(cl{i});
                cl{i} = tor_extract_rois([], clu, clu);

            end


            cluster_orthviews(cl{i}, 'handle', i, 'add');

        end

        switch disptype
            case 'pa|t'
                %spm_orthviews_change_colormap([.2 .2 .4], [1 1 0], [.9 .6 .1]);  %slate to orange to yellow
                spm_orthviews_change_colormap([.5 0 0], [1 1 0], [1 .5 0]);  %blue to yellow

            case  'lr'
                spm_orthviews_change_colormap(colors2, colors1);  %blue to white to yellow

            case 'pt|a'
                spm_orthviews_change_colormap(colors2, colors1);

        end

    case 'lr surface'

        poscm = colormap_tor([.2 .2 .4], [1 1 0], [.9 .6 .1]);  %slate to orange to yellow
        negcm = colormap_tor([0 0 1], [0 .3 1]);  % light blue to dark blue (shouldn't be needed!)


        for i = 1:bayes_model.nclasses

            cl{i} = mask2clusters(bayes_model.image_output.likerationame{i});

            fprintf('Displayed: regions with positive likelihood ratio for each task.\n');

            % Left medial and lateral surface
            % ----------------------------------------------------------
            create_figure('Brain Surface'); scn_export_papersetup(800);

            cluster_surf(cl{i}, 4, 'heatmap', 'colormaps', poscm, negcm, 'hires left');
            sh = cluster_surf(cl{i}, 4, 'heatmap', 'colormaps', poscm, negcm, 'brainstem');
            saveas(gcf, [bayes_model.image_output.likerationame{i}(1:end-4) '_left_medial'], 'png');

            view(270,0); lightRestoreSingle;
            saveas(gcf, [bayes_model.image_output.likerationame{i}(1:end-4) '_left_lateral'], 'png');

            % Right medial and lateral surface
            % ----------------------------------------------------------
            create_figure('Brain Surface'); cluster_surf(cl{i}, 4, 'heatmap', 'colormaps', poscm, negcm, 'hires right');
            sh = cluster_surf(cl{i}, 4, 'heatmap', 'colormaps', poscm, negcm, 'brainstem');
            scn_export_papersetup(800); saveas(gcf, [bayes_model.image_output.likerationame{i}(1:end-4) '_right_lateral'], 'png');

            view(270,0); lightRestoreSingle;
            saveas(gcf, [bayes_model.image_output.likerationame{i}(1:end-4) '_right_medial'], 'png');

            % Limbic
            % ---------------------------------------------------------
            create_figure('Brain Surface'); cluster_surf(cl{i}, 4, 'heatmap', 'colormaps', poscm, negcm, 'limbic');
            saveas(gcf, [bayes_model.image_output.likerationame{i}(1:end-4) '_limbic'], 'png');
            saveas(gcf, [bayes_model.image_output.likerationame{i}(1:end-4) '_limbic'], 'fig');
        end


    case 'map plot'

        % MAP plot

        colors = {'b' 'r' 'k' 'g' 'm'};
        create_figure('MAP Plot across observations');
        for i = 1:bayes_model.nclasses

            trueclassi = find(bayes_model.Xi(:,i) == 1);
            myclassmap = bayes_model.apparent.map(trueclassi, i);

            plot(bayes_model.apparent.map(:, i), 'color', colors{i});
            plot(trueclassi, bayes_model.apparent.map(find(bayes_model.Xi(:,i) == 1), i), 'o', 'Color', colors{i},'MarkerFaceColor',colors{i});

            incorrect_obs = ~(myclassmap == max(bayes_model.apparent.map(trueclassi, :), [], 2));
            incorrect_obs_vec = zeros(bayes_model.nobs, 1);
            incorrect_obs_vec(trueclassi(incorrect_obs)) = 1;

            plot(find(incorrect_obs_vec), myclassmap(incorrect_obs), 'kx', 'MarkerSize', 12, 'LineWidth', 2);
        end
        xlabel('Observations'); ylabel('log(MAP) (up to scaling constant)');



    case 'class plot'

        pr = bayes_model.priors;
        fr = [linspace(1, 0, 10)' linspace(0, 1, 10)'];   % frequencies of guesses, from always A to always B

        xfr = linspace(1, 0, 10);

        %hr_byclass = fr .* repmat(pr, size(fr, 1), 1)
        colors = {'b' 'r' 'k' 'g' 'm'};

        hr = fr * pr';

        create_figure('Summary Line Plot');
        %plot(xfr, hr_byclass);
        plot(xfr, hr, 'k', 'LineWidth', 2);

        phighestfreq = sum(bayes_model.xval.class_est == wh) ./ length(bayes_model.xval.class_est);  % p chosen for highest-frequency class
        plot(phighestfreq, bayes_model.xval.prop_correct, 'ko', 'LineWidth', 2, 'MarkerFaceColor', [.5 .5 .5], 'MarkerSize', 12);

        xx = repmat(phighestfreq, 2, bayes_model.nclasses);
        yy = [bayes_model.xval.prop_correct_by_class'; bayes_model.xval.prop_correct_by_class'];

        for i = 1:size(yy, 2)
            plot(xx(:,i), yy(:,i), 'o', 'Color', colors{i}, 'LineWidth', 2, 'MarkerFaceColor', colors{i}, 'MarkerSize', 12);
        end

        % Standard error of hit rate
        % var(hits for each category) = p*(1-p)*N for that category
        % where p is the prior p(correct) for that category, and N is the number of
        % choices in that category given the frequency of sampling that category.

        se_hr = sqrt( hr .* (1 - hr) ./ bayes_model.nobs );
        ci95_hr = 1.96 * se_hr;

        plot(xfr, hr + ci95_hr, 'k--', 'LineWidth', 1);
        plot(xfr, hr - ci95_hr, 'k--', 'LineWidth', 1);
        fill_around_line(hr, ci95_hr, 'k', xfr);

        ylabel('Correct classification rate');
        xlabel('Prop. choices from highest class');


        % N = fr .* bayes_model.nobs;     % for each of k categories, in a row

        % pq = repmat(pr .* (1 - pr), size(fr,1), 1);
        % varsum = sum(pq .* N, 2)
        %
        %
        % se_hr = sqrt( sum(repmat(pr .* (1 - pr), size(fr,1), 1) ./ (bayes_model.nobs .* fr), 2) )
        % se_hr(isinf(se_hr)) = 0;


    otherwise
        error('Unknown display option.')
end

end
