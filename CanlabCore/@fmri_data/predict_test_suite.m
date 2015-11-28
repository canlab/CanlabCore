function [allcverr, allyhat] = predict_test_suite(dat, varargin)
% Run a set of cross-validated prediction algorithms on an fmri_data object
% and plot the outcome.
%
% :Usage:
% ::
%
%     [allcverr, allyhat] = predict_test_suite(dat, [optional inputs])
% 
% :Functionality:
%   - Requires matlab 2012a or later for full functionality
%   - Handles categorical or continuous outcomes automatically
%
% :Inputs:
%
%   **dat:**
%        an fMRI data object. 
%        dat.Y must be assigned, and must have continuous or binary outcomes assigned.
%
% :Optional:
%
% **quick:**
%        Skip extended output
%
% **nfolds:**
%        Followed by number of folds or custom holdout vector (default = 5-fold balanced)
%
% :Examples:
% ::
%
%    predict_test_suite(dat, 'nfolds', subjid);
%
% ..
%    Tor Wager, copyright 2012. Initial version: Nov 2012
%
%    Programmers' notes:
%
%    Here are possible extensions to the functionality:
%    1) Test different data scaling methods in preprocess(dat) and/or
%       rescale(dat)
% 
%    2) Feature selection: Test effects of thresholding weight maps
%       (and others)
%
%    3) Test "best case" feature selection effects - circular analysis
%
%    Lasso trace plots?  Optimal alpha/lambda in elastic net?
% ..

create_figure('predicted-actual corr', 2, 4);

nlevels = length(unique(dat.Y));
nfolds = 5;

% SET UP ALGORITHMS
% ------------------------------------------------------------------
switch nlevels
    case {1, 0}
        error('NO VARIANCE IN obj.Y')
        
    case 2
        methodnames = {'cv_svm' 'cv_lassopcrmatlab'};
        errtype = 'mcr';
        outputsummaryfield = 'cverr';
        
    otherwise
        methodnames = {'cv_svr' 'cv_pcr' 'cv_lassopcrmatlab' 'cv_lassopcr'}; %'cv_univregress' 
        errtype = 'mse';
        outputsummaryfield = 'pred_outcome_r';

end

% SET UP OPTIONAL INPUTS
% Pass all valid arguments specified here on to predict.m
% ------------------------------------------------------------------
optargs = {};
for i = 1:length(varargin)
        if ischar(varargin{i})
            switch varargin{i}
                    % functional commands for predict
                case {'nfolds', 'useparallel'} 
                    optargs{end+1} = varargin{i};
                    optargs{end+1} = varargin{i+1};
            end
        end
end

if ~any(strcmp(optargs, 'nfolds'))
    optargs{end+1} = 'nfolds';
    optargs{end+1} = nfolds;
end

    

allyhat = zeros(size(dat.Y, 1), length(methodnames));

allcverr = zeros(1, length(methodnames));

% RUN ALL ALGORITHMS
% ------------------------------------------------------------------
for i = 1:length(methodnames)
    
    try
        
        [cverr, stats, optout] = predict(dat, 'algorithm_name', methodnames{i}, 'error_type', errtype, optargs{:});
        
        allyhat(:, i) = stats.yfit;
        
        allcverr(i) = cverr;
        
        subplot(2, 4, i);
        
        plot_correlation_samefig(stats.yfit, stats.Y);
        xlabel('Predicted (5-fold crossval)');
        ylabel('Observed')
        title(sprintf('%s', methodnames{i}));
        
        plot_folds(stats)
        drawnow
        
    catch
        fprintf('ERROR running %s\n', methodnames{i})
        keyboard
        %rethrow(lasterr)
    end
    
end % loop through methods

% SUMMARY OUTPUT
% ------------------------------------------------------------------
rmtx = corrcoef(allyhat);
disp('Correlations among predicted values across methods')
print_matrix(rmtx, methodnames, methodnames);

% Quit now if we don't want extended output

if any(strcmp(varargin, 'quick')), return, end

% IMPACT OF COMPONENT SELECTION
% ------------------------------------------------------------------
numc = round(linspace(1, max(2, size(dat.dat, 2) - 2), 10));
numc = unique(numc);

if nlevels == 2
    mymethodname = 'cv_lassopcrmatlab';
else
    mymethodname = 'cv_pcr';
end

for i = 1:length(numc)
    
    [cverr_comp(i), stats] = predict(dat, 'numcomponents', numc(i), 'algorithm_name', mymethodname, 'error_type', errtype, optargs{:});
    
    r_comp(i) = stats.(outputsummaryfield);
    
    if nlevels == 2
        % reverse so that higher = better
        r_comp(i) = 1 - r_comp(i);
    end
    
end

whplot = length(methodnames) + 1;
subplot(2, 4, whplot)
axh(1) = gca;

plot(numc, r_comp, 'o-', 'Color', [.3 .3 .3], 'MarkerFaceColor', [.5 .5 1], 'LineWidth', 2, 'MarkerSize', 10);
xlabel('Components');
if nlevels == 2
    ylabel('Classification accuracy');
else
    ylabel('Pred-outcome correlation');
end
title(sprintf('Component selection: %s', mymethodname));


% IMPACT OF PENALIZATION 
% ------------------------------------------------------------------
mymethodname = 'cv_lassopcrmatlab';

for i = 1:length(numc)
    
    [cverr_comp(i), stats] = predict(dat, 'lasso_num', numc(i), 'algorithm_name', mymethodname, 'error_type', errtype, optargs{:});
    
    r_comp(i) = stats.(outputsummaryfield);
    
    if nlevels == 2
        % reverse so that higher = better
        r_comp(i) = 1 - r_comp(i);
    end
end

whplot = length(methodnames) + 2;
subplot(2, 4, whplot)
axh(2) = gca;

plot(numc, r_comp, 'o-', 'Color', [.3 .3 .3], 'MarkerFaceColor', [.5 .5 1], 'LineWidth', 2, 'MarkerSize', 10);
xlabel('Lasso: components retained');
if nlevels == 2
    ylabel('Classification accuracy');
else
    ylabel('Pred-outcome correlation');
end
title(sprintf('LASSO penalization: %s', mymethodname));

subplot(2, 4, 8); axis off

equalize_axes(axh);

end % main function


function plot_folds(stats)


k = length(stats.teIdx);

colors = scn_standard_colors(k);

for i = 1:length(stats.teIdx)

    wh = find(stats.teIdx{i});
    
    if length(wh) < 3 % don't do it for 1-2 points only
        continue
    end
    
    yfit = stats.yfit(wh);
    y = stats.Y(wh);
    
    plot(yfit, y, 'o', 'MarkerFaceColor', colors{i}); 
    
    
    % ref line
    b = glmfit(yfit, y);
    
    plot([min(yfit) max(yfit)], [b(1)+b(2)*min(yfit) b(1)+b(2)*max(yfit)], 'Color', colors{i}, 'LineWidth', 2);
    
end

axis tight

end


