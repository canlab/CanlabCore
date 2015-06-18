function wh_holdout_final = xval_select_holdout_set(Y, covs, nfolds, holdout_proportion, varargin)
    %
    % wh_holdout_final = xval_select_holdout_set(Y, covs, nfolds, holdout_proportion, [doplot])
    %
    % Selects one or more holdout sets balancing on an outcome and set of
    % covariates.  The idea is inspired by Rubin's propensity score method, and
    % the logic is that selecting holdout sets that are balanced on the outcome
    % and other variables of no interest will reduce the variance in the
    % estimated predictive accuracy in a cross-validated parameter optimization
    % scheme.  The goal of this algorithm is thus to choose holdout sets where
    % the probability of holdout assignment is balanced on (independent of) the covariates.
    %
    % For assessing accuracy with little bias, a leave-one-out holdout method
    % can be a good choice, as it has minimum bias.  However, the accuracy
    % estimates have high variance, which often motivates k-fold schemes with
    % larger holdout sets.  For model/parameter selection (e.g., in inner
    % cross-validation schemes, or when model selection is the primary goal),
    % larger holdout sets and small training sets may be desirable.  See Shao,
    % 1998, JASA.
    %
    % It is also advantageous to have training/holdout sets that are independent
    % of one another, but this is not optimized here.
    %
    % Input variables:
    % Y = outcome variable
    % covs = covariates to balance on
    % nfolds = number of holdout sets desired
    % holdout_proportion = how much of the sample is reserved for holdout, [0 1]
    %
    % Output variables:
    % wh_holdout = a matrix of [obs x nfolds], containing optimized holdout sets
    %
    % Example inputs:
    % Y = randn(30, 1);
    % covs = randn(N, 1);
    % nfolds = 10;
    % holdout_proportion = .5;
    %
    % Example function call:
    % wh_holdout = xval_select_holdout_set(Y, covs, 20, .5);
    %
    % Notes:
    % There is a plot option hard-coded in the script for demonstration
    % purposes.

    doplot = 0;
    if ~isempty(varargin), doplot = varargin{1}; end

    anyzero = 1;
    wh_holdout_final = [];

    wh_holdout_final = generate_balanced_sets(Y, covs, nfolds, holdout_proportion, doplot);

    % find any missing observations
    % Make sure all observations are selected equally often, to the
    % degree possible
    % Find the holdout sets that select the most frequently used observations
    obssum = sum(wh_holdout_final, 2);
    ismissing = double(obssum == 0);

    while any(ismissing)

        wh_holdout_new = generate_balanced_sets(Y, covs, nfolds, holdout_proportion, doplot); % obs x sets
        
        new_to_save = wh_holdout_new' * ismissing; % these holdout sets have the missing obs from before
    
        %frequent_obs_score = wh_holdout_new' * obssum; % higher values mean the holdout set uses more frequently selected observations (bad)
        
        % Replace the ones with the most frequent scores with the new ones
        % to keep
        frequent_obs_score = wh_holdout_final' * obssum; % higher values mean the holdout set uses more frequently selected observations (bad)
        [f2, indx] = sort(frequent_obs_score, 1, 'descend');  
        
        if any(new_to_save)
            
            new_to_save = find(new_to_save);
            wh_to_replace = indx(1:length(new_to_save));
            wh_holdout_final(:, wh_to_replace) = wh_holdout_new(:, new_to_save);

            % it is possible that we've replaced one that causes another obs to
            % be missing, so we need the while loop
            obssum = sum(wh_holdout_final, 2);
            ismissing = double(obssum == 0);
        else
            
            %% do nothing   ismissing = 1;
            
        end

    end
    


end % function



function wh_holdout_final = generate_balanced_sets(Y, covs, nfolds, holdout_proportion, doplot)

    N = length(Y);

    niter = nfolds * 10;
    [betavals, tvals] = deal(zeros(1, niter));
    wh_holdout = false(N, niter);

    warning off % iteration limit...
    
    for i = 1:niter

        wh_holdout(:, i) = zeros(N, 1);
        nh = round(N .* holdout_proportion);

        wh = randperm(N);
        wh = wh(1:nh);
        wh_holdout(wh, i) = 1;

        % For testing purposes
        %covs = mvnrnd([1 3], [1 .4; .4 1], N); % random outcome and covs

        % p-value is linear with respect to betas near the origin (50%)
        % so we can select subsets based on beta alone.
        %
        b = glmfit([Y covs], wh_holdout(:, i), 'binomial', 'link', 'logit');

        % This is used for testing purposes
        if doplot
            [b, dev, stat] = glmfit([Y covs], wh_holdout(:, i), 'binomial', 'link', 'logit');
            
            %
            %         % plot partial fit wrt Y
            %         yfit = [min(Y):.1:max(Y)];
            %         xfit = glmval(b(1:2), yfit, 'logit');
            %
            %         create_figure; plot(Y, wh_holdout(:, i), 'ko','MarkerFaceColor', 'k');
            %         set(gca, 'XLim', [-.5 1.5]);
            %         plot(yfit, xfit, 'r');
            %         ylabel('P(holdout)'); xlabel('outcome');

            tvals(i) = stat.t(2);
        end

        betavals(i) = b(2);
        %

    end

    warning on
    
    % Save the most balanced holdout sets

    cutoff = prctile(abs(betavals), 10);  % save top 10% -- that ensures us the right size holdout set
    wh_to_keep = abs(betavals) <= cutoff;
    wh_holdout_final = wh_holdout(:, wh_to_keep);

    if doplot
        create_figure('plot'); plot_correlation_samefig(betavals, tvals);
        plot(betavals(wh_to_keep), tvals(wh_to_keep), 'go');
        xlabel('Logistic reg beta'); ylabel('T-value');
        title('Selection set estimates for logistic reg, outcome')

        % plot partial fit wrt Y for one example holdout set
        [b, dev, stat] = glmfit([Y covs], wh_holdout(:, 1), 'binomial', 'link', 'logit');

        yfit = [min(Y):.1:max(Y)];
        xfit = glmval(b(1:2), yfit, 'logit');

        create_figure('Example holdout set'); plot(Y, wh_holdout(:, 1), 'ko','MarkerFaceColor', 'k');
        set(gca, 'XLim', [-.5 1.5]);
        plot(yfit, xfit, 'r');
        ylabel('P(holdout)'); xlabel('Outcome');
        axis auto
        title('Example holdout set: P(holdout|Y)')
        drawnow
    end

end

