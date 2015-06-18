function [trialw, allw, wh_bad, ortho_trialX] = single_trial_weights(trialX, trial_beta_indices, bf)
    % [trialw, allw, wh_bad, ortho_trialX] = single_trial_weights(trialX, trial_beta_indices, [basis_set_name or basis set vectors])
    %
    % Get weights for single trial model based on estimability of the beta
    % coefficients for each trial.
    %
    % These weights can be used in a 2nd level analysis (single subject,
    % condition effects) to weight trials by the quality of their
    % estimation at the 1st (single trial) level.
    %
    % Because multiple basis functions are combined in a nonlinear way,
    % this function uses the minimum efficiency for the regressors used to
    % model the trial.
    %
    % Note: uses inverse, so you should not have completely redundant columns in the
    % design matrix.
    %
    % Tor Wager, April 1, 2008
    %
    % Examples:
    % For NSF1 or Expectancy1:
    % load exp1227/exp1227_design.mat
    % trial_beta_indices = 1:size(trialX_of_interest, 2);  % Trial est are first xx betas  1:192;
    % basis_set_name = trialmodels.basistype               % or 'pain10';
    % [trialw, allw] = single_trial_weights(trialmodels.trialX, trial_beta_indices, trialmodels.basistype);
    %
    % OR
    % load exp1227/exp1227_design.mat
    % [trialw, allw, wh_bad] = single_trial_weights(trialmodels.trialX, 1:trialmodels.ntrials * size(trialmodels.bf, 2), trialmodels.bf);
    %
    % Example for running through all subjects in an experiment:
    % % cd /Volumes/Expectancy/Imaging/Analyses/new_analysis
    % % load EXPT
    % %
    % % for i = 1:length(EXPT.subjects)
    % %     fname = [EXPT.subjects{i} filesep EXPT.subjects{i} '_design.mat'];
    % %     clear trialmodels
    % %
    % %     if exist(fname, 'file')
    % %
    % %         load(fname)
    % %
    % %         [trialw, allw, wh_bad, ortho_trialX] = single_trial_weights(trialmodels.trialX, 1:trialmodels.ntrials * size(trialmodels.bf, 2), trialmodels.bf);
    % %
    % %         title(fname)
    % %         orient(gcf, 'portrait')
    % %         print('-dpsc2', '-append', 'trial_weight_images');
    % %
    % %
    % %         trial_weights = struct('name', fname, 'trialw', trialw, 'allw', allw, 'wh_bad', wh_bad, 'ortho_trialX', ortho_trialX);
    % %
    % %         savename = [EXPT.subjects{i} filesep EXPT.subjects{i} '_trial_weights'];
    % %         disp(['Saving: ' savename]);
    % %         save(savename, 'trial_weights');
    % %
    % %     else
    % %         disp([fname ' is not a valid file!']);
    % %     end
    % %
    % % end

    % Notes:
    % betas, t, p values are invariant to scale of weights.
    % sigma is not.
    % Therefore, in calculating weights for next level, we can ignore sigma,
    % which is a scaling factor, and we'll still get the right betas, t-values,
    % p-values

    % Two equivalent ways of getting stes for betas:
    % sigma = 3.7;
    % v1 = (  sigma^2 * diag(inv(trialX' * trialX)) ) .^ .5
    % v2 =  sigma * ( diag(inv(trialX' * trialX)) ) .^ .5
    % v3 =  3.7 * sqrt ( diag(inv(trialX' * trialX)) )

    % W = 1  ./ (sigma * sqrt ( diag(inv(trialX' * trialX)) )
    % So we can use 1 ./ sqrt ( diag(inv(trialX' * trialX)) ), which equals the
    % original W up to a scaling factor of 1/sigma
    % This is convenient because we can ignore sigma, and the weights depend
    % purely on the design matrix (X) and are the same across the brain.


    if ischar(bf)
        basis_set_name = bf;
        [dummy, dummy, bf] = trial_level_beta3('onsets', [], 'output', 'design', 'basistype', basis_set_name);
    end

    % Remove pairwise redundant columns
    % These do not pose a problem for model fitting with pinv if they are
    % nuisance regressors, but they do mess up weights.
    % ----------------------------------------------------
    rk = rank(trialX);
    if rk < size(trialX, 2)
        fprintf('Warning: Design matrix is rank deficient. Only OK if redundancy is in nuisance covs and using pinv to fit.\n');

        c = corrcoef(trialX);
        c = c - eye(size(c));
        [rows, cols] = find(triu(abs(c) > (1 - 100*eps)));
        for i = 1:length(cols)
            fprintf('Column %3.0f is redundant with %3.0f.  Removing it.\n', cols(i), rows(i));
        end
        trialX(:, cols) = [];

        rk = rank(trialX);
        if rk < size(trialX, 2)
            fprintf('Warning: Design matrix is still rank deficient. I could not fix it.\n');
        end

    end

    allw =  1 ./ (sqrt ( diag(inv(trialX' * trialX)) ));

    nbf = size(bf, 2);
    windx = 1;

    % Save weights for the trial estimates of interest
    allw = allw(trial_beta_indices);
    lastb = length(allw);

    for i = 1:nbf:lastb

        trialw(windx) = min(allw(i : i + nbf - 1));

        windx = windx + 1;

    end

    trialw = trialw ./ sum(trialw);
    create_figure('Trial Weights', 2, 1) ;
    plot(trialw, 'ko', 'MarkerFaceColor', [.3 0 .7], 'MarkerSize', 12, 'LineWidth', 2)
    ylabel('Trial weights (scaled to sum to 1)');
    xlabel('Trial Number');
    title('Trial Weights');
    drawnow

    % Stuff for getting Variance Inflation Factors
    % Must make sure trial basis sets are orthogonal, because they will
    % otherwise enhance VIFS.  in modeling, BFs combine to produce best-fitting
    % estimate, so it is not necessary that we orthogonalize-- the fits, and
    % thus the trial amplitudes, will be the same either way.  But the VIFs
    % will be large if trial BFs are not orthogonalized.

    disp('Orthogonalizing trial basis sets');
    for i = 1:nbf:lastb
        fprintf('%3.0f ', i);

        tx = trialX(:, i : i + nbf - 1);

        % Orthogonalize each BF to preceding ones
        for j = 2:size(tx, 2)

            prevcols = tx(:, 1 : j - 1);
            r = tx(:, j) - prevcols * pinv(prevcols) * tx(:, j);

            tx(:, j) = r;

        end

        trialX(:, i : i + nbf - 1) = tx;

    end

    fprintf('\n');



    disp('Getting VIFs...wait a minute...');
    vifs = getvif(trialX, 1);
    vifs = vifs(trial_beta_indices);  % need to use ALL to compute, because they depend on ALL regressors in model

    % Save max VIF for the trial estimates of interest
    windx = 1;
    for i = 1:nbf:lastb

        trialvifs(windx) = max(vifs(i : i + nbf - 1));

        windx = windx + 1;

    end

    wh_bad = trialvifs > 2.5;

    subplot(2, 1, 2)
    plot(trialvifs, 'ko', 'MarkerFaceColor', [.3 0 .9], 'MarkerSize', 12, 'LineWidth', 2)
    ylabel('Trial Max VIF');
    title('Trial Variance Inflation Factors')

    if any(wh_bad)
        plot(find(wh_bad), trialvifs(wh_bad), 'ko', 'MarkerFaceColor', [.8 0 .3], 'MarkerSize', 12, 'LineWidth', 2)
        legend({'VIF < 2.5' 'VIF > 2.5'});
    else
        legend({'VIF < 2.5'});
    end

    han = plot_horizontal_line(1, 'k--');
    han = plot_horizontal_line(2.5, 'r--');

    ortho_trialX = trialX;

end

