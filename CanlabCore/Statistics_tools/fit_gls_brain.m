function fit_gls_brain(imgs, X, arorder, conditionnames, maskimg, varargin)
% Run Generalized Least Squares model with AR-model autoregression (optional)
% at each voxel in a set of images.
% 
% :Usage:
% ::
%
%     fit_gls_brain(imgs, X, arorder, conditionnames, maskimg, ['contrasts', contrast_mtx, 'contrastnames', contrastnames], ['weights', weights])
%
% Single trial model to get trial amplitude, width, AUC, etc.
% Take one of those (e.g., AUC) and get a list of trial images
% Pass that into this function.
%
% fit_gls_brain(trial_amp_imgs, eventdesign{1}, contrasts
%
% You must add the intercept to X yourself if you want one!
%
% contrasts can be empty, and if it is, this function will write images
% for statistics on individual predictors.  If you enter contrast values,
% then it will write images for contrasts instead.
%
% conditionnames = column or contrast image names, in cell array
%
% e.g. eventnames{1} = {'high'    'medium'    'low'    'warm'}
%
% :Examples:
%
% This one looks @ significance for betas in X model:
% ::
%
%    load ../Multilev_mediation-try4(resliced)_10k/mediation_SETUP.mat
%    imgs = SETUP.data.M{1};
%    X = eventdesign{1};
%    arorder = 1;
%    contrasts = [];
%    conditionnames = eventnames{1};
%    maskimg = spm_select(1); % try gray matter mask...
%    fit_gls_brain(imgs, X, arorder, contrasts, conditionnames, maskimg)
%
% Now define contrasts and re-run on contrast values:
% ::
%
%    contrasts = [3 1 -1 -3; .25 .25 .25 .25; 1 0 0 -1]'
%    conditionnames = {'Linearpain' 'Average_resp' 'High-Low'};
%    fit_gls_brain(imgs, X, arorder, contrasts, conditionnames, maskimg)


    % ..
    %    Optional inputs
    % ..
    nobs = size(imgs, 1);
    weights = ones(nobs, 1);        % default weights
    contrasts = [];
    contrastnames = {};

    for i = 1:length(varargin)
        if ischar(varargin{i})
            switch varargin{i}
                % functional commands
                case 'contrasts', contrasts = varargin{i+1};

                case 'weights', weights = varargin{i+1};

                case 'contrastnames', contrastnames = varargin{i+1};

                otherwise, warning(['Unknown input string option:' varargin{i}]);
            end
        end
    end
    
    % Check weights
    wtol = max(sqrt(weights))*eps^(1/3);
    if any(weights == 0) || any(weights < wtol)
        fprintf('Some weights are too low or badly scaled. Setting zero weights to %3.4f\n.', wtol);
        weights(weights < wtol) = wtol;
    end

    W = diag(weights);

    % [t, df, beta/contrast, Phi, sigma,stebeta, F] =
    % fit_gls(y,X,c,p,[PX, equal to pinv(X), for speed])


    % Map mask image into space of imgs, if not already, and write
    % 'mask.img'
    scn_map_image(maskimg, imgs(1, :), 'write', 'mask.img');
    maskimg = 'mask.img';
    
    % Setup inputs and function handle
    % -------------------------------------------------------------
    wh = intercept(X, 'which');
    if isempty(wh), disp('Warning! No intercept found in model.'); end

    %PX = pinv(X);
    PX = inv(X' * W * X) * X' * W;  % General weighted
    
    % Check if X is OK
    condx = cond(X);
    fprintf('Condition number of design matrix: %3.2f\n', condx);
    if condx > 10000, error('Design matrix is not estimable!'); end
    
    fhandle = @(y) fit_gls(y, X, contrasts, arorder, PX, weights);


    % Setup output image names
    % -------------------------------------------------------------
    % [beta, t, pvals, convals, con_t, con_pvals, sigma, Phi, df, stebeta, conste, F]
    
    %if isempty(contrasts), betaconname = 'beta'; else betaconname = 'contrast'; end

    names = cell(1, 12);
    names{7} = 'sigma.img';
    names{8} = 'phi.img';
    names{9} = 'df.img';
    names{12} = 'omnibus_F.img';

    for i = [1:6 10:11]  % Potential Multi-image outputs
        for j = 1:length(conditionnames)
            switch i
                case 1
                    names{i} = char(names{i}, ['beta_' conditionnames{j} '.img']);
                case 2
                    names{i} = char(names{i}, ['t_' conditionnames{j} '.img']);
                case 3
                    names{i} = char(names{i}, ['p_' conditionnames{j} '.img']);



                case 10
                    names{i} = char(names{i}, ['ste_beta_' conditionnames{j} '.img']);

            end
        end

        for j = 1:length(contrastnames)
            switch i
                case 4
                    names{i} = char(names{i}, ['con_beta_' contrastnames{j} '.img']);
                case 5
                    names{i} = char(names{i}, ['con_t_' contrastnames{j} '.img']);
                case 6
                    names{i} = char(names{i}, ['con_p_' contrastnames{j} '.img']);
                case 11
                    names{i} = char(names{i}, ['ste_con_' contrastnames{j} '.img']);
            end
        end

    end
    
    for i = [1:6 10:11]
        names{i} = names{i}(2:end, :);
    end
    
    save fit_gls_brain_SETUP imgs X arorder contrasts names conditionnames contrastnames maskimg

    % Run
    % -------------------------------------------------------------
    [beta, t, pvals, convals, con_t, con_pvals, sigma, Phi, df, stebeta, conste, F] = ...
        image_eval_function(imgs, fhandle, 'mask', maskimg, 'outimagelabels' , names);
    
    % other optional args: 'outimagelabels' , 'connames'

end




% % %
% % %
% % %             [t, df, beta, Phi, sigma, stebeta] = fit_gls(h, trial2ndX, [], 1, trial2ndpx);
% % %
% % %             conditiont(i(k), j(k), :) = t;
% % %             conditionste(i(k), j(k), :) = stebeta;
% % %             conditionmean(i(k), j(k), :) = beta;
% % %             trialphi(i(k), j(k), :) = Phi;
% % %
% % %             conditionp(i(k), j(k), :) = 2 .* (1 - tcdf(abs(t), df));
% % %
% % %             if ~isempty(c)
% % %                 stecons = diag(c' * diag(stebeta).^2 * c).^.5;
% % %                 con_est = c' * beta;
% % %                 tcons = con_est ./ stecons;
% % %
% % %                 contrastest(i(k), j(k), :) = con_est;
% % %                 contrastt(i(k), j(k), :) = tcons;
% % %                 contrastste(i(k), j(k), :) = stecons;
% % %                 contrastp(i(k), j(k), :) = 2 .* (1 - tcdf(abs(tcons), df));
% % %             end
% % %
% % %
% % %             % 2nd-level estimation of ratings model with AR(1) model
% % %             % use first param, ignore intercept
% % %             % -------------------------------------------------------------
% % %             if exist('ratings','var')
% % %                 [t, df, beta, Phi, sigma, stebeta] = fit_gls(h, ratings, [], 1, ratingspx);
% % %
% % %                 % 1st col. is rating effect
% % %                 ratingsest(i(k), j(k), 1) = beta(1);
% % %                 ratingst(i(k), j(k), 1) = t(1);
% % %                 ratingsste(i(k), j(k), 1) = stebeta(1);
% % %                 ratingsp(i(k), j(k), 1) = 2 .* (1 - tcdf(abs(t(1)), df));
% % %
% % %                 % 2nd column is intercept
% % %                 intcptest(i(k), j(k), 1) = beta(2);
% % %                 intcptt(i(k), j(k), 1) = t(2);
% % %                 intcptste(i(k), j(k), 1) = stebeta(2);
% % %                 intcptp(i(k), j(k), 1) = 2 .* (1 - tcdf(abs(t(2)), df));
% % %
% % %             end
% % %         end
% % %
% % %         if k == 1000, fprintf(1, '%3.0f s per 1000 vox.', etime(clock, et)), end
