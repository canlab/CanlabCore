function nuisance_cov_estimates(X, images, SETUP, varargin)
% Estimates F-map for model parameters of interest
% Writes to disk:
% 'F_cols_of_interest.img' 'p_cols_of_interest.img' (3-D images)
% 'resid_full_model.img' (a 4-D image)
%
% Locates voxels whose activity is unrelated to the model
%
% Extracts principal components from these voxels for use as covariates
% in subsequent models
%
% Note: invalidates statistical inference in subsequent models based on
% these data, but may improve predictive accuracy and single-trial model.
%
% :Usage:
% ::
%
%     nuisance_cov_estimates(X, images, SETUP, varargin)
%
% :Defining the SETUP structure with inputs:
%
% SETUP.(fields)
%
%   **.wh_of_interest:**
%        vector of which columns of X matrix are of interest
%                       % columns of interest in X matrix; tests var explained by these with F-test
%
%   **.mask:**
%        name of mask image
%
%   **.TR:**
%        repetition time of volume (image) acquisition
%
%   **.HPlength:**
%        high-pass filter length, in s
%
%   **.scans_per_session:**
%        vector of # volumes in each run, e.g., [128 128 128 128 128]
%
%   **.dummyscans:**
%        indices of images in each run that will be modeled
%        with separate dummy variables
%
%   **.startslice:**
%        start at slice #...
%
% :SETUP Optional Inputs:
%   **nopreproc:**
%        to skip preprocessing (i.e., for trial-level inputs)
%
% ..
%    Tor Wager, Nov 07
% ..

    max_to_save = 10;

    % ---------------------------------------------------------------------
    % Set up preprocessing
    % To skip, enter 'nopreproc' as var. arg.
    % ---------------------------------------------------------------------

    [preprochandle, SETUP] = filter_setup(SETUP, X, varargin{:});

    % ---------------------------------------------------------------------
    % Set up F-test
    % F-test: set up anonymous function for one voxel
    % ---------------------------------------------------------------------

    [fhandle, SETUP] = f_test_setup(X, SETUP, images, preprochandle);

    save nuisance_SETUP SETUP

    % ---------------------------------------------------------------------
    % Run preprocessing and analysis
    % ---------------------------------------------------------------------
    print_banner('Finding voxels that are not task-related.');

    if ~isfield(SETUP, 'startslice') || isempty(SETUP.startslice), SETUP.startslice = 1; end

    [v1, v2, v3] = image_eval_function(images, fhandle, 'mask', SETUP.mask, 'preprochandle', preprochandle, 'outnames', SETUP.names, 'start', SETUP.startslice);
    v1; v2; v3;

    % save results
    fprintf('Printing results in spm .ps file');
    spm_image('init', 'F_cols_of_interest.img');
    spm_orthviews_name_axis('F values for X of interest', 1);
    spm_print;

    % ---------------------------------------------------------------------
    % Get voxels that are not task-responsive
    % ---------------------------------------------------------------------
    [yy, xyz, wh_not_modeled, volInfo] = get_non_task_responsive_voxels(images);


    % ---------------------------------------------------------------------
    % Get principal components of these voxels
    % and save them.
    % ---------------------------------------------------------------------

    [noise_components, pc, latent, num_to_save] = get_pcs(yy, max_to_save);

    NOISE = struct('xyz', xyz, 'images', images, 'pc', pc, 'latent', latent, 'wh_not_modeled', wh_not_modeled, 'noise_components', noise_components, 'num_to_save', num_to_save);

    % rsquare for each predictor
    % --------------------------------------------------------------------
    for i = 1:size(X, 2)
        [b, bint, r, rint, stats] = regress(X(:,i), [noise_components ones(size(X, 1), 1)]);
        NOISE.descrip = 'Regression of noise components on columns of design matrix:';
        NOISE.noise_rsquare_for_each_predictor(i) = stats(1);
        NOISE.noise_pval_for_each_predictor(i) = stats(3);
    end

    save nuisance_SETUP -append NOISE

    % spatial maps of noise components
    % --------------------------------------------------------------------
    spatial_map_names = [];
    for i = 1:num_to_save

        dat = zeros(volInfo.n_inmask, 1);
        dat(wh_not_modeled) = pc(:,i);
        name = sprintf('noise_maps.img, %s', num2str(i));
        iimg_reconstruct_3dvol(dat, volInfo, 'outname', name, 'descrip', 'Created by nuisance_cov_estimates.m' );
        fprintf('Writing: %s\n', name);
        spatial_map_names = strvcat(spatial_map_names, name);
    end

    NOISE.spatial_map_names = spatial_map_names;
    save nuisance_SETUP -append NOISE

    % write adjusted (denoised) images: residuals from nuisance params
    SETUP = remove_noise_components(SETUP, NOISE);
    save nuisance_SETUP -append SETUP


    % show figures of results
    % --------------------------------------------------------------------

    X = SETUP.data.X;

    create_figure('Eigenvariate Plot', NOISE.num_to_save, 1)
    for i = 1:NOISE.num_to_save, subplot(NOISE.num_to_save, 1, i); plot(NOISE.noise_components(:,i)); title(['Comp' num2str(i)]); axis auto; axis off; end

    % save results
    scn_export_papersetup(600);
    spm_print('Eigenvariate Plot');

    create_figure('Rsquare'); bar(NOISE.noise_rsquare_for_each_predictor); set(gca, 'XTick', 0:size(X, 2) + 1)
    title('Var explained by nuisance covs in each column of X');

    % save results
    scn_export_papersetup(600);
    spm_print('Rsquare');

    spm_check_registration(NOISE.spatial_map_names);

    global st

    for i = 1:NOISE.num_to_save
        spm_orthviews_name_axis(['Comp' num2str(i)], i);

        vv = spm_read_vols(spm_vol(NOISE.spatial_map_names(i,:)));
        vv = vv(:); mn = mean(vv) - 3 * std(vv); mx = mean(vv) + 3 * std(vv);
        spm_orthviews('Window', i, [mn mx])

        axes(st.vols{i}.ax{2}.ax)
        title(sprintf('%3.4f %3.4f', mn, mx), 'FontSize', 12);

    end

    % save results
    spm_print;
end  % main function





% --------------------------------------------------------------------
% --------------------------------------------------------------------

% sub-functions

% --------------------------------------------------------------------
% --------------------------------------------------------------------


% ---------------------------------------------------------------------
% STEP 1
% ---------------------------------------------------------------------
function [fhandle, SETUP] = f_test_setup(X, SETUP, images, preprochandle)
    % reduced model
    Xred = X; Xred(:,SETUP.wh_of_interest) = [];

    % computational work
    px = pinv(X); pxred = pinv(Xred);

    fhandle = @(y) F_test_full_vs_red(y, X, Xred, px, pxred);

    try
        disp('Testing analysis function.');
        fhandle(randn(size(images, 1), 1));
    catch
        disp('Error in function eval handle: Check function, inputs, and sizes.');
        rethrow(lasterr);
    end


    SETUP.names = {'F_cols_of_interest' 'p_cols_of_interest' 'resid_full_model'};

    SETUP.preprochandle = preprochandle;
    SETUP.fhandle = fhandle;

    SETUP.data.descrip = 'Data in noise component estimation function';
    SETUP.data.X = X;
    SETUP.data.Xred = Xred;
    SETUP.data.images = images;
end


% ---------------------------------------------------------------------
% STEP 2
% ---------------------------------------------------------------------
function [yy, xyz, wh_not_modeled, volInfo] = get_non_task_responsive_voxels(images)
    name = 'p_cols_of_interest.img';
    [volInfo, pvals] = iimg_read_img(name, 2);
    pvals = pvals(volInfo.image_indx);
    wh_not_modeled = pvals > .05;
    n = sum(wh_not_modeled);

    iimg_reconstruct_3dvol(wh_not_modeled, volInfo, 'outname', 'mask_of_non_task_related.img');
    spm_image('init', 'mask_of_non_task_related.img');
    spm_orthviews_name_axis('Non-task related voxels', 1);
    spm_print;

    fprintf('Non-zero p-values: %3.0f\n Voxels with little variability explained by model: %3.0f \n', volInfo.n_inmask, n);
    if n < 1
        disp('Not enough voxels. Exiting.');
    end

    % get voxel coordinates
    xyz = volInfo.xyzlist(wh_not_modeled,:);
    xyz(:,end+1) = 1; xyz = xyz';

    % get data from residual images: residualized wrt task and other nuisance
    % covariates.  this makes pcs non-redundant with other filtering and
    % nuisance covs.

    fprintf('Counting residual images...');
    tic
    n = scn_num_volumes('resid_full_model.img');
    resid_imgs = expand_4d_filenames('resid_full_model.img', n);
    fprintf('Found %3.0f, done in %3.0f sec\n', n, toc);

    fprintf('Extracting data for non-responsive voxels...');
    tic
    yy = spm_get_data(resid_imgs, xyz);
    fprintf('%3.0f sec\n', toc);
end


% ---------------------------------------------------------------------
% STEP 3
% ---------------------------------------------------------------------
function [noise_components, pc, latent, num_to_save] = get_pcs(yy, max_to_save)
    print_banner('Getting principal components from non-task-related voxels.');

    [pc, score, latent] = princomp(yy, 'econ');

    % Eigenvalue plot
    create_figure('Eigenvalue Plot'), bar(latent);
    xlabel('Components'), ylabel('Variance')
    scn_export_papersetup(500);
    set(gca, 'XLim', [0 50]);
    drawnow
    spm_print('Eigenvalue Plot')

    num_to_save = min(max_to_save, sum(latent > 1));
    if num_to_save == 0
        disp('Uh-oh! No components to save.  data problem.');
        keyboard
    end

    num_to_save = max(num_to_save, 1);

    fprintf('Saving %3.0f components in NOISE struct in nuisance_SETUP.mat\n', num_to_save);
    noise_components = score(:,1:num_to_save);
end


% ---------------------------------------------------------------------
% STEP 4
% ---------------------------------------------------------------------
function SETUP = remove_noise_components(SETUP, NOISE)
    % write adjusted (denoised) images: residuals from nuisance params
    print_banner('Creating denoised_images.img');

    X = NOISE.noise_components;

    % add covariates from original X
    covs = SETUP.data.X; covs(:,SETUP.wh_of_interest) = [];
    X = [X covs];
    X(:,end+1) = 1; % may not be necessary, if intercepts already in model, but shouldn't hurt.
    Xred = X(:,end);
    px = pinv(X); pxred = pinv(Xred);

    fhandle = @(y) F_test_full_vs_red(y, X, Xred, px, pxred);
    SETUP.denoised_names = {'F_nuisance_covs' 'p_nuisance_covs' 'denoised_images'};
    [v1, v2, v3] = image_eval_function(SETUP.data.images, fhandle, 'mask', SETUP.mask, 'preprochandle', SETUP.preprochandle, 'outnames', SETUP.denoised_names, 'start', SETUP.startslice);

    SETUP.denoising_fhandle = fhandle;
end



% ---------------------------------------------------------------------
% Subfcn: set up
% ---------------------------------------------------------------------
% Set up preprocessing
function [preprochandle, SETUP] = filter_setup(SETUP, X, varargin)

    preprochandle = [];
    wh_elim = [];
    hpflag = 1;  % only does it if requested, though

    for i = 1:length(varargin)
        if ischar(varargin{i})
            switch varargin{i}
                % reserved keywords
                case 'custompreproc'
                    preprochandle = varargin{i + 1};  % e.g., 'custompreproc', @(data) scale(data) for z=scores;

                    hpflag = 0;
                    SETUP.TR = NaN;
                    SETUP.HPlength = [];
                    SETUP.dummyscans = [];
                    wh_elim = i;

                case {'nopreproc'}
                    hpflag = 0;
                    SETUP.preproc = 0;
                    SETUP.TR = NaN;
                    SETUP.HPlength = [];
                    SETUP.dummyscans = [];
                    wh_elim = i;

                    % We need to allow mediation SETUPions here, so eliminate this from list and do not error check here.
                    %otherwise, warning(['Unknown input string SETUPion:' varargin{i}]);
            end
        end
    end

    varargin(wh_elim) = [];

    % required args
    N = {'TR', 'mask', 'scans_per_session', 'preproc', 'HPlength', 'dummyscans', 'startslice'};
    %N = fieldnames(SETUP);
    for i = 1:length(N)
        if ~isfield(SETUP, N{i}) || isempty(SETUP.(N{i}))
            switch N{i}
                case {'TR', 'mask', 'scans_per_session', 'preproc'}
                    error(['Enter SETUP.' N{i}]);

                case 'HPlength'
                    SETUP.(N{i}) = [];

                case 'dummyscans'
                    SETUP.(N{i}) = 1:2;

                case 'startslice'
                    SETUP.startslice = 1;

                otherwise
                    disp('Warning! Unrecognized field in SETUPions structure SETUP.');
            end
        end
    end

    if SETUP.preproc && hpflag
        [tmp, I, S] = hpfilter(X(:,1), SETUP.TR, SETUP.HPlength, SETUP.scans_per_session, SETUP.dummyscans); % creates intercept and smoothing matrices

        preprochandle = @(Y) hpfilter(Y, [], S, SETUP.scans_per_session, I);  % function handle with embedded fixed inputs

        % test preprocessing
        try
            disp('Testing preprocessing function.');
            test = preprochandle(X(:,1));
        catch
            disp('Testing of preprocessing failed!  User should diagnose this error.');
            rethrow(lasterr);
        end

        if any(isnan(test))
            disp('Warning! Some preprocessed data values in test are NaN.  Problem with model or preproc function?');
        end
    end
end



function print_banner(str)
    fprintf('\n---------------------------------------------------------------------\n')
    fprintf('%s', str)
    fprintf('\n---------------------------------------------------------------------\n')
end

