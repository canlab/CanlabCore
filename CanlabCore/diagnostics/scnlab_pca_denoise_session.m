% ..
%    scnlab_pca_denoise_session(image_names, basename, nuisX, designX, [mask name])
%
%    Tor Wager, Feb 2008
% ..

function scnlab_pca_denoise_session(image_names, basename, nuisX, designX, varargin)
    diary('qc_report.txt')

    fprintf('---------------------------------------------------------\n')
    fprintf('PCA denoising and visualization\n');
    fprintf(['Basename: ' basename '\n']);
    fprintf('---------------------------------------------------------\n')

    img_ext = get_image_ext(image_names);

    k = 10; % number of components

    % Create mask and get info
    % -----------------------------------------------------------------
    if nargin < 5 || isempty(varargin)
        disp('Creating mask from data.')
        mask_thresh = fmri_mask_thresh(image_names(1,:));
        mask_thresh = mask_thresh .* .8;
        [maskInfo, dat] = iimg_read_img(expand_4d_filenames(image_names(1,:), 1));
        dat = dat > mask_thresh;

        maskname = ['mask' img_ext];
        iimg_reconstruct_3dvol(dat, maskInfo, 'outname', maskname);
        fprintf('Written %s\n', maskname);
    else
        disp('Reading mask from file.')
        maskname = varargin{1};
        [maskInfo, dat] = iimg_read_img(maskname);
    end


    % read maskInfo for in-mask voxels
    maskInfo = iimg_read_img(maskname, 2);

    create_figure('PCA results', 2, 3);
    montage_image_Worsley(maskname, 'noscale');
    title([basename ' Mask']);

    axis tight
    axis off
    drawnow

    %% Get data and principal components
    % -----------------------------------------------------------------
    fprintf('Reading data. ');
    [dat, volInfo] = iimg_get_data(maskInfo, image_names);

    fprintf('Components. ');
    % We're de-noising so no need to use PCA on correlations.  Covs is probably
    % good...
    %[eigvec, compscore, eigval] = princomp(dat, 'econ');
    [eigvec, compscore, eigval] = princomp_largedata(dat, 'econ');

    subplot(2, 3, 2);  plot(eigval, 'ko-');
    title('Eigenvalues');
    set(gca, 'XLim', [0 min(k + 3, size(compscore, 2))] );
    plot_vertical_line(k + .5);
    drawnow

    eigvec = eigvec(:,1:k);
    compscore = compscore(:,1:k);

    fprintf('Writing pca_spatial%s 4-D images.\n', img_ext);
    if ~strcmpi(spm('Ver'), 'spm5'), disp('Warning! Unless you use SPM5 or above, all components will not be saved/viewable.'); end

    % reconstruction and writing of spatial maps
    iimg_reconstruct_vols(eigvec, maskInfo, 'outname', [basename '_pca_spatial' img_ext]);


    %[eigvec2, compscore2, eigval2] = princomp(zscore(dat), 'econ');

    %% View spatio-temporal components

    create_figure('Spatial and temporal components', 2, 1);
    plot_matrix_cols(compscore);


    subplot(2, 1, 2);
    montage_image_Worsley([basename '_pca_spatial' img_ext], 'pcacov');
    axis tight


    spatialnames = expand_4d_filenames([basename '_pca_spatial' img_ext]);
    for i = 1:size(spatialnames, 1)
        create_figure(['Component ' num2str(i)]);
        montage_image_Worsley(deblank(spatialnames(i,:)), 'pcacov');
    end

    drawnow

    nk = size(compscore, 2);
    if ~isempty(designX)
        create_figure('Components With Fits', nk, 1);
        
        wh_real_regressors = (range(designX) ~= 0); %exclude intercepts
        
        % FITS
        xx = [designX nuisX];
        px = pinv(xx);
        f = zeros(size(compscore));

        for i = 1:nk
            yy = compscore(:,i);

            b = px * yy;
            % f(:,i) = designX * b(1:size(designX, 2));
            f(:,i) = designX(:,wh_real_regressors) * b(wh_real_regressors);

            subplot(nk, 1, i);

            plot(compscore(:,i), 'k');
            plot(f(:,i), 'r');
            axis tight
        end
    end

    subplot(nk, 1, 1);
    title('Components with task fits');
    print('-dpsc2', '-append', 'qc_report');

    %% Must select components to exclude here. If we have nuisance (movement,
    % physio) and task-related stuff, show this now.

    [nuisance_ratio, rsquare_design, rsquare_nuis] = scn_component_rsquare(compscore, nuisX, designX);

    wh_rem = input('Enter components to remove from data in [ ], or return to remove none: ');

    %
    %% reconstruction and writing of residual images (Selected components removed)
    n = size(dat, 1);

    %reconstructed = repmat(mean(dat, 1), n, 1) + compscore(:,1:ndim) * eigvec(:,1:ndim)';

    %leave mean in
    reconstructed = compscore(:,wh_rem) * eigvec(:,wh_rem)';
    residuals = dat - reconstructed;

    create_figure('PCA results', 2, 3, 1);
    subplot(2, 3, 3);
    text(.1, .9, {'Task R^2'})
    text(.1, .5, num2str(rsquare_design(wh_rem)))
    text(.5, .9, {'Nuis R^2'})
    text(.5, .5, num2str(rsquare_nuis(wh_rem)))
    text(.8, .9, {'Ratio R^2'})
    text(.8, .5, num2str(nuisance_ratio(wh_rem)))
    
    create_figure('PCA results', 2, 3, 1);
    subplot(2, 3, 4);
    imagesc(dat')
    colorbar
    axis tight
    xlabel('Time'); ylabel('Voxels');
    title('Data without denoising');
    
    subplot(2, 3, 5);
    set(gca, 'YDir', 'Reverse')
    imagesc(reconstructed')
    colorbar
    axis tight
    xlabel('Time'); ylabel('Voxels')
    title('Removed from data');

    subplot(2, 3, 6);
    imagesc(residuals')
    colorbar
    axis tight
    xlabel('Time'); ylabel('Voxels');
    title('Adjusted data');

    print('-dpsc2', '-append', 'qc_report');
    scn_export_papersetup(600);
    saveas(gcf, [basename 'pca_session'], 'png');

    % Fig of removed components
    create_figure('Removed Components', 3, 1);
    plot_matrix_cols(compscore(:,wh_rem));
    xlabel('Time');
    title([basename ' Removed Components']);
    subplot(3, 1, 2);
    if(~isempty(wh_rem))
        montage_image_Worsley(spatialnames(wh_rem,:), 'pcacov');
    end
    drawnow
    subplot(3, 1, 3);
    plot_matrix_cols(compscore);
    hold on
    plot_matrix_cols(compscore(:,wh_rem),'denoising');
    xlabel('Time');
        axis tight
    title([basename ' Removed Components in Red']);
    print('-dpsc2', '-append', 'qc_report');
    scn_export_papersetup(600);
    saveas(gcf, [basename 'pca_removed'], 'png');

    disp(['Saved variables in .mat file: ' basename '_pca_denoise_output']);
    save([basename '_pca_denoise_output'], 'eig*', 'comp*', 'wh_rem', '*X');

    % write images
    outname = [basename '_denoised' img_ext];
    fprintf('Writing %3.0f volumes in %s\n', n, outname);

    iimg_reconstruct_vols(residuals', maskInfo, 'outname', outname);

    %%
    diary off
end

function img_ext = get_image_ext(image_names)
    image_names = cellstr(image_names);
    [d f img_ext] = fileparts(image_names{1});
end
