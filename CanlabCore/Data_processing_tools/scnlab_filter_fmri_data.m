function names = scnlab_filter_fmri_data(imgs, mvmt, mask, tr, spersess, hp)
% Outlier and artifact removal for one subject
% Writes new output images for timeseries
%
% :Usage:
% ::
%
%     names = scnlab_filter_fmri_data(imgs, mvmt, mask, tr, spersess, hp)
%
% :Examples:
% ::
%
%    OUT = mean_image(imgs, 'mean_ravol.img',ones(size(imgs,1),1));
%    spm_imcalc_ui('mean_ravol.img', 'graymatter.img', 'i1 > 0');
%    spm_image('init', 'graymatter.img');
%    names = scnlab_filter_fmri_data(imgs, mvmt, 'graymatter.img', 2, repmat(184, 1, 6), 80);

    tic
    imgs = check_valid_imagename(imgs);

    disp('Loading data')
    [dat, volInfo] = iimg_get_data(mask, imgs);


    % -------------------------------------------------------

    disp('Setting up filtering, and logging to scn_filter_output.txt')
    diary scn_filter_output.txt
    OPT = scnlab_outlier_id('setup', 'tr', tr, 'spersess', spersess, 'dummy', 1:2, 'hp', hp, 'mad', 4, 'niter', 5, 'mvmt', mvmt);
    OPT.volInfo = volInfo;
    diary off

    % -------------------------------------------------------

    disp('Test run...')
    fhandle = @(y) scnlab_outlier_id('data', y, 'options', OPT);
    fhandle(dat(:,1));

    disp('Saving OPT structure in scnlab_filter_setup.mat')
    save scnlab_filter_setup OPT

    disp('Printing graphics image to scnlab_filter_output.ps');
    h = findobj('Tag', 'Data Detail');
    if ~isempty(h) && ishandle(h)
        print(h, '-dpsc','-r200','scnlab_filter_output.ps','-append')
    end

    % -------------------------------------------------------

    disp('Running on all data.')
    OPT.doplot = 0;
    OPT.verbose = 0;
    fhandle = @(y) scnlab_outlier_id('data', y, 'options', OPT);
    [filt_dat, outliers, num_outliers, mvmt_rsquare, mvmt_baseline_rsquare, ypercent_change, raw_equalvar_p, ...
        raw_equalvar_F, perc_equalvar_p, perc_equalvar_F ] = matrix_eval_function(dat, fhandle);
    clear ybaseline


    ypercent_change = ypercent_change';

    filt_dat = filt_dat';
    avg = mean(filt_dat);
    mabsdev = mad(filt_dat);


    % -------------------------------------------------------

    disp('writing mean, MAD, num_outlier, and mvmt_rsquare, etc. images in current directory.')


    iimg_reconstruct_3dvol(mabsdev', volInfo, 'outname', 'scn_filter_MAD.img', 'descrip', 'Created with scnlab_filter_fmri_data.m');
    iimg_reconstruct_3dvol(avg', volInfo, 'outname', 'scn_filter_mean.img', 'descrip', 'Created with scnlab_filter_fmri_data.m');
    iimg_reconstruct_3dvol(num_outliers, volInfo, 'outname', 'scn_filter_num_outliers.img', 'descrip', 'Created with scnlab_filter_fmri_data.m');
    iimg_reconstruct_3dvol(mvmt_rsquare, volInfo, 'outname', 'scn_filter_mvmt_rsquare.img', 'descrip', 'Created with scnlab_filter_fmri_data.m');
    iimg_reconstruct_3dvol(mvmt_baseline_rsquare, volInfo, 'outname', 'scn_filter_mvmt_baseline_rsquare.img', 'descrip', 'Created with scnlab_filter_fmri_data.m');

    iimg_reconstruct_3dvol(raw_equalvar_p, volInfo, 'outname', 'scn_filter_raw_equalvar_p.img', 'descrip', 'Created with scnlab_filter_fmri_data.m');
    iimg_reconstruct_3dvol(raw_equalvar_F, volInfo, 'outname', 'scn_filter_raw_equalvar_F.img', 'descrip', 'Created with scnlab_filter_fmri_data.m');
    iimg_reconstruct_3dvol(perc_equalvar_p, volInfo, 'outname', 'scn_filter_perc_equalvar_p.img', 'descrip', 'Created with scnlab_filter_fmri_data.m');
    iimg_reconstruct_3dvol(perc_equalvar_F, volInfo, 'outname', 'scn_filter_perc_equalvar_F.img', 'descrip', 'Created with scnlab_filter_fmri_data.m');

    rawvsperc = perc_equalvar_F - raw_equalvar_F;
    iimg_reconstruct_3dvol(rawvsperc, volInfo, 'outname', 'scn_filter_perc_vs_raw_equalvar_F.img', 'descrip', 'Created with scnlab_filter_fmri_data.m');

    spm_check_registration(char('scn_filter_mean.img', 'scn_filter_MAD.img', 'scn_filter_mvmt_rsquare.img', 'scn_filter_mvmt_baseline_rsquare.img', ...
        'scn_filter_num_outliers.img', ...
        'scn_filter_raw_equalvar_F.img', 'scn_filter_perc_equalvar_F.img', 'scn_filter_perc_vs_raw_equalvar_F.img'));
    scale_windows



    disp('Printing graphics image to scnlab_filter_output.ps');
    h = findobj('Tag', 'Graphics');
    if ~isempty(h) && ishandle(h)
        print(h, '-dpsc','-r200','scnlab_filter_output.ps','-append')
    end

    % -------------------------------------------------------
    dowrite = 1;
    names = [];

    if dowrite
        disp('Writing f* time series images in original image directory.')
        for i = 1:size(filt_dat, 1)
            if i == 1
                names = write_img(filt_dat(i, :)', imgs(i, :), volInfo, 'f');
            else

                names = char(names, write_img(filt_dat(i, :)', imgs(i, :), volInfo, 'f'));
            end
        end

        if OPT.dopercent
            disp('Writing perc* time series images in original image directory.')
            for i = 1:size(filt_dat, 1)
                if i == 1
                    names = write_img(ypercent_change(i, :)', imgs(i, :), volInfo, 'perc');
                else

                    names = char(names, write_img(ypercent_change(i, :)', imgs(i, :), volInfo, 'perc'));
                end
            end
        end

        disp('Saving filtered image names as ''names'' in scnlab_filter_setup.mat');
        save scnlab_filter_setup -append names
    else
        disp('Not writing any images: flag is off.')
    end

    % -------------------------------------------------------

    disp('All done.')
    toc

end





% -------------------------------------------------------
% -------------------------------------------------------
% -------------------------------------------------------


function name = write_img(y, name, volInfo, prefix)

    name = make_img_filename(name, prefix);
    iimg_reconstruct_3dvol(y, volInfo, 'outname', name, 'descrip', 'Created with scnlab_filter_fmri_data.m');

end


function name = make_img_filename(name, prefix)
    [d, name] = fileparts(name);
    name(name == '.') = '_';
    name = fullfile(d, [prefix name '.img']);
end


function scale_windows

    global st

    warning off
    
    win = 5;
    v = spm_read_vols(spm_vol(st.vols{win}.fname));
    v = v(:);
    v(v == 0) = [];
    lb = min(v);
    ub = prctile(v, 95);
    spm_orthviews('Window', 3, [lb ub])
    title(st.vols{win}.ax{2}.ax, sprintf('%s\nDisplay Rng: %3.0f, %3.0f',st.vols{win}.fname,lb, ub),'FontSize', 12)
    title(st.vols{win}.ax{3}.ax, sprintf('Mean: %3.0f\nRange: %3.0f to %3.0f', mean(v), min(v),max(v)),'FontSize', 10)


    for win = [ 2 3 4 6 7 ]
        v = spm_read_vols(spm_vol(st.vols{win}.fname));
        v = v(:);
        v(v == 0) = [];
        lb = 0;
        ub = prctile(v, 95);
        spm_orthviews('Window', win, [lb ub])
        title(st.vols{win}.ax{2}.ax, sprintf('%s\nDspRng: %3.4f, %3.4f',st.vols{win}.fname,lb, ub),'FontSize', 12)
        title(st.vols{win}.ax{3}.ax, sprintf('Mean: %3.4f\nRng: %3.4f to %3.4f', mean(v), min(v),max(v)),'FontSize', 10)
    end

    for win = [ 1 8 ]
        v = spm_read_vols(spm_vol(st.vols{win}.fname));
        v = v(:);
        v(v == 0) = [];
        lb = mean(v)-3*std(v);
        ub = mean(v)+3*std(v);
        spm_orthviews('Window', win, [lb ub])
        title(st.vols{win}.ax{2}.ax, sprintf('%s\nDspRng: %3.4f, %3.4f',st.vols{win}.fname,lb, ub),'FontSize', 12)
        title(st.vols{win}.ax{3}.ax, sprintf('Mean: %3.4f\nRng: %3.4f to %3.4f', mean(v), min(v),max(v)),'FontSize', 10)
    end

    spm_orthviews('interp',0)
    colormap jet

    warning on
end
