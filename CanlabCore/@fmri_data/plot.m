function [wh_outlier_uncorr, wh_outlier_corr, X] = plot(fmridat, varargin)
% Custom plots for fmri_data object. Includes mulitple plots and outlier detection.
%
% :Usage:
% ::
%
% plot(fmri_data_object, [optional inputs])
% Ke edit
% :Inputs:
%
%   **plotmethod:**
%        followed by plot method:
%        'data' [default]
%        'means_for_unique_Y' : Average within Y levels before plot
%
%   **noorthviews:**
%        Skip orthviews plot
%
%   **nooutliers:**
%        Skip outlier ID and plot
%
%   **'montages'**
%       Montage of mean only
%       Shows run-by-run means if .images_per_session has values in object
%
%   **'orthviews'**
%      Orthviews only 
% 
% **'means_for_unique_Y'**
%      Show montages of means for unique values of .Y 
%
% :Outputs:
%
% The 'data' [default] plot method returns a composite of several plots,
% including an outlier analysis (which can be optionally skipped).
% 
% 6 plots and an SPM orthviews presentation of the data.  In the below
% and elsewhere, "image" connotes a 3D brain volume captured every TR.
%
%   **subplot 1:**
%         Data matrix (top left). Also called a "carpet plot" in neuroimaging.
%         The color reflects the intensity of signal in each voxel (column)
%         for each image (row).
%         Here you can look for images that are bright or dark across
%         the image, which would indicate a global shift in values or difference
%         in scale across the images (rows).
%         These can be produced by artifacts that are broadly spatially distributed,
%         or by scanner drift. Most datasets have some of these.
%         This plot can also show that the pattern across voxels in one row (image)
%         may be different than another, highlighting a difference in the spatial pattern.
%         If you are plotting a time series dataset from one participant, unusual images
%         could result from physiological noise, head movement, or other
%         scanner artifacts. They could reflect after-effects of non-linear
%         interactions across these, or a change in the images after a
%         participant has moved their head to a new position.
%         If you are plotting a group dataset with one contrast image per
%         participant (i.e., what you would subject to a group analysis),
%         standard statistical assumptions include that observations
%         (participants) are all on the same scale, with the same variance.
%         It is common for these assumptions to be violated, and you can
%         sometimes see these violations here.
%         Lastly, the range shown in the color bar on the right side of this plot
%         will be quite large if there are large outliers in the data.
%         The units are in contrast unit values, but it is important to check if a
%         few extreme values are forcing all other values to be in the middle of the scale.
%         In this case, there will be very little color variation in the plot, which is a
%         "red flag" indicating outliers or extreme values.
%         Other plots show you different representations of this dataset,
%         in ways that make it easier to see some of its properties.
%
%   **subplots 2 and 3:**
%         Covariance and Correlation Matrices (Top Middle/Top Left):
%         These plots both show similarity across images. Both should show bright main
%         diagonals and off-diagonals that are non zero or generally positive, depending
%         on the dataset.
%
%         Covariance Matrix:   $\frac{{X}'X}{n-1}$, X is mean-centered, n is number of rows in X
%         The diagonal reflects the image variance, shows whether the variances
%         of each image are equal. If the variances are all of approximately the same scale,
%         the diagonal will be a single color.
%
%        Correlation Matrix $\frac{{X}'X}{n-1}$, colums of X are Z-scored, n is number of rows in X
%        This plot shows the correlation between each image and the others.
%        Now, the main diagonal will always be one because each image is
%        perfectly correlated with itself. The off-diagonals should be positive
%        if the images are similar to one another.
%
%   **subplot 4:**
%         Histogram (Bottom Left): A histogram of values across all images and voxels.
%         Depending on your input images, low values could reflect out-of-brain voxels,
%         as there is no signal there. Values of exactly 0 are excluded as missing data
%         in all image operations.
%         For a group dataset (e.g., contrast images), it is expected that
%         the distribution of contrast values will have roughly mean 0.
%         There are other tools for looking at this distribution for each individual
%         image in the dataset, and each image x tissue type (assuming
%         MNI-space images). See help fmri_data.histogram.
%
%   **subplot 5:**
%         The Global Mean Values (Bottom Middle) for each image.
%         The means ought to be similar to one another. Ideally, these means
%         will all be in a similar range. The error bars show one standard deviation.
%         each point is an image.  The point's X value is the mean
%         intensity of every voxel in that image, and the Y value is the
%         stdev of intensities for all voxels in that image.
%
%   **subplot 6:**
%         Mahalanobis Distance (Bottom Right) is a measure of how far away each
%         image is from the rest in the sample. This a standard measure of multivariate
%         distance for each of a set of multivariate observations (here, images).
%         The larger the distance, the more dissimilar it is from other images.
%         High values generally indicate extreme values/potential outliers, but
%         in any normally distributed dataset, there are going to be some values
%         that lie farther out. Also note that what to do about extreme values
%         can be complex, and there is much discussion about how to handle
%         them.
%
%         Potential outliers are identified using fmri_data.mahal.
%         To identify outliers, we assume that the points are distributed according
%         to a chi-square distribution. The expected distance is based on multivariate
%         normally distributed data for the percentile of the dataset that corresponds
%         to each image. even the most extreme values may not be greater than what
%         one expects by chance. The analysis produces p-values at uncorrected and
%         Bonferroni-corrected levels, and any image that is marked as an outlier
%         is one that exceeds the expected value with a p-value of less than 0.05.
%         Those will be marked in darker red. Anything considered to be an outlier with
%         p < 0.05 Bonferroni corrected will be marked with an even darker red.
%         The outliers from this plot are returned to the workspace as output.
%         Note: Outlier identification here does not use global values, only Mahalanobis distance.
%
%   **Orthviews:**
%        SPM Orthviews: These use the spm_orthviews function from SPM software.
%        The X,Y, and Z coordinates correspond to the distance in millimeters
%        from the set origin. Three images are shown:
%        - The mean image is the voxel-wise average across images in the dataset
%        - The STD image is the voxel-wise standard deviation across images in the dataset
%        - The mean/STD image is a simple estimate of effect size (Cohen's D) for every voxel.
%          If a contrast with one image per person is passed to plot(), this plot gives you
%          the effect size in a simple group analysis (e.g., one-sample t-test) across
%          the brain. It is reasonable to use this plot to look for distortions such as high values for means and mean/STDs in the ventricles.
%
%   **X:**
%        X is a struct that contains the following fields:
%        - X.r: The correlation matrix
%        - X.p: The p-value matrix that corresponds to the correlation
%           matrix.
%        - X.sig: The logical sognificance matrix that corresponds to a
%        p-values less than .05.
%       X can be passed into such methods such as plot_correlation_matrix()
%       and other methods that use such structs.


% Programmers' notes: tor - 4/6/2018, changed mahalanobis distance method
% to 'corr', using correlation matrix, and changing how many/how prin comps
% are retained before mahal.  Added wh_outlier_uncorr, wh_outlier_corr
% output.
%
% 12/2/2021  Tor : Added some input options and documentation, and
% run-by-run montages of means.
%
% 05/18/2022 Lukas: maximized windows for figures systematically, to
% improve image quality when publishing
%
% 07/02/2023 Michael Sun: Added the X struct as an output.


[wh_outlier_uncorr, wh_outlier_corr] = deal([]);

% -------------------------------------------------------------------------
% OPTIONAL INPUTS
% -------------------------------------------------------------------------

plotmethod = 'data';  % or 'means_for_unique_Y'
doorthviews = true;
dooutliers = true;
dorunmeanmontages = true;

allowable_inputs = {'plotmethod' 'doorthviews' 'dooutliers' 'dorunmeanmontages'};

keyword_inputs = {'noorthviews' 'nooutliers' 'norunmontages' 'nomontages' 'montage' 'montages' 'means_for_unique_Y'};

% optional inputs with default values - each keyword entered will create a variable of the same name

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}

            case allowable_inputs

                eval([varargin{i} ' = varargin{i+1}; varargin{i+1} = [];']);

            case keyword_inputs
                % Skip deal with this below

            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end

% This pattern will flexibly assign arguments based on keywords.
% The names of the input keyword and variable created do not need to match.
% Multiple input keywords can be mapped to the same variable.

% optional inputs with default values
for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}

            case {'noorthviews'}, doorthviews = false;

            case 'nooutliers', dooutliers = false;

            case {'norunmontages', 'nomontages'}, dorunmeanmontages = false;

            case {'montage' 'montages'}, plotmethod = 'montages'; % montages only

            case 'orthviews', plotmethod = 'orthviews'; % orthviews only

            case 'means_for_unique_Y', plotmethod = 'means_for_unique_Y'; 

            case allowable_inputs
                % skip - handled above

            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end



switch plotmethod
    %  ==============================================================
    case 'data'
    %  ==============================================================

        if isempty(fmridat.dat)
            warning('No data in .dat field.');
            return
        end

        create_figure('fmri data matrix', 2, 3);
        set(gcf,'WindowState','maximized'); %@lukasvo76 added to improve figures when publishing

        % center voxels - 9/9/18 Tor
        [~, nimg] = size(fmridat.dat);
        mm = mean(fmridat.dat, 2, 'omitnan'); %  nanmean(fmridat.dat, 2);
        dat_vox_centered = fmridat.dat - repmat(mm, 1, nimg);
        imagesc(dat_vox_centered');
        colorbar;

        axis tight; set(gca, 'YDir', 'Reverse')
        title('fmri data .dat Data matrix');
        xlabel('Voxels'); ylabel('Images');
        drawnow;

        tmp = mean(dat_vox_centered(:), 'omitnan');  % nanmean(dat_vox_centered(:));
        stmp = std(dat_vox_centered(:), 'omitnan');

        % if std == 0, won't work...
        if stmp < 1000*eps, stmp = 1000*eps; end
        myrange = [tmp - 3*stmp tmp + 3*stmp];
        set(gca, 'CLim', myrange);
        drawnow;


        if ~isempty(fmridat.Y)
            p = get(gca, 'Position'); ystart = p(2); ylen = p(4);

            axes('Position', [.05 ystart .03 ylen]);   % axh = axes('Position', [.05 ystart .03 ylen]);
            imagesc(fmridat.Y);
            title('Y');
            axis tight;
        end
        drawnow;

        % ---------------------------------------------------------------
        % Covariance
        % ---------------------------------------------------------------

        subplot(2, 3, 2);
        covmtx = cov(fmridat.dat);
        imagesc(covmtx);
        axis tight; set(gca, 'YDir', 'Reverse');
        title('Spatial covariance across images');
        colorbar;
        
        cm = colormap_tor([0 0 1], [1 0 0], [1 1 1]);
        colormap(gca, cm)
        climval = max(abs(get(gca, 'CLim')));
        set(gca, 'CLim', [-climval climval]);

        drawnow;

        if ~isempty(fmridat.Y)
            p = get(gca, 'Position'); ystart = p(2); ylen = p(4);

            axes('Position', [.05 ystart .03 ylen]); % axh = axes('Position', [.05 ystart .03 ylen]);
            imagesc(fmridat.Y);
            title('Y');
            axis tight;

        end
        drawnow

        % ---------------------------------------------------------------
        % Histogram
        % ---------------------------------------------------------------

        subplot(2, 3, 4);

        if size(fmridat.dat,2) < 100
            % Figure with individual histograms for each image
            % This is slow with many images

            histogram(fmridat, 'nofigure', 'by_image', 'singleaxis');

        else
            histogram(fmridat, 'nofigure');
        end

        drawnow

        clear dattmp

        globalmean = mean(fmridat.dat, 1, 'omitnan');  % global mean of each obs
        globalstd = std(fmridat.dat, 0, 1, 'omitnan');  % global mean of each obs
        nobs = length(globalmean);
        sz = rescale_range(globalstd, [1 6]); % marker size related to global std
        sz(sz == 0) = 1;
        %sz(sz < .5) = .5;

        %         % ---------------------------------------------------------------
        %         % Global mean vs. std
        %         % ---------------------------------------------------------------
        %         subplot(2, 3, 3); hold on;
        %         plot(globalmean, globalstd, 'k.');
        %         title('Image mean vs. std across voxels');
        %         xlabel('Image mean');
        %         ylabel('Image std');

        % Global mean vs. std
        % ---------------------------------------------------------------

        if size(fmridat.dat,2) > 1

            r = corr(double(globalmean'), double(globalstd'));

            mystr = sprintf('Corr between image mean and spatial std: %3.2f', r);
            subplot(2, 3, 2);
            xlabel(mystr)

        end

        % ---------------------------------------------------------------
        % Correlation matrix
        % ---------------------------------------------------------------
        subplot(2, 3, 3); hold on;
        [cormtx p] = corr(fmridat.dat);
        imagesc(cormtx);
        sig=p<.05;

        X.r=cormtx;
        X.p=p;
        X.sig=sig;

        axis tight; set(gca, 'YDir', 'Reverse');
        title('Spatial correlation across images');
        colorbar;
        
        cm = colormap_tor([0 0 1], [1 0 0], [1 1 1]);
        colormap(gca, cm)
        set(gca, 'CLim', [-1 1]);

        drawnow;

        if ~isempty(fmridat.Y)
            p = get(gca, 'Position'); ystart = p(2); ylen = p(4);

            axh = axes('Position', [.05 ystart .03 ylen]);
            imagesc(fmridat.Y);
            title('Y');
            axis tight;

        end
        drawnow

        % ---------------------------------------------------------------
        % Global mean vs. image series (i.e., time)
        % ---------------------------------------------------------------
        if size(fmridat.dat,2) > 1

            subplot(2, 3, 5);  hold on;

            plot(globalmean, '-');
            axis tight

            Y = 1:nobs;
            Yname = 'Case number';

            plot_horizontal_line(mean(globalmean), 'k--');

            for i = 1:nobs
                plot(Y(i), globalmean(i), 'ko', 'MarkerSize', sz(i), 'LineWidth', 1);
            end

            if nobs <= 500
                errorbar(Y, globalmean, globalstd, 'Color', [.1 .3 .6]);
            else
                upperline = globalmean + globalstd;
                lowerline = globalmean - globalstd;
                xdata = [Y fliplr(Y) Y(1)];
                ydata = [upperline fliplr(lowerline) upperline(1)];
                patch(xdata, ydata, 'y', 'linestyle', 'none', 'FaceColor', [.1 .3 .6], 'faceAlpha', .3);
            end

            ylabel('Global mean');
            xlabel([Yname ' (err bars = 1 sd)']);
            title('Global mean values (size = spatial std)');
            axis tight;
            drawnow;


            % ---------------------------------------------------------------
            % Mahalanobis distance
            % ---------------------------------------------------------------

            subplot(2, 3, 6);

            if ~dooutliers
                axis off

            else

                % do without time series - more generic
                [wh_outlier_uncorr, wh_outlier_corr] = outliers(fmridat, 'notimeseries');

            end % dooutliers


        end % of > 1 image

        % [coeff, score, latent] = princomp(fmridat.dat, 'econ');
        % %d2 = mahal(score, score);
        % plot(latent)


        % ---------------------------------------------------------------
        % Orthviews
        % ---------------------------------------------------------------
        if doorthviews

            plot_orthviews(fmridat)

        end % if doorthviews



        % ---------------------------------------------------------------
        % Montages
        % ---------------------------------------------------------------

        if dorunmeanmontages

            montage_run_means(fmridat);

        end


        %  ==============================================================
    case 'orthviews' % Orthviews only
        %  ==============================================================

        plot_orthviews(fmridat)

        %  ==============================================================
    case 'montages' % Montages only
        %  ==============================================================

        montage_run_means(fmridat);

        %  ==============================================================
    case 'means_for_unique_Y'
        %  ==============================================================

        u = unique(fmridat.Y);

        [v, n] = size(fmridat.dat);
        nu = length(u);

        if nu > 20
            error('More than 20 unique values of Y.  For means_by_condition, Y should be discrete integer-valued.');
        end

        [means, stds] = deal(zeros(nu, v));

        for i = 1:nu
            means(i, :) = nanmean(fmridat.dat(:, fmridat.Y == u(i))');
            stds(i, :) = nanstd(fmridat.dat(:, fmridat.Y == u(i))');
        end

        create_figure('means by condition (unique Y values)', 2, 1);
        set(gcf,'WindowState','maximized'); %@lukasvo76 added to improve figures when publishing
        imagesc(means);
        colorbar;
        axis tight; set(gca, 'YDir', 'Reverse');
        title('Means by condition');
        xlabel('Voxels');
        if iscell(fmridat.Y_names) && ~isempty(fmridat.Y_names)
            set(gca, 'YTick', u, 'YTickLabel', fmridat.Y_names);
        else
            ylabel('Unique Y values');
        end

        drawnow;

        subplot(2, 1, 2)
        imagesc(stds);
        colorbar;
        axis tight; set(gca, 'YDir', 'Reverse');
        title('Standard deviations by condition');
        xlabel('Voxels');
        if iscell(fmridat.Y_names) && ~isempty(fmridat.Y_names)
            set(gca, 'YTick', u, 'YTickLabel', fmridat.Y_names);
        else
            ylabel('Unique Y values');
        end
        drawnow

        % ---------------------------------------------------------------
        % Orthviews
        % ---------------------------------------------------------------
        if ~isempty(fmridat.volInfo)

            vecs_to_reconstruct = means';
            m = fmridat;
            m.dat = vecs_to_reconstruct;
            orthviews(m);

            %create_orthviews(vecs_to_reconstruct, fmridat);
            n = size(vecs_to_reconstruct, 2);

            if iscell(fmridat.Y_names) && ~isempty(fmridat.Y_names) && length(fmridat.Y_names) == n
                axnames = fmridat.Y_names;
            else
                for i = 1:n, axnames{i} = sprintf('Y = %3.3f', i); end
            end

            for i = 1:n
                spm_orthviews_name_axis(axnames{i}, i);
            end
            set(gcf, 'Name', 'Orthviews_means_by_unique_Y');


            % ---------------------------------------------------------------
            % Montage: mean across conditions
            % ---------------------------------------------------------------
            vecs_to_reconstruct = mean(means)';
            vecs_to_reconstruct(vecs_to_reconstruct < prctile(vecs_to_reconstruct, 70)) = 0;
            fig_handle = create_montage(vecs_to_reconstruct, fmridat);
            set(fig_handle, 'Name', 'Montage_mean_across_conditions');
            

            vecs_to_reconstruct = std(means)' ./ mean(means)';
            vecs_to_reconstruct(vecs_to_reconstruct < prctile(vecs_to_reconstruct, 70)) = 0;
            fig_handle = create_montage(vecs_to_reconstruct, fmridat);
            set(fig_handle, 'Name', 'Montage_coeff_of_var_across_conditions');


        end


    otherwise
        error('Unknown plot method');
end

end




function fig_handle = create_montage(vecs_to_reconstruct, fmridat)

n = size(vecs_to_reconstruct, 2);
overlay = which('SPM8_colin27T1_seg.img');

for i = 1:n

    dat = vecs_to_reconstruct(:, i);
    % top and bottom 10%
    dat(dat > prctile(dat, 10) & dat < prctile(dat, 90)) = 0;

    cl{i} = iimg_indx2clusters(dat, fmridat.volInfo);

    fig_handle(i) = montage_clusters(overlay, cl{i}, [2 2]);

    set(fig_handle, 'Name', sprintf('Montage %3.0f', i), 'Tag', sprintf('Montage %3.0f', i));

end

end


function rx = rescale_range(x, y)
% re-scale x to range of y
m = range(y)./range(x);

if isinf(m)
    % no range/do not rescale
    rx = x;
else
    x = x - min(x);
    rx = y(1) + x * ((y(2) - y(1)) ./ max(x));
end

end


% ---------------------------------------------------------------
% Montage of run means
% ---------------------------------------------------------------

function montage_run_means(dat)

use_generic_labels = false;

if isempty(dat.images_per_session)
    % They may not really be runs...use generic labels
    use_generic_labels = true;

    dat.images_per_session = size(dat.dat, 2);
end

n_runs = length(dat.images_per_session);
Xi = intercept_model(dat.images_per_session);

create_figure('Slice_montage', n_runs, 1);

% Set figure position
set_figure_position(n_runs / 1, 6)

% Create run means
% ----------------------------------------------------------------------
m = cell(1, n_runs);

for i = 1:n_runs

    m{i} = mean(get_wh_image(dat, find(Xi(:, i))));

end

% Set color limits for all runs together
% ----------------------------------------------------------------------
% colormap('default'); % Ensure a default colormap

mm = cat(m{:});
alldat = mm.dat(:);

% clim = [mean(alldat) - 3*std(alldat) mean(alldat) + 3*std(alldat)];
% have to omitnan
clim = [mean(alldat, 'omitnan') - 3*std(alldat, 'omitnan') mean(alldat, 'omitnan') + 3*std(alldat, 'omitnan')];

% indices = reshape(1:n_runs^2, n_runs, [])';
% indices_except_first=indices;
% indices_except_last=indices;
% 
% indices_except_first(:, 1) = [];
% indices_except_last(:, n_runs) = [];

% last_column_indices = n_runs:n_runs:n_runs^2;

for i = 1:n_runs

    subplot(n_runs, 1, i);

    han = display_slices(m{i}, 'axial', 'slices_per_row', 10, 'spacing', 8, 'startslice', -30, 'endslice', 40, 'clim', clim);

    if ~use_generic_labels
        title(sprintf('Run %3.0f', i))
    end

end

% fix axis limits
myax = subplot(n_runs, 1, 1);
mypos = get(myax, 'Position');

for i = 1:n_runs

    myax = subplot(n_runs, 1, i);
    myaxpos = get(myax, 'Position');
    myaxpos(4) = mypos(4);
    set(myax, 'Position', myaxpos);
end

colormap gray
drawnow



dataMatrix = [];
imageNames= {};
for i = 1:numel(m)
    structure = m{i};
    datField = structure.dat;
    dataMatrix = [dataMatrix, datField];
    imgNames{i} = ['Run ', num2str(i)];
end

% subplot(n_runs, n_runs, sort(indices_except_last(:)'));

if n_runs==2
    figure; plot(m{1}.dat, m{2}.dat, 'k.'); refline
    % figure; corr(m{1}.dat, m{2}.dat, 'k.'); refline
    corr(m{1}.dat, m{2}.dat)
    title(sprintf('r = %3.2f',  corr(m{1}.dat, m{2}.dat)))
    xlabel('Run-1')
    ylabel('Run-2')
    drawnow
end

if n_runs > 2
    plot_correlation_matrix(dataMatrix, 'var_names', imgNames, 'circles', false);
    drawnow
end





end % montage_run_means



% ---------------------------------------------------------------
% Figure position (montages)
% ---------------------------------------------------------------

function  set_figure_position(i, j)
figaspect = j / i;  % aspect ratio, width / height

screensz = get(0, 'screensize');
%maxwh = min(screensz(3), screensz(4)) .* .8;

maxw = screensz(3) .* .8; % max dimension, 80% of full size
maxh = screensz(4) .* .8;

if j > i
    w = maxw;
    h = min(maxh, maxw ./ figaspect);
else
    h = maxh;
    w = min(maxw, maxh .* figaspect);
end

pos = [30 screensz(4) - 30 w h]; % set to top corner, with width and height

set(gcf, 'Position', pos);

end


% ---------------------------------------------------------------
% Orthviews
% ---------------------------------------------------------------

function plot_orthviews(fmridat)

% check to be sure:
fmridat.dat(isnan(fmridat.dat)) = 0;

m = mean(fmridat.dat',1)'; %mean values of each voxel
s = std(fmridat.dat',1)'; %std of each voxel
d = m./s;
d(m == 0 | s == 0) = 0;

if size(fmridat.dat,2) > 1 % if there is more than one image, show std and snr too
    vecs_to_reconstruct = [m s d];
else
    vecs_to_reconstruct = [m];% else just show mean image
end

m = fmridat;
m.dat = vecs_to_reconstruct;

if isempty(fmridat.volInfo)
    disp('.volInfo is empty. Skipping orthviews and other brain plots.');
else
    %create_orthviews(vecs_to_reconstruct, fmridat);
    orthviews(m);

    spm_orthviews_name_axis('Mean data', 1);
    if size(fmridat.dat,2) > 1
        spm_orthviews_name_axis('STD of data', 2);
        spm_orthviews_name_axis('Mean / STD', 3);
    end
    set(gcf, 'Name', 'Orthviews_fmri_data_mean_and_std');
end

end


