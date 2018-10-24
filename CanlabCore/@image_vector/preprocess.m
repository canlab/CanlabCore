function [obj, varargout] = preprocess(obj, meth, varargin)
% Preprocesses data in an image_vector (e.g., fmri_data) object; many options for filtering and outlier id
%
% Data is observations (i.e., voxels, subjects) x images, so operating on the columns operates on
% images, and operating on the rows operates on voxels (or variables more
% generally) across images.
%
% :Inputs:
%
%   **meth:** Options
%
%   **resid:**
%        Residualize voxels with respect to covariates
%        Uses obj.covariates, obj.dat.
%        Adds intercept automatically. You can tell it to add the mean response per voxel back in:
%        obj = preprocess(obj, 'resid', [add mean back in flag])
%
%   **hpfilter:**
%        High-pass filter and remove run intercepts and first two
%        images per run. Uses obj.dat, obj.images_per_session
%        obj = preprocess(obj, 'hpfilter', HPlen in s, TR)
%
%   **windsorize:**
%        Windsorize entire data matrix to 3 STD
%
%   **windsorizevoxels:**
%        Windsorize each time series in data matrix to 3 STD
%
%   **session_outliers:**
%        Identify session-wise (run-wise) outliers with significant
%        based on mahalanobis distance with FDR-corrected P-values in chi-square test.
%        Impute session grand mean outliers.
%
%   **outliers:**
%        Identify outlier time points for each session based on
%        mahalanobis distance (see above) across global mean for slices and
%        spatial STD for slices, as in scn_session_spike_id.
%        Outliers at 3 SD based on timeseries added to obj.covariates.
%
%   **outliers_rmssd:**
%        Identify outlier time points for each session based on
%        root-mean-square successive differences between images (across voxels.)
%        this is the std (across voxels) of the successive diffs across images.
%        Outliers at 3 SD based on timeseries added to obj.covariates.
%
%        [obj, rmssd, wh_outliers_rmssd] = preprocess(obj, 'outliers_rmssd');
%
%   **smooth:**
%         Smoothed images with Gaussian filter
%           - obj = preprocess(obj, 'smooth', FWHM in mm)
%         *NOTE* SMOOTHING KERNEL MAY BE IN VOX, AS VOL INFO IS NOT PASSED IN
%
%   **interp_images:**
%        Interpolate all voxels in a series of images specified
%        by logical vector whout.
%          - obj = preprocess(obj, 'interp_images', whout);
%
%   **remove_white_csf:**
%        Extract data values for gray, white, CSF
%        Regress gray matter mean on first 5 components of white-matter and CSF
%        Remove the fitted values from the images image-wise
%        This adjusts for variations in overall image intensity that are explainable by variations in white-matter and CSF
%        By default, uses a highly eroded standard mask for gray/white/CSF,
%        that avoids mixing signal components coming from gray matter.
%        Requires that images be registered in MNI space to work
%        appropriately.
%        This effectively removes a scalar multiple of each white/CSF regressor from
%        all voxels in each image set.
%        Estimating parameters for each voxel independently would add a lot of
%        variability.  this way, we estimate the overall location of the image
%        (shift up/down from zero in gray matter) that is predictable from
%        gray/white variables, and remove that.
%        we can apply this to any data - time series, contrast images, beta
%        images, signature response values.
%
%   **remove_csf:**
%        Same as remove_white_csf, but use only CSF predictors
%        This is because work, e.g., that of John Gore, suggests there is likely
%        real signal in WM, so removing it can cause real signal to be
%        removed
%
%   **rescale_by_csf:**
%       This attempts to model and remove scale inhomogeneity across
%       images.  It estimates the relationship between the spatial median abs.
%       deviation (MAD) for CSF voxels and for GM voxels. The proportion of the GM deviation 
%       for each image fitted by CSF is divided out. 

% :Examples:
% ::
%
%   % two complementary ways to get and plot outliers:
% ---------------------------------------------------------------------
%    dat = preprocess(dat, 'outliers', 'plot');
%    subplot(5, 1, 5); % go to new panel...
%    dat = preprocess(dat, 'outliers_rmssd', 'plot');
%
%    Concatenate a set of image objects and then regress out white/CSF components
% ---------------------------------------------------------------------
%    DATA_CAT = cat(DATA_OBJ{:});
%    for i = 1:size(DATA_OBJ, 2), sz(i) = size(DATA_OBJ{i}.dat, 2); end
%    DATA_CAT.images_per_session = sz;
%    DATA_CAT.removed_images = 0;
%
%    DATA_CAT = preprocess(DATA_CAT, 'remove_csf');
%
%    For a 2nd-level image dataset, Regress out and rescale by CSF values
%    Plot before and after
% ---------------------------------------------------------------------
% obj = load_image_set('emotionreg');
% histogram(obj, 'by_tissue_type', 'byimage');
% obj = preprocess(obj, 'remove_csf');
% obj = preprocess(obj, 'rescale_by_csf');
% history(obj)

switch meth
    
    % ---------------------------------------------------------------------
    case {'resid', 'residuals'}
        % ---------------------------------------------------------------------
        
        add_mean_in = 0;
        if ~isempty(varargin) > 0 && varargin{1}
            add_mean_in = 1;
        end
        
        obj.dat = resid(obj.covariates, obj.dat', add_mean_in)';
        
        obj.history{end + 1} = 'Residualized voxels with respect to covariates.';
        
        % ---------------------------------------------------------------------
    case 'hpfilter'
        % ---------------------------------------------------------------------
        
        if length(varargin) == 0
            error('Enter HP filter cutoff in sec as 3rd argument.');
        end
        
        if length(varargin) == 1
            error('Enter TR in sec as 4th argument.');
        end
        
        hpcutoff = varargin{1};
        TR = varargin{2};
        
        whgood = ~all(obj.dat == 0, 2);
        
        fprintf('removing intercepts, 1st 2 images/run, high-pass filtering %3.0f voxels at %3.0f sec', sum(whgood), hpcutoff);
        
        imageseriesmean = mean(obj.dat(whgood, :), 2);
        
        if isempty(obj.images_per_session)
            obj.images_per_session = size(obj.dat, 2); % n images
        end
        
        obj.dat(whgood, :) = hpfilter(obj.dat(whgood, :)', TR, hpcutoff, obj.images_per_session, [], 1:2)';
        
        obj.dat(whgood, :) = obj.dat(whgood, :) + repmat(imageseriesmean, 1, size(obj.dat, 2));
        
        obj.history{end+1} = 'HP filtered and residualized with respect to session intercepts and first 2 images of each run, mean added back in';
        
        fprintf('\n');
        
        % ---------------------------------------------------------------------
    case {'windsor', 'windsorize'}
        % ---------------------------------------------------------------------
        
        % Windsorize entire data matrix
        %m = mean(obj.dat(:));
        % note: the line above was giving different values for operation on
        % single and double values in Matlab 2012a.  Very disturbing! The
        % correct value is given by double() or the line below.
        m = sum(sum(obj.dat)) ./ numel(obj.dat); % save memory space
        
        s = 3 .* std(obj.dat(:));
        
        whbad = obj.dat < m - s;
        obj.dat(whbad) = m - s;
        whbad2 = obj.dat > m + s;
        obj.dat(whbad2) = m + s;
        
        nbad = sum(whbad(:)) + sum(whbad2(:));
        percbad = 100 .* nbad ./ numel(obj.dat);
        
        obj.history{end+1} = sprintf('Windsorized data matrix to 3 STD; adjusted %3.0f values, %3.1f%% of values', nbad, percbad);
        disp(obj.history{end});
        
        % ---------------------------------------------------------------------
    case 'windsorizevoxels'
        % ---------------------------------------------------------------------
        whbad = all(obj.dat == 0, 2);
        nok = sum(~whbad);
        
        fprintf('windsorizing %3.0f voxels to 3 std: %05d', nok, 0);
        for i = 1:nok
            if mod(i, 100) == 0, fprintf('\b\b\b\b\b%05d', i); end
            obj.dat(i, :) = trimts(obj.dat(i, :)', 3, [])';
        end
        
        obj.history{end+1} = 'Windsorized data voxel-wise to 3 STD';
        
        % ---------------------------------------------------------------------
    case 'session_outliers'
        % ---------------------------------------------------------------------
        if ~isa(obj, 'fmri_data'), error('method only defined for fmri_data objects.'), end
        if isempty(obj.images_per_session), error('method only works if images_per_session field is filled.'); end
        
        nscan = obj.images_per_session;  % num images per session
        I = intercept_model(nscan);
        
        obj.history{end+1} = 'Imputed mean for session outliers in GM/GSTD space: ';
        
        create_figure('session_outliers', 3, size(I, 2));
        
        for i = 1:size(I, 2)
            
            wh = find(I(:, i));
            y = obj.dat(:, wh);
            
            mv = mean(y')';  % mean across time, for each voxel
            m = mean(y)';    % mean across voxels, for each time point
            s = std(y)';     % spatial std
            
            % mahalanobis: strange patterns across slices
            d2 = mahal([m s], [m s]);
            
            % threshold
            p = 1 - chi2cdf(d2, 2);
            p(p == 0) = eps;
            pthr = max(.05 / length(d2), FDR(p, .05));  % max of FDR or Bonferroni
            if isempty(pthr), pthr = .05 / length(d2); end
            whb = p < pthr; % which images are outliers
            
            % impute mean image for this session
            y(:, whb) = repmat(mv, 1, sum(whb));
            obj.dat(:, wh) = y;
            
            % report
            obj.history{end} = [obj.history{end} sprintf(' S%02d:%2.0f imgs', i, sum(whb))];
            
            subplot(3, size(I, 2), i)
            plot(m, s, 'k.')
            plot(m(whb), s(whb), 'ro', 'LineWidth', 2);
            xlabel('Image mean');
            if i == 1, ylabel('Image std'); end
            title(['Session' num2str(i)])
            
            subplot(3, size(I, 2), size(I, 2)+i)
            plot(m, 'k.-')
            plot(find(whb), m(whb), 'ro', 'LineWidth', 2);
            xlabel('Time');
            if i == 1, ylabel('Image mean');  end
            
            m = mean(y)';
            subplot(3, size(I, 2), 2*size(I, 2)+i)
            plot(m, 'k.-')
            plot(find(whb), m(whb), 'ro', 'LineWidth', 2);
            xlabel('Time');
            if i == 1, ylabel('Mean after adjustment');    end
            
            drawnow
            
        end % session
        
        disp(obj.history{end});
        
        % ---------------------------------------------------------------------
    case 'outliers'
        % ---------------------------------------------------------------------
        
        % regress out outliers
        stdev = 3;
        
        doplot = any(strcmp(varargin, 'plot'));
        
        nscan = obj.images_per_session;
        if isempty(nscan), error('Enter obj.images_per_session to use this method.'); end
        
        I = intercept_model(nscan);
        for i = 1:size(I, 2)
            
            clear gslice stdslice; % Wani added to fix a bug
            
            fprintf('Session %3.0f: ', i)
            
            wh = find(I(:, i));
            y = obj.dat(:, wh);
            
            g = mean(y)';
            
            if isempty(obj.volInfo), error('obj.volInfo must be complete!'); end
            
            z = obj.volInfo.xyzlist(:, 3);
            
            if size(z, 1) == length(obj.removed_voxels)
                z(obj.removed_voxels) = [];
            end
            
            u = unique(z);
            
            for t = 1:size(y, 2) % time within session
                
                for j = 1:length(u)
                    gslice(j, t) = mean( obj.dat(z == u(j), t) );
                    stdslice(j, t) = std( obj.dat(z == u(j), t) );
                end
                
            end
            
            wh_no_data = ~any(gslice')' | ~any(stdslice')';
            gslice(wh_no_data,:) = [];
            stdslice(wh_no_data,:) = [];
            
            % mahalanobis: strange patterns across slices
            d2 = mahal([stdslice' gslice'], [stdslice' gslice']);
            
            [gtrim, dummy, gspikes] = trimts(g, stdev, [], 1);
            
            [dummy, dummy, mahalspikes] = trimts(d2, stdev, [], 1);
            spikes = unique([gspikes; mahalspikes]);
            
            spikesperimg = 100*length(spikes) ./ length(g);
            snr = mean(g) ./ std(g);
            
            fprintf('%3.0f Potential outliers\t%%Spikes: %3.2f\tGlobal SNR (Mean/STD): %3.2f\n', length(spikes), spikesperimg , snr)
            
            nuisance_covs{i} = intercept_model(length(gtrim), spikes);
            % get rid of intercept, because we may want to add this
            % separately
            nuisance_covs{i} = nuisance_covs{i}(:, 2:end);
            
            if doplot
                gall{i} = g';
                gspikesall{i} = gspikes';
                gsliceall{i} = scale(gslice', 1)';
                stdsliceall{i} = scale(stdslice', 1)';
                d2all{i} = d2';
            end
            
        end % session
        
        covs =  blkdiag(nuisance_covs{:});
        
        % Plot
        
        if doplot
            create_figure('sessions', 5, 1);
            
            subplot(5, 1, 1);
            imagesc(cat(2, gsliceall{:})); axis tight
            ylabel('slice');
            title('Global mean signal')
            
            subplot(5, 1, 2);
            imagesc(cat(2, stdsliceall{:})); axis tight
            ylabel('slice');
            title('Standard deviation across slice voxels')
            
            subplot(5, 1, 3);
            toplot = cat(2, gall{:});
            draw_session_boxes(obj.images_per_session, toplot)
            plot(toplot); axis tight
            set(gca, 'YLim', [min(toplot) max(toplot)])
            title('Global in-brain signal')
            
            subplot(5, 1, 4);
            toplot = cat(2, d2all{:});
            draw_session_boxes(obj.images_per_session, toplot)
            plot(toplot); axis tight
            set(gca, 'YLim', [min(toplot) max(toplot)])
            ylabel('d^2');
            title('Mahalanobis distance')
            xlabel('Time (images)')
            for i = 1:size(covs, 2)
                lineh = plot_vertical_line(find(covs(:, i)));
                set(lineh, 'Color', 'r');
            end
            
            subplot(5, 1, 5); axis off; % for rmssd plot...
            
        end
        
        obj.covariates = [obj.covariates covs];
        for i = 1:size(covs, 2), obj.covariate_names{end+1} = ['outlier_cov']; end
        
        obj.history{end+1} = sprintf('Added %3.0f global/mahal outlier covariates to covariates field.', size(covs, 2));
        disp(obj.history{end});
        
        % ---------------------------------------------------------------------
    case 'outliers_rmssd'
        % ---------------------------------------------------------------------
        % outliers from root mean square successive differences across
        % images
        
        % remove intercepts first...
        nscan = obj.images_per_session;
        if ~isempty(nscan) %, error('Enter obj.images_per_session to use this method.'); end
            
            I = intercept_model(nscan);
            datadj = obj.dat - (I* (pinv(I) * obj.dat'))';
            
        else
            % skip intercept model
            datadj = obj.dat;
        end
        
        sdlim = 3;
        sdiffs = diff(datadj')'; % successive differences
        sdiffs = [mean(sdiffs, 2) sdiffs]; % keep in image order; shift over one
        mysd = std(sdiffs(:));
        
        rmssd = (mean(sdiffs .^ 2)) .^ .5; %spatsd = std(sdiffs);  % this is not quite the same thing!
        
        % avoid first time point being very different and influencing distribution and plots.
        rmssd(1) = mean(rmssd);
        
        out_rmssd = [rmssd > mean(rmssd) + sdlim*std(rmssd)];
        
        doplot = any(strcmp(varargin, 'plot'));
        if doplot
            %create_figure('outliers_rmssd', 2, 1);
            draw_session_boxes(obj.images_per_session, rmssd)
            plot(rmssd, 'b');
            title('Root mean square successive diffs (rmssd) across images');
            axis tight;
            set(gca, 'YLim', [min(rmssd) max(rmssd)])
            ylabel('RMSSD');
            xlabel('Time (images)');
        end
        
        
        % add covariates
        wh = find(out_rmssd);
        covs = zeros(length(out_rmssd), length(wh));
        for i = 1:length(wh)
            covs(wh(i), i) = 1;
            
            if doplot, h = plot_vertical_line(wh(i)); set(h, 'Color', 'r'); end
        end
        
        obj.covariates = [obj.covariates covs];
        for i = 1:size(covs, 2), obj.covariate_names{end+1} = ['outlier_rmssd_cov']; end
        
        obj.history{end} = [obj.history{end} sprintf('\nOutliers in RMSSD images: %3.0f%%, %2.0f imgs.\n', sum(out_rmssd)./length(out_rmssd), sum(out_rmssd))];
        disp(obj.history{end});
        
        varargout{1} = rmssd';
        varargout{2} = out_rmssd';
        
        % ---------------------------------------------------------------------
    case 'smooth'
        % ---------------------------------------------------------------------
        
        if length(varargin) == 0
            error('Enter smoothing FWHM in mm as 3rd argument.');
        end
        
        obj = replace_empty(obj);
        obj.dat = iimg_smooth_3d(obj.dat, obj.volInfo, varargin{1});
        
        obj.history{end+1} = sprintf('Smoothed images with %3.0f FWHM filter', varargin{1});
        disp(obj.history{end});
        
        % ---------------------------------------------------------------------
    case 'interp_images'
        % ---------------------------------------------------------------------
        
        if length(varargin) == 0
            error('whout argument missing. Must be logical vector of 1s and 0s for images to interpolate.');
        else
            whout = varargin{1};
        end
        
        if any(whout ~= 1 & whout ~= 0)
            error('whout must be logical vector of 1s and 0s for images to interpolate.');
        end
        whout = logical(whout);
        wh = find(whout);
        whin = ~whout;                      % do just once
        
        
        fprintf(1, 'Interpolating outlier images. Setup: ')
        tic
        obj = remove_empty(obj);
        
        t = 1:size(obj.dat, 2);
        t(wh) = [];
        v = size(obj.dat, 1);
        
        fprintf('%3.0f sec. Running: ', toc);
        
        tic
        for i = 1:v
            Vq = interp1(t, obj.dat(i, whin), wh, 'linear','extrap');
            % figure; plot(t, obj.dat(i, whin)); hold on; plot(wh, Vq, 'ro');
            obj.dat(i, wh) = Vq;
        end
        
        fprintf('Done in %3.0f sec.\n', toc);
        
        obj.history{end + 1} = sprintf('Interpolated %3.0f images with 1-D linear interp.', length(wh));
        
        % ---------------------------------------------------------------------
    case 'remove_white_csf'
        % ---------------------------------------------------------------------
        
        doplot = any(strcmp(varargin, 'plot'));
        
        [means, components] = extract_gray_white_csf(obj);
        
        x = cat(2, components{2:3});        % predictors: white and csf only
        
        wh_gray_white = 1:size(x, 2);
        
        meangray = means(:, 1);             % mean gray matter for each image
        
        if isa(obj, 'fmri_data') && ~isempty(obj.images_per_session)
            
            xi = intercept_model(obj.images_per_session);
            x = [x xi];                       % model
            gray_tmp = resid(xi, meangray);   % for partial correlation
            
        else
            gray_tmp = meangray;
        end
        
        % regress out signal explainable by gray/white (fits)
        % estimate regression coefficients for
        % gray matter average predicted by white-matter and CSF covariates
        
        b = pinv(x) * meangray;
        
        fit = x(:, wh_gray_white) * b(wh_gray_white);
        
        % get correlation value
        
        r = corr(gray_tmp, fit);
        
        if doplot
            create_figure('gray predicted from white/CSF components')
            plot_correlation_samefig(fit, meangray);
        end
        
        sz = size(obj.dat);
        to_subtract = repmat(fit', sz(1), 1);
        
        obj.dat = obj.dat - to_subtract;
        
        obj.history{end + 1} = sprintf('Regressed out white/CSF components image-wise. \nCorrelation between predicted and actual mean gray matter before removal: r = %3.4f', r);
        
        % ---------------------------------------------------------------------
    case 'remove_csf'
        % ---------------------------------------------------------------------
        
        obj = remove_empty(obj);
        
        doplot = any(strcmp(varargin, 'plot'));
        
        [means, components] = extract_gray_white_csf(obj);
        
        x = cat(2, components{3});
        % predictors: csf only. This is because work, e.g., that of John Gore, suggests there may be real signal in WM
        
        wh_gray_white = 1:size(x, 2);       % here, CSF, not really gray.white
        
        meangray = means(:, 1);             % mean gray matter for each image
        
        if isa(obj, 'fmri_data') && ~isempty(obj.images_per_session)
            
            % create orthogonal contrast set with session/set intercepts
            % doesn't actually need to be orthogonalized if we are not
            % removing the session intercepts

            xi = intercept_model(obj.images_per_session);

            c = create_orthogonal_contrast_set(size(xi, 2));

            xi = xi * c';                     % orthogonal set of nsessions - 1
            xi(:, end) = 1;                   % now we can add the intercept

            x = [x xi];                       % model
            gray_tmp = resid(xi, meangray);   % for partial correlation
            
        else
            gray_tmp = meangray;
        end
        
        % regress out signal explainable by gray/white (fits)
        % estimate regression coefficients for
        % gray matter average predicted by white-matter and CSF covariates
        
        b = pinv(x) * meangray;
        
        % partial fit of nuisance
        fit = x(:, wh_gray_white) * b(wh_gray_white);
        
        % get correlation value
        
        r = corr(gray_tmp, fit);
        
        if doplot
            create_figure('gray predicted from CSF components')
            plot_correlation_samefig(fit, meangray);
            xlabel('Predicted mean from CSF');
            ylabel('Observed gray-matter mean');
        end
        
        % Subtract fitted value for overall GM as predicted from CSF for
        % each image.
        % Then add back in the intercept (overall mean)
        sz = size(obj.dat);
        to_subtract = repmat(fit', sz(1), 1);
        
%         obj.dat = remove_empty(obj.dat);  % components in reduced space
        obj.dat = obj.dat - to_subtract;
        
        obj.history{end + 1} = sprintf('Regressed out CSF components image-wise. \nCorrelation between predicted and actual mean gray matter before removal: r = %3.4f', r);
        
        % ---------------------------------------------------------------------
    case 'rescale_by_csf'
        % ---------------------------------------------------------------------
        
        obj = remove_empty(obj);
        
        doplot = any(strcmp(varargin, 'plot'));
        
        [means, components, gwcsfobj] = extract_gray_white_csf(obj);
        
        madgray = mad(gwcsfobj{1}.dat)';    % median abs deviation gray matter for each image
        madcsf = mad(gwcsfobj{3}.dat)';    % median abs deviation CSF for each image
        
        x = madcsf;
        % predictors: csf only. This is because work, e.g., that of John Gore, suggests there may be real signal in WM
        
        wh_gray_white = 1:size(x, 2);
        
        % for correlation plot
        if isa(obj, 'fmri_data') && ~isempty(obj.images_per_session)
            
            xi = intercept_model(obj.images_per_session);
            x = [x xi];                       % model
            gray_tmp = resid(xi, madgray);   % for partial correlation
            
        else
            gray_tmp = madgray;
        end
        
        % regress out signal explainable by gray/white (fits)
        % estimate regression coefficients for
        % gray matter average predicted by white-matter and CSF covariates
        
        b = pinv(x) * madgray;
        
        fit = x(:, wh_gray_white) * b(wh_gray_white);
        
        % get correlation value
        
        r = corr(gray_tmp, fit);
        
        if doplot
            create_figure('gray predicted from CSF components')
            plot_correlation_samefig(fit, meangray);
            xlabel('Predicted MAD from CSF');
            ylabel('Observed gray-matter MAD');
        end
        
        % Scale by fitted value for overall GM as predicted from CSF for
        % each image.
        sz = size(obj.dat);
        to_divide_by = repmat(fit', sz(1), 1);
        
        obj.dat = obj.dat ./ to_divide_by;
                
        obj.history{end + 1} = sprintf('Scaled images by fitted relationship between MAD GM and MAD CSF. \nCorrelation between CSF and GM scale before removal: r = %3.4f', r);
        
       
    otherwise
        error('Unknown preprocessing method.')
        
        
end


end % function





function draw_session_boxes(images_per_session, y)


nsess = length(images_per_session);
colors = {'k' 'k' 'k' 'k' 'k' 'k' 'k' 'k'};  %{'ro' 'go' 'bo' 'yo' 'co' 'mo' 'ko' 'r^' 'g^' 'b^' 'y^' 'c^' 'm^' 'k^'};
while nsess > length(colors), colors = [colors colors]; end

yval = min(y(:)) - .05 * min(y(:));
yheight = (max(y(:)) + .05 * max(y(:))) - yval;

c = cumsum(images_per_session);
st = [0 c(1:end-1)];

for jj = 1:2:nsess
    h1 = drawbox(st(jj), images_per_session(jj), yval, yheight, colors{jj}(1));
    set(h1, 'FaceAlpha', .10, 'EdgeColor', 'none');
end

end



function [r,X] = resid(X,y, varargin)
% [r,X] = resid(X,y, [add mean back in flag])
%
% tor wager
% residuals from model fit
%
% adds intercept, if missing; adds mean response for each variable if
% requested
% Last edited: 6/2013

add_int = 0;
if ~isempty(varargin) > 0 && varargin{1}
    add_int = 1;
end

% find intercept
wint = all(X == repmat(X(1,:), size(X,1), 1));
if ~any(wint)
    X(:,end+1) = 1;
    wint = zeros(1, size(X, 2));
    wint(end) = 1;
end

y = double(y);

% break up for efficiency - otherwise Matlab chokes with large
% matrices...(why?)
px = pinv(X);
pxy = px * y;
xpxy = X * pxy;


r = y - xpxy; % X * pinv(X) * y;

if add_int
    m = mean(y);
    ym = repmat(m, size(X, 1), 1);
    r = r + ym;
end

end




