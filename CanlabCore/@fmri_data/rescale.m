function fmridat = rescale(fmridat, meth, varargin)
% Rescales data in an fmri_data object
% Data is observations x images, so operating on the columns operates on
% images, and operating on the rows operates on voxels (or variables more
% generally) across images.
%
% :Usage:
% ::
%
%    fmridat = rescale(fmridat, meth)
%
% :Inputs:
%
%   **Methods:**
%     Rescaling of voxels
%     - centervoxels        subtract voxel means
%     - zscorevoxels        subtract voxel means, divide by voxel std. dev
%     - rankvoxels          replace values in each voxel with rank across images
%
%                           *Note: these methods must exclude invalid (0 or
%                           NaN) voxels image-wise. Some images (but not others) in an
%                           object may be missing some voxels.
%
%     Rescaling of images
%     - centerimages        subtract image means
%     - zscoreimages        subtract image means, divide by image std. dev
%
%                           *Note: these methods must exclude invalid (0 or
%                           NaN) voxels image-wise. Some images (but not others) in an
%                           object may be missing some voxels.
%     - prctileimages       - Convert image values to percentile scores. This
%                           normalizes the range and distribution of all
%                           images. It is like ranking but MUCH faster for
%                           images with many elements.
%                           - Normalize each image by number of valid values, so images in set
%                           will be on the same scale if some are missing voxels
%
%     - l2norm_images       divide each image by its l2 norm, multiply by sqrt(n valid voxels)
%     - divide_by_csf_l2norm  divide each image by CSF l2 norm. requires MNI space images for ok results!
%     - rankimages          rank voxels within image;
%     - csf_mean_var        Subtract mean CSF signal from image and divide by CSF variance. (Useful for PE
%                               maps where CSF mean and var should be 0).
%
%     Other procedures
%     - windsorizevoxels    Winsorize extreme values to 3 SD for each image (separately) using trimts
%     - percentchange       Scale each voxel (column) to percent signal change with a mean of 100
%                           based on smoothed mean values across space (cols), using iimg_smooth_3d
%     - tanh                Rescale variables by hyperbolic tangent function
%                           this will shrink the tails (and outliers) towards the mean,
%                           making subsequent algorithms more robust to
%                           outliers. Affects normality of distribution.
%     - doublecenter        Double-center voxels and images, so the mean of each voxel and each image is 0
%
% Whole run-level scaling:
%     - 'session_grand_mean_scaling_spm_style'
%
% Appropriate for multi-session (time series):
%     - session_global_percent_change
%     - session_global_z
%     - session_multiplicative
%     - percentchange
%
% See also fmri_data.preprocess

switch meth
    
    case 'centervoxels'
        
        % Consider missing voxels, and exclude case-wise (image-wise)
        ismissing = fmridat.dat == 0 | isnan(fmridat.dat);
        fmridat.dat(ismissing) = NaN;
        
        m = nanmean(fmridat.dat');
        fmridat.dat = (fmridat.dat' - m)';  % Note: Sizes are wrong without repmat, but Matlab 2020b at least figures this out.
        
        % Old - this does not handle different missing voxels in different images
        fmridat.dat = scale(fmridat.dat', 1)';
        
        fmridat.dat(ismissing) = 0;  % Replace with 0 for compatibility with image format
        
        fmridat.history{end+1} = 'Centered voxels (rows) across images';
        
    case 'zscorevoxels'
        
        % Consider missing voxels, and exclude case-wise (image-wise)
        ismissing = fmridat.dat == 0 | isnan(fmridat.dat);
        fmridat.dat(ismissing) = NaN;
        
        m = nanmean(fmridat.dat');
        s = nanstd(fmridat.dat');
        fmridat.dat = ((fmridat.dat' - m) ./ s)';  % Note: Sizes are wrong without repmat, but Matlab 2020b at least figures this out.
        
        fmridat.dat(ismissing) = 0;  % Replace with 0 for compatibility with image format
        
        % Old - this does not handle different missing voxels in different images
        %         fmridat.dat = scale(fmridat.dat')';
        
        fmridat.history{end+1} = 'Z-scored voxels (rows) across images';
        
    case 'rankvoxels'
        for i = 1:size(fmridat.dat, 1) % for each voxel
            
            d = fmridat.dat(i, :)';
            
            % Consider missing voxels, and exclude case-wise (image-wise)
            ismissing = d == 0 | isnan(d);
            d(ismissing) = NaN;
            
            if ~all(d == 0)
                fmridat.dat(i, ~ismissing) = rankdata(d(~ismissing))';
            end
            
        end
        
        ismissing = fmridat.dat == 0 | isnan(fmridat.dat);
        fmridat.dat(ismissing) = 0; % Replace with 0 for compatibility with image format
        
        fmridat.history{end+1} = 'Ranked voxels (rows) across images';
        
    case 'rankimages'
        dat = zeros(size(fmridat.dat));
        
        parfor i = 1:size(fmridat.dat, 2)
            
            d = fmridat.dat(:,i);
            
            % Consider missing voxels, and exclude case-wise (image-wise)
            ismissing = d == 0 | isnan(d);
            d(ismissing) = NaN;
            
            if ~all(d == 0)
                d(~ismissing, i) = rankdata(d(~ismissing));
            end
            
            dat(:, i) = d;
            
        end
        
        dat(isnan(dat)) = 0; % Replace with 0 for compatibility with image format
        
        fmridat.dat = dat;
        
        fmridat.history{end+1} = 'Ranked images (columns) across voxels';
        

    case 'prctileimages'

        res = 0.1;  
        
        newobj = fmridat;

        for i = 1:size(fmridat.dat, 2)

            d = fmridat.dat(:,i);

            % Consider missing voxels, and exclude case-wise (image-wise)
            ismissing = d == 0 | isnan(d);
            d(ismissing) = NaN;

            if ~all(d == 0 | isnan(d))

                x = [-Inf prctile(d, 0:res:100)];

                for j = 1:length(x) - 1

                    wh = d <= x(j + 1) & d > x(j);
                    newobj.dat(wh, i) = j;

                end % bins

                % normalize by number of valid values, so images in set
                % will be on the same scale if some are missing voxels
                d = d ./ length(d(~ismissing));

            end % if image has values

        end % images

        newobj.dat(isnan(newobj.dat)) = 0; % Replace with 0 for compatibility with image format

        fmridat = newobj;

        fmridat.history{end+1} = 'Converted each image to nearest of 1000 bins';


    case 'centerimages'
        
        % center images (observations)
        
        % Consider missing voxels, and exclude case-wise (image-wise)
        ismissing = fmridat.dat == 0 | isnan(fmridat.dat);
        fmridat.dat(ismissing) = NaN;
        
        m = nanmean(fmridat.dat);
        fmridat.dat = (fmridat.dat - m);  % Note: Sizes are wrong without repmat, but Matlab 2020b at least figures this out.
        
        fmridat.dat(ismissing) = 0;  % Replace with 0 for compatibility with image format
        
        % Old - this does not handle different missing voxels in different images
        % fmridat.dat = scale(fmridat.dat, 1);
        
        fmridat.history{end+1} = 'Centered images (columns) across voxels';
        
    case 'zscoreimages'
        
        % Consider missing voxels, and exclude case-wise (image-wise)
        ismissing = fmridat.dat == 0 | isnan(fmridat.dat);
        fmridat.dat(ismissing) = NaN;
        
        m = nanmean(fmridat.dat);
        s = nanstd(fmridat.dat);
        fmridat.dat = (fmridat.dat - m) ./ s;  % Note: Sizes are wrong without repmat, but Matlab 2020b at least figures this out.
        
        fmridat.dat(ismissing) = 0;  % Replace with 0 for compatibility with image format
        
        % Old - this does not handle different missing voxels in different images
        %         fmridat.dat = scale(fmridat.dat);
        
        fmridat.history{end+1} = 'Z-scored each image (columns) across voxels, excluding missing values (0 or NaN) for the image';
        
        
    case 'doublecenter'
        
        % Consider missing voxels, and exclude case-wise (image-wise)
        ismissing = fmridat.dat == 0 | isnan(fmridat.dat);
        fmridat.dat(ismissing) = NaN;
        
        imagemeans = nanmean(fmridat.dat);
        [v, n] = size(fmridat.dat);
        imagemeanmatrix = repmat(imagemeans, v, 1);
        dat_doublecent = fmridat.dat - imagemeanmatrix;
        
        voxelmeans = nanmean(dat_doublecent, 2);
        voxelmeanmatrix = repmat(voxelmeans, 1, n);
        dat_doublecent = dat_doublecent - voxelmeanmatrix;
        fmridat.dat = dat_doublecent;
        
        fmridat.dat(ismissing) = 0;  % Replace with 0 for compatibility with image format
        
        fmridat.history{end+1} = 'Double-centered data matrix across images and voxels';
        
    case 'l1norm_images'
        
        % Vector L1 norm of vector
        % divide by this value to normalize image
        
%         wh = x ~= 0 & ~isnan(x);  % valid values
%         nvalid = sum(wh);
%         l1norm = double(sum(abs(x(wh))))
%         normalized_x = x .* nvalid ./ l1norm;  % divides each voxel x in image by the mean valid voxel value.


        normfun = @(x) sum(abs(x));
        
        x = fmridat.dat;
        
        % Consider missing voxels, and exclude case-wise (image-wise)
        ismissing = fmridat.dat == 0 | isnan(fmridat.dat);
        x(ismissing) = NaN;
        
        for i = 1:size(x, 2)
            % remove nans, 0s
            xx = x(~ismissing(:, i), i);
            % divides each voxel xx in image by the mean valid voxel value.
            n(i) = normfun(xx);
            xx = xx.* length(xx) ./ n(i);
            x(~ismissing(:, i), i) = xx;
        end
        
        x(ismissing) = 0;  % Replace with 0 for compatibility with image format
                
        fmridat.dat = x;
    
    case 'l2norm_images'
        
        % Vector L2 norm / sqrt(length) of vector
        % divide by this value to normalize image
        
        normfun = @(x) sum(x .^ 2) .^ .5;
        
        x = fmridat.dat;
        
        % Consider missing voxels, and exclude case-wise (image-wise)
        ismissing = fmridat.dat == 0 | isnan(fmridat.dat);
        x(ismissing) = NaN;
        
        for i = 1:size(x, 2)
            
            % remove nans, 0s
            xx = x(~ismissing(:, i), i);
            %             isbad = xx == 0 | isnan(xx);
            %             xx(isbad) = [];
            
            % divide by sqrt(length) so number of elements will not change scaling
            n(i) = normfun(xx) ./ sqrt(length(xx));
            
            xx = xx ./ n(i);
            
            %             x(:, i) = zeroinsert(isbad, xx);
            
            x(~ismissing(:, i), i) = xx;
        end
        
        x(ismissing) = 0;  % Replace with 0 for compatibility with image format
                
        fmridat.dat = x;
        
        
    case 'divide_by_csf_l2norm'
        
        [~, ~, ~, l2norms] = extract_gray_white_csf(fmridat);
        
        % divide each column image by its respective ventricle l2norm
        fmridat.dat = bsxfun(@rdivide, fmridat.dat, l2norms(:, 3)') ;
        
        fmridat.history{end+1} = 'Divided each image by its L2 norm';
        
    case 'session_global_percent_change'
        
        nscan = fmridat.images_per_session;  % num images per session
        I = intercept_model(nscan);
        for i = 1:size(I, 2)
            
            wh = find(I(:, i));
            y = fmridat.dat(:, wh)';   % y is images x voxels
            gm = mean(y); % mean at each voxel, 1 x voxels
            
            % subtract mean at each vox, divide by global session mean
            y = (y - repmat(gm, size(y, 1), 1)) ./ std(y(:));
            fmridat.dat(:, wh) = y;
        end
        
        fmridat.history{end+1} = 'Rescaled voxels to percent change relative to the grand mean of each run';
        
    case 'session_grand_mean_scaling_spm_style'
        % SPM's default method of global mean scaling
        % useful for replicating SPM analyses or comparing SPM's scaling to
        % other methods.
        % not implemented yet because tor decided to use spm_global on images for comparison; this could be done though...
        % see help spm_global
        nscan = fmridat.images_per_session;  % num images per session

        if isempty(nscan)
            nscan = size(fmridat.dat, 2);  % single run with all images
        end

        I = intercept_model(nscan);
        for i = 1:size(I, 2)
            
            wh = find(I(:, i));
            y = fmridat.dat(:, wh);   % y is vox  x images
            
            y_vec = y(:);
            y_vec(y_vec == 0 | isnan(y_vec)) = [];  % exclude invalid values from scaling
            
            gm = mean(y_vec); % grand mean for session
            
            % subtract mean at each vox, divide by global session mean
            y = 100 .* y ./ gm;
            fmridat.dat(:, wh) = y;
        end
        
        fmridat.history{end+1} = 'Scaled each run to grand mean of 100';
        
    case 'session_global_z'
        
        % scale each session so that global brain mean and global brain std
        % across time are the same for each session
        
        % underlying model:
        % session-specific shifts in mean signal and scaling exist and are
        % independent
        
        % minus global (whole-brain) mean / global (whole-brain) std.
        % across time
        
        nscan = fmridat.images_per_session;  % num images per session
        I = intercept_model(nscan);
        for i = 1:size(I, 2)
            
            wh = find(I(:, i));
            y = fmridat.dat(:, wh);
            gm = mean(y); % mean at each time point
            
            % subtract mean at each vox, divide by session global std
            y = (y - repmat(gm, size(y, 1), 1)) ./ std(y(:));
            fmridat.dat(:, wh) = y;
        end
        
        fmridat.history{end+1} = 'Z-scored voxels using pooled std across the brain';

    case 'session_multiplicative'
        
        % scale - multiplicative
        % underlying model:
        % a is a process constant across time, but different for each voxel
        % (T2 contrast as a function of voxel properties)
        % b is a process constant across voxels, but different at each time
        % point within a session
        % (overall scaling, which varies across time)
        % these interact multiplicatively to create variation in both
        % global mean and std deviation jointly, which is why global mean
        % and global std are intercorrelated correlated.
        %
        % image std scales with b, and so does image mean.
        % model assumes a linear relationship between mean and std with
        % slope = 1
        
        gm = mean(fmridat.dat, 1); % mean across brain, obs series
        gs = std(fmridat.dat, 1);
        
        create_figure('global mean vs std', 1, 2);
        plot(gm, gs, 'k.')
        xlabel('Image mean');
        ylabel('Image std');
        title('Mean vs. std, before scaling, each dot is 1 image')
        drawnow
        
        %         if length(varargin) == 0 || isempty(varargin{1})
        %             error('Must enter number of images in each session as input argument');
        %         end
        
        if isempty(fmridat.images_per_session)
            fmridat.images_per_session = size(fmridat.dat, 2);
        end
        
        nscan = fmridat.images_per_session;  % num images per session
        I = intercept_model(nscan);
        
        for i = 1:size(I, 2)
            
            wh = find(I(:, i));
            y = fmridat.dat(:, wh);
            y(y < 0) = 0;  % images should be all positive-valued
            
            a = mean(y')';  % mean across time, for each voxel
            b = mean(y);    % mean across voxels, for each time point
            
            % ystar is reconstruction based on marginal means
            ystar = a*b ./ (mean(a));
            
            fmridat.dat(:, wh) = y ./ (ystar + .05*mean(ystar(:)));  % like m-estimator; avoid dividing by zero by adding a constant
            
            % fmridat.dat(:, wh) = y ./ repmat(b, size(y, 1), 1);  % intensity normalization
        end
        
        gm = mean(fmridat.dat, 1); % mean across brain, obs series
        gs = std(fmridat.dat, 1);
        
        subplot(1, 2, 2)
        plot(gm, gs, 'k.')
        xlabel('Image mean');
        ylabel('Image std');
        title('After scaling')
        drawnow
        
        fmridat.history{end+1} = 'Intensity-normalized each run with multiplicative scaling';

    case 'windsorizevoxels'
        
        whbad = all(fmridat.dat == 0, 2);
        nok = sum(~whbad);
        
        fprintf('windsorizing %3.0f voxels to 3 std: %05d', 0);
        for i = 1:nok
            if mod(i, 100) == 0, fprintf('\b\b\b\b\b%05d', i); end
            fmridat.dat(i, :) = trimts(fmridat.dat(i, :)', 3, [])';
        end
        
        fmridat.history{end+1} = 'Windsorized each voxel data series to 3 sd';
        
    case 'tanh'
        % rescale variables by hyperbolic tangent function
        % this will shrink the tails (and outliers) towards the mean,
        % making subsequent algorithms more robust to outliers.
        % However, it also truncates and flattens the distribution
        % (non-normal)
        
        fmridat.dat = tanh(zscore(fmridat.dat'))';
        
        fmridat.history{end+1} = 'Transformed by tanh(zscore(values))';

    case 'percentchange'
        % scale each voxel (column) to percent signal change with a mean of 100
        % based on smoothed mean values across space (cols), using iimg_smooth_3d
        
        m = mean(fmridat.dat',1)'; % mean at each voxel, voxels x 1
        
        sfwhm = 16;
        ms = iimg_smooth_3d(m, fmridat.volInfo, sfwhm, fmridat.removed_voxels);
        
        % subtract mean at each vox, divide by global session mean
        fmridat.dat = 100 + 100 .* (fmridat.dat - repmat(m, 1, size(fmridat.dat, 2))) ./ repmat(ms, 1, size(fmridat.dat, 2));
        
        fmridat.history{end+1} = 'Rescaled to mean 100, voxelwise % signal change with 16 mm fwhm smoothing of divisor mean';
        
    case 'correly'
        
        % correl with y
        % provides implicit feature selection when using algorithms that
        % are scale-dependent (e.g., SVM, PCA)
        
    case 'csf_mean_var'
        
        [~,~,tissues] = extract_gray_white_csf(fmridat);
        csfStd = nanstd(tissues{3}.dat);
        csfMean = nanmean(tissues{3}.dat);
        
        fmridat = fmridat.remove_empty;
        dat = fmridat.dat;
        
        dat = bsxfun(@minus,dat,csfMean);
        dat = bsxfun(@rdivide,dat,csfStd);
        
        fmridat.dat = dat;
        fmridat = fmridat.replace_empty;
        
    otherwise
        error('Unknown scaling method.')
        
end

end
