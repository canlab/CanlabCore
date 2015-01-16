function fmridat = rescale(fmridat, meth, varargin)
% fmridat = rescale(fmridat, meth)
%
% Rescales data in an fmri_data object
% Data is observations x images, so operating on the columns operates on
% images, and operating on the rows operates on voxels (or variables more
% generally) across images.
%
% Methods:
% 'centervoxels'
% 'zscorevoxels'
% 'centerimages'
% 'zscoreimages'
% 'rankvoxels'
%
% 'windsorizevoxels'
% 'percentchange'
% 'tanh'
%
% Appropriate for multi-session (time series) only:
% 'session_global_percent_change'
% 'session_global_z'
% 'session_multiplicative'
%
% see also fmri_data.preprocess

switch meth
    
    case 'centervoxels'
        
        fmridat.dat = scale(fmridat.dat', 1)';
        
        fmridat.history{end+1} = 'Centered voxels (rows) across images';
        
    case 'zscorevoxels'
        
        fmridat.dat = scale(fmridat.dat')';
        
        fmridat.history{end+1} = 'Z-scored voxels (rows) across images';
        
    case 'rankvoxels'
        for i = 1:size(fmridat.dat, 1) % for each voxel
            
            d = fmridat.dat(i, :)';
            if ~all(d == 0)
                fmridat.dat(i, :) = rankdata(d)';
            end
            
        end
        
    case 'centerimages'
        
        % center images (observations)
        fmridat.dat = scale(fmridat.dat, 1);
        
        fmridat.history{end+1} = 'Centered images (columns) across voxels';
        
    case 'zscoreimages'
        
        fmridat.dat = scale(fmridat.dat);
        
        fmridat.history{end+1} = 'Z-scored imagesc(columns) across voxels';
        
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
        
    case 'session_spm_style'
        % SPM's default method of global mean scaling
        % useful for replicating SPM analyses or comparing SPM's scaling to
        % other methods.
        % not implemented yet because tor decided to use spm_global on images for comparison; this could be done though...
        % see help spm_global
        %         nscan = fmridat.images_per_session;  % num images per session
        %         I = intercept_model(nscan);
        %         for i = 1:size(I, 2)
        %
        %             wh = find(I(:, i));
        %             y = fmridat.dat(:, wh)';   % y is images x voxels
        %             gm = mean(y); % mean at each voxel, 1 x voxels
        %
        %             % subtract mean at each vox, divide by global session mean
        %             y = (y - repmat(gm, size(y, 1), 1)) ./ std(y(:));
        %             fmridat.dat(:, wh) = y;
        %         end
        
        
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
        
    case 'percentchange'
        % scale each voxel (column) to percent signal change
        % with a mean of 100
        % based on smoothed mean values across space (cols)
        
        m = mean(fmridat.dat')'; % mean at each voxel, voxels x 1
        
        sfwhm = 16;
        ms = iimg_smooth_3d(m, fmridat.volInfo, sfwhm, fmridat.removed_voxels);
        
        % subtract mean at each vox, divide by global session mean
        fmridat.dat = 100 + 100 .* (fmridat.dat - repmat(m, 1, size(fmridat.dat, 2))) ./ repmat(ms, 1, size(fmridat.dat, 2));
        
        fmridat.history{end+1} = 'Rescaled to mean 100, voxelwise % signal change with 16 mm fwhm smoothing of divisor mean';
        
    case 'correly'
        
        % correl with y
        % provides implicit feature selection when using algorithms that
        % are scale-dependent (e.g., SVM, PCA)
        
    otherwise
        error('Unknown scaling method.')
        
end

end