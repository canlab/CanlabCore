function plot(fmridat, plotmethod)
% plot(fmridat, [plotmethod])
%
% Plot methods:
% ----------------------------------------
% Plot data matrix
% plot(fmri_data_object)
%
% Plot means by condition
% plot(fmri_data_object, 'means_for_unique_Y')
%
%
%
% Output:
% ---------------------------------------
% 5 plots and an SPM orthviews presentation of the data.  In the below and elsewhere, "image" connotes a 3D brain volume
% captured every TR.
%
%   subplot 1:  the fMRI data itself.  Color is intensity of signal
%   subplot 2:  presented as a histogram of values for every voxel collected.
%   The low values are typically out-of-brain voxels, as there is no signal
%   there
%   subplot 3:  each point is an image.  The point's X value is the mean
%   intensity of every voxel in that image, and the Y value is the stdev of
%   intensities for all voxels in that image
%   subplot 4:  covariance between images
%   subplot 5:  each point is an image (case = image).  X value is image
%   number in the run, Y is image mean intensity, and the size of the
%   circular marker represents stdev for that image
%
%   Orthviews: mean and STD for a given voxel averaged over time.  Note that the values for mean and STD here are higher than in
%   the plots above.  That is because mean and STD are calculated here by
%   voxel, but in the plots above they are calculated by image.  Images also include out-of-brain areas.


if nargin < 2
    plotmethod = 'data';
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
        imagesc(fmridat.dat');
        colorbar

        axis tight; set(gca, 'YDir', 'Reverse')
        title('fmri data .dat Data matrix');
        xlabel('Voxels'); ylabel('Images');
        drawnow
        
        tmp = mean(fmridat.dat(:));
        stmp = nanstd(fmridat.dat(:));
        myrange = [tmp - 3*stmp tmp + 3*stmp];
        set(gca, 'CLim', myrange);
        drawnow

        
        if ~isempty(fmridat.Y)
            p = get(gca, 'Position'); ystart = p(2); ylen = p(4);
            
            axh = axes('Position', [.05 ystart .03 ylen]);
            imagesc(fmridat.Y);
            title('Y');
            axis tight;
        end
        drawnow

        % ---------------------------------------------------------------
        % Covariance
        % ---------------------------------------------------------------
        
        subplot(2, 3, 4);
        covmtx = cov(fmridat.dat);
        imagesc(covmtx);
        axis tight; set(gca, 'YDir', 'Reverse')
        title('cov(images)');
        colorbar
        drawnow
        
        if ~isempty(fmridat.Y)
            p = get(gca, 'Position'); ystart = p(2); ylen = p(4);
            
            axh = axes('Position', [.05 ystart .03 ylen]);
            imagesc(fmridat.Y);
            title('Y');
            axis tight;
            
        end
        drawnow
         
        % ---------------------------------------------------------------
        % Histogram
        % ---------------------------------------------------------------
        dattmp = fmridat.dat(:);
        subplot(2, 3, 2);
        [h, x] = hist(dattmp, 100);
        han = bar(x, h);
        set(han, 'FaceColor', [.3 .3 .3], 'EdgeColor', 'none');
        axis tight;
        xlabel('Values'); ylabel('Frequency');
        title('Histogram of values');
        drawnow

        clear dattmp
        
        globalmean = nanmean(fmridat.dat);  % global mean of each obs
        globalstd = nanstd(fmridat.dat);  % global mean of each obs
        nobs = length(globalmean);
        sz = rescale_range(globalstd, [1 6]); % marker size related to global std
        sz(sz == 0) = 1;
        %sz(sz < .5) = .5;

        % ---------------------------------------------------------------
        % Global mean vs. std
        % ---------------------------------------------------------------
        subplot(2, 3, 3); hold on;
        plot(globalmean, globalstd, 'k.');
        title('Image mean vs. std across voxels'); 
        xlabel('Image mean'); 
        ylabel('Image std');

        % ---------------------------------------------------------------
        % Global mean vs. time
        % ---------------------------------------------------------------
        if size(fmridat.dat,2) > 1
            subplot(2, 3, 5);  hold on;

            plot(globalmean, '-');
            axis tight

            Y = 1:nobs;
            Yname = 'Case number';

            for i = 1:nobs
                plot(Y(i), globalmean(i), 'ko', 'MarkerSize', sz(i), 'LineWidth', 1);
            end
            ylabel('Global mean');
            xlabel(Yname);
            title('Globals for each case (size = spatial std)')
            axis tight
            drawnow

            if ~isempty(fmridat.Y) && size(fmridat.Y, 1) == nobs

                subplot(2, 3, 6);
                Y = fmridat.Y;
                Yname = 'Y values in fmri data obj';
                for i = 1:nobs
                    plot(Y(i), globalmean(i), 'ko', 'MarkerSize', sz(i), 'LineWidth', 1);
                end
                ylabel('Global mean');
                xlabel(Yname);
                title('Globals for each case (size = spatial std)')
                axis tight
                drawnow
            end
        end
        
        
        % [coeff, score, latent] = princomp(fmridat.dat, 'econ');
        % %d2 = mahal(score, score);
        % plot(latent)
        
        
        % ---------------------------------------------------------------
        % Orthviews
        % ---------------------------------------------------------------
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
        
        if isempty(fmridat.volInfo)
            disp('.volInfo is empty. Skipping orthviews and other brain plots.');
        else
            create_orthviews(vecs_to_reconstruct, fmridat);
            spm_orthviews_name_axis('Mean data', 1);
            if size(fmridat.dat,2) > 1
                spm_orthviews_name_axis('STD of data', 2);
                spm_orthviews_name_axis('Mean / STD', 3);
            end
            set(gcf, 'Name', 'Orthviews_fmri_data_mean_and_std');
        end
        
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
        imagesc(means);
        colorbar
        axis tight; set(gca, 'YDir', 'Reverse')
        title('Means by condition');
        xlabel('Voxels');
        if iscell(fmridat.Y_names) && ~isempty(fmridat.Y_names)
            set(gca, 'YTick', u, 'YTickLabel', fmridat.Y_names);
        else
            ylabel('Unique Y values');
        end
        
        drawnow
        
        subplot(2, 1, 2)
        imagesc(stds);
        colorbar
        axis tight; set(gca, 'YDir', 'Reverse')
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
            create_orthviews(vecs_to_reconstruct, fmridat);
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
            set(fig_handle, 'Name', 'Montage_mean_across_conditions')
            
            vecs_to_reconstruct = std(means)' ./ mean(means)';
            vecs_to_reconstruct(vecs_to_reconstruct < prctile(vecs_to_reconstruct, 70)) = 0;
            fig_handle = create_montage(vecs_to_reconstruct, fmridat);
            set(fig_handle, 'Name', 'Montage_coeff_of_var_across_conditions')
            

        end
        
        
    otherwise
        error('Unknown plot method');
end

end


function create_orthviews(vecs_to_reconstruct, fmridat)

vecs_to_reconstruct = zeroinsert(fmridat.removed_voxels, vecs_to_reconstruct);

n = size(vecs_to_reconstruct, 2);
overlay = which('SPM8_colin27T1_seg.img');
spm_check_registration(repmat(overlay, n, 1));

for i = 1:n
    
    cl{i} = iimg_indx2clusters(vecs_to_reconstruct(:, i), fmridat.volInfo);
    cluster_orthviews(cl{i}, 'add', 'handle', i);
    
    spm_orthviews_change_colormap([0 0 1], [1 1 0], [0 1 1], [.5 .5 .5], [1 .5 0]);
    
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
    
    set(fig_handle, 'Name', sprintf('Montage %3.0f', i), 'Tag', sprintf('Montage %3.0f', i))
    
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


