function [g, spikes, gtrim, nuisance_covs, spikesperimg, snr] = scn_session_spike_id(imgs, varargin)
% Gets global image values for a session, and uses trimts.m to find
% outliers. The optional input MADs allows one to lower or raise the
% threshold for identifying scans as spikes (default = 10).
%
% :Usage:
% ::
%
%     [g, spikes, gtrim, nuisance_covs, spikesperimg, snr] = scn_session_spike_id(imgs,'mask',[mask name],'MADs',[MADs],'doplot',[0/1])
%
% Multi-session mode returns much more output and more images, and
% takes in a cell array with images (preferably 4-D) for each session
% (run).
% 
% :Inputs:
%
%   **'mask',[pathtomaskfile]:**
%        mask images using the mask in pathtomaskfile, default: implicit mask
%
%   **'MADs',[scalar]:**
%        change Mahalanobis distance, default: 10
%
%   **'doplot',[0 / 1]:**
%        plot result figures, default: true
%
% Returns:
%
%   **g:**
%        global values
%
%   **spikes:**
%        identified spikes
%
%   **gtrim:**
%        trimmed/adjusted global values, can be used as covariate in GLM
%
%
%   **nuisance_covs:**
%        a matrix of 1)gtrim and 2) dummy regressors that can be used to minimize
%        spike influence in GLM
%
% We may want to save norms on the number of outliers found.
%
% :Examples:
% ::
%
%    % Get image names
%    for i = 1:6, sess_images{i} = filenames(sprintf('run%02d/vol0*img', i), 'char', 'absolute'); end
%
%    % Run
%    [g, spikes, gtrim, nuisance_covs, snr] = scn_session_spike_id(sess_images);
%
% ..
%    Tor Wager
%    figure options, new input format Stephan
% ..

basedir = pwd;
yamlfilename = 'qc_results.yaml'; % for database integration
global doplot


% Cell mode

if iscell(imgs)
    % ---------------------------------------------------------------------
    % ---------------------------------------------------------------------
    % Multi-session
    % ---------------------------------------------------------------------
    % ---------------------------------------------------------------------
    
    nsessions = length(imgs);
    
    g = cell(1, nsessions);
    spikes = cell(1, nsessions);
    gtrim = cell(1, nsessions);
    nuisance_covs = cell(1, nsessions);
    spikesperimg = cell(1, nsessions);
    snr = cell(1, nsessions);
    
    fprintf('\nMulti-session spike ID output will be written to:\n%s\n', basedir);
    
    for i = 1:nsessions
        fprintf('\n-------------------------------------\nSession %3.0f of %3.0f\n', i, nsessions)
        [g{i}, spikes{i}, gtrim{i}, nuisance_covs{i}, spikesperimg{i}, snr{i}] = scn_session_spike_id(imgs{i},varargin{:});
        
    end
    
    summarize_multisession_output(doplot);
    
    % write a .yaml file with info to upload to QC database
    append_to_yaml_file;
    
    
    
else
    % ---------------------------------------------------------------------
    % ---------------------------------------------------------------------
    % Single-session
    % ---------------------------------------------------------------------
    % ---------------------------------------------------------------------
    
    g = [];
    spikes = [];
    gtrim = [];
    nuisance_covs = [];
    spikesperimg = [];
    snr = [];
    
    if isempty(imgs), disp('No Images.'); return, end
    
    % defaults
    stdev = 3;
    useimplicitmask = 1;
    doplot = 1;
    MADs = 10;
 
    % parse inputs
    for k = 1:length(varargin)
        if ischar(varargin{k})
            switch varargin{k}
                case 'mask',
                    useimplicitmask = 0;
                    maskvalorname = varargin{k + 1};
                    fprintf(1,'\nUsing input mask: %s\n',maskvalorname);
                    
                case 'MADs', MADs = varargin{k + 1};
                    
                case 'doplot', doplot = varargin{k + 1};
            end
        end
    end
    
    % create implicit mask (default)
    if useimplicitmask == 1
        maskvalorname = 'implicit_mask.img';
        [dummy, dummy, inmaskvox] = fmri_mask_thresh_canlab(imgs, maskvalorname,'mean',doplot);
    end
    if ~exist(maskvalorname, 'file')
        fprintf('Mask file does not exist:\n%s\n', maskvalorname)
        error('scn_session_spike_id: Quitting')
    end
    
    
    if doplot
        % Stuff to save file (png)
        % ---------------------------------------------------------------------
        qcdir = fullfile(basedir, 'qc_images');
        if ~exist('qc_images', 'dir'), mkdir('qc_images'); end
        
        i = 1;
        fname = ['qc_images' filesep 'scn_session_spike_s' num2str(i) '.png'];
        while exist(fname, 'file')
            i = i + 1;
            fname = ['qc_images' filesep 'scn_session_spike_s' num2str(i) '.png'];
        end
        fprintf('Will save image file: %s\n', fname);
    end
    
    % Do the work
    % ---------------------------------------------------------------------
    [g, gslice, stdslice] = tor_global(imgs, maskvalorname);
    wh_no_data = ~any(gslice')' | ~any(stdslice')';
    gslice(wh_no_data,:) = [];
    stdslice(wh_no_data,:) = [];
    
    if doplot
        % save tor_global images
        drawnow
        h = findobj('Tag', 'Implicit Mask');
        if ~isempty(h)
            fname2 = fullfile('qc_images', ['implicit_mask_histogram_s' num2str(i) '.png']);
            figure(h)
            scn_export_papersetup(450);
            saveas(h, fname2);
            %close(h);
        end
        
        h = findobj('Tag', 'montage_axial'); %'SCNlab_Montage');
        if ~isempty(h)
            fname2 = fullfile('qc_images', ['implicit_mask_montage_s' num2str(i) '.png']);
            figure(h)
            %scn_export_papersetup(450);
            scn_export_papersetup(300); % probs with invalid drawable otherwise
            saveas(h, fname2);
            %close(h);
        end
    end
    
    
    % mahalanobis: strange patterns across slices
    d2 = mahal([stdslice' gslice'], [stdslice' gslice']);
    
    
    [gtrim, dummy, gspikes] = trimts(g, stdev, [], 1, 3, MADs);
    
    
    [dummy, dummy, mahalspikes] = trimts(d2, stdev, [], 1, 3, MADs);
    spikes = unique([gspikes; mahalspikes]);
    
    spikesperimg = 100*length(spikes) ./ length(g);
    snr = mean(g) ./ std(g);
    
    fprintf('%3.0f Potential outliers\t%%Spikes: %3.2f\tGlobal SNR (Mean/STD): %3.2f\n', length(spikes), spikesperimg , snr)
    
    nuisance_covs = intercept_model(length(gtrim), spikes);
    nuisance_covs(:, 1) = gtrim;
    
    if doplot, 
        plot_session_figure();
    end
    
    
end

% ---------------------------------------------------------------------
% ---------------------------------------------------------------------
% Inline functions
% ---------------------------------------------------------------------
% ---------------------------------------------------------------------

    function summarize_multisession_output(plotfigs)
        %         fprintf('\n')
        %         for i = 1:length(spikes), fprintf('%3.0f Potential outliers\n', length(spikes{i})), end
        fprintf('\nTotal spikes:\t%3.0f\t%%Spikes/Image:\t%3.2f\tAvg SNR:\t%3.2f\n', ...
            length(cat(1, spikes{:})), mean(cat(2, spikesperimg{:})), mean(cat(2, snr{:})))
        
        save scn_session_spike_id_output g spikes gtrim nuisance_covs spikesperimg snr
        
        % count images per run and add lines for sessions
        for i = 1:nsessions
            % If 4-D images, then print size of file and number of volumes (3-D
            % acquisitions)
            if size(imgs{i},1) == 1
                % assume 4-D
                str4d{i} = '4-D images';
                n_images_per_run(i) = scn_num_volumes(imgs{i});
            else
                str4d{i} = '3-D images';
                n_images_per_run(i) = size(imgs{i},1);
            end
        end
        fprintf('\nimages_per_run:\t');
        fprintf('%3.0f\t', n_images_per_run);
        fprintf('\n');
        
        gg = cat(1, g{:});

        if plotfigs
            create_figure('All globals');
            plot(gg, 'k', 'LineWidth', 2);
            title('All global values. Green = Sessions, Red = outliers');
            hh = plot_vertical_line(cumsum(n_images_per_run) + .5, 'g');
            set(hh, 'LineWidth', 2);
        end
        
        cc = [0 cumsum(n_images_per_run)];
        for i = 1:nsessions
            whspikes_sess = spikes{i} + cc(i);
            if plotfigs
                hh = plot(whspikes_sess, gg(whspikes_sess), 'ro', 'LineWidth', 2);
            end
        end
        
        
        % build model for variance explained
        cumulative_spikes = [];
        for i = 1:nsessions
            cumulative_spikes = [cumulative_spikes; spikes{i} + cc(i)];
        end
        
        allspikes = zeros(sum(n_images_per_run), length(cumulative_spikes));
        for i = 1:length(cumulative_spikes)
            allspikes(cumulative_spikes(i), i) = 1;
        end
        
        Xint = intercept_model(n_images_per_run);
        Xint(:, end) = []; % no need for rank deficient mtx
        
        for i = 1:nsessions, lineffect{i} = (1:n_images_per_run(i))' - mean((1:n_images_per_run(i))); end
        lineffect = blkdiag(lineffect{:});
        
        % get numbers for % var due to session, drift within session, etc.
        % allspikes + intercepts + linear drift + error
        g_std = std(gg);
        
        X = [Xint lineffect allspikes]; X(:, end+1) = 1;
        r = gg - X * pinv(X) * gg;
        res_var = var(r);
        res_std = std(r);
        
        X = [lineffect allspikes]; X(:, end+1) = 1;
        r = gg - X * pinv(X) * gg;
        blk_std = sqrt(var(r) - res_var);
        
        X = [Xint allspikes]; X(:, end+1) = 1;
        r = gg - X * pinv(X) * gg;
        lin_std = sqrt(var(r) - res_var);
        
        X = [Xint lineffect]; X(:, end+1) = 1;
        r = gg - X * pinv(X) * gg;
        spike_std = sqrt(var(r) - res_var);
        
        nuisanceX_sess_lin_spike = [Xint lineffect allspikes];
        
        fprintf('\nSources of variation in global signal (std, raw units)\n');
        fprintf('Total variation:\t%3.2f\n', g_std);
        fprintf('Session effects:\t%3.2f\n', blk_std);
        fprintf('Linear trends:\t%3.2f\n', lin_std);
        fprintf('Spikes:\t%3.2f\n', spike_std);
        fprintf('Within-session error:\t%3.2f\n', res_std);
        
        
        disp('Saved in scn_session_spike_id_output.mat :');
        disp('-----------------------------------------------');
        disp('Images embedded in html files and qc_report subdirectory');
        disp('g: global in-brain mean for each session');
        disp('spikes: indices of which images are global outliers');
        disp('spikesperimg: 100*spikes / num images');
        disp('snr: global mean / global stdev');
        disp('gtrim: trimmed global in-brain mean for each session');
        disp('nuisance_covs: spike nuisance covariates (per session)');
        disp(' ');
        disp('cumulative_spikes, allspikes: indices of outlier images across all runs');
        disp('nuisanceX_sess_lin_spike: nuisance covariate design matrix across all runs');
        disp('res_std: variation in global signal related to within-session variation');
        disp('blk_std: variation in global signal related to session effects');
        disp('lin_std spike_std: linear variation in global signal across time (global drift)');
        disp('spike_std: variation in global signal related to estimated outliers');
        disp('-----------------------------------------------');
        disp(' ');
        
        save scn_session_spike_id_output -append cumulative_spikes allspikes nuisanceX_sess_lin_spike res_std blk_std lin_std spike_std
        
        if plotfigs
            create_figure('sources of global variation'); pie([blk_std lin_std spike_std res_std].^2, {'Session', 'Linear', 'Spikes', 'Residual'});
            colormap winter
            hh = findobj(gcf,'Type', 'text'); set(hh, 'FontSize', 24);
            axis off
            scn_export_papersetup;
            fname2 = fullfile('qc_images', 'Sources_of_variance.png');
            saveas(gcf, fname2);
            
            X = [Xint lineffect allspikes]; X(:, end+1) = 1;
            r = gg - X * pinv(X) * gg;
            create_figure('FFT of unexplained global variation');
            fft_plot_scnlab(r, 1, 'samefig');
            xlabel('Frequency (1/images)');
            scn_export_papersetup;
            fname2 = fullfile('qc_images', 'FFT_of_unexplained_global_signal.png');
            saveas(gcf, fname2);
        end
    end



    function plot_session_figure()
        fht = create_figure('scn_session_spike_id', 4, 1);
        % adjust figure size to proper width (SG. 8/25/2017);
        screensz = get(0, 'screensize');
        set(fht,'Position',[30 5 screensz(3)*.85 screensz(4)*.95]);
        
        if isempty(gslice)
            warning('scn_session_spike_id:emptyData', 'Global slice data is empty!! No voxels in mask??');
            return
        end
        
        imagesc(scale(gslice')')
        colorbar('horiz'), title('Global mean by slice x time (zscore within slice)');
        xlabel('Time'); ylabel('Slice Number');
        axis auto, axis tight
        
        subplot(4, 1, 2);
        imagesc(scale(stdslice')')
        colorbar('horiz'), title('Spatial standard deviation by slice x time (zscore within slice)');
        xlabel('Time'); ylabel('Slice Number');
        axis auto, axis tight
        
        subplot(4, 1, 3);
        plot(g, 'k'); hold on; plot(gtrim,'b');
        legend({'Global' 'Adjusted'})
        plot_vertical_line(gspikes, 'r');
        title(['Global val: Spikes Identified: ' num2str(length(gspikes)) ' Global SNR: ' num2str(snr)]);
        set(gca, 'XLim', [0 length(g)]);
        
        subplot(4, 1, 4);
        plot(d2,'b');
        plot_vertical_line(mahalspikes, 'r');
        title(['Mahalanobis dist: Spikes Identified: ' num2str(length(mahalspikes)) ' Global val: ' num2str(mean(d2))]);
        set(gca, 'XLim', [0 length(g)]);
        
        % this *may* be needed for publishing to html
        try
            drawnow;
            snapnow;
        catch
        end
        
        drawnow;
        orient(gcf, 'portrait')
        %%print('-dpsc2', '-append', 'qc_report');
        scn_export_papersetup(800);
        saveas(gcf,fname);
        
    end




    function append_to_yaml_file
        DB.subject_dir = basedir;
        
        DB.mean_global = mean(cat(1, gtrim{:})); % inmaskvox
        DB.mean_spikes_per_image = mean(cat(2, spikesperimg{:}));
        DB.mean_temporal_snr = mean(cat(2, snr{:}));
        DB.std_temporal_snr_across_runs = std(cat(2, snr{:}));
        
        for i = 1:nsessions
            gmean(i) = mean(gtrim{i});
            gstd(i) = std(gtrim{i});
        end
        
        DB.std_global_across_runs = std(gmean);
        
        r = corrcoef(gmean, gstd);
        
        if length(r) > 1 % only for multiple runs
            DB.corrcoef_global_mean_std_across_runs = r(1,2); % possibly useful for testing whether % signal change reduces inter-session variance
        end
        
        load scn_session_spike_id_output res_std blk_std lin_std spike_std
        DB.global_std_description = 'Sources of nuisance variance, in standard deviation units.';
        DB.global_std_withinsession = res_std;
        DB.global_std_sessioneffect = blk_std;
        DB.global_std_lineardrift = lin_std;
        DB.global_std_spikes = spike_std;
        
        struct2yaml(yamlfilename, DB, 'add', 'replace');
        
    end

end  % main function

