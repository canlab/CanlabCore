function outlier_tables = slice_movie(dat, varargin)
% Movie of slice timeseries (sagittal slice)
% Enter an image_vector or fmri_data object (usually with time series)
%
% :Usage:
% ::
%
%  [outlier_table, outlier_regressor_matrix_uncorr, outlier_regressor_matrix_corr] = slice_movie(dat, [full_path_of_movie_output_file,image_skip_interval])
%
% Images usually change slowly over time, and sudden changes in intensity can also often be a sign of bad things
% -- head movement artifact or gradient misfires, interacting with the magnetic field to create distortion
% across the brain.
% RMSSD tracks large changes across successive images, regardless of what the sign of the changes is or where they are.
% In addition, images with unusually high spatial standard deviation across voxels may be outliers with image
% intensity distortions in some areas of the image but not others (e.g., bottom half of brain vs. top half,
% or odd vs. even slices).
% The CANlab method slice_movie( ), for fmri_data objects, creates a visual movie so you can see what the
% image-to-images changes are. It pauses where they're unusual, as defined by the method fmri_data.outliers( ). 
% It also returns a table, outlier table,  and matrices you can use as covariates in design
% outlier_regressor_matrix_uncorr and outlier_regressor_matrix_corr.
% These have an indicator regressor (1 or 0 values) for every image that is quite different from the preceding ones
% (the pause point in the movie). This is based on two things: (1) rmssd and others being > a cutoff number of standard
% deviations from the mean, (2) spatial standard deviation of the images being > a cutoff number of
% standard deviations from the mean. The cutoff is 3 median absolute deviations. by default.
% This matrix can be added to your design matrix as a set of nuisance covariates of no interest.
%
% *Optional Inputs:
%
%   **'movieoutfile', <filepath>:**
%        Followed by a char array detailing the full path to save the
%           movie file
%
%   **'image_interval', n:**
%        Followed by an integer value describing the interval
%        between images in each subsequent frame of the movie
%       (default = 1). Higher values will skip, showing every n images
%
% :Examples:
% ::
%
%   obj = load_image_set('kragel18_alldata');
%   obj2 = rescale(obj, 'l2norm_images');
%   slice_movie(obj2);
%
%   Show a movie of RMSSD and write a movie file to disk in the qc_images subdirectory:
%   rmssd_movie(fmri_dat, 'movie_output_file', '/Subj1/qc_images/rmssd_movie', 'image_interval', 5)
%
%   This would save an movie based on the images in fmri_dat to the
%   above directory, with an interval of 5 images between each
%   frame (so, the movie would show image 1, 6, 11, 16, etc)
%

% Programmers Notes:
% ..
%
% ..

% -------------------------------------------------------------------------
% DEFAULT ARGUMENT VALUES
% -------------------------------------------------------------------------
madlim = 3;     % median absolute deviation limits

writetofile = false;
movieoutfile = [];

showmovie = true;
image_interval = 1;

doverbose = true;

domontage = false;

% %deal out varargin
% i = 1;
% while i <= length(varargin)
%     if ischar(varargin{i})
%         movieoutfile = varargin{i};
%         writetofile = 1;
%     elseif isnumeric(varargin{i})
%         image_interval = varargin{i};
%     end
%     i=i+1;
% end
%

% -------------------------------------------------------------------------
% OPTIONAL INPUTS
% -------------------------------------------------------------------------

% This is a compact way to assign multiple variables. The input argument
% names and variable names must match, however:

allowable_inputs = {'madlim' 'movieoutfile' 'showmovie' 'image_interval'};

keyword_inputs = {'writetofile' 'nomovie' 'nodisplay', 'montage'};

% optional inputs with default values - each keyword entered will create a variable of the same name

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            
            case allowable_inputs
                
                eval([varargin{i} ' = varargin{i+1}; varargin{i+1} = [];']);
                
            case keyword_inputs
                % Skip, deal with these below
                
            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end

% 2nd pass: Keyword inputs. These supersede earlier inputs
for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            
            case 'writetofile'
                writetofile = true;
                if isempty(movieoutfile)
                    movieoutfile = fullfile(pwd, 'slice_movie_file.tiff');
                    fprintf('Writing file with default name:\n%s\n', movieoutfile);
                end
                
            case {'nomovie' 'nodisplay'}
                showmovie = false;
            
            case ('montage')
                domontage = true;
                
        end
    end
end

% -------------------------------------------------------------------------
% MAIN FUNCTION
% -------------------------------------------------------------------------

% Get potential outliers to pause at

[slow, ~, outlier_tables] = outliers(dat, 'madlim', madlim, 'doverbose', doverbose);

mm = mean(dat);


% -------------------------------------------------------------------------
% DISPLAY MOVIE
% -------------------------------------------------------------------------

if showmovie
    
    fh = create_figure('image_view', 2, 2); % set fig size using 2,2
    for i = 1:4, myax = subplot(2, 2, i); delete(myax); end
    
    ax1 = axes('Position', [.1 .75 .8 .2], 'FontSize', 16); 
    title('Root mean square successive diffs');
    hold on
    
    ax2 = axes('Position', [.1 .1 .4 .55], 'FontSize', 16);
    title('Sagittal slices');
    
    % slice axes
    ax3 = axes('Position', [.55 .1 .4 .55], 'FontSize', 16);
    title('Axial slices'); 
    
    vdat = reconstruct_image(mm);
    wh = round(size(vdat, 1)./2);     % Middle Sagittal Slice
    wh_z = round(size(vdat, 3)./2);   % Middle Axial Slice
    
    sdiffs = diff(dat.dat')';
    sdiffs = [mean(sdiffs, 2) sdiffs]; % keep in image order
    mysd = std(sdiffs(:));
    mylim = [mean(vdat(:)) - madlim*mysd mean(vdat(:)) + madlim*mysd];

    plot(ax1, outlier_tables.score_table.rmssd_dvars, 'k'); axis tight
    axes(ax1), hold on;
    cutoff_val = mean(outlier_tables.score_table.rmssd_dvars) + madlim * std(outlier_tables.score_table.rmssd_dvars);
    hh = plot_horizontal_line(cutoff_val);
    set(hh, 'LineStyle', '--');
    
    plot(find(slow), cutoff_val * ones(sum(slow, 1)), 'o', 'Color', [1 .5 0], 'MarkerFaceColor', [.7 .3 0]);
    
    axes(ax2); 
    if domontage
        display_slices(get_wh_image(dat, 1), 'sagittal');
    else
        get_sagittal_slice(vdat, wh, mylim)
    end
    
    drawnow
    hold off
    colormap gray
    
    axes(ax3)
    if domontage
        display_slices(get_wh_image(dat, 1), 'axial');
    else
        get_axial_slice(vdat, wh_z, mylim)
    end

    % Get color limits
    clim = prctile(dat.dat(:), [.05 99.5]);
    
    for i = 1:image_interval:size(dat.dat, 2)  % changed from sdiffs
        
        vh = plot(ax1, i, outlier_tables.score_table.rmssd_dvars(i), 'ro', 'MarkerFaceColor', 'r');
        
        axes(ax2)
        if domontage
            display_slices(get_wh_image(dat, 1), 'sagittal');
        else
            get_sagittal_slice(vdat, wh, mylim)
        end
        
        axes(ax3)
        if domontage
            display_slices(get_wh_image(dat, 1), 'axial');
        else
            get_axial_slice(vdat, wh_z, mylim)
        end
        
        set([ax2 ax3], 'CLim', clim)
        
        drawnow
        
        if writetofile
            F = getframe(fh);
            if i == 1
                imwrite(F.cdata, movieoutfile,'tiff', 'Description', dat.fullpath(1,:), 'Resolution', 30);
            else
                imwrite(F.cdata, movieoutfile,'tiff', 'WriteMode', 'append', 'Resolution', 30);
            end
            
        elseif slow(i)
            pause(1);
        end
        
        delete(vh);
%         delete(vh2);
    end
    
end % showmovie

end % main function


function get_sagittal_slice(vdat, wh, clim)
    imagesc(rot90(squeeze(vdat(wh,:,:))), clim);
    axis image
    title('Sagittal Slice')
end


function get_axial_slice(vdat, wh, clim)
    imagesc(squeeze(vdat(:,:,wh)), clim)
    axis image
    title('Axial Slice')
end