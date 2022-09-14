function [rmssd, rmssd_outlier_regressor_matrix] = rmssd_movie(dat, varargin)
% Movie of successive differences (sagittal slice)
% Enter an image_vector or fmri_data object (usually with time series)
%
% :Usage:
% ::
%
% [rmssd, rmssd_outlier_regressor_matrix] = rmssd_movie(dat, [full_path_of_movie_output_file,image_skip_interval])
%
% Images usually change slowly over time, and sudden changes in intensity can also often be a sign of bad things
% -- head movement artifact or gradient misfires, interacting with the magnetic field to create distortion
% across the brain.
% RMSSD tracks large changes across successive images, regardless of what the sign of the changes is or where they are.
% In addition, images with unusually high spatial standard deviation across voxels may be outliers with image
% intensity distortions in some areas of the image but not others (e.g., bottom half of brain vs. top half,
% or odd vs. even slices).
% The CANlab method rmssd_movie( ), for fmri_data objects, creates a visual movie so you can see what the
% image-to-images changes are. It pauses where they're unusual. It also returns a matrix rmssd_outlier_regressor_matrix,
% which has an indicator regressor (1 or 0 values) for every image that is quite different from the preceding ones
% (the pause point in the movie). This is based on two things: (1) rmssd being > a cutoff number of standard
% deviations from the mean, (2) spatial standard deviation of the images being > a cutoff number of
% standard deviations from the mean. The cutoff is 3 s.d. by default.
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
%    Show a movie of RMSSD and write a movie file to disk in the qc_images subdirectory:
%    rmssd_movie(fmri_dat, 'movie_output_file', '/Subj1/qc_images/rmssd_movie', 'image_interval', 5)
%
%   This would save an movie based on the images in fmri_dat to the
%   above directory, with an interval of 5 images between each
%   frame (so, the movie would show image 1, 6, 11, 16, etc)
%

% Programmers Notes:
% ..
%    Edited 8/7/14 by Scott
%      - added skip interval
%      - updated help
%
%    8/7/14 Programmer Note: if more varargin options are desired in the
%    future, the function call will likely need to be re-written. The
%    current form exists for backwards compatibility - obviously changing
%    the function call will mean that other functions that use this
%    (preproc) will need to be modified
%
%    Modified 5/2021 by Tor Wager, to return rmssd values and covariates
%    based on estimated outliers at 3 standard deviations.
% ..

% -------------------------------------------------------------------------
% DEFAULT ARGUMENT VALUES
% -------------------------------------------------------------------------
sdlim = 3;

writetofile = false;
movieoutfile = [];

showmovie = true;
image_interval = 1;

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

allowable_inputs = {'sdlim' 'movieoutfile' 'showmovie' 'image_interval'};

keyword_inputs = {'writetofile' 'nomovie' 'nodisplay'};

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
                    movieoutfile = fullfile(pwd, 'rmssd_movie_file');
                    fprint('Writing file with default name:\n%s\n', movieoutfile);
                end
                
            case {'nomovie' 'nodisplay'}
                showmovie = false;
                
        end
    end
end

% -------------------------------------------------------------------------
% MAIN FUNCTION
% -------------------------------------------------------------------------

mm = mean(dat);
sdiffs = diff(dat.dat')';
sdiffs = [mean(sdiffs, 2) sdiffs]; % keep in image order

%mymean = mean(dat.dat); % global mean
spatsd = std(sdiffs);  % variation in successive diffs - like rmssd but without abs mean diff. artifacts in some parts of image...
rmssd = ( mean(sdiffs .^ 2) ) .^ .5; % rmssd - root mean square successive diffs

% avoid first time point being very different and influencing distribution and plots.
spatsd(1) = median(spatsd);
rmssd(1) = median(rmssd);

% z-scores of rmssd
% rmssd_z = scale(rmssd);

% used later, for plot axes
mysd = std(sdiffs(:));
mylim = [mean(sdiffs(:)) - sdlim*mysd mean(sdiffs(:)) + sdlim*mysd];

% Get potential outliers to pause at
slow1 = abs(rmssd) > mean(rmssd) + sdlim*std(rmssd);
slow2 = abs(spatsd) > mean(spatsd) + sdlim*std(spatsd);
slow = slow1 | slow2;

%rmssd_outlier_regressor_matrix = full(ind2vec(find(slow)', size(dat.dat, 2)));

rmssd_outlier_regressor_matrix = intercept_model(size(dat.dat, 2), find(slow));
rmssd_outlier_regressor_matrix = rmssd_outlier_regressor_matrix(:, 2:end);  % remove initial intercept column

% -------------------------------------------------------------------------
% DISPLAY MOVIE
% -------------------------------------------------------------------------

if showmovie
    
    fh = create_figure('succ_diffs', 2, 1);
    ax1 = subplot(2, 1, 1);
    title('Root mean square successive diffs');
    
    ax2 = subplot(2, 1, 2);
    set(ax1, 'Position', [.13 .75 .77 .2]);
    
    ax3 = axes('Position', [.13 .48 .77 .2]);
    set(ax3, 'FontSize', 16);
    
    axes(ax3);
    title('STD of successive diffs');
    hold on;
    axes(ax2);
    
    vdat = reconstruct_image(mm);
    wh = round(size(vdat, 1)./2);
    
    plot(ax1, rmssd, 'k'); axis tight
    axes(ax1), hold on;
    hh = plot_horizontal_line(mean(rmssd) + sdlim*std(rmssd));
    set(hh, 'LineStyle', '--');
    
    plot(ax3, spatsd, 'k'); axis tight
    axes(ax3), hold on;
    hh = plot_horizontal_line(mean(spatsd) + sdlim*std(spatsd));
    set(hh, 'LineStyle', '--');
    
    axes(ax2);
    imagesc(squeeze(vdat(wh, :, :))', mylim);
    drawnow
    hold off
    colormap gray
    
    for i = 1:image_interval:size(sdiffs, 2)
        
        vh = plot(ax1, i, rmssd(i), 'ro', 'MarkerFaceColor', 'r');
        vh2 = plot(ax3, i, spatsd(i), 'ro', 'MarkerFaceColor', 'r');
        
        mm.dat = sdiffs(:, i);
        vdat = reconstruct_image(mm);
        imagesc(squeeze(vdat(wh, :, :))', mylim);
        axis image;
        set(ax2, 'YDir', 'Normal')
        xlabel('Successive differences (sagittal slice)');
        
        drawnow
        
        if writetofile
            F = getframe(fh);
            if i == 1
                imwrite(F.cdata, movieoutfile,'tiff', 'Description', dat.fullpath(1,:), 'Resolution', 30); % Edit description so that imwrite doesn't crash if dat is a cell-array instead of of text - Michael Sun 2/17/2022
            else
                imwrite(F.cdata, movieoutfile,'tiff', 'WriteMode', 'append', 'Resolution', 30);
            end
            
        elseif slow(i)
            pause(1);
        end
        
        delete(vh);
        delete(vh2);
    end
    
end % showmovie

end % main function


