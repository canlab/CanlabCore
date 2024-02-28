function animate(data, tr, varargin)
    % Helper script to create condition vectors for hrf_fit_one_voxel()
    % Michael Sun, Ph.D.
    % - Takes fmri_data() object or the number of TRs in a 4D object.
    % - conditions: cell-vector of cellstrings for each condition e.g., {'hot', 'warm', 'imagine'}
    % - onsets: SPM-style onsets cell array in seconds from first-level model, e.g., {{1,2,3}, {4, 5, 6}, {7, 8, 9}}
    % - duration: SPM-style duration cell array in seconds from first-level model, e.g., {{12, 12 ,12}, {12, 12, 12}, {12, 12, 12}} 
    %
    % *Usage:
    % ::
    %    Condition = generateConditionTS(image_obj, {'hot','warm','imagine'}, onsets, durations})
    %
    % Note 1: Preset a SPIKES or SPIKETRAINS variable in order to toggle the
    % generation of Spikes (single 1s) or Spiketrains (a train of 1s) to
    % represent each event.
    %
    % Note 2: If there was slice-timing correction performed, then you will
    % want to correct for it by adding 0.5TRs to all of your times.

    % Default to modeling single spike instead of duration of events.

    % Validate input arguments

%     if nargin < 4
%         error('You must provide fmri_data, tr, outfile name, and outtype "avi" or "gif"');
%     end
    
    if numel(varargin)==1
        outfile=varargin{1};
        outtype='avi';
    elseif numel(varargin)==2
        outfile=varargin{1};
        outtype=varargin{2};
    else
        outfile='animated_file.avi';
        outtype='avi';
    end

    % Use fileparts to split the path
    [parentFolder, ~, ~] = fileparts(outfile);

    % Check if the folder exists
    if ~exist(parentFolder, 'dir')
        % Folder does not exist, so create it
        mkdir(parentFolder);
        disp('Folder was created.');
    end


    if strcmp(outtype,'avi')
        makeFmriAvi(data, tr, outfile);
    elseif strcmp(outtype, 'gif')
        makeFmriGif(data, tr, outfile);
    else
        disp([outtype, ' is not supported yet.'])
    end

end

% function makeFmriAvi(data, tr, outfile)
%     % data: fmri_data object
%     % varargin: 
%     %   'TR', numeric
%     %   'outfile', cellstr
%     %   '', 
%     
%     
%     % Define the animation parameters
%     
%     frame_rate = 1/tr; % e.g., 2.17; for a TR of 0.46
%     
% %     % Create a movie object
% %     writerObj = VideoWriter(outfile);
% %     writerObj.FrameRate = frame_rate;
% %     open(writerObj);
% %     
% %     % Create a loop that iterates over each time point in the fMRI data
% %     for i = 1:size(data.dat,2)
% %         % Replace with preferred plotting logic using 'data'.
% %         get_wh_image(data, i).montage('full', 'cmaprange', [-7 7]);
% %     
% %         % Capture the image as a frame in the movie object
% %         frame = getframe(gcf);
% %         writeVideo(writerObj,frame);
% %         close;
% %     end
% 
%     % Create a movie object
%     writerObj = VideoWriter(outfile);
%     writerObj.FrameRate = frame_rate;
%     open(writerObj);
% 
%     % Preallocate figure and axes
% %     fig = figure;
% %     ax = axes('Parent', fig);
% % 
% %     % Process each frame
% %     parfor i = 1:size(data.dat,2)
% %         % Update plot in the preallocated figure
% %         cla(ax); % Clear axes
% %         get_wh_image(data, i).montage('full', 'cmaprange', [-7 7], 'Parent', ax);
% % 
% %         % Capture the image as a frame in the movie object
% %         frame = getframe(fig);
% %         writeVideo(writerObj, frame);
% %     end
% 
% 
%     % Parallel generation of frames with index
%     % Parallel Approach
%     nFrames = size(data.dat,2);
%     frames(nFrames) = struct('cdata',[],'colormap',[]); % Preallocate structure array
%     
%     parfor i = 1:nFrames
%         frames(i).data = frameCapture(data, i);
%         frames(i).index = i; % Store the index
%     end
%     
%     % Sort frames based on index to ensure correct order
%     [~, order] = sort([frames.index]);
%     sortedFrames = frames(order);
%     
%     % Sequentially write sorted frames to video
%     for i = 1:nFrames
%         writeVideo(writerObj, sortedFrames(i).data);
%     end
% 
% 
%     
%     % Close the movie object
%     close(writerObj);
% %     close(fig);
% 
% end

function makeFmriAvi(data, tr, outfile)




    frame_rate = 1/tr;
    writerObj = VideoWriter(outfile);
    writerObj.FrameRate = frame_rate;
    open(writerObj);

    nFrames = size(data.dat,2);
    tempDir = tempname; % Create a temporary directory name
    mkdir(tempDir); % Make the temporary directory

    % Parallel generation of frames
    parfor i = 1:nFrames
        frameFilename = fullfile(tempDir, sprintf('frame_%06d.png', i));
        frame = frameCapture(data, i);
        imwrite(frame.cdata, frameFilename); % Write frame to file
    end

    % Sequentially read frames from files and write to video
    for i = 1:nFrames
        frameFilename = fullfile(tempDir, sprintf('frame_%06d.png', i));
        if exist(frameFilename, 'file')
            frame = imread(frameFilename);
            writeVideo(writerObj, frame);
        else
            error(['Frame file missing: ', frameFilename]);
        end
    end

    % Close the movie object
    close(writerObj);

    % Cleanup: Delete the temporary frames
    rmdir(tempDir, 's'); % Remove the directory and its contents
end




function makeFmriGif(data, fr, outfile)
    % data: fmri_data object or appropriate data for the animation
    % fr: frame rate
    % outfile: output file name
    
    % Ensure that 'data' contains valid information for the animation.
    % Ensure that 'outfile' is a string and ends with '.gif'.
    
    % Number of frames based on the data
    numFrames = size(data.dat,2);
    mov(numFrames) = struct('cdata',[],'colormap',[]);

    % Loop through each frame of the data
    for i = 1:numFrames
        % Replace with preferred plotting logic using 'data'.
        get_wh_image(data, i).montage('full', 'cmaprange', [-7 7]);
        
        drawnow
        % Capture the current plot as a movie frame
        mov(i) = getframe(gcf);
        close;
    end
    
    % Convert the movie frames to indexed images
    [imind, cm] = rgb2ind(mov(1).cdata, 256, 'nodither');
    imind(1,1,1,numFrames) = 0;
    for i = 1:numFrames
        imind(:,:,1,i) = rgb2ind(mov(i).cdata, cm, 'nodither');
    end
    
    % Save the indexed images as an animated gif
    delay = 1/fr; % delay between frames in seconds
    loopcount = inf; % number of times to repeat animation (inf = indefinitely)
    imwrite(imind, cm, outfile, 'gif', 'DelayTime', delay, 'LoopCount', loopcount);
end


function frame = frameCapture(data, i)

    % Call montage, which creates its own figure
    get_wh_image(data, i).montage('full', 'cmaprange', [-7 7], 'noverbose');

    % Find the most recently created figure
    figs = findall(groot, 'Type', 'figure');
    latestFig = figs(end);

    % Make the figure visible and set the size
    set(latestFig, 'Position', [100, 100, 2560, 1080], 'Visible', 'on');

    drawnow; snapnow;

    % Dynamic waiting for figure to update
%     waitTime = 0;
%     maxWaitTime = 400; % Slightly longer than your longest frame processing time
%     while waitTime < maxWaitTime
%         drawnow; snapnow; % Update figure window
%         pause(5); % Check every 5 seconds
%         waitTime = waitTime + 5;
%         % Add any additional checks here if possible to confirm rendering is complete
%     end

    % Capture the frame from the latest figure
    frame = getframe(latestFig);

    % Close the latest figure
    close(latestFig);


end










