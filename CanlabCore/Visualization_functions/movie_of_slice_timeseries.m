function movie_of_slice_timeseries(imgs, slicenumber, moviename, orientation)
% :Inputs:
%
%   **imgs:**
%        Image names
%
%   **slicenumber:**
%        which slice to view (in volume)
%
%   **moviename:**
%        something like: 'slice10_timeseries.avi'
%
%   **orientation:**
%        'axial' or 'sagittal'
%
% :Examples:
% ::
%
%    slicenumber = 10;
%    moviename = 'slice10_timeseries.avi';

    if(~exist('orientation', 'var') || isempty(orientation))
        orientation = 'axial';
    end

    disp('Mapping volumes.')
    V = spm_vol(imgs);

    disp('Loading data.')

    sl = timeseries_extract_slice(V, slicenumber, orientation);

    disp('Setting limits.')
    mymad = mad(sl(:));
    mymean = mean(sl(:));
    limits = [0 mymean + 3 * mymad];


    slice_fig = figure();

    nimgs = size(sl, 3);

    disp('Starting movie.')
    movlength = 60; % in seconds

    [mov, nframes, axh] = setup_movie([], movlength, moviename);

    for i = 1:nimgs
        H = gca;
        imagesc(sl(:, :, i));
        set(H, 'CLim', limits)
        colorbar
        axis image

        [dd, ff, ee] = fileparts(V(i).fname);
        title(sprintf('Image %3.0f: %s', i, ff), 'FontSize', 24)

        mov = add_a_frame(mov, H);
    end

    tmp = close(mov); %#ok;
    close(slice_fig);
    
    disp('Success!')
end


function [mov, nframes, axh] = setup_movie(mov, movlength, moviename)
    axh = gca;

    fps = 5;
    nframes = movlength .* fps;

    if isempty(mov)
        mov = avifile(moviename, 'Quality', 75, 'Compression', 'None', 'Fps', fps);
    end

    % add to existing
    %O = struct('add2movie', [], 'zoom', 1, 'azOffset', [], 'elOffset', [], 'timecourse', [], 'fps', 5, 'length', 6);

    axis image
end



function mov = add_a_frame(mov, H)
    drawnow

    lightRestoreSingle(H);
    try
        mov = addframe(mov, H);
    catch
        disp('Cannot write frame.  Failed to set stream format??')
        mov = close(mov);
    end
end

