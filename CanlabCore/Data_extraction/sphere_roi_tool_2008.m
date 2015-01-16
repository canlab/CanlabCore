% [cl, all_data] = sphere_roi_tool_2008(imgs, [radius], [n x 3 coord list to use], ['useexisting'])
%
% This function extracts data from a sphere defined around a coordinate you choose.
% You enter a set of images. It will bring up the first image in SPM
% orthviews.  You can click on a voxel in the image, and it will extract
% values from the timeseries.  
%
% It returns an output structure, cl, that has the sphere coordinates,
% extracted data, and other information. It also makes a plot of the average values.
% 
% all_data is all the data from each voxel, [images x voxels]
%
% Example:
%-----------------------------------------------
% imgs = spm_select();
% spm_defaults
% 
% [cl, all_data] = sphere_roi_tool_2008(imgs);
% 
% Now that you have the cl sphere structure defined, you can visualize it
% in other ways:
% Let's make a surface plot (if your ROI is near the left medial surface):
%
% cluster_surf(cl, 5, 'left');
%
% Let's make a montage plot:
% montage_clusters([], cl, {'r'});
%
%
% You can also use this just to get a sphere, of any
% arbitrary radius, with no images to get data from.
% The command below uses a default template image and an 8 mm radius.
% cl = sphere_roi_tool_2008([], 8);
%
% If you enter an n x 3 list of coordinates, it will make spheres around
% them without you entering anything interactively:
% 
% note: the snr in the sphere is mean signal for all voxels over time,
% divided by the *spatial* sd, ignoring variance over time
%
% Tor Wager, 2008; Update 2011 for automatic coords

function [cl, all_data] = sphere_roi_tool_2008(imgs, varargin)

radius = 10;
useexisting = 0;
coords = [];
if length(varargin) > 0, radius = varargin{1}; end

for i = 1:length(varargin)
    if ~ischar(varargin{i}) && size(varargin{i}, 2) == 3
        % coords
        coords = varargin{i};
    end
end

if any(strcmp(varargin, 'useexisting'))
    useexisting = 1; 
end

if isempty(imgs), imgs = which('avg152T1.nii'); end

if ~useexisting
    spm_image('init', imgs(1,:))
end

create_figure('Average data values', 1, 1);

quitnow = [];
cl = [];

while isempty(quitnow)

    if isempty(coords)
        % graphical selection of coordinates
        
        quitnow = input('Click on sphere center and press return, or 1+return to quit');
        
    else
        % automatically decide if done
        if length(cl) == size(coords, 1)
            quitnow = 1;
        end
    end
    
    if quitnow > 0
        disp('Exiting')
        continue

    else
        if isempty(coords)
            pos = spm_orthviews('Pos');
        else
            pos = coords(length(cl) + 1, :);
        end
        
        cl = [cl get_sphere(pos, imgs, radius)];

        create_figure('Average data values', 1, 1, 1);
        plot(cl(end).timeseries, 'ko-');
        xlabel('Images'); ylabel('Values');
    end

end

all_data = cl(end).all_data;

end



function cl = get_sphere(pos, imgs, radius)

    spm_orthviews('Reposition', pos);

    disp('Extracting data');
    [data, XYZvoxSphere, XYZmmSphere] = iimg_sphere_timeseries(imgs, pos, radius);

    datamean = mean(data, 2);
    glev = mean(datamean);
    snr = glev ./ std(datamean);

    % Viz sphere
    V = spm_vol(deblank(imgs(1,:)));

    cl(1).XYZ = XYZvoxSphere;
    cl(1).XYZmm = XYZmmSphere;
    cl(1).Z = ones(1, size(XYZvoxSphere, 2));
    cl(1).title = 'Sphere';
    cl(1).voxSize = diag(V(1).mat(1:3, 1:3))';
    cl(1).M = V(1).mat; % use first frame only if 4-D image
    cl(1).timeseries = datamean;
    cl(1).snr = snr;
    cl(1).global_lev = glev;
    cl(1).all_data = data;
    
    cl(1).mm_center = mean(cl(1).XYZmm, 2)';
    cl(1).numVox = size(cl(1).XYZ, 2);
    
    cl(1).dim = V.dim;
    
    cl(1).title = 'Sphere';
    
    cluster_orthviews(cl,  {[1 0 0]}, 'add', 'solid');
    drawnow

end

