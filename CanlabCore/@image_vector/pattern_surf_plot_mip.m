function [mip, x, y, voldata] = pattern_surf_plot_mip(m, varargin)
% pattern_surf_plot_mip Axial maximum intensity projection pattern surface plot.
%
% Reconstructs the image_vector data into a 3-D volume, optionally smooths
% it with a 3-D Gaussian, takes the axial (Z) maximum intensity
% projection, trims empty borders, and renders the result as both a 2-D
% image and a 3-D surface plot.
%
% :Usage:
% ::
%
%     [mip, x, y, voldata] = pattern_surf_plot_mip(m, [optional inputs])
%
% :Inputs:
%
%   **m:**
%        image_vector (e.g., fmri_data or statistic_image) object.
%
% :Optional Inputs:
%
%   **'nofigure':**
%        Skip the figure creation / plotting step.
%
%   **'nosmooth':**
%        Skip 3-D Gaussian smoothing of the volume (default: smoothing on).
%
%   **'smoothbox':**
%        Smoothing kernel size (passed to smooth3). Default: 5.
%
%   **'sd':**
%        Gaussian SD for smoothing (passed to smooth3). Default: 2.
%
%   **'xlim' / 'ylim' / 'zlim':**
%        Spatial limits (currently parsed; reserved for future use).
%
%   **'noverbose':**
%        Suppress verbose output.
%
% :Outputs:
%
%   **mip:**
%        The maximum intensity projection (Y-by-X) with empty borders
%        trimmed and out-of-mask voxels set to NaN.
%
%   **x:**
%        Vector of x mm coordinates corresponding to mip columns.
%
%   **y:**
%        Vector of y mm coordinates corresponding to mip rows.
%
%   **voldata:**
%        The trimmed (and optionally smoothed) 3-D volume reconstructed
%        from m. Not used by this function but returned for downstream
%        operations such as rigid-body transformation.
%
% :See also:
%   - reconstruct_image
%   - smooth3
%   - voxel2mm
%
% ..
%    Tor Wager, 2016
%    Update: Aug 2016
% ..

% ------------------------------------------------------
% parse inputs
% ------------------------------------------------------
dofigure = true;
dosmooth = true;
mysmoothbox = 5;
sd = 2;

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            % reserved keywords
            
            case 'nofigure', dofigure = false;
                
            case 'nosmooth', dosmooth = false;
            case 'smoothbox', mysmoothbox = varargin{i + 1}; varargin{i + 1} = [];
            case 'sd', mygaussstd = varargin{i + 1}; varargin{i + 1} = [];
               
            case 'xlim', xlim = varargin{i + 1}; varargin{i + 1} = [];
            case 'ylim', ylim = varargin{i + 1}; varargin{i + 1} = [];
            case 'zlim', zlim = varargin{i + 1}; varargin{i + 1} = [];
                
                
            case 'noverbose'
                doverbose = false;
                
                %otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end

dat3d = reconstruct_image(m);

if dosmooth
    dat3d = smooth3(dat3d, 'gaussian', mysmoothbox, sd);
    dat3d_for_exclude = dat3d;
else
    % just use to determine exclude area with zero values
    dat3d_for_exclude = smooth3(dat3d, 'gaussian', mysmoothbox, sd);
end

% get x and y coords in mm
% ---------------------------

XX = (1:size(dat3d, 1))'; % x
XX(:, 2) = 0; XX(:, 3) = 0;
XYZmm = voxel2mm(XX', m.volInfo.mat);
x = XYZmm(1, :)';

XX = (1:size(dat3d, 2))';  % y
XX = [zeros(length(XX), 1) XX]; XX(:, 3) = 0;
XYZmm = voxel2mm(XX', m.volInfo.mat);
y = XYZmm(2, :)';

% get maximum intensity projection (mip)
% ---------------------------
mip = max(dat3d, [], 3);
voldata = dat3d;

mip_for_exclude = max(dat3d_for_exclude, [], 3); 

% cut down to size - remove empty areas
whzero =  all(mip_for_exclude == 0, 1); % all x is 0

mip(:, whzero) = [];
mip_for_exclude(:, whzero) = [];

y(whzero) = [];
voldata(:, whzero, :) = [];

whzero =  all(mip_for_exclude == 0, 2); % y
mip(whzero, :) = [];
mip_for_exclude(whzero, :) = [];

x(whzero) = [];
voldata(whzero, :, :) = [];

whzero =  squeeze(all((all(voldata == 0)))); % z
voldata(:, :, whzero) = [];

% Exclude zero values
mip(mip_for_exclude == 0) = NaN;

mip = mip'; % put y on y-axis and x on x-axis

% smooth


% Plots
if dofigure
    
    create_figure('surf', 1, 2);
    
    han = imagesc(x, y, mip);
    ylabel('y'); xlabel('x');
    axis tight image
    colorbar
    
    set(han, 'AlphaDataMapping', 'scaled', 'AlphaData', abs(mip) .^ .5);
    
    
    subplot(1, 2, 2);
    han = surf(x, y, mip);
    axis tight vis3d
    ylabel('y'); xlabel('x');
    view(36, 46);
    camzoom(.8)
    
    set(han, 'AlphaDataMapping', 'scaled', 'AlphaData', abs(mip) .^ .5, 'FaceColor', 'interp', 'FaceAlpha', 'interp', 'EdgeColor', 'interp');
    
    %axis off
end


end
