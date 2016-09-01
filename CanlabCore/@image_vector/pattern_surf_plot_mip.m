function [mip, x, y, voldata] = pattern_surf_plot_mip(m, varargin)
%
% axial maximum intensity projection pattern surface plot
% needs documentation
%
% voldata is not used, just an output for other processes, e.g., rigid-body
% transformation
%
% m : image_vector (e.g., fmri_data or statistic_image) object
%
% Optional inputs:
% 'nofigure', dofigure = false;
% 'nosmooth', dosmooth = false;
% 'smoothbox', mysmoothbox = varargin{i + 1}; varargin{i + 1} = [];
% 'sd', mygaussstd = varargin{i + 1}; varargin{i + 1} = [];
%
% Tor Wager, 2016
% Update: Aug 2016

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
