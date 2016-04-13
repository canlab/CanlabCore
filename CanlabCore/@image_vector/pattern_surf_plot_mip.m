function [mip, x, y, voldata] = pattern_surf_plot_mip(m)
%
% axial maximum intensity projection pattern surface plot
% needs documentation
%
% voldata is not used, just an output for other processes, e.g., rigid-body
% transformation

dat3d = reconstruct_image(m);

dat3d = smooth3(dat3d, 'gaussian', 5, 2);

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

% cut down to size - remove empty areas
whzero =  all(mip == 0, 1); % all x is 0

mip(:, whzero) = [];
y(whzero) = [];
voldata(:, whzero, :) = [];

whzero =  all(mip == 0, 2); % y
mip(whzero, :) = [];
x(whzero) = [];
voldata(whzero, :, :) = [];

whzero =  squeeze(all((all(voldata == 0)))); % z
voldata(:, :, whzero) = [];


mip = mip'; % put y on y-axis and x on x-axis

% smooth


% Plots
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
