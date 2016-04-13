% function hewma_timeseries_plot
%
% graphic display of significant voxels in whole-brain hewma analysis
% and plotting of timeseries points you click on in the image.
%
% ..
%    tor wager
% ..
if ~exist('overlayimg') == 1, overlayimg = 'hewma_sig.img';,end
disp(['Overlay image is: ' overlayimg])

% used in button-up fcn callback
E2 = EXPT;
clear EXPT

global VOL
global f
global f2
global EXPT
EXPT = E2;

    % -------------------------------------------------------------------
    % * essential stuff for the viewing
    % -------------------------------------------------------------------
    
    cl = mask2clusters(overlayimg);
    
% view clusters
cluster_orthviews(cl,'unique');
set(gcf,'WindowButtonUpFcn','[dat,files,stats,mycov] = hewma_plot_coord_btnupfcn;')


% get coordinate mapping matrix
VOL = struct('M',cl(1).M);

% prepare figure
f1 = figure('Color','w','Name','Hewma plots');
f = f1;     % button callback uses figure f


cd ..

%return






