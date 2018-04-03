function han = make_surface_figure(varargin)

% This is a function to create a figure with surface-brain plots from 
% gifti-files (.gii) as available from HCP. The returned handles are passed 
% to the function plot_surface_map, which overlays functional data on the
% brain surfaces. plot_surface_map calls this function, so you don't need
% to call it directly.
%
% You need the gifti-toolbox installed
% (http://www.artefact.tk/software/matlab/gifti/) and HCP S1200 release
% surface files on your MATLAB path 
% (e.g., S1200.L.inflated_MSMAll.32k_fs_LR.surf.gii, 
%        S1200.R.inflated_MSMAll.32k_fs_LR.surf.gii), unless you specify 
%  your own surface files (see 'surfacefile' input options below).
% 
%
%
% :Usage:
% ::
%     h = make_surface_figure(varargin);
%
% See plot_surface_map for plotting functional data.
%
% 
% 
% :Optional inputs:
%
%   **'surface'**
%       followed by the name of the HCP surface type to use. Will be 
%       ignored if you specify 'surfacefile' (see below).
%       options are: inflated [default], midthickness, pial, very_inflated
%
%   **'surfacefiles'**
%       followed by 2x1 cell array with the fullpath to your individually
%       selected surface files, e.g. {path_to_left_hem; path_to_right_hem}.
%       If specified, the 'surface' option above will be ignored. If not 
%       specified, the function looks for the HCP files specified with the 
%       'surface' option on the matlab path.
%
%   **'facecol'**
%       followed by a [1 x 3] RGB vector specifying the color of the brain 
%       surfaces. Default is [0.6 0.6 0.6].
% 
%   **'facealpha'**
%       followed by a scalar value for the datamap facealpha value
%       (default=1).
% 
%   **'bgcolor'**
%       followed by a [1 x 3] RGB vector setting the figure background color
%
%   **'newfig'**
%       will create a new figure instead of overwriting any existing
%       surface figure (default).
%
% 
% 
% ..
%     Author and copyright information:
%     -------------------------------------------------------------------------
%     Copyright (C) 2018 Stephan Geuter
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
% ..
% 

%
% ..
%    Programmers' notes:
%    List dates and changes here, and author of changes
%
%   1/15/2018 - created
%   Stephan Geuter, sgeuter@jhmi.edu
%
%   2/21/2018 - updated surface file options
%   Fred S Barrett & Stephan Geuter
%


newfig = 0;
surftype = 'inflated';
facecol  = [0.6 0.6 0.6];
bgcol    = [1 1 1];
txtcol   = [0 0 0];
facealpha = 1;

if numel(varargin)>0
    
    % loop input arguments
    for j=1:numel(varargin)
        
        switch varargin{j}
            case {'surface'}, surftype = varargin{j+1}; varargin{j+1} = '';
            % valid arguments: inflated [default], midthickness, pial, very_inflated
            
             case {'surfacefile','surfacefiles'}
                 fileL = varargin{j+1}{1};
                 fileR = varargin{j+1}{2}; 
                 varargin{j+1} = '';
            % if not provided, will default to S1200.[L/R].[surface]_MSMAll.32k_fs_LR.surf.gii
            
            case {'facecol'}, facecol = varargin{j+1}; varargin{j+1} = '';
            % color for surfaces. default = [0.6 0.6 0.6]
            
            case {'facealpha'}, facealpha = varargin{j+1}; varargin{j+1} = '';    
            % facealpha, default = 1
            
            case {'bgcolor'}, bgcol = varargin{j+1}; varargin{j+1} = '';
                            if mean(bgcol)<0.5, txtcol = [1 1 1]; end
            
            case {'newfig','newfigure'}, newfig = 1;
            % make a new figure to plot, default = 0
            
            otherwise
        end
    end
end
    
if exist('gifti','file')==0
    warning('''gifti'' function not found on path. Trying to add.');
    addpath(genpath('/Users/sgeuter/Documents/MATLAB/toolbox/gifti'));
end

% check for surface files
if ~exist('fileL','var'), fileL = which(sprintf('S1200.L.%s_MSMAll.32k_fs_LR.surf.gii',surftype)); end
if ~exist('fileR','var'), fileR = which(sprintf('S1200.R.%s_MSMAll.32k_fs_LR.surf.gii',surftype)); end

if isempty(fileL) || isempty(fileR)
    error('Surface files not found on MATLAB path. Check for: \n%s\n%s',...
        sprintf('S1200.L.%s_MSMAll.32k_fs_LR.surf.gii',surftype),...
        sprintf('S1200.R.%s_MSMAll.32k_fs_LR.surf.gii',surftype));
else
    % load surfaces
    hemL = gifti(fileL);
    hemR = gifti(fileR);
end


%% start drawing 

f = findobj('Type','figure','-and','Tag','SurfaceFig');
if newfig || isempty(f)
    f = figure('Tag','SurfaceFig');
end
figure(f(1)); clf; set(f,'color',bgcol,'DefaultTextColor',txtcol,'position',[1 600 780 350]);



%%% left lateral
ax(1)=axes('position',[.025 .31 .39 .66]);
p(1) = patch('Faces',hemL.faces,'Vertices',hemL.vertices,'FaceColor',facecol, ...
        'EdgeColor','none','SpecularStrength',.2,'FaceAlpha',facealpha,'SpecularExponent',200);

view(-90,0);
axis off; axis image; material dull;

h(1) = light;
lightangle(h(1),-90,5); % lateral light
lighting gouraud;
        
%%% left medial
ax(2)=axes('position',[.245 .02 .245 .4]);
p(2) = patch('Faces',hemL.faces,'Vertices',hemL.vertices,'FaceColor',facecol, ...
        'EdgeColor','none','SpecularStrength',.2,'FaceAlpha',facealpha,'SpecularExponent',200);

view(90,0);
axis off; axis image; material dull;

h(2) = light;
lightangle(h(2),90,5); % medial light
lighting gouraud;
        

%%% right lateral
ax(3)=axes('position',[.585 .31 .39 .66]);
p(3) = patch('Faces',hemR.faces,'Vertices',hemR.vertices,'FaceColor',facecol, ...
        'EdgeColor','none','SpecularStrength',.2,'FaceAlpha',facealpha,'SpecularExponent',200);

view(90,0);
axis off; axis image; material dull;

h(3) = light;
lightangle(h(3),90,5); % lateral light
lighting gouraud;


%%% right medial
ax(4)=axes('position',[.51 .02 .245 .4]);
p(4) = patch('Faces',hemR.faces,'Vertices',hemR.vertices,'FaceColor',facecol, ...
        'EdgeColor','none','SpecularStrength',.2,'FaceAlpha',facealpha,'SpecularExponent',200);

view(-90,0);
axis off; axis image; material dull;

h(4) = light;
lightangle(h(4),-90,5); % medial light
lighting gouraud;

%%%
drawnow;
figure(f);

% output
han.fig = f;
han.ax = ax;
han.obj = p;
han.light = h;
han.label = {'left lateral','left medial','right lateral','right medial'};
han.map = {};

end
