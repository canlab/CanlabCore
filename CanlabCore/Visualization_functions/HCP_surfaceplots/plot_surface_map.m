function h = plot_surface_map(dat,varargin)

% This is a method to plot cifti (grayordinate) data on cortical surfaces
% (e.g., from the HCP project). Read data with fieldtrip into matlab and
% call this function. Right now, only dense scalar files (*.dscalar.nii)
% with one map per file are supported.
% 
% :Usage:
% ::
%     d = ft_read_cifti('yourdata.dscalar.nii');
%     h = plot_surface_map(d, varargin);
%
% See make_surface_figure for the underlaying rendering of the brain surfaces.
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
% 
% :Inputs:
%
%   **dat:**
%        an cifti-data object from ft_read_cifti with fields 'dscalar'
%        holding the functional information and 'brainstructure' coding the
%        respective brainstructure for each functional datapoint in
%        'dscalar'. 'brainstructure' assumes HCP grayordinate conventions 
%        with 1=left hemisphere and 2=right hemisphere.
%
%
% :Optional inputs:
%
%   **'surface'**
%       followed by the name of the HCP surface type to use. options are: inflated [default], 
%       midthickness, pial, very_inflated
%
%   **'colmap'**
%       followed by a colormap matrix to use for the data (m x 3) or a 
%       string with name of a MATLAB build-in colormap (eg. 'parula'). If
%       you have the cbrewer2 package installed, you can call this, 
%       eg. plot_surface_map(dat,'colormap',cbrewer2('RdBu',64));
%
%   **'color'**
%       followed by a [1 x 3] RGB vector for a single color to use, eg.
%       when plotting a mask. if 'color' is entered, 'colmap' is ignored.
%       Colorbar plotting is also suppressed.
%
%   **'nocolorbar'**
%       do not plot a colorbar in the lower left corner
%
%   **'figure'**
%       followed by a structure with handles as returned by 
%       make_surface_figure(). Will plot into an existing figure with 
%       existing surface renderings provided in the struct.
% 
%   **'facealpha'**
%       followed by a scalar value for the datamap facealpha value
%       (default=0.8).
%
%   **'title'**
%       followed by a string with a figure title. Will be plotted in the
%       top center, between the lateral surfaces.
%
% 
% 
% 
% :Examples:
% ::
%
%  % Example 1: plot two binary masks on top of each other
%  % ------------------------------------------------------------------------- 
%  d1 = ft_read_cifti('mask1.dscalar.nii');
%  d2 = ft_read_cifti('mask2.dscalar.nii');
%  h = plot_surface_map(d1,'color',[1 .5 0],'facealpha',0.8);
%  h = plot_surface_map(d2,'figure',h,'title','mask 1 and mask 2','color',[0 .5 1],'facealpha',0.5);
%
%


%
% ..
%    Programmers' notes:
%    List dates and changes here, and author of changes
%
%   1/15/2018 - created
%   Stephan Geuter, sgeuter@jhmi.edu
%
%


%%% defaults %%%
surftype = 'inflated'; % surface type
docolorbar = 1; % colorbar
newfig = 1; % make a new figure
drange = []; % limits for colormap/colorbar
onecolor = []; % plot a mask in a single color
facealpha = 0.7; % face alpha value for data
figtitle = []; % title for figure to pring
colmap = [[113 220 247]; % default colormap
           [113 220 247];
           [0 150 255];
           [0 150 255];
           [0 99 255];
           [0 99 255];
           [0 50 255];
           [0 50 255];
%           [153 153 153];
           [255 50 0];
           [255 50 0];
           [255 105 0];
           [255 105 0];
           [255 205 0];
           [255 205 0];
           [255 255 100]
           [255 255 100]]/255;

       
% process input arguments
if numel(varargin)>0
    
    % loop input arguments
    for j=1:numel(varargin)
        
        switch lower(varargin{j})
            case {'surface'}, surftype = varargin{j+1}; varargin{j+1} = '';
            % valid arguments: inflated [default], midthickness, pial, very_inflated
            
            case {'colmap','colormap'}, colmap = varargin{j+1}; varargin{j+1} = '';
                             if ischar(colmap), colmap = eval(colmap); end
                             
            case {'nocolorbar','nolegend'}, docolorbar = 0;
             
            case {'figure','fighandle'}, newfig = 0; h = varargin{j+1}; varargin{j+1} = '';
                
            case {'cmaplim','cmaprange'}, drange = varargin{j+1}; varargin{j+1} = '';
                
            case {'color'}, onecolor = varargin{j+1}; varargin{j+1} = '';
             
            case {'facealpha'}, facealpha = varargin{j+1}; varargin{j+1} = '';
                
            case {'title'}, figtitle = varargin{j+1}; varargin{j+1} = '';
                
            case ''
                
            otherwise
                disp('unknown input option');
        end
    end
end



%%% get to work
% process input data
ldat = single(dat.dscalar(dat.brainstructure==1,:));
ldat(ldat==0) = NaN;
rdat = single(dat.dscalar(dat.brainstructure==2,:));
rdat(rdat==0) = NaN;

% no color map limits entered
if isempty(drange)
    % colormap
    if isempty(onecolor)
        drange = [nanmin([ldat;rdat]) nanmax([ldat;rdat])];
        drange = [-max(abs(drange)),max(abs(drange))];
        
    else
        % flat color
        drange = [0 nanmax(abs([ldat;rdat]))];
    end
end

if any(isnan(drange))
    drange = [0 1];
end

% ldat = rescale(ldat,drange(1),drange(2));
% rdat = rescale(rdat,drange(1),drange(2));
% % number of colors in colmap is odd, remove data for the mid/neutral color
% if mod(size(colmap,1),2)==1
%     drmv = size(colmap,1) * 2;
% else
%     drmv = 64;
% end
% ldat(ldat>(drange(1)/drmv) & ldat<(drange(2)/drmv)) = NaN;
% rdat(rdat>(drange(1)/drmv) & rdat<(drange(2)/drmv)) = NaN;
    
% make background surface figure
if newfig == 1
    h = make_surface_figure('surface',surftype);
end
h.n_maps = numel(h.map);


% set colormap
if isempty(onecolor)
    set(h.fig,'colormap',colmap);
else
    docolorbar = 0; % no colorbar for masks
end

% plot data on top of the surface
for j=1:length(h.obj)
    if strfind(h.label{j},'left')
        plotdat = ldat;
    elseif strfind(h.label{j},'right')
        plotdat = rdat;
    end

    % map onto colors
    if isempty(onecolor)
        % draw the map using a colormap
        set(h.fig,'colormap',colmap);
        h.map{h.n_maps+1}(j) = patch(h.ax(j),'Faces',h.obj(j).Faces,'Vertices',h.obj(j).Vertices,'FaceVertexCData',plotdat,...,
        'FaceColor','interp','EdgeColor','none','SpecularStrength',.2,'FaceAlpha',facealpha,'SpecularExponent',200);
    
    else
        % one color
        pidx = isfinite(plotdat);
        cDatCol = (pidx * onecolor);
        cDatCol(~pidx,:) = NaN;
        
        % draw the map
        h.map{h.n_maps+1}(j) = patch(h.ax(j),'Faces',h.obj(j).Faces,'Vertices',h.obj(j).Vertices,...
            'CData',plotdat,'FaceVertexCData',cDatCol,'FaceVertexAlphaData',pidx*facealpha,'FaceAlpha','flat',...
            'FaceColor','flat','EdgeColor','none','SpecularStrength',.2,'SpecularExponent',200);
    end
    
    material(h.map{h.n_maps+1}(j),'dull');
    lighting(h.ax(j),'gouraud');
    caxis(h.ax(j),double(drange));
end

h.n_maps = numel(h.map);


% add colobar
if docolorbar
   h.legendax = axes('position',[0.03 0.1 0.15 0.1],'visible','off','xlim',[0 1],'ylim',[0 1]);
   colormap(h.legendax,colmap);
   caxis(drange);
   h.colorbar=colorbar(h.legendax,'Location','South','AxisLocation','out','Ticks',[drange(1) 0 drange(2)],'Ticklabels',{num2str(drange(1),3),'0',num2str(drange(2),3)});
   h.colorbar.FontName = 'Helvetica Neue';
   drawnow;
end

% add title
if ischar(figtitle)
    % title exists, just replace text
    if isfield(h,'title')
        h.title.String = figtitle; drawnow;
    else
        h.titleax = axes('position',[0 0.9 1 0.1],'visible','off','xlim',[0 1],'ylim',[0 1]);
        h.title = text(0.5,0.5,figtitle,'FontName','Helvetica Neue','FontSize',12,'FontWeight','b','HorizontalALignment','center','VerticalAlignment','middle');
    end
end



end
