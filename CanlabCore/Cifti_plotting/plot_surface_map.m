function h = plot_surface_map(dat,varargin)

% This is a method to plot cifti (grayordinate) data on cortical surfaces
% (e.g., from the HCP project). Read data with fieldtrip into matlab and
% call this function. Right now, only dense scalar files (*.dscalar.nii)
% with one map per file are supported.
% 
% :Usage:
% ::
%     d = ft_read_cifti('yourdata.dscalar.nii');
%               -OR-
%     d = 'path_to_your_data.dscalar.nii';
%               -OR-
%     d = {lh_data,rh_data};
%
%
%     h = plot_surface_map(d, varargin);
%
%
% See make_surface_figure for the underlaying rendering of the brain surfaces.
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
%           -OR-
%
%        a string with the full path to a cifti-file (*.dscalar.nii)
%
%           -OR-
%
%        a 1x2 cell array with the functional data to plot for left and
%        right hemispheres in each cell, respectively. Here, dat{1} 
%        contains functional information for the left hemisphere, and 
%        dat{2} contains functional information for the right hemisphere. 
%        Data for each hemispere are v x 1 vectors with one
%        datapoint for each surface vertex. Default for HCP standard space
%        is 32,492 vertices per hemisphere. If you provide individual
%        hemisphere surface files via the 'surfacefiles' option, the number
%        of vertices may differ. In this case, the number of data points in
%        data{1:2} must match the vertices of the 'surfacefiles'
%
%
% :Optional inputs:
%
%   **'surface'**
%       followed by the name of the HCP surface type to use. options are: 
%       inflated [default], midthickness, pial, very_inflated
%       This option will be ignored if you specify 'surfacefiles'
%       explicitly.
%
%   **'surfacefiles'**
%       followed by 2x1 cell array with the fullpath to your individually
%       selected surface files, e.g. {path_to_left_hem; path_to_right_hem}.
%       If specified, the 'surface' option above will be ignored. If not 
%       specified, the function looks for the HCP files specified with the 
%       'surface' option on the matlab path.
%
%   **'wh_map'**
%       followed by a scalar value with the map index to plot for
%       dscalar.nii-files with multiple maps. default=1.
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
%   **outline**
%       will only draw the borders of clusters or ROIs. Recommended for
%       marking ROIs from an atlas with an .dlabel-cifti. E.g., create a
%       cifti with a few regions of interest, each coded by a different
%       integer value and plot using the 'outline' option.
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
% % Example 1: plot the first effect size map from HCP S1200
% % ------------------------------------------------------------------------- 
% d = ft_read_cifti(which('HCP_S1200_997_tfMRI_ALLTASKS_level2_cohensd_hp200_s4_MSMSulc.dscalar.nii'),'mapname','array');
% h = plot_surface_map(d,'facealpha',0.8);
%
%
% % Example 2: plot two binary masks on top of each other
% % ------------------------------------------------------------------------- 
% d1 = ft_read_cifti('mask1.dscalar.nii');
% d2 = ft_read_cifti('mask2.dscalar.nii');
% h = plot_surface_map(d1,'color',[1 .5 0],'facealpha',0.8);
% h = plot_surface_map(d2,'figure',h,'title','mask 1 and mask 2','color',[0 .5 1],'facealpha',0.5);
%
%
% % Example 3: plot the Glasser et al (2016) Nature atlas parcellation
% % ------------------------------------------------------------------------- 
% glassermap = which('Q1-Q6_RelatedParcellation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR.dlabel');
% h = plot_surface_map(glassermap);
%
%
% Example 4: plot only the borders of the Glasser atlas regions on pial surface 
% ------------------------------------------------------------------------- 
% glassermap = which('Q1-Q6_RelatedParcellation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR.dlabel');
% h = plot_surface_map(glassermap,'surface','pial','outline');
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
%   2/21/2018 - updated surface file options, added dlabel.nii option
%   Fred S Barrett & Stephan Geuter
%
%   3/30/2018 - added outline and background color options
%   Stephan Geuter
%


%%% defaults %%%
surftype = 'inflated'; % surface type
docolorbar = 1; % colorbar
newfig = 1; % make a new figure
drange = []; % limits for colormap/colorbar
onecolor = []; % plot a mask in a single color
facealpha = 0.7; % face alpha value for data
bgcol = [1 1 1]; % background color
txtcol= [0 0 0]; % text color
figtitle = []; % title for figure to pring
dooutline = 0; % only plot the outlines of clusters 
wh_map = 1;
datafield = 'dscalar'; % default fieldname of dscalar-structure
colmap = [[113 220 247]; % default colormap
           [113 220 247];
           [0 150 255];
           [0 150 255];
           [0 99 255];
           [0 99 255];
           [0 50 255];
           [0 50 255];
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
            
            case {'surfacefile','surfacefiles'}
                 fileL = varargin{j+1}{1};
                 fileR = varargin{j+1}{2}; 
                 varargin{j+1} = '';
            % if not provided, will default to S1200.[L/R].[surface]_MSMAll.32k_fs_LR.surf.gii
            
            case {'wh_map','wh'}, wh_map = varargin{j+1}; varargin{j+1} = '';
                
            case {'outline'}, dooutline = 1;    
            
            case {'colmap','colormap'}, colmap = varargin{j+1}; varargin{j+1} = '';
                             if ischar(colmap), colmap = eval(colmap); end
                             
            case {'nocolorbar','nolegend'}, docolorbar = 0;
             
            case {'figure','fighandle'}, newfig = 0; h = varargin{j+1}; varargin{j+1} = '';
                
            case {'cmaplim','cmaprange'}, drange = varargin{j+1}; varargin{j+1} = '';
                
            case {'color'}, onecolor = varargin{j+1}; varargin{j+1} = '';
             
            case {'facealpha'}, facealpha = varargin{j+1}; varargin{j+1} = '';
            
            case {'bgcolor'}, bgcol = varargin{j+1}; varargin{j+1} = '';
                            if mean(bgcol)<0.5, txtcol = [1 1 1]; end
                
            case {'title'}, figtitle = varargin{j+1}; varargin{j+1} = '';
                
            case ''
                
            otherwise
                disp('unknown input option');
        end
    end
end

% get surface files
if ~exist('fileL','var'), fileL = which(sprintf('S1200.L.%s_MSMAll.32k_fs_LR.surf.gii',surftype)); end
if ~exist('fileR','var'), fileR = which(sprintf('S1200.R.%s_MSMAll.32k_fs_LR.surf.gii',surftype)); end
surffiles = {fileL; fileR};


%%% get to work
% process input data
if iscell(dat) % data in cell array
  ldat = dat{1};
  rdat = dat{2};
  
else % data in cifti-struct or string with filepath
    
    if ischar(dat)  % read from file
    
        if contains(dat,'dscalar.nii')
           
            % dscalar file
            dat = ft_read_cifti(dat,'mapname','array');
            
        elseif contains(dat,'dlabel.nii')
            
            % dlabel file
            dat = read_dlabel_cifti(dat);
        end
    end
    
    % dcheck if is dscalar or dlabel
    if isstruct(dat) && ~isfield(dat,'dscalar')
        
        % dlabel data
        datafield = 'indexmax';
        
        % use the colors specified in the dlabel file, no colorbar here
        colmap = dat.color;
        docolorbar = 0;
    end
    
    % separate left and right hem data
    ldat = single(dat.(datafield)(dat.brainstructure==1,wh_map));
    ldat(ldat==0) = NaN;
    rdat = single(dat.(datafield)(dat.brainstructure==2,wh_map));
    rdat(rdat==0) = NaN;
end





% colormap specs and colorbar info
% no color map limits entered
if isempty(drange)
    % colormap
    if isempty(onecolor)
        drange = [nanmin([ldat;rdat]) nanmax([ldat;rdat])];
        if min(drange)<0
            drange = [-max(abs(drange)),max(abs(drange))];
        end
            
    else
        % flat color
        drange = [0 nanmax(abs([ldat;rdat]))];
    end
end

if any(isnan(drange))
    drange = [0 1];
end

   
% make background surface figure
if newfig == 1
    h = make_surface_figure('surfacefiles',surffiles,'bgcolor',bgcol);
end
h.n_maps = numel(h.map);


% set colormap
if isempty(onecolor)
    set(h.fig,'colormap',colmap);
else
    docolorbar = 0; % no colorbar for masks
end
set(h.fig,'color',bgcol);


% convert to outlines
if dooutline
    if any(onecolor)
        outlinecol = onecolor;
    else
        outlinecol = colmap;
    end
    [outlineF, outlineCData] = surface_outlines({ldat,rdat},{h.obj(1),h.obj(2)},outlinecol);
    ldat = outlineCData{1};
    rdat = outlineCData{2};
end



% plot data on top of the brain surfaces
for j=1:length(h.obj)
    
    % get hemisphere data
    if any(strfind(h.label{j},'left'))
        plotdat = ldat;
    elseif any(strfind(h.label{j},'right'))
        plotdat = rdat;
    end
    
    % get faces from hemisphere patch or outline faces
    if dooutline && any(strfind(h.label{j},'left'))
        F = outlineF{1};
    elseif dooutline && any(strfind(h.label{j},'right'))
        F = outlineF{2};
    else
       F = h.obj(j).Faces;
    end
    
    if dooutline || any(onecolor)
        fcolormode = 'flat';
    else
        fcolormode = 'interp';
    end
 
    
    % map onto colors
    if isempty(onecolor)
        % draw the map using a colormap
        set(h.fig,'colormap',colmap);
        h.map{h.n_maps+1}(j) = patch(h.ax(j),'Faces',F,'Vertices',h.obj(j).Vertices,'FaceVertexCData',plotdat,...,
        'FaceColor',fcolormode,'EdgeColor','none','SpecularStrength',.2,'FaceAlpha',facealpha,'SpecularExponent',200);
    
    else
        % one color
        pidx = isfinite(plotdat(:,1));
        cDatCol = (pidx * onecolor);
        cDatCol(~pidx,:) = NaN;
        
        % draw the map
        h.map{h.n_maps+1}(j) = patch(h.ax(j),'Faces',F,'Vertices',h.obj(j).Vertices,...
            'CData',plotdat,'FaceVertexCData',cDatCol,'FaceVertexAlphaData',pidx*facealpha,'FaceAlpha','flat',...
            'FaceColor',fcolormode,'EdgeColor','none','SpecularStrength',.2,'SpecularExponent',200);
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
   % colorbar ticks
   if min(drange)<0
       cticks=[drange(1) 0 drange(2)];
       cticklabel = {num2str(drange(1),3),'0',num2str(drange(2),3)};
   elseif size(colmap,1)>6
       cticks = linspace(drange(1),drange(2),4);
       cticklabel = split(num2str(cticks))';
   else
       cticks=linspace(drange(1),drange(2),size(colmap,1));
       cticklabel = split(num2str(cticks))';
   end
   
   h.colorbar=colorbar(h.legendax,'Location','South','AxisLocation','out',...
        'Ticks',cticks,'Ticklabels',cticklabel);
   h.colorbar.FontName = 'Helvetica Neue';
   h.colorbar.Color = txtcol;
   drawnow;
end


% add title
if ischar(figtitle)
    % title exists, just replace text
    if isfield(h,'title')
        h.title.String = figtitle;
    else
        h.titleax = axes('position',[0 0.9 1 0.1],'visible','off','xlim',[0 1],'ylim',[0 1]);
        h.title = text(0.5,0.5,figtitle,'FontName','Helvetica Neue','FontSize',12,'FontWeight','b','HorizontalALignment','center','VerticalAlignment','middle');
    end
    h.title.Color = txtcol;
    drawnow;
end



end
