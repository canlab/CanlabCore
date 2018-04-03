function [newF,cData] = surface_outlines(plotdat,patchobj,colmap)

% helper function to compute the borders of ROIs on surface data
%
% input:
%
% plotdat   - cell array (1x2) with the data to plot for a single 
%             hemisphere in each cell. 
%             plotdat{1} = left hemisphere, plotdat{2} = right hemisphere
% 
% patchobj  - cell array (1x2) containing the patch objects for left {1}
%             and right {2} hemispheres. Only need the Face and Vertex info
%             here, so could be structs in the cells. 
%             e.g. patchobj{hem}.Faces = yourFaces;  
%             patchobj{hem}.Vertices = yourVertices;
% 
% colmap    - colormap to color the ROI outlines with. Can also be a single
%             color for all borders (e.g. [0.2 0.2 0.2]).
% 
% output:
% 
% newF      - the new Faces containing only the outline faces around each ROI
% 
% cData     - the color information for the new outline faces. use with 
%             patch as patch(...,'FaceVertexCdata',cData) 
%                 
%         
%  Examples:
%  % make surface figure
%  h = make_surface_figure('surfacefiles',surffiles);
%  [outlineF, outlineCData] = surface_outlines({ldat,rdat},{h.obj(1),h.obj(2)}, [0 0 0]);
%  % left lateral hem
%  patch(h.ax(1),'Faces',outlineF{1},'Vertices',h.obj(1).Vertices,'FaceVertexCData',outlineCData{1},...,
%         'FaceColor',fcolormode,'EdgeColor','none','SpecularStrength',.2,'FaceAlpha',facealpha,'SpecularExponent',200);
%  % right lateral hem 
%  patch(h.ax(3),'Faces',outlineF{2},'Vertices',h.obj(3).Vertices,'FaceVertexCData',outlineCData{2},...,
%         'FaceColor',fcolormode,'EdgeColor','none','SpecularStrength',.2,'FaceAlpha',facealpha,'SpecularExponent',200);
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
%   3/30/2018 - created
%   Stephan Geuter, sgeuter@jhmi.edu
%



% regionIDs = [199 214 280 360  20 100 170];
% col = lines(numel(regionIDs));
% 



% loop hemispheres
for hem=1:2
    
    % get data for this hemisphere
    atlasDat = plotdat{hem};
    atlasDat(isnan(atlasDat)) = 0;
    
    % get Faces, Vertices of the brain surface
    F = patchobj{hem}.Faces;
    V = patchobj{hem}.Vertices;
    
    % unique ROI/blob ID's
    hemRegions = unique(atlasDat(isfinite(atlasDat)));
    
    newF{hem}  = double.empty(0,3);
    cData{hem} = double.empty(0,3);
    
    % loop ROIs
    for r=1:numel(hemRegions)
        
        % find the vertices that belong to the current ROI
        idx = atlasDat .* single(ismember(atlasDat,hemRegions(r)));
        
        if numel(unique(idx))>1 && any(unique(idx)~=0)
            % compute the new Faces
            Fclasses = idx(F);
            roiFace  = F(any(diff(Fclasses,1,2),2) , :);
            
            % get color for current ROI outline
            if size(colmap,1) == 1
                roiCol   = repmat(colmap,size(roiFace,1),1);
            else
                roiCol   = repmat(colmap(hemRegions(r),:),size(roiFace,1),1);
            end
            
            % cat data across ROIs
            newF{hem} = vertcat(newF{hem},roiFace);
            cData{hem}= vertcat(cData{hem},roiCol);
        end
    end
    
    
%     for j=find(cellregexp(h.label,hemstr))
%         hp(j) = patch(h.ax(j),'Faces',newF{hem},'Vertices',V,'FaceColor','flat',...
%             'Edgecolor','none','FaceVertexCdata',cData{hem},...
%             'SpecularStrength',.2,'SpecularExponent',200);
%         material(hp(j),'dull');
%         lighting(h.ax(j),'gouraud');
%     end
end





