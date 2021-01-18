%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Provided by the Herzenberg Lab at Stanford University 
%   License: BSD 3 clause
%

function [H, edgeBinIdxs]=gridEdge(fgOrDensityObject, allBorders, ...
    clusterIds, color, ax, markerSize, marker, lineStyle, lineWidth)
H=0;
if nargin<7
    lineStyle='-';
    lineWidth=.5;
end
if isa(fgOrDensityObject, 'Density')
    density=fgOrDensityObject;
else
    if ischar(clusterIds)
        clusterIds=fgOrDensityObject.toClust(clusterIds);
    end
    density=fgOrDensityObject.density;
    if isempty(density.pointers)
        fgOrDensityObject.dbm(false, false);
        density=fgOrDensityObject.density;
    end
end
binIdxs=[];
for i=1:length(clusterIds)
    binIdxs=[binIdxs find(density.pointers==clusterIds(i))];
end

gce=edu.stanford.facs.swing.GridClusterEdge(density.M);
if ~allBorders
    gce.compute(binIdxs, nargout>1, density.mins, density.deltas)
    [xx, yy]=clockwise([gce.x gce.y]);
else
    gce.computeAll(binIdxs, density.mins, density.deltas)
    if strcmp('none', lineStyle)
        xx=gce.x;
        yy=gce.y;
    else
        [xx, yy]=clockwise([gce.x gce.y]);
    end
    
end
if nargout>1
    edgeBinIdxs=gce.edgeBins;
end
if ~isa(fgOrDensityObject, 'Density')
    if length(clusterIds)==1
        fgOrDensityObject.edgeX{clusterIds}=xx;
        fgOrDensityObject.edgeY{clusterIds}=yy;
    end
end
if nargin>3
    if isempty(color)
        color=[.6 .6 .6];
    end
    if nargin<7 || isempty(marker)
        marker='d';
    end
    if nargin<6 || isempty(markerSize)
        markerSize=3;
    end
    if nargin<4 || isempty(ax)
        H=plot(xx, yy, 'marker', marker, 'MarkerSize', markerSize, ...
            'Color', color, 'LineStyle', lineStyle, 'LineWidth', lineWidth);
    else
        H=plot(ax, xx, yy, 'marker', marker, 'MarkerSize', markerSize, ...
            'Color', color, 'LineStyle', lineStyle, 'LineWidth', lineWidth);
    end
end
end
