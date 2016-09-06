function [hPatch,outP,FV, cl, myLight] = mask2surface(P,varargin)
% ::
%
%    [hPatch,outP,FV, cl, myLight] = mask2surface(P)
%
% :Input:
%
%   **P:**
%        a file name for a mask file, OR a clusters structure
%
% :Inputs:
%
%   **hPatch:**
%        handle to surface patch
%
%   **outP:**
%        file name, with path, of mat file containing faces and vertices
%
%   **FV:**
%        isosurface
%
%   **cl:**
%        clusters structure with coords and meshgrid info
%
%   **varargin:**
%        suppress lighting (0) or do lighting (1)
%
%   **varargin{2}:**
%        color for patch objects
%
% ..
%    tor wager 9/23/02
% ..
%
% uses get_cluster_volume, adapted from Sergey Pakhomov
%
% use with cluster_surf.m to map activations onto surfaces:
% cluster_surf(clusters1,clusters2,outP,10,{[0 1 0] [1 0 0]})
%
% or use getvertexcolors.m to map colors onto surface using hPatch.

hPatch = []; outP = [];  cl = []; myLight = [];
mycolor = [.5 .5 .5];
myalpha = 1;

P = which(P);
if isempty(P), error('Cannot find image.');, end

if isstr(P), cl = mask2clusters(P);, else, cl = P;, P=cl(1).P;, end

height_threshold = .9;
dolight = 1;

if length(varargin) > 0, dolight = varargin{1};, end
if length(varargin) > 1, mycolor = varargin{2};, end

set(0,'CurrentFigure', gcf)

omit = cat(1,cl.numVox);
cl(find(omit==1)) = [];

for i = 1:length(cl)
    
    % added for SPM2 compatibility
    cl(i).voxSize = abs(cl(i).voxSize);
    
    cl2(i) = get_cluster_volume(cl(i));
    V = smooth3(cl2(i).vTal,'gaussian',[5 5 5]);
    
	FV(i) = isosurface(cl2(i).xTal, cl2(i).yTal, cl2(i).zTal, V, height_threshold);
	hPatch(i) = patch(FV);
	isonormals(cl2(i).xTal, cl2(i).yTal, cl2(i).zTal, cl2(i).vTal, hPatch(i));
    
    try
		set(hPatch, 'Tag', P, 'FaceColor', mycolor, 'EdgeColor', 'none', 'FaceAlpha', myalpha);
	catch
		set(hPatch, 'Tag', P, 'FaceColor', mycolor, 'EdgeColor', 'none');
	end
   
end

cl = cl2;  clear cl2

if dolight
    view(135,30);    
    myLight = standardMRIlighting('full',[hPatch hPatch]);
    material dull
end

[d f e] = fileparts(P);
outP = fullfile(d, filesep,['surf_' f '.mat']);
%if outP(1) == filesep, outP = outP(2:end);end
save(outP, 'FV');

return
