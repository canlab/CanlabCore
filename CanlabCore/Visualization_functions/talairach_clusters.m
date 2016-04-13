function cl = talairach_clusters(xyz,L,str,varargin)
% :Usage:
% ::
%
%    cl = talairach_clusters(xyz,L,str,[color, e.g., 'r'])
%
% Saves and/or visualizes contiguous clusters given coordinates from the Carmack Talairach atlas and a 
% textlabel (e.g., "Amygdala") for coordinates to retrieve 
%
% :Examples:
% ::
%
%    load talairach_info L5 x y z; xyz = [x y z]; 
%    cl = talairach_clusters(xyz,L5,'Amygdala','y'); save cl_amy cl
%
%    % OR
%
%    load 
%    cl = talairach_clusters(cl,[],[],'y');
%
% get clusters structure given an XYZ mm list, a list of labels L in a cell array, and a
% string to match labels to
%
% Complete example of making a limbic figure (some parts shown)
% ::
%
%    load talairach_info L5 L3 x y z; xyz = [x y z];
%    cl = talairach_clusters(xyz,L5,'Amygdala','y');
%    cl = talairach_clusters(xyz,L3,'Caudate','b');
%    cl = talairach_clusters(xyz,L5,'Putamen','g'); save cl_putamen cl
%    cl = talairach_clusters(xyz,L5,'Lateral Globus Pallidus','c'); save cl_glo1 cl
%    cl = talairach_clusters(xyz,L5,'Medial Globus Pallidus','c'); save cl_glo2 cl
%    axis off
%    set(gcf,'Color','w')
%    h = lightangle(0,30);
%    material dull
%    h2 = lightangle(0,-30);
%    h3 = lightangle(90,-30);
%
% ..
%    Tor Wager
% ..

doplot = 1;
color = 'r';
if length(varargin) > 0, color = varargin{1}; ,end


if isstruct(xyz)
    % already have clusters, do nothing
    cl = xyz;
else
    fprintf(1,'Finding xyz -> ');
    P = which('Carmack_single_subj_T1.img');
    P = which('single_subj_T1_gray.img');   % just needs 2 x 2 x 2 mm voxels -- or res desired

    wh = strmatch(str,L); 
    fprintf(1,'Found %3.0f coords -> ', size(wh,1));
    xyz = xyz(wh,:);
    
    if isempty(xyz), fprintf(1,'Exiting. \n');,return, end
    cl = xyz2clusters(xyz,P);
end

if doplot
    % plot
    for i = 1:length(cl)
        [out] = imageCluster('cluster',cl(i),'color',color);
    end
end

return
