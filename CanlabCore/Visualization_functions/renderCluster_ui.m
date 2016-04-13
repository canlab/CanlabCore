function O = renderCluster_ui(varargin)
% :Usage:
% ::
%
%    O = renderCluster_ui([opt] O)
%
% This ui is set to render on the default single_subj_T1, in default colors
% More flexibility is available if you use the functions.
%
% main functions used:
%  - tor_3d.m - images head with cutaway views
%  - imageCluster.m - images a cluster isosurface
%  - mni_TSU.m and tor_ihb_TalSpace.m - to get clusters
%    these and related functions are part of Talairach Space Utility
%    written by Sergey Pakhomov, 2001
%    modified very slightly by Tor Wager to not convert to Talairach Space
%    and use MNI coordinates instead.
%
%   Use of the TSU functions for getting clusters require SPM99.
%
% :Output to workspace:
%
%   Isosurface handles are in p, for head isosurfaces, and cH, for the cluster isosurface
%
%   D = image data, Ds = smoothed data, hdr = img header
%
% Add more clusters by using:
% ::
%
%    cH(2) = imageCluster('cluster',clusters(i));
%
% :O: = option structure with fields
%
%   **dohead:**
%        y/n add head surface
%
%   **head:**
%        filename or 'default' for default canonical brain
%
%   **dobrain:**
%        y/n/filename add (transparent) brain surface, y for default brain or enter filename
%
%   **get_from:**
%        workspace/file/TSU/TSUfigure/none get clusters from here
%        if workspace, enter clusters in O.clusters
%
%   **which_cl:**
%        vector of clusters to image (from list)
%
%   **clcol:**
%        cluster colors - 3 el. vector, single letter, or string of letters
%        enter letter or letter string in single quotes. e.g., O.clcol = 'yrgb';
%
%   **whichc:**
%        letter string (no quotes) - which axes to cut along - xyzw are choices
%
%   **addtext:**
%        y/n add text to clusters
%
%   **textfield:**
%        field in cluster structure containing text
%
%   **textcol:**
%        character code (r, b, g, etc.) for color of text
%
%   **bestCoords:**
%        coordinates in mm to define x,y, and z cuts
%
%   **revx:**
%        text string to reverse x cut direction, enter 1 to do it.
% 
% for solid brain rendering without the scalp, where dohead gives you a
% brain image, use brain_render_T1.img
%
% see cluster_cutaways.m for an easy-to-use version.
%
% for transparent brain rendering of a set of clusters, try:
% ::
%
%    figure('Color','w');
%    O = struct('dohead','n','dobrain','y','get_from','workspace', ...
%                'clusters',clusters, ...
%                'which_cl',1:length(clusters),'whichc','y','bestCoords',[0 0 0],'clcol','y','addtext','n', ...
%                'head','single_subj_T1');
%    renderCluster_ui(O)
%
% ..
%    By Tor Wager, 10/3/2001, last edit 1/31/01
% ..

if nargin > 0, O = varargin{1};, else, O.dummy = [];, end

N = fieldnames(O);
for n = 1:length(N)
    eval(['fval = O.' N{n} ';']) 
    if isempty(fval),O = rmfield(O,N{n});,end
end

p = [];

%------------------------------------------------------------------------------
% Get head, if not default
%------------------------------------------------------------------------------
if ~(isfield(O,'dohead')),O.dohead = input('Add head surface? (y/n) ','s');,end
if strcmp(O.dohead,'y'), dohead = 1;, elseif strcmp(O.dohead,'n'), dohead = 0;, else error('choose y or n'), end

fileName = 'avg152T1';

if dohead
    cflag = 0;
    while cflag == 0
        if ~(isfield(O,'head')),
            O.head = input('Use default head avg152T1...?  (type default or file name, no .img extension)','s');
        end
        switch O.head
        case 'default', fileName = 'avg152T1';
        otherwise fileName = O.head;
        end
        if exist([fileName '.img']) == 2, 
	    cflag = 1;
        else
	        disp(['Can''t locate file: ' fileName '.img'])
            pause(3)
        end
    end
end

%------------------------------------------------------------------------------
% Add surface?
%------------------------------------------------------------------------------
cflag = 0;
while cflag == 0
    if ~(isfield(O,'dobrain')),O.dobrain = input('Add brain surface (y\n\filename,no .img extension)...?','s');,end
    switch O.dobrain
    case 'y', fileName2 = ['surf_' fileName];, dosurf = 1;
    case 'n', fileName2 = [];, cflag = 1;, dosurf = 0;
    otherwise fileName2 = O.dobrain;, dosurf = 1;
    end
    if dosurf & exist([fileName2 '.mat']) == 2, 
	    cflag = 1;
    elseif dosurf
	    disp(['Can''t locate file: ' fileName2 '.mat'])
    end
end


%------------------------------------------------------------------------------
% Get clusters, if necessary
%------------------------------------------------------------------------------
cflag = 0;
while cflag == 0
    if ~(isfield(O,'get_from')),O.get_from = input('Get clusters from...?  (workspace,file,TSU,figure,none)','s');,end
    switch O.get_from
    case 'workspace', cflag = 1;
        clusters = O.clusters;
    case 'file', [filename,pathname] = uigetfile('*.mat','pick a mat file containing ''clusters'' variable.');
        try, load(fullfile(pathname,filename)), cflag = 1;, catch, cflag = 0; disp('Problem with file'), end 
        if ~(exist('clusters')) == 1, disp('No ''clusters'' var found in file.'), cflag = 0;, end
            
    case 'TSU', clusters = imageCluster('getclusters');, cflag = 1;
    case 'figure', clusters = imageCluster('getfigclusters');, cflag = 1;
    case 'none',clusters = [];,cflag = 1;
    otherwise disp(['Not a valid choice...type entry string with no quotes.'])
    end
end

omit = find(cat(1,clusters.numVox) < 2);
if ~isempty(omit)
    clusters(omit) = [];
    disp('OMITTED clusters with 1 or 0 voxels:')
    omit
end

for i = 1:length(clusters), clusters(i).cluster_num = i;, end

%------------------------------------------------------------------------------
% Exclude clusters that have only one voxel; cannot image
%------------------------------------------------------------------------------
%if ~isempty(clusters), omit = cat(1,clusters.numVox); omit = omit == 1; clusters(omit) = [];, end

%------------------------------------------------------------------------------
% Get number of cluster to image 
%------------------------------------------------------------------------------
if ~isempty(clusters)
    cluster_names = str2mat(clusters.name);
    clNums = num2str((1:length(clusters))');
    blankSpcs = repmat(' ',length(clusters),1);
    cluster_names_wNums = [clNums blankSpcs cluster_names]
    clear clNums,clear blankSpcs

    cflag =0;
    while cflag == 0
        if ~(isfield(O,'which_cl')),O.which_cl = input('Enter number of cluster to image, or vector of clusters, ex. [2 4] : ');,end
        if any(O.which_cl > length(clusters)), disp(['Invalid: no cluster with number ' num2str(O.which_cl(O.which_cl>length(clusters))) '.']),
            O = rmfield(O,'which_cl');
        else cflag = 1;     % ok to go ahead
        end
    end
else
    clNums = [];
    O.which_cl = [];
end

%------------------------------------------------------------------------------
% Get whether to add text, and which field has the text
%------------------------------------------------------------------------------
if ~(isfield(O,'addtext')),O.addtext = input('Add text to image? (y/n) ','s');, end
cflag = 0;
while cflag == 0
    switch O.addtext
    case 'y',
        addtext = 1; cflag = 1;
        if ~(isfield(O,'textfield')),O.textfield = input('Type name of field containing text or numbers (cluster_num for cluster number) ','s');, end
        if ~(isfield(O,'textcol')),O.textcol = input('Enter character code for text color no quotes ','s');, end
    case 'n',
        addtext = 0; cflag = 1;
    otherwise
        disp('Not a valid choice - choose y or n')
    end
end


%------------------------------------------------------------------------------
% If clusters is missing voxSize (voxel size), then prompt (should be col vector)
%------------------------------------------------------------------------------
if ~(isfield(clusters,'voxSize'))
    if ~(isfield(O,'voxSize'))
        voxSize = input('Enter x y z voxel sizes in mm: ')';
    else
        voxSize = O.voxSize;
    end
    for i = 1:length(clusters), clusters(i).voxSize = voxSize;,end
end

%------------------------------------------------------------------------------
% Define cluster color 
%------------------------------------------------------------------------------
cflag =0;
while cflag == 0
    disp('Colors: enter vector, single letter (all clusters) or string of letters for different clusters')
    if ~(isfield(O,'clcol')),O.clcol = input('Enter color(s) of these clusters, strings in quotes : ');, end
    cflag = 1;
end

for i = 1:length(O.which_cl), 
    if ischar(O.clcol) & length(O.clcol) > 1,
        myColors{i} = O.clcol(i);
    else
        myColors{i} = O.clcol;
    end
end


%------------------------------------------------------------------------------
% Define cluster and cluster center (-5mm for easy viewing)
%------------------------------------------------------------------------------
clOrig = clusters(O.which_cl);
cl = [];

for i = 1:length(O.which_cl)

    %------------------------------------------------------------------------------
    % get xMin xMax and vTal if not already defined (if clusters made with RoiUtility)
    %------------------------------------------------------------------------------
    if ~(isfield(clOrig(i),'xMin')) & clOrig(i).numVox > 1
	    mycl = get_cluster_volume(clOrig(i));
    elseif ~(isfield(clOrig(i),'xMin')) & clOrig(i).numVox == 1
        mycl = clOrig(i);
        disp(['Note: Cluster ' num2str(i) ' has one voxel.'])
        mycl.xMin = clOrig(i).XYZmm(1); mycl.xMax = mycl.xMin;
        mycl.yMin = clOrig(i).XYZmm(2); mycl.yMax = mycl.yMin;
        mycl.zMin = clOrig(i).XYZmm(3); mycl.zMax = mycl.zMin;
        
    else
	    mycl = clOrig(i);
    end
    
    if i == 1,
	    cl = mycl;
    else
        cl(i) = mycl;
    end 

        coords(i,1) = mean([cl(i).xMin cl(i).xMax]) - 5;
        coords(i,2) = mean([cl(i).yMin cl(i).yMax]) - 5;
        coords(i,3) = mean([cl(i).zMin cl(i).zMax]) - 5;
end



if isempty(coords), 
    if isfield(O,'bestCoords')
        bestCoords = O.bestCoords;
    else
        bestCoords = input('Enter [x y z] coords to define cuts. ');, 
    end
    coords = bestCoords;
else
    if size(coords,1) > 1, 
        bestCoords = min(coords);, 
        %tmp = abs(coords(coords(:,1) <= 0,1)); tmp3 = find(tmp == min(tmp)); % closest to mid on L
        %tmp2 = abs(coords(coords(:,1) > 0,1)); tmp4 = find(tmp2 == min(tmp)); % closest on R
        %tmp5 = min([tmp(tmp==max(tmp)) tmp2(tmp2==max(tmp2))]);           % whichever is closer to midline
        %tmp6 = find(tmp5 == min(tmp5)); tmp6 = tmp6(1);
        %tmp7 = [tmp3 tmp4]; tmp7 = tmp7(tmp6);
        
        %bestCoords(1) = coords(tmp7,1);    % replace x coordinate with closest abs to midline.

    else, bestCoords = coords;, 
    end
    
    %bestCoords = [0 0 0];
end

if isfield(O,'bestCoords')
    bestCoords = O.bestCoords;
end
        
disp(['Coordinates to use are: ' num2str(bestCoords)])

%------------------------------------------------------------------------------
% Get views to cut away
%------------------------------------------------------------------------------
if ~(isfield(O,'whichc')),O.whichc = input('Enter axes to cut along.  Valid choices are xyzw. Example yz: ','s');,end

%------------------------------------------------------------------------------
% Image the clusters (first) and add the head surface if specified
%------------------------------------------------------------------------------
% set clusters first or make3davi won't work. ??? or maybe it's the figure size.
cH = [];
for i = 1:length(cl)
    if cl(i).numVox > 1,
        cH(i) = imageCluster('cluster',cl(i),'color',myColors{i},'alpha',.7); drawnow
        set(cH(i),'Tag',['cluster' num2str(i)])
    else
        lineH(i) = plot3(cl(1).XYZmm(1,:),cl(1).XYZmm(2,:),cl(1).XYZmm(3,:),[myColors{i}(1) 'o'],'MarkerFaceColor',myColors{i},'MarkerSize',10,'LineWidth',2);
        set(lineH(i),'Tag',['cluster' num2str(i)])    
    end
    
    if addtext
        eval(['mytext = cl(i).' O.textfield])
        if ~ischar(mytext), mytext = round(mytext * 100) / 100;, mytext = num2str(mytext);,end
        textpos = mean(cl(i).XYZmm,2);
        textp(i) = text(textpos(1),textpos(2),textpos(3)+5,mytext,'Color',O.textcol,'FontSize',14);
    end
    set(cH,'SpecularStrength',.4);
end


if ~isfield(O,'revx'), O.revx = 0;, end

if dohead
    if O.revx,
        disp('Reversing x-cut direction')
        [D,Ds,hdr,O.headp,bestCoords] = tor_3d('whichcuts',O.whichc,'coords',bestCoords,'filename',fileName,'revx');
    else
        [D,Ds,hdr,O.headp,bestCoords] = tor_3d('whichcuts',O.whichc,'coords',bestCoords,'filename',fileName);
    end
end
    
        
end
colormap(gray(100))

O.H = cH;
if addtext, O.textp = textp;,end

if dosurf
    eval(['load ' fileName2])
    try
        p(end+1) = patch('Faces',faces,'Vertices',vertices,'FaceColor',[.5 .5 .5],'EdgeColor','none','SpecularStrength',.2,'FaceAlpha',.4,'SpecularExponent',200);
    catch   % transparency doesn't work?  try without
        p(end+1) = patch('Faces',faces,'Vertices',vertices,'FaceColor',[.5 .5 .5],'EdgeColor','none','SpecularStrength',.2,'SpecularExponent',200);
    end
    O.p = p;
    if ~dohead,myLight = standardMRIlighting('full',[p(end) p(end)]);,end
end  


saveas(gcf,'temp','fig')
save temp_handles p cH

%set(p(1),'SpecularExponent',200)
%set(p(3),'SpecularExponent',200)
%set(p(5),'SpecularExponent',200)
axis off
set(gcf,'Color','k')
%set(p(1),'FaceColor',[.8 .5 .4]);
%set(p(3),'FaceColor',[.8 .5 .4]);
%set(p(5),'FaceColor',[.8 .5 .4]);
%set(p(1),'SpecularStrength',.2)
%set(p(3),'SpecularStrength',.2)
%set(p(5),'SpecularStrength',.2)

% make3davi(O)
return
