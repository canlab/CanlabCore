function O = renderCluster_ui4(O)
% :Usage:
% ::
%
%    O = renderCluster_ui4(O)
%
% :Inputs:
%
%   **O**
%
%   **O.head:**
%        name of head to use, no .img extension
%
%   **O.surf:**
%        'y' or 'n': image surface
%
%   **O.sets:**
%        number of sets
%
%   **O.getfrom:**
%        cell array of where to get clusters from: 'workspace' 'file' etc.
%
%   **O.varname:**
%        cell array of variable names
%
%   **O.color:**
%        cell array of colors - strings or vectors
%
%   **O.numbers:**
%        'y' or 'n': add numbers to plot
%
%   also need individual fields named contents of varname, which contain clusters

if ~isfield(O,'clusterAlpha'),O.clusterAlpha = .7;, end
if ~isfield(O,'coords'),O.coords = [30 30 -5];, end
if ~isfield(O,'whichcuts'),O.whichcuts = 'z';, end
if ~isfield(O,'dohead'),O.dohead = 1;, end

%------------------------------------------------------------------------------
% Get head, if not default
%------------------------------------------------------------------------------
cflag = 0;
while cflag == 0
    if ~isfield(O,'head')
        a = input('Use default head single_subj_T1...?  (type default or file name, no .img extension)','s');
    else a = O.head;
    end
    
    switch a
    case 'default', fileName = 'single_subj_T1';
    otherwise fileName = a;
    end
    if exist([fileName '.img']) == 2, 
	cflag = 1;
    else
	disp(['Can''t locate file: ' fileName '.img'])
    O = rmfield(O,'head');
    end
end

%------------------------------------------------------------------------------
% Add surface?
%------------------------------------------------------------------------------
cflag = 0;
while cflag == 0
    if ~isfield(O,'surf')
        a = input('Add brain surface (y\n\filename,no .img extension)...?','s');
    else a = O.surf;
    end
    
    switch a
    case 'y', fileName2 = ['surf_' fileName];, dosurf = 1;
    case 'n', fileName2 = [];, cflag = 1;, dosurf = 0;
    otherwise fileName2 = a;, dosurf = 1;
    end
    if dosurf & exist([fileName2 '.mat']) == 2, 
	    cflag = 1;
    elseif dosurf
	    disp(['Can''t locate file: ' fileName '.mat'])
        O = rmfield(O,'surf');
    end
end

if ~isfield(O,'sets')
    numSets = input('Enter number of cluster sets: ');
else numSets = O.sets;
end

%------------------------------------------------------------------------------
% Add numbers?
%------------------------------------------------------------------------------
cflag = 0;
while cflag == 0
    if ~isfield(O,'numbers')
        addNums = input('Add numbers to plot? ','s');
    else addNums = O.numbers;
    end
    cflag = 1;
end

for i = 1:numSets
    
    disp(['Cluster Set ' num2str(i)])
    disp('_______________________________________________')
    
    %------------------------------------------------------------------------------
    % Get clusters, if necessary
    %------------------------------------------------------------------------------
    cflag = 0;
    while cflag == 0
       if ~isfield(O,'getfrom')
            a = input('Get clusters from...?  (workspace,file,figure)','s');
        else 
            a = O.getfrom{i};
        end
        
       switch a
       case 'workspace', 
            cflag = 1;
            if ~isfield(O,'varname')
                varname = input('Enter name of workspace variable: ','s');
                eval(['clusters = ' varname])
            else 
                varname = O.varname{i};
                eval(['clusters = O.' varname ';'])
            end
            
            
            if isempty(clusters),disp('No clusters or not a valid variable name.'),
                cflag=0;
                O = rmfield(O,'varname');    
            end
       case 'file', clusters = imageCluster('getclusters');, cflag = 1;
       case 'figure', clusters = imageCluster('getfigclusters');, cflag = 1;
       otherwise disp(['Not a valid choice...type entry string with no quotes.'])
       end
    end

    if ~isfield(O,'color')
        mycolor = input('Enter color string in single quotes or color vector for clusters: ');
    else 
        mycolor = O.color{i};
    end
    
    
    cH = [];
    for j = 1:length(clusters)
        cH(j) = imageCluster('cluster',clusters(j),'color',mycolor,'alpha',O.clusterAlpha); drawnow
        set(cH(j),'Tag',['set' num2str(i) 'cluster' num2str(j)])
    end
    
    clusterHandles{i} = cH;
    
    if strcmp(addNums,'y')
          numberHandles{i} = cl_line_plots(clusters);
    end
    
end

%------------------------------------------------------------------------------
% * lighting stuff
%------------------------------------------------------------------------------
lighting gouraud;
[az,el] = view;
myLight = lightangle(az,el);
set(myLight,'Tag','myLight')
rotate3d on

% set callback to light follow the camera
%------------------------------------------------------------------------------
if exist('lightFollowView') == 2
     set(gcf, 'WindowButtonUpFcn', 'lightFollowView');
else
     warning('Cannot find lightFollowView.m to set light position.')
end
axis vis3d
axis image
drawnow

%------------------------------------------------------------------------------
% * head surface
%------------------------------------------------------------------------------
if dosurf
    p = [];
    eval(['load ' fileName2])
    try
        p(end+1) = patch('Faces',faces,'Vertices',vertices,'FaceColor',[.5 .5 .5],'EdgeColor','none','SpecularStrength',.2,'FaceAlpha',.5,'SpecularExponent',200);
    catch
        p(end+1) = patch('Faces',faces,'Vertices',vertices,'FaceColor',[.5 .5 .5],'EdgeColor','none','SpecularStrength',.2,'SpecularExponent',200);
    end
    axis vis3d
    axis image
    drawnow
    O.p = p;
end  

O.clusterHandles = clusterHandles;
O.numberHandles = numberHandles;

if O.dohead
    [D,Ds,hdr,p,coords,X,Y,Z] = tor_3d('coords',O.coords,'whichcuts',O.whichcuts);
    set(gcf,'Color','k');axis off;colormap(gray)
end
