clear mycolor
clear all_clusters
clear all_clusterHandles

%------------------------------------------------------------------------------
% Get head, if not default
%------------------------------------------------------------------------------
cflag = 0;
while cflag == 0
    a = input('Use default head single_subj_T1...?  (type default or file name, no .img extension)','s');
    switch a
    case 'default', fileName = 'single_subj_T1';
    otherwise fileName = a;
    end
    if exist([fileName '.img']) == 2, 
	cflag = 1;
    else
	disp(['Can''t locate file: ' fileName '.img'])
    end
end

%------------------------------------------------------------------------------
% Add surface?
%------------------------------------------------------------------------------
cflag = 0;
while cflag == 0
    a = input('Add brain surface (y\n\filename,no .img extension)...?','s');
    switch a
    case 'y', fileName2 = ['surf_' fileName];, dosurf = 1;
    case 'n', fileName2 = [];, cflag = 1;, dosurf = 0;
    otherwise fileName2 = a;, dosurf = 1;
    end
    if dosurf & exist([fileName2 '.mat']) == 2, 
	    cflag = 1;
    elseif dosurf
	    disp(['Can''t locate file: ' fileName '.mat'])
    end
end

numSets = input('Enter number of cluster sets: ');

for i = 1:numSets
    
    disp(['Cluster Set ' num2str(i)])
    disp('_______________________________________________')
    
    %------------------------------------------------------------------------------
    % Get clusters, if necessary
    %------------------------------------------------------------------------------
    cflag = 0;
    while cflag == 0
       a = input('Get clusters from...?  (workspace,file,figure)','s');
       switch a
       case 'workspace', 
            cflag = 1;
            varname = input('Enter name of workspace variable: ','s');
            eval(['clusters = ' varname ';'])
            if isempty(clusters),disp('No clusters or not a valid variable name.'),cflag=0;,end
       case 'file', clusters = imageCluster('getclusters');, cflag = 1;
       case 'figure', clusters = imageCluster('getfigclusters');, cflag = 1;
       otherwise disp(['Not a valid choice...type entry string with no quotes.'])
       end
    end

    mycolor{i} = input('Enter color string in single quotes or color vector for clusters: ');
    all_clusters{i} = clusters;
end


%------------------------------------------------------------------------------
% Make 3 views; loop through cluster sets
%------------------------------------------------------------------------------

figure; 
myView = {[0 90] [270 0] [180 0]}

for myAxis = 1:3
    subplot(1,3,myAxis); hold on
    
    
    for i = 1:numSets
        cH = [];
        for j = 1:length(all_clusters{i})
            cH(j) = imageCluster('cluster',all_clusters{i}(j),'color',mycolor{i},'alpha',.6); drawnow
            set(cH(j),'Tag',['set' num2str(i) 'cluster' num2str(j)])
        end
    
        clusterHandles{i} = cH;
    
    end

    all_clusterHandles{myAxis} = clusterHandles;
    
    %------------------------------------------------------------------------------
    % * lighting stuff
    %------------------------------------------------------------------------------
    view(myView{myAxis}(1),myView{myAxis}(2))
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
            p(end+1) = patch('Faces',faces,'Vertices',vertices,'FaceColor',[.5 .5 .5],'EdgeColor','none','SpecularStrength',.2,'FaceAlpha',.2,'SpecularExponent',200);
        catch
            p(end+1) = patch('Faces',faces,'Vertices',vertices,'FaceColor',[.5 .5 .5],'EdgeColor','none','SpecularStrength',.2,'SpecularExponent',200);
        end
        axis vis3d
        axis image
        drawnow
    end  

end     % loop through axes