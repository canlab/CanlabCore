%clear xyz; 
%clear L5
%clear L5names

% if these exist already, use them
global xyz
global L5
global L5names

doL3 = 1;

if doL3
    global L3
    global L3names
    
    if isempty(xyz) | isempty(L3) | isempty(L3names) 

        disp('Loading Carmack brain talairach_info.  This may take a while.');
        load talairach_info L3 x y z; xyz = [x y z]; 
        L3names = unique(L3); %L3names = unique(L3);
    end

    disp('Region names to choose from:');

    disp(L3names)

    regname = input('Enter name of region to get ROIs for: ','s');


    figure;
    cl = talairach_clusters(xyz,L3,regname,'y');
    addbrain;
    drawnow;
    
else

    if isempty(xyz) | isempty(L5) | isempty(L5names) 

        disp('Loading talairach_info.  This may take a while.');
        load talairach_info L5 x y z; xyz = [x y z]; 
        L5names = unique(L5); %L3names = unique(L3);

    end

    disp('Region names to choose from:');

    disp(L5names)

    regname = input('Enter name of region to get ROIs for: ','s');


    figure;
    cl = talairach_clusters(xyz,L5,regname,'y');
    addbrain;
    drawnow;

end
