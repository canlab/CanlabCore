function [cl,cl2] = hewma_plot_bivariate(varargin)
% Visualize a change_point map stored in hewma_cp.img
% and a run-length map stored in hewma_runlen.img
% (output from hewma2)
%
% :Usage:
% ::
%
%     [cl,cl2] = hewma_plot_bivariate([Method],[optional overlay image])
%
% and classify voxels into groupings based on locations in the bivariate
% space

Method = 'group';
if length(varargin), Method = varargin{1};, end
if length(varargin) > 1, ovl = varargin{2}; else, ovl = which('scalped_single_subj_T1.img');, end


doplot = 1;


% ----------------------------------------------------
% *
% *
% * get data, x, and list of all voxels to go with each data point, CLU
% *
% *
% ----------------------------------------------------
switch Method
    
    case 'group'

        mask = 'hewma_sig.img';
        t = 'hewma_t.img';

        name1 = 'hewma_cp.img';
        name2 = 'hewma_runlen.img';
        name2 = 'gausslongest.img';

    case 'difference'
        mask = 'hewma_sigdiff.img';
        t = 'hewma_tdiff.img';


        name1 = 'hewma_cpdiff.img';
        name2 = 'hewma_runlendiff.img';

    otherwise
        error('Method (first arg) must be group or difference.');
end

% ----------------------------------------------------
% initial viewing of sig map
% ----------------------------------------------------

cl = mask2clusters(mask,t);

if doplot, cluster_orthviews(cl,'overlay',ovl,'bivalent'); %,{[1 0 0]});,
    %cluster_orthviews(cl,{[1 0 0]},'add');,
    %cluster_orthviews(cl,{[1 0 0]},'add');,
end

CLU = clusters2CLU(cl);

% ----------------------------------------------------
% get mask for pos and neg activation
% ----------------------------------------------------
v = spm_read_vols(spm_vol(t));


% ----------------------------------------------------
% load data, get data for each voxel in CLU
% ----------------------------------------------------

v1 = spm_read_vols(spm_vol(name1));
v2 = spm_read_vols(spm_vol(name2));

for i = 1:size(CLU.XYZ,2)
    % data in voxel list
    x(i,1) = v1(CLU.XYZ(1,i),CLU.XYZ(2,i),CLU.XYZ(3,i));
    x(i,2) = v2(CLU.XYZ(1,i),CLU.XYZ(2,i),CLU.XYZ(3,i));
end

x = round(x);



go = input('Make tables?');
if go, mincs = input('Min cluster size?  (in voxels): ');, end

% ----------------------------------------------------
% POSITIVE
% ----------------------------------------------------

% select positive points and CLU
wh = find(CLU.Z > 0);
dat = x; dat = dat(wh,:); 
CLU2 = CLU; CLU2.XYZ = CLU2.XYZ(:,wh); CLU2.XYZmm = CLU2.XYZmm(:,wh); CLU2.Z = CLU2.Z(:,wh);



[cl2,nclasses,colors] = cluster_kmeans_parcel(dat,CLU2,1);


% -------------------------------------------------------------------
% * tables
% -------------------------------------------------------------------

if go
    
    fprintf(1,'\nPositive deviations\n-----------------------------------\n');
    
    for i = 1:length(cl2);
        
        whomit = cat(1,cl2{i}.numVox); whomit(whomit >= mincs) = 0; cl2{i}(find(whomit)) = [];
        
        fprintf(1,'Class %3.0f, color = %s\n-----------------------------------\n',i,num2str(colors{i}));
        if isempty(cl2{i})
            fprintf(1,'No clusters large enough.\n'); 
        else
            cluster_table(cl2{i});
        end
        fprintf(1,'\n-----------------------------------\n')
    end
end

cl = cl2;

% ----------------------------------------------------
% NEGATIVE
% ----------------------------------------------------

dodecrease = input('Look at deactivations?');

if dodecrease

    % select negative points and CLU
    wh = find(CLU.Z < 0);
    dat = x; dat = dat(wh,:); 
    CLU2 = CLU; CLU2.XYZ = CLU2.XYZ(:,wh); CLU2.XYZmm = CLU2.XYZmm(:,wh); CLU2.Z = CLU2.Z(:,wh);



    [cl2,nclasses,colors] = cluster_kmeans_parcel(dat,CLU2,1);


    % -------------------------------------------------------------------
    % * tables
    % -------------------------------------------------------------------

    if go

        fprintf(1,'\nNegative deviations\n-----------------------------------\n');

        for i = 1:length(cl2);

           whomit = cat(1,cl2{i}.numVox); whomit(whomit >= mincs) = 0; cl2{i}(find(whomit)) = [];

            fprintf(1,'Class %3.0f, color = %s\n-----------------------------------\n',i,num2str(colors{i}));
            if isempty(cl2{i})
                fprintf(1,'No clusters large enough.\n'); 
            else
                cluster_table(cl2{i});
            end
            fprintf(1,'\n-----------------------------------\n')
        end
    end


    cl = [cl cl2];

end

cl2 = cl; % cl2 has output in cell array, sep. by class
cl = cat(2,cl{:}); % cl is single vector

% -------------------------------------------------------------------
% * cluster extraction
% -------------------------------------------------------------------
        
go = input('Save timeseries in clusters?');

if go
    nm= 'hewma_timeseries.mat'; ctr = 1;
    
    while exist(nm) == 2
        nm = [nm(1:end-4) num2str(ctr) '.mat'];
        ctr = ctr+1;
    end
    fprintf(1,'Extracting cluster data and saving in %s\n',nm);
    mincs = input('Min cluster size?  (in voxels): ');
    whomit = cat(1,cl.numVox); whomit(whomit >= mincs) = 0; cl(find(whomit>0)) = [];
    
    cl = hewma_save_timeseries(cl,mincs);
    eval(['save ' nm ' cl']);
    
end


% -------------------------------------------------------------------
% * timeseries plotting
% -------------------------------------------------------------------
        
go = input('Plot timeseries interactively?');

if go
    if ~(exist('EXPT') == 1)
        file = spm_get(1,'*mat','Select EXPT.mat',pwd);
        [dd,ff,ee] = fileparts(file);
        cd(dd)
        load(file)
    end
    
    % used in button-up fcn callback
    E2 = EXPT;
    clear EXPT

    global VOL
    global f
    global f2
    global EXPT
    EXPT = E2;


    set(gcf,'WindowButtonUpFcn','[dat,files,stats,mycov] = hewma_plot_coord_btnupfcn;')

    % get coordinate mapping matrix
    VOL = struct('M',cl(1).M);

    % prepare figure
    f1 = figure('Color','w','Name','Hewma plots');
    f = f1;     % button callback uses figure f


    stop = input('Click on figure to plot.  Press return to exit.');
    
end


return



