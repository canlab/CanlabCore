function montage_clusters_maxslice(ovl,clusters,varargin)
% :varargin: (in any order) = 
%   a) additional clusters structures 
%   b) cell array of colors (text format), must be ROW vector {'r' 'g'} etc...
%       if length of color string is longer than number of clusters inputs,
%       additional colors are interpreted as '2-intersection' and 'all-intersection' 
%       colors, in that order.  This overrides single string argument (c, below) for
%       color input 
%   c) single string color argument for overlaps (intersections)
%       plots intersections of ANY TWO clusters right now.
%       also color for plotting points, if entered.
%       use 'nooverlap' as an input argument to suppress this.
%   d) [n x 3] matrix of points to plot on map 
%   e) text labels for points, must be cell array in COLUMN vector
%   f) single number, 1/0 for whether to plot overlapping coordinates in overlap colors
%       default is 1.
%   g) color limit vector [min max] indicates color mapping for blobs rather than
%       solid colors.  This will do hot/cool mapping
%
% Intersections of 2 colors are magenta, and ALL colors are yellow
% unless otherwise specified
%
% This maxslice version does 1) Not create a new figure, and 2) finds max
% slice through clusters to create single slice image.
%
% Useful for displaying next to timecourse plots (for example)
% try:
% ::
%
%    CLU = clusters2clu(clusters);
%    spm_orthviews('AddColouredBlobs',1,CLU.XYZ,CLU.Z,CLU.M,[1 0 0])
%
% ..
%    by Tor Wager  last edit 12/11/02
% ..

% ..
%    defaults
% ..

if isempty(ovl), ovl = which('scalped_single_subj_T1.img');, end
myc = {'b' 'r' 'g' 'c' 'm' 'w'};
customcolors = 0;
bcolor = 'm';
acolor = 'y';
XYZmm_both = [];
clindx = 1;
XYZpts = [];
myplab = [];
cl = [];
plotovl = 1;
doverb = 0;
colorlimits = [];

% recursive operation if clusters is a cell array
% indicating multiple clusters structures
%if iscell(clusters),
%    for i = 1:length(clusters),
%        montage_clusters(ovl,clusters{i},varargin)
%    end
%    return
%end

XYZmm = cat(2,clusters.XYZmm);

% Single slice stuff: Find most frequent XYZmm
tmp = unique(XYZmm(3,:)); for i = 1:length(tmp), m(i) = sum(XYZmm(3,:) == tmp(i));,end
wh = tmp(find(m==max(m))); wh = wh(1);
XYZmm = XYZmm(:,find(XYZmm(3,:) == wh));


XYZmm_all = XYZmm;


% ----------------------------------------
% * process input arguments
% ----------------------------------------

if length(varargin) > 0
    for i = 1:length(varargin)
        if isstruct(varargin{i})
            % struct inputs interpreted as clusters variables
            cl{clindx} = varargin{i};
            clXYZmm = cat(2,cl{clindx}.XYZmm);
            clindx = clindx + 1;
            a = XYZmm'; b = clXYZmm';
            XYZmm_both = [XYZmm_both intersect(a,b,'rows')'];
            if ~isempty(XYZmm_both), XYZmm_all = intersect(XYZmm_all',b,'rows')';, else, XYZmm_all = [];, end
            XYZmm = [XYZmm clXYZmm];
           
        elseif isstr(varargin{i})
            % string arguments interpreted as overlap colors
            if strcmp(varargin{i},'nooverlap'),plotovl = 0;,
            else
                bcolor = varargin{i};
            end
        elseif iscell(varargin{i})
            % cell array vectors with one row interpreted as color inputs
            if size(varargin{i},1) == 1
                customcolors = 1;
                myc = varargin{i};
            elseif size(varargin{i},2) == 1
                % cell array vectors with one column interpreted as point coordinates to plot
                myplab = varargin{i};
            else
                error('cell array must be row (for colors) or column (for text labels) vector.')
            end
        elseif prod(size(varargin{i})) == 1
            % single integers interpreted as 'do overlap plot' flag, can be 1 or 0
            % default is 1
            plotovl = varargin{i};
        elseif any(size(varargin{i}) == 3)
            % any 3-vector matrix interpreted as list of points to plot on 
            XYZpts = varargin{i};
            %XYZmm = [XYZmm XYZpts];
            % then it's a list of points to plot
        elseif size(varargin{i},1) == 1 & size(varargin{i},2) == 2
            % 1 x 2 double vector, means we want height-mapped colors instead of solid
            colorlimits = varargin{i};
        else
            error('Unknown input argument type.')
        end
    end
end

if length(cl) < 1, plotovl = 0;, end

if customcolors
    if length(myc) > length(cl)+2, acolor = myc{length(cl) + 3};, end
    if length(myc) > length(cl)+1, bcolor = myc{length(cl) + 2};, end
end
  
if doverb
    disp([num2str(length(cl) + 1) ' Clusters found.'])
    if plotovl
        disp(['Two cluster overlap: ' bcolor ', found ' num2str(length(XYZmm_both)) ' voxels.'])
        disp(['All cluster overlap: ' acolor ', found ' num2str(length(XYZmm_all)) ' voxels.'])
    else
        disp(['No overlap plotting.'])
    end
end

% ----------------------------------------
% * overlay image
% ----------------------------------------
try
    V = spm_vol(ovl);
    oimg = spm_read_vols(V);
    V.M = V.mat;
catch
    disp('Cannot use spm_vol: No SPM?')
    [oimg,hdr,h] = readim2(ovl);
    V.mat = diag([hdr.xsize hdr.ysize hdr.zsize 1]); V.mat(:,4) = [hdr.origin(1:3); hdr.SPM_scale];
    orig = (hdr.origin(1:3) - [hdr.xdim hdr.ydim hdr.zdim]') .* [hdr.xsize hdr.ysize hdr.zsize]';
    V.mat(1:3,4) = orig;
    V.mat = [ 
     2     0     0   -92
     0     2     0  -128
     0     0     2   -74
     0     0     0     1];
     
    V.M = V.mat;
end
    
textx = size(oimg,1) - 50;
texty = 6; %size(oimg,2) - 6;

%[array,hdr,h,whichslices,rows,cols,figh] = readim2(ovl,'p');

% how many slices, which ones

XYZ = mm2voxel(XYZmm,V,1)';    % 1 no re-ordering, allows repeats             %2 THIS RE-ORDERS VOXELS, BUT IS FAST, AND ORDER SHOULDN'T MATTER HERE.
whsl = unique(XYZ(3,:));
nsl = length(whsl) + 1;
rc = ceil(sqrt(nsl));
h = [];

colormap gray; set(gcf,'Color','w')
cc = colormap(gray); cc(1:10,:) = repmat([1 1 1],10,1);    % [0 0 .3] for dark blue
colormap(cc)


% ----------------------------------------
% * plot first cluster structure
% ----------------------------------------

index = plot_cluster(clusters,oimg,rc,V,whsl,myc{1},textx,texty,1,size(oimg),colorlimits,1);
if ~isempty(XYZpts), ph = plot_points(XYZpts,rc,V,whsl,bcolor,myplab);, end
% plot points for single-voxel clusters
if isempty(colorlimits)
    for i = 1:length(clusters), 
        try,if clusters(i).numVox == 1, plot_points(clusters(i).XYZmm,rc,V,whsl,myc{1},myplab);, end,catch,end
    end
end

% ----------------------------------------
% * plot additional cluster structures
% ----------------------------------------

if length(cl) > 0
    for i = 1:length(cl)
        index = plot_cluster(cl{i},[],rc,V,whsl,myc{i+1},textx,texty,i+1,size(oimg),colorlimits,1);
        if isempty(colorlimits)
            for j = 1:length(cl{i}), 
                try,if cl{i}(j).numVox == 1, plot_points(cl{i}(j).XYZmm,rc,V,whsl,myc{i+1},myplab);,end,catch,end
            end
        end
    end
end


% ----------------------------------------
% * plot overlap areas
% ----------------------------------------

if plotovl
    
    if ~isempty(XYZmm_both)
        bcl.XYZmm = XYZmm_both; bcl.M = cl{1}(1).M;
        plot_cluster(bcl,[],rc,V,whsl,bcolor,textx,texty,i+1,size(oimg),colorlimits,1);
    end

    if ~isempty(XYZmm_all)
        bcl.XYZmm = XYZmm_all; bcl.M = cl{1}(1).M;
        plot_cluster(bcl,[],rc,V,whsl,acolor,textx,texty,i+1,size(oimg),colorlimits,1);
    end

end

% fix bug at end - replot XYZ slice text
%for z = whsl
%    subplot(rc,rc,index); 
%    zmm = voxel2mm([1 1 z]',V.mat);
%    text(textx,texty,['z = ' num2str(zmm(3))],'Color','w')
%    zi(z) = zmm(3);
%end
%zi(zi > 0)

%try
%    enlarge_axes(gcf)
%catch
%    disp('Error enlarging axes. Skipping.')
%end

return




% ----------------------------------------
%
% * Sub-functions
%
% ----------------------------------------

function index = plot_cluster(clusters,oimg,rc,V,whsl,myc,textx,texty,clind,odims,colorlimits,varargin)
% varargin suppresses end text

% take only first row of Z scores
for i = 1:length(clusters), clusters(i).Z = clusters(i).Z(1,:);,end


% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
% Set up
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

XYZmm = cat(2,clusters.XYZmm);
XYZ = mm2voxel(XYZmm,V,1)';     % no unique suppression

% get relative voxel sizes for fill
xs = diag(clusters(1).M(1:3,1:3));
xs = xs ./ diag(V.mat(1:3,1:3));
% add 25% to make sure no gaps in plot
xs = xs + .25*xs;

if ~isempty(colorlimits)
    
    if colorlimits(1) < 1 & colorlimits(2) < 1
       
        % -------------------------------------------------------------
	    % define 4 split color maps - red/yellow, green/blue
        % -------------------------------------------------------------
    
    end
    
    
    % -------------------------------------------------------------
	% define color maps - biscale hot/cool
    % -------------------------------------------------------------

	% color map - hot
	% --------------------------------------------
	h1 = (0:1/99:1)';
	h2 = ones(size(h1)); 
	h3 = zeros(size(h1));
	h = [h1 h3 h3; h2 h1 h3; h2 h2 h1];
	h(1:50,:) = []; % take only red values to start
	% in new matlab: h = colormap(hot(300));

	% color map - winter
	% --------------------------------------------
	h1 = (0:1/249:1)';
	h2 = (1:-1/(249*2):.5)';
	h3 = zeros(size(h1));
	hc = [h3 h1 h2];

    % -------------------------------------------------------------
	% determine overall z-score range
    % -------------------------------------------------------------
    
	allz = cat(2,clusters.Z); 
    zrange = allz;
	tmp = zrange(zrange > 0);
	tmpc = zrange(zrange < 0);

	if ~isempty(tmp)
		zrange = [min(tmp) max(tmp)];
		zh = zrange(1):(zrange(2)-zrange(1))./249:zrange(2);
		zh = round(zh*100);
    else
        zh = [];
    end

	if ~isempty(tmpc)
		zrangec = [min(tmpc) max(tmpc)];
		zhc = zrangec(1):(zrangec(2)-zrangec(1))./249:zrangec(2);
        if isempty(zhc), zhc = 1:length(hc);,end
		zhc = round(zhc*100);
    else, 
        zhc = [];
    end
    
    if isempty(zh) & isempty(zhc), error('No Z-scores in cluster'),end

else

    % surface patch method
    % ----------------------------------------------------------------------------------------
    vol = voxel2mask(XYZ',odims);
    vol = smooth3(vol);
    
end

% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
% Loop through slices
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

index = 1;
for z = whsl
    
    %subplot(rc,rc,index); 
    
    
    if ~isempty(oimg)
        imagesc(oimg(:,:,z)')
        set(gca,'YDir','normal');
        hold on; axis image; axis off
        zmm = voxel2mm([1 1 z]',V.mat);
        text(textx,texty,['z = ' num2str(zmm(3))],'Color','k','FontSize',18)
    else
        hold on
    end
    
    doimgpatch = 0;
    
    if doimgpatch   %isempty(colorlimits)
        
        % --------------------------------------------------------------------------------
        % solid color patches
        % --------------------------------------------------------------------------------
        if z>1,
            mvol = vol(:,:,z-1:z); for i = 1:size(mvol,3),myvol(:,:,i) = mvol(:,:,i)';,end
            %FVC = isocaps(vol(:,:,z-1:z)',0,'zmax');
        else 
            mvol = vol(:,:,z:z+1); for i = 1:size(mvol,3),myvol(:,:,i) = mvol(:,:,i)';,end
            %FVC = isocaps(vol(:,:,z:z+1)',0,'zmax');
        end
        FVC = isocaps(myvol,0,'zmax');
    
        try
	        patch(FVC,'EdgeColor','none','FaceColor',myc,'FaceAlpha',1)
        catch
	        patch(FVC,'EdgeColor','none','FaceColor',myc)
        end
    
    else
        % --------------------------------------------------------------------------------
        % color-mapped points
        % --------------------------------------------------------------------------------
    
        myxyz = XYZ(1:2,XYZ(3,:) == z);
        
        if ~isempty(colorlimits)
            myz = allz(XYZ(3,:) == z);
            clear h2, clear wh
        end
        
        %plot(myXYZ(2,:),myXYZ(1,:),[myc 's'],'MarkerFaceColor',myc,'MarkerSize',3)
        hold on
        for j = 1:size(myxyz,2)
            
            % added to replace image patch, solid color point fill for pos Z values only
            if isempty(colorlimits)
                h2(j) = fill([myxyz(1,j) myxyz(1,j) myxyz(1,j)+xs(1) myxyz(1,j)+xs(1)], ...
                    [myxyz(2,j) myxyz(2,j)+xs(2) myxyz(2,j)+xs(2) myxyz(2,j)],myc, ...
                    'EdgeColor','none');
                
            else
                    
            if myz(j) < 0
                tmp = find((zhc-round(myz(j)*100)).^2 == min((zhc-round(myz(j)*100)).^2));
                wh(j) = tmp(1);
                %h2(j) = plot(myxyz(1,j),myxyz(2,j),'Color',hc(wh(j),:),'MarkerSize',3,'MarkerFaceColor',hc(wh(j),:));
                

                h2(j) = fill([myxyz(1,j) myxyz(1,j) myxyz(1,j)+xs(1) myxyz(1,j)+xs(1)],[myxyz(2,j) myxyz(2,j)+xs(2) myxyz(2,j)+xs(2) myxyz(2,j)],hc(wh(j),:),'EdgeColor','none');
            else
                tmp = find((zh-round(myz(j)*100)).^2 == min((zh-round(myz(j)*100)).^2));
                wh(j) = tmp(1);
                %h2(j) = plot(myxyz(1,j),myxyz(2,j),'Color',h(wh(j),:),'MarkerSize',3,'MarkerFaceColor',h(wh(j),:));
                h2(j) = fill([myxyz(1,j) myxyz(1,j) myxyz(1,j)+xs(1) myxyz(1,j)+xs(1)],[myxyz(2,j) myxyz(2,j)+xs(2) myxyz(2,j)+xs(2) myxyz(2,j)],h(wh(j),:),'EdgeColor','none');   
            end
            %if exist('h2') == 1, set(h2,'Marker','square'),end
            
        end %if isempty colorlimits
        end
    
    end % if solid or color-mapped
    
    index = index + 1;
    drawnow
    
end % slice loop


% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
% Text and color bars
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

if length(varargin) == 0
    subplot(rc,rc,index)
    a = pwd; a = a(end-6:end);
    b = num2str(clusters(1).threshold);
    axis off
    mfig = gcf;
    c = num2str(length(clusters));
    text(0,clind-1,[myc ': ' clusters(1).title ' ' a ' u = ' b ', ' c ' clusters'])
    axis([0 1 -1 clind])
    
    if ~isempty(colorlimits)
        if ~isempty(tmp)
    		bfig = figure('Color','w');
            hold on;
            
            if any(allz > 0)
		        zh2 = zh./100;
		        for i = 2:size(h,1), fill([zh2(i-1) zh2(i-1) zh2(i) zh2(i)],[0 1 1 0],h(i,:),'EdgeColor','none');, end
                plot([min(zh2) max(zh2)],[0 0],'k-')
		        set(gca,'YTickLabel',''); 
                xlabel('Z-score','FontSize',14)
            end
    	end
        
        if any(allz < 0)
            %if ~exist('bfig'), 
                %bfig = figure('Color','w'); hold on;, 
            %end
            zh2 = zhc./100;
		    for i = 2:size(hc,1), fill([zh2(i-1) zh2(i-1) zh2(i) zh2(i)],[0 1 1 0],hc(i,:),'EdgeColor','none');, end
            plot([min(zh2) max(zh2)],[0 0],'k-')
		    set(gca,'YTickLabel',''); 
            xlabel('Z-score','FontSize',14)
        end
            
    end
    
    %set(gca,'Position',[0.1300    0.1100    0.7750    0.8150])
    set(gca,'FontSize',18)
    xlabel('Z-score','FontSize',18)
    %set(gca,'Position',[0.1300    0.1100    0.7750    0.8150])
end


return




function ph = plot_points(XYZmm,rc,V,whsl,myc,myplab)

XYZ = mm2voxel(XYZmm,V,1)';     % suppress unique voxel output
index = 1;
phind = 1;

for z = whsl
    
    subplot(rc,rc,index);
    hold on
    myXYZ = XYZ(:,XYZ(3,:) == z);
    if ~isempty(myplab), myplz = myplab(XYZ(3,:) == z);, end
    
    for i = 1:size(myXYZ,2)
        ph(phind) = plot3(myXYZ(1,i),myXYZ(2,i),100,[myc(1) '.'],'MarkerFaceColor',myc(1),'MarkerSize',8);
        
        if ~isempty(myplab)
            text(myXYZ(1,i),myXYZ(2,i),100,myplz{i},'Color',myc(1))
        end
        
        phind = phind + 1;    
    end
    
    index = index + 1;   
    %view(0,90)
    %plot(0,0,'kd')
end
