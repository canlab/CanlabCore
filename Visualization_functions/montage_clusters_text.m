function montage_clusters_text(ovl,clusters,varargin)
% montage_clusters_text(ovl,clusters,varargin)
% 
% Tor Wager
% varargin = additional clusters structures
% this function puts text cluster numbers on cluster centers

% color = cell array of text strings indicating colors {'r' 'g'} etc...


if isempty(ovl), ovl = which('scalped_single_subj_T1.img');, end

myc = {'b' 'r' 'y' 'g' 'm' 'c'};

XYZmm = cat(2,clusters.XYZmm);

if length(varargin) > 0
    for i = 1:length(varargin)
        cl = varargin{i};
        clXYZmm = cat(2,cl.XYZmm);
        XYZmm = [XYZmm clXYZmm];
    end
end

% ----------------------------------------
% * overlay image
% ----------------------------------------
V = spm_vol(ovl);
oimg = spm_read_vols(V);
V.M = V.mat;

textx = size(oimg,2) - 50;
texty = size(oimg,1) - 6;

%[array,hdr,h,whichslices,rows,cols,figh] = readim2(ovl,'p');

% how many slices, which ones

XYZ = mm2voxel(XYZmm,V,2)';
whsl = unique(XYZ(3,:));
nsl = length(whsl) + 1;
rc = ceil(sqrt(nsl));
h = [];

% get centers for text plotting
cen = getcenters(clusters);

figure; colormap gray; set(gcf,'Color','w')


% plot first cluster structure
for i = 1:length(clusters)
    index = plot_cluster(clusters(i),oimg,rc,V,whsl,myc{1},textx,texty,1,size(oimg),cen(i,:),num2str(i));
    drawnow
end

% plot additional cluster structures

if length(varargin) > 0
    for i = 1:length(varargin)
        cl = varargin{i};
        cen = getcenters(cl);
        index = plot_cluster(cl,[],rc,V,whsl,myc{i+1},textx,texty,i+1,size(oimg),cen);
        
    end
end



return



% sub-functions
%

function index = plot_cluster(clusters,oimg,rc,V,whsl,myc,textx,texty,clind,odims,varargin)
% var argument is list of centers, x y z are rows, cols are coords

XYZmm = cat(2,clusters.XYZmm);
XYZ = mm2voxel(XYZmm,V)';

% text numbers, if plotting text at xyz coords of centers (varargin)
if length(varargin) > 0
    cencoo = varargin{1};
    cencoo = mm2voxel(cencoo,V)';
    centext = varargin{2};
    %cenind = 1;
end

% surface patch method
% ----------------------------------------------------------------------------------------
vol = voxel2mask(XYZ',odims);
vol = smooth3(vol);


index = 1;
for z = whsl
    
    subplot(rc,rc,index); 
    
    if ~isempty(oimg)
        set(gca,'YDir','reverse');
        imagesc(oimg(:,:,z))
        hold on; axis image; axis off
        zmm = voxel2mm([1 1 z]',V.mat);
        text(textx,texty,['z = ' num2str(zmm(3))],'Color','w')
    else
        hold on
    end
    

    if z>1,FVC = isocaps(vol(:,:,z-1:z),0,'zmax');
    else FVC = isocaps(vol(:,:,z:z+1),0,'zmax');
    end
    try
	patch(FVC,'EdgeColor','none','FaceColor',myc,'FaceAlpha',.7)
    catch
	patch(FVC,'EdgeColor','none','FaceColor',myc)
    end
    
    % text numbers, if plotting text at xyz coords of centers (varargin)
    if length(varargin) > 0
        cencooz = cencoo(:,cencoo(3,:) == z);
        for i = 1:size(cencooz,2)
            text(cencooz(2,i),cencooz(1,i),centext,'Color','k'); %num2str(cenind),'Color','k')
            %cenind = cenind + 1;
        end
    end
    
    % plot method
    
    %myXYZ = XYZ(1:2,XYZ(3,:) == z);
    %plot(myXYZ(2,:),myXYZ(1,:),[myc 's'],'MarkerFaceColor',myc,'MarkerSize',3)
    
    index = index + 1;
    drawnow
    
end

subplot(rc,rc,index)
a = pwd; a = a(end-6:end);
b = num2str(clusters(1).threshold);

c = num2str(length(clusters));
text(0,clind-1,[myc ': ' clusters(1).title ' ' a ' u = ' b ', ' c ' clusters'])
axis off
axis([0 1 -1 clind])

return


function cen = getcenters(clusters)
    for i = 1:length(clusters)
        cen(i,:) = clusters(i).mm_center;
    end
    cen = cen';
return
    
