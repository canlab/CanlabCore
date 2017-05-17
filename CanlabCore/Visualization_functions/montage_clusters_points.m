function ph = montage_clusters_points(ovl,clusters,XYZpts,varargin)
% ::
%
%    ph = montage_clusters_points(ovl,clusters,XYZpts,varargin)
%
% :Inputs:
%
%   **varargin:**
%        additional clusters structures
%
%   **XYZpts:**
%        XYZ mm coordinates of points to plot
%
%   **ph:**
%        point handles
%
% ..
%    Tor Wager
% ..

if isempty(ovl), ovl = which('scalped_single_subj_T1.img');, end
myc = {'b' 'r' 'y' 'g' 'm' 'c'};

if size(XYZpts,1) ~= 3, XYZpts = XYZpts';,end

XYZmm = cat(2,clusters.XYZmm);
XYZmm = [XYZmm XYZpts];

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

XYZ = mm2voxel(XYZmm,V)';
whsl = unique(XYZ(3,:));
nsl = length(whsl) + 1;
rc = ceil(sqrt(nsl));
h = [];


figure; colormap gray; set(gcf,'Color','w')


% plot first cluster structure

index = plot_cluster(clusters,oimg,rc,V,whsl,myc{1},textx,texty,1,size(oimg));
ph = plot_points(XYZpts,rc,V,whsl,myc{2});

% plot additional cluster structures

if length(varargin) > 0
    for i = 1:length(varargin)
        cl = varargin{i};
        index = plot_cluster(cl,[],rc,V,whsl,myc{i+1},textx,texty,i+1,size(oimg));
        
    end
end



return



% sub-functions
%

function index = plot_cluster(clusters,oimg,rc,V,whsl,myc,textx,texty,clind,odims)

XYZmm = cat(2,clusters.XYZmm);
XYZ = mm2voxel(XYZmm,V)';

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



function ph = plot_points(XYZmm,rc,V,whsl,myc)

XYZ = mm2voxel(XYZmm,V)';
index = 1;
phind = 1;
for z = whsl
    
    subplot(rc,rc,index);
    hold on
    myXYZ = XYZ(:,XYZ(3,:) == z);
    
    for i = 1:size(myXYZ,2)
        ph(phind) = plot3(myXYZ(2,i),myXYZ(1,i),myXYZ(3,i),[myc(1) '.'],'MarkerFaceColor',myc(1),'MarkerSize',8);
        phind = phind + 1;    
    end
    
    index = index + 1;
        
end

return
