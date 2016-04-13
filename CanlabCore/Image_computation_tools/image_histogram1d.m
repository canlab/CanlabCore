function [cl2,classes] = image_histogram1d(varargin)
% Visualize a change_point map stored in hewma_cp.img
% (output from hewma2)
% and classify voxels into groupings based on CP
%
% :Usage:
% ::
%
%     [cl2,classes] = image_histogram1d([image name],[overlay])
%

name = 'hewma_cp.img';
%name = 'hewma_runlen.img';
overlay = [];

if length(varargin) > 0, name = varargin{1};,end
if length(varargin) > 1, overlay = varargin{2};,end

% ----------------------------------------------------
% initial viewing of CP map
% ----------------------------------------------------

[CLU.XYZ,CLU.XYZmm,CLU.Z,CLU.V] = img2voxel(name);
CLU.XYZ = CLU.XYZ(1:3,:);
%cl = mask2clusters(name);

%cluster_orthviews(cl,'overlay',overlay);
v2 = CLU.Z;

% ----------------------------------------------------
% load data
% ----------------------------------------------------

%v = spm_read_vols(spm_vol(name));
%v2 = v(:);

%discretize for convenience
v2 = round(v2.*100);
%v2(abs(v2) < eps | isnan(v2)) = [];  % shoul


% ----------------------------------------------------
% make histogram and get clusters (classifications) of CPs
% ----------------------------------------------------

nbins = input('Number of bins: '); %unique(v2);
tor_fig; f1 = gcf; [h,x] = hist(v2./100,nbins); hh = bar(x,h); set(hh,'FaceColor',[.5 .5 .5])

nclasses = input('How many classes?');
str = sprintf('%3.0f Classes: color image overlay starting at which class (e.g., 1 for all classes): ',nclasses)
startclass = input(str);

err = 1; indx = 1;
while err
    try
        classes = kmeans(v2, nclasses);   % ,'start','uniform');
        err = 0;
    catch
    end
    indx = indx + 1;
    if indx == 11, disp('kmeans: tried 10 times.  No solution.'); err = 0;, return, end
end



% For each x, find the modal class
% this is used to color the histogram
x100 = x*100;
binwid = diff(x100); binwid = binwid(1)./2;

for i =1:length(x100)
    xrange = [x100(i)-binwid x100(i)+binwid];
    xclass = classes(v2>xrange(1) & v2<=xrange(2));
    binclass(i) = mode(xclass);
end

%v2 = round(v2.*100);
% % classmap's elements are the range of values in v2
% mincp = min(v2);
% maxcp = max(v2);
% for i = mincp:maxcp, tmp = unique(classes(find(v2==i)));
%     if isempty(tmp), classmap(i) = 0;,
%     else,classmap(i) = tmp(1);, end
% end

% CLU = clusters2CLU(cl);
% CLU.Z = v2';
% CLU.Z(abs(CLU.Z)<eps) = NaN;

CLU.cp = v2';
v2(v2<=0) = 1;   %make sure no negative numbers or zeros for some reason
CLU.Z = classes';  %classmap(v2);



% ----------------------------------------------------
% define colors and sort by class size
% ----------------------------------------------------

%colors = {[1 0 0] [0 1 0] [0 0 1] [1 1 0] [1 0 1] [0 1 1]};
%while length(colors) < nclasses, colors = [colors colors];,end

colors = hot(nclasses+2 - startclass + 1);
colors = colors(2:end-1,:);

%un-plotted classes
colors = [repmat(.5,startclass-1,3); colors];

for i = 1:length(colors), col{i} = colors(i,:);,end
colors = col;
%colors = col(2:end);    %(floor(length(col)./3):end);

% sort by change point, ascending
for i = 1:nclasses,
    meancp(i) = mean(v2(classes==i));,
    nvox(i) = sum(classes==i);
    indx(i) = i;
end
[meancp,i] = sort(meancp,2,'ascend');
%nvox(i) = nvox(1:length(i));
%indx(i) = indx(1:length(i));
nvox  = nvox(i);
indx = indx(i);

meancp = meancp ./ 100; % convert back to original units.  resolution = 2 dec. places.

%table
fprintf(1,'\nClass\tNum. Voxels\tMean Value\tColor\t\t\n');
for j =1:nclasses
    fprintf(1,'%3.0f\t%3.0f\t%3.2f\t%3.2f\t%3.2f\t%3.2f\t\n',indx(j),nvox(j),meancp(j),colors{j}(1),colors{j}(2),colors{j}(3));
end
fprintf(1,'\n');

% sort by class size
%for i = 1:nclasses, nvox(i) = sum(classes==i);,end
%[nvox,i] = sort(nvox,2,'descend');

%colors(i) = colors(1:length(i));


% ----------------------------------------------------
% re-plot histogram with color codes
% ----------------------------------------------------
figure(f1);

clear i
for i = startclass:nclasses
    %     wh = find(classmap == indx(i)); range = [min(wh) max(wh)+1];
    %     wh = find(x*100 <= range(2) & x*100 >= range(1));
    %     hh = bar(x(wh),h(wh)); set(hh,'FaceColor',colors{i});

    wh = find(binclass == indx(i));
    hh = bar(x(wh),h(wh)); set(hh,'FaceColor',colors{i});

end
xlabel('Image value')
ylabel('Number of voxels')


% ----------------------------------------------------
% re-make separate clusters for each class
% and plot on brain
% ----------------------------------------------------

clear cl2
for i = startclass:nclasses
    CLUtmp = CLU;
    wh = find(CLUtmp.Z == indx(i));
    CLUtmp.XYZmm = CLUtmp.XYZmm(:,wh);
    CLUtmp.XYZ = CLUtmp.XYZ(:,wh);
    CLUtmp.Z = CLUtmp.Z(:,wh);
    CLUtmp.cp = CLUtmp.cp(:,wh);
    cl2{i} = tor_extract_rois([],CLUtmp,CLUtmp);

    if i == startclass
        cluster_orthviews(cl2{i},colors(i),'overlay',overlay);
    else
        cluster_orthviews(cl2{i},colors(i),'add','overlay',overlay);
    end
end

% write mask of sig. results
domask = input('Write mask of colormapped voxels? (1 or 0): ');
if domask
    xyz = [];
    for i=1:length(cl2)
        if ~isempty(cl2{i}), xyz = [xyz cl2{i}.XYZ];, end
    end
    V = spm_vol(name);
    V.fname = 'img_hist_mask.img';
    mask = voxel2mask(xyz', V.dim(1:3));
    spm_write_vol(V,mask);
end


return



