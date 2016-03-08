function cluster_orthviews_overlap(mask1, mask2, varargin)
% Plot two clusters on the orthviews, and their intersections for
% activations and deactivations in intermediate colors
%
% :Usage:
% ::
%
%    cluster_orthviews_overlap(mask1, mask2, [colors cell])
%
% ..
%  
% :Inputs:
% 
%     **mask1**
%     path to nifit image
% 
%     **mask2**
%     path to nifit image
%     
% :Optional Inputs:
%     
%     **colors**
%     cell array of RGB colors like: {[1 0 0] [0 0 1] [1 1 0] [0 1 1]}
% 
% ..
%    tor wager, june 2010
%    A good function, but not polished yet.  needs debugging for various
%    minor issues.
% ..

    colors = {[1 0 0] [0 0 1] [1 1 0] [0 1 1]};
    
    if length(varargin) > 0
        colors = varargin{1};
    end
    
% split into pos and neg
[maskInfo, dat] = iimg_read_img(mask1, 2);
pos = dat; pos(dat < 0) = 0;
neg = dat; neg(dat > 0) = 0;
posmask1 = pos;
negmask1 = neg;

[maskInfo2] = iimg_read_img(mask2, 2);
dat2 = scn_map_image(mask2, mask1);
dat2 = dat2(:);
% t2 = scn_map_image(mask2, tdat2);
% t2 = t2(:);
% dat2 = dat2 .* t2;

pos = dat2; pos(dat2 < 0) = 0;
neg = dat2; neg(dat2 > 0) = 0;
posmask2 = pos;
negmask2 = neg;

% get intersections

int_dat = iimg_intersection(mask1, mask2, 'posneg'); % intersections
pos12 = int_dat(:, 1);
neg12 = int_dat(:, 4);

%remove intersections from main mask dat
posmask1(pos12) = 0;
posmask2(pos12) = 0;
negmask1(neg12) = 0;
negmask2(neg12) = 0;

% get clusters

poscl1 =  iimg_indx2clusters(posmask1, maskInfo);
poscl2 =  iimg_indx2clusters(posmask2, maskInfo2);

negcl1 =  iimg_indx2clusters(negmask1, maskInfo);
negcl2 =  iimg_indx2clusters(negmask2, maskInfo2);

poscl12 = iimg_indx2clusters(pos12, maskInfo);
negcl12 = iimg_indx2clusters(neg12, maskInfo);

% montage

% pos1 neg1 pos2 neg2
%colors = {[1 0 0] [0 0 1] [1 1 0] [0 1 1]};
% reverse
%colors = {[0 0 1] [1 0 0] [0 1 1] [1 1 0]};

cluster_orthviews('whitebg');

if ~isempty(poscl1), cluster_orthviews(poscl1, colors(1), 'solid', 'add'); end
if ~isempty(negcl1), cluster_orthviews(negcl1, colors(2), 'solid', 'add'); end

if ~isempty(poscl2), cluster_orthviews(poscl2, colors(3), 'solid', 'add'); end
if ~isempty(negcl2), cluster_orthviews(negcl2, colors(4), 'solid', 'add'); end

if ~isempty(poscl12), cluster_orthviews(poscl12, {mean([colors{1}; colors{3}])}, 'solid', 'add'); end
if ~isempty(negcl12), cluster_orthviews(negcl12, {mean([colors{2}; colors{4}])}, 'solid', 'add'); end

cluster_orthviews_montage(6, 'axial');
h = findobj(gcf,'Type', 'Text'); delete(h);

cluster_orthviews_montage(6, 'sagittal');
h = findobj(gcf,'Type', 'Text'); delete(h);

disp('Positive overlap (blue)');
disp(' ')
cluster_table(poscl12, 0, 0);
disp(' ')
disp('Negative overlap (yellow)');
cluster_table(negcl12, 0, 0);
disp(' ')

