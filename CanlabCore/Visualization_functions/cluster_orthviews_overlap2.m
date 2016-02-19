function cluster_orthviews_overlap2(masks, varargin)
% Plot blobs on the orthviews, and their intersections for
% activations and deactivations in intermediate colors
%
% :Usage:
% ::
%
%    cluster_orthviews_overlap2(masks, ['colors', colors cell], ['surface'], ['negative'] )
%
% Positive effects only!
%
% :Inputs
% 
%     **masks**
%     image_vector  or fmri_data object
%     
% :Optional Inputs
% 
% 	**colors**
%     'colors',{[RGB],[RGB],..} cell array of RGB triplets
%     
%     **surface**
%     'surface',[0/1]
%     default =1 
%     
%     **negative**
%     'negative' string to only use negative clusters, otherwise only positive clusters
%     
%     **nomontage**
%     'nomontage',[0/1]
%     don't do montage, default do a montage
%     
% ..
%    tor wager, june 2010
%    A good function, but not polished yet.  needs debugging for various
%    minor issues.
% ..

if isa(masks, 'image_vector')
    loadfromfile = 0;
    n = size(masks.dat, 2);
    
elseif isstruct(masks)
    % cl structure
    loadfromfile = 0;
    n = length(masks);
    
elseif ischar(masks)
    % file names
    loadfromfile = 1;
    n = size(masks, 1);
    
else
    error('Inputs must be all cl structs or mask filenames');
end

colors = scn_standard_colors(n);

dosurface = 0;
domontage = 1;
signstr = 'positive';
    
for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            
            % functional commands
            case 'colors', colors = varargin{i+1};
            case 'surface', dosurface = 1;
            case 'negative', signstr = 'negative';
            case 'nomontage', domontage = 0;
                
            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end

if loadfromfile
    
    error('Not implemented yet.  Use object-oriented tools to get fmri_data object');
    
%     [maskInfo1, dat1] = iimg_read_img(mask1, 2);
%     [maskInfo2, dat2] = iimg_read_img(mask2, 2);
%     [maskInfo3, dat3] = iimg_read_img(mask3, 2);
%     
%     % positive only
%     dat1 =  dat1 > 0;
%     dat2 =  dat2 > 0;
%     dat3 =  dat3 > 0;
%     
%     int_dat = iimg_intersection(mask1, mask2, 'posneg'); % intersections
%     pos12 = int_dat(:, 1);
%     
%     int_dat = iimg_intersection(mask1, mask3, 'posneg'); % intersections
%     pos13 = int_dat(:, 1);
%     
%     int_dat = iimg_intersection(mask2, mask3, 'posneg'); % intersections
%     pos23 = int_dat(:, 1);
%     
%     pos123 = pos12 & pos13;
%     
%     intersect_all = any([pos12 pos13 pos23 pos123], 2);
%     
%     %remove higher-order intersections from  mask dat
%     dat1(intersect_all) = 0;
%     dat2(intersect_all) = 0;
%     dat3(intersect_all) = 0;
%     
%     pos12(pos123) = 0;
%     pos23(pos123) = 0;
%     pos13(pos123) = 0;
%     
%     % make clusters
%     clear cl
%     
%     cl{1} =  iimg_indx2clusters(dat1, maskInfo1);
%     cl{2} =  iimg_indx2clusters(dat2, maskInfo2);
%     cl{3} =  iimg_indx2clusters(dat3, maskInfo3);
%     
%     cl{4} =  iimg_indx2clusters(pos12, maskInfo1);
%     cl{5} =  iimg_indx2clusters(pos13, maskInfo2);
%     cl{6} =  iimg_indx2clusters(pos23, maskInfo3);
%     
%     cl{7} =  iimg_indx2clusters(pos123, maskInfo1);
    
else
    % get data from indicators
    
    dat = masks.dat;
    if isa(masks, 'statistic_image') && size(masks.sig, 1) == size(masks.dat, 1)
        dat = double(dat) .* double(masks.sig);
    end
    
    if strcmp(signstr, 'positive')
        dat = dat > 0; % positive
    elseif strcmp(signstr, 'negative')
        dat = dat < 0; % negative
    end
    
end

disp('Getting overlap colors');

u = unique(dat, 'rows');
u(sum(u, 2) == 0, :) = [];

if isempty(u)
    disp('No valid voxels');
    return
end

dat_final = dat;

for i = 1:size(u, 1)
    % each unique combo of colors
    
    mycolors = colors(u(i, :));
    mycolors = cat(1, mycolors{:});
    
    [C,IA] = intersect(dat, u(i, :), 'rows');
    
    to_add = false(size(dat, 1), 1);
    to_add(IA) = 1;
    
    dat_final(:, end+1) = to_add;
    
    colors{end+1} = mean(mycolors, 1);
    
end

disp('Getting clusters');

for i = 1:size(dat_final, 2)
   
    cl{i} = iimg_indx2clusters(double(dat_final(:, i)), masks.volInfo);
    
end

disp('Displaying clusters');


cluster_orthviews();
for i = 1:size(dat_final, 2)
    if ~isempty(cl{i}), cluster_orthviews(cl{i}, colors(i), 'add', 'solid'); end
end

if domontage
    
    cluster_orthviews_montage(8, 'axial', [], 'onerow');
    h = findobj(gcf,'Type', 'Text'); delete(h);
    
    cluster_orthviews_montage(8, 'sagittal', [], 'onerow');
    h = findobj(gcf,'Type', 'Text'); delete(h);
    
end


if dosurface
    
    create_figure;
    s = addbrain('hires left'); set(s, 'FaceColor', [.5 .5 .5], 'FaceAlpha', 1);
    
    for i = 1:size(dat_final, 2)
        cluster_surf(cl{i}, colors(i), 2, s);
    end
    
    scn_export_papersetup; saveas(gcf, 'cluster_overlap2_Lmedial.png');
    view(270, 0); lightFollowView
    scn_export_papersetup; saveas(gcf, 'cluster_overlap2_Lateral.png');
    
    create_figure('s2');
    s = addbrain('hires right'); set(s, 'FaceColor', [.5 .5 .5], 'FaceAlpha', 1);
    
    for i = 1:size(dat_final, 2)
        cluster_surf(cl{i}, colors(i), 2, s);
    end
    
    lightFollowView
    scn_export_papersetup; saveas(gcf, 'cluster_overlap_3colors_Rmedial.png');
    view(90, 0); lightFollowView
    scn_export_papersetup; saveas(gcf, 'cluster_overlap_3colors_Rateral.png');
    
    create_figure('s3');
    s = addbrain('limbic');
    s = [s addbrain('brainstem')];
    axis auto
    set(s, 'FaceColor', [.5 .5 .5]);
    set(s(end-1), 'FaceAlpha', .15);
    view(135, 10)

    for i = 1:size(dat_final, 2)
        cluster_surf(cl{i}, colors(i), 2, s);
    end
    
    scn_export_papersetup; saveas(gcf, 'cluster_overlap_3colors_Lmedial.png');
    
end

end
