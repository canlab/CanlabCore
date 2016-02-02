function cluster_orthviews_overlap_3colors(mask1, mask2, mask3, varargin)
% Plot three clusters on the orthviews, and their intersections for
% activations and deactivations in intermediate colors
%
% :Usage:
% ::
%
%    cluster_orthviews_overlap_3colors(mask1, mask2, mask3, ['colors', colors cell], ['surface'] )
%
% Positive effects only!
%
% :Inputs:
%
%   mask1 mask2 mask3 are either image names (preferred!) or clusters
%   structures.  for clusters structures, need to add .dim field
%   and you need the 2010 object-oriented code in the canlab repository.
%
% ..
%    tor wager, june 2010
%    A good function, but not polished yet.  needs debugging for various
%    minor issues.
% ..

colors = {[0 0 1] [1 0 0] [1 1 0]};
dosurface = 0;

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            % reserved keywords
            case 'betas'
            case 'design'
                
                % functional commands
            case 'colors', colors = varargin{i+1};
            case 'surface', dosurface = 1;
                
            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end

if isstruct(mask1) && isstruct(mask2) && isstruct(mask3)
    loadfromfile = 0;
elseif ischar(mask1) && ischar(mask2) && ischar(mask3)
    loadfromfile = 1;
else
    error('Inputs must be all cl structs or mask filenames');
end

if loadfromfile
    
    [maskInfo1, dat1] = iimg_read_img(mask1, 2);
    [maskInfo2, dat2] = iimg_read_img(mask2, 2);
    [maskInfo3, dat3] = iimg_read_img(mask3, 2);
    
    % positive only
    dat1 =  dat1 > 0;
    dat2 =  dat2 > 0;
    dat3 =  dat3 > 0;
    
    int_dat = iimg_intersection(mask1, mask2, 'posneg'); % intersections
    pos12 = int_dat(:, 1);
    
    int_dat = iimg_intersection(mask1, mask3, 'posneg'); % intersections
    pos13 = int_dat(:, 1);
    
    int_dat = iimg_intersection(mask2, mask3, 'posneg'); % intersections
    pos23 = int_dat(:, 1);
    
    pos123 = pos12 & pos13;
    
    intersect_all = any([pos12 pos13 pos23 pos123], 2);
    
    %remove higher-order intersections from  mask dat
    dat1(intersect_all) = 0;
    dat2(intersect_all) = 0;
    dat3(intersect_all) = 0;
    
    pos12(pos123) = 0;
    pos23(pos123) = 0;
    pos13(pos123) = 0;
    
    % make clusters
    clear cl
    
    cl{1} =  iimg_indx2clusters(dat1, maskInfo1);
    cl{2} =  iimg_indx2clusters(dat2, maskInfo2);
    cl{3} =  iimg_indx2clusters(dat3, maskInfo3);
    
    cl{4} =  iimg_indx2clusters(pos12, maskInfo1);
    cl{5} =  iimg_indx2clusters(pos13, maskInfo2);
    cl{6} =  iimg_indx2clusters(pos23, maskInfo3);
    
    cl{7} =  iimg_indx2clusters(pos123, maskInfo1);
    
else
    % Work from clusters directly
    
    % convert to masks
    
    if ~isfield(mask1, 'dim') || isempty(mask1(1).dim) || ...
            ~isfield(mask2, 'dim') || isempty(mask2(1).dim) || ...
            ~isfield(mask3, 'dim') || isempty(mask3(1).dim)
        
        error('Enter .dim (image dimension in voxels) field in all cl structs. e.g., [91 109 91]');
    end
    
    r = cluster2region(mask1);
    m1 = region2imagevec(r);
    
    r = cluster2region(mask2);
    m2 = region2imagevec(r);
    
    r = cluster2region(mask3);
    m3 = region2imagevec(r);
    
    m1 = replace_empty(m1);
    m2 = replace_empty(m2);
    m3 = replace_empty(m3);
    
    dat1 =  m1.dat > 0;
    dat2 =  m2.dat > 0;
    dat3 =  m3.dat > 0;
    
    pos123 = all([dat1 dat2 dat3], 2);
    pos12 = all([dat1 dat2], 2);
    pos23 = all([dat2 dat3], 2);
    pos13 = all([dat1 dat3], 2);
    
    intersect_all = any([pos12 pos13 pos23 pos123], 2);
    
    %remove higher-order intersections from  mask dat
    dat1(intersect_all) = 0;
    dat2(intersect_all) = 0;
    dat3(intersect_all) = 0;
    
    pos12(pos123) = 0;
    pos23(pos123) = 0;
    pos13(pos123) = 0;
    
    % make clusters
    clear cl
    
    cl{1} =  iimg_indx2clusters(dat1, m1.volInfo);
    cl{2} =  iimg_indx2clusters(dat2, m2.volInfo);
    cl{3} =  iimg_indx2clusters(dat3, m3.volInfo);
    
    cl{4} =  iimg_indx2clusters(pos12, m1.volInfo);
    cl{5} =  iimg_indx2clusters(pos13, m1.volInfo);
    cl{6} =  iimg_indx2clusters(pos23, m1.volInfo);
    
    cl{7} =  iimg_indx2clusters(pos123, m1.volInfo);
    
    
    
end


% set colors for overlap
% colors{end+1} = {mean([colors{1}; colors{2}])};
% colors{end+1} = {mean([colors{1}; colors{3}])};
% colors{end+1} = {mean([colors{2}; colors{3}])};
% 
% colors{end+1} = {mean([colors{1}; colors{2}; colors{3}])};

colors(end+1) = {2*mean([colors{1}; colors{2}])}; colors{end}(colors{end} > 1) = 1;
colors(end+1) = {2*mean([colors{1}; colors{3}])}; colors{end}(colors{end} > 1) = 1;
colors(end+1) = {2*mean([colors{2}; colors{3}])}; colors{end}(colors{end} > 1) = 1;

colors(end+1) = {3*mean([colors{1}; colors{2}; colors{3}])};


cluster_orthviews();
for i = 1:7
    if ~isempty(cl{i}), cluster_orthviews(cl{i}, colors(i), 'add', 'solid'); end
end

cluster_orthviews_montage(10, 'axial', [], 'onerow', 'range', [-35 55]);
h = findobj(gcf,'Type', 'Text'); set(h, 'FontSize', 18) %delete(h);

cluster_orthviews_montage(10, 'sagittal', [], 'onerow', 'range', [-40 40]);
h = findobj(gcf,'Type', 'Text'); set(h, 'FontSize', 18) %delete(h);


if dosurface
    
    create_figure;
    s = addbrain('hires left'); set(s, 'FaceColor', [.5 .5 .5], 'FaceAlpha', 1);
    for i = 1:7
        cluster_surf(cl{i}, colors(i), 2, s);
    end
    for i = 4:7
        cluster_surf(cl{i}, colors{i}(1), 2, s);
    end
    
    scn_export_papersetup; saveas(gcf, 'cluster_overlap_3colors_Lmedial.png');
    view(270, 0); lightFollowView
    scn_export_papersetup; saveas(gcf, 'cluster_overlap_3colors_Lateral.png');
    
    create_figure('s2');
    s = addbrain('hires right'); set(s, 'FaceColor', [.5 .5 .5], 'FaceAlpha', 1);
    for i = 1:3
        cluster_surf(cl{i}, colors(i), 2, s);
    end
    for i = 4:7
        cluster_surf(cl{i}, colors{i}(1), 2, s);
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
    for i = 1:3
        cluster_surf(cl{i}, colors(i), 2, s);
    end
    for i = 4:7
        cluster_surf(cl{i}, colors{i}(1), 2, s);
    end
    scn_export_papersetup; saveas(gcf, 'cluster_overlap_3colors_Lmedial.png');
    
    
    
end

end
