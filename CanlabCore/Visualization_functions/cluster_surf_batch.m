function surf_handles = cluster_surf_batch(varargin)
% :Usage:
% ::
%
%    surf_handles = cluster_surf_batch(varargin)
%
% :Examples:
% ::
%
%    % Single-map visualization
%    P2 = threshold_imgs('rob_tmap_0001.img',tinv(1-.005,12),15,'pos');
%    cluster_surf_batch(P2);
%    surf_handles = cluster_surf_batch(cl,{[1 0 0]},cl2);
%
%    % Two maps with overlap
%    surf_handles = cluster_surf_batch(cl,{[1 0 0] [0 1 0] [1 1 0]},cl2);


surf_handles = [];

% Input arguments
% cl, then colors, then cl2, in that order

if length(varargin) == 0
    disp('Choose thresholded image to get clusters from.')
    imgname = spm_get(1);
    cl = mask2clusters(imgname);

elseif isstruct(varargin{1})
    cl = varargin{1};

elseif isstr(varargin{1})
    cl = mask2clusters(varargin{1});

else
    error('You must input either nothing, an image name, or a cluster structure.');
end

morecl = {};
ovlstr = [];
mydistance = 5;
color = {[1 0 0]};

for i = 2:length(varargin)
    v = varargin{i};
    if iscell(v)
        color = v;
    elseif isstruct(v)
        morecl{end+1} = v;
    elseif ischar(v)
        ovlstr = v;     % could enter any string that goes into cluster_orthviews
    elseif length(v) == 1
        mydistance = v;
    end
end


% if length(varargin) > 1
%     color = varargin{2};
% end
%
% cl2 = [];
% if length(varargin) > 2
%     cl2 = varargin{3};
% end
% cl3 = [];
% if length(varargin) > 3
%     cl3 = varargin{4};
% end
% cl4 = [];
% if length(varargin) > 4
%     cl4 = varargin{5};
% end
while length(color) < length(morecl)+2, color{end+1} = rand(1,3); end

if strcmp(ovlstr,'nooverlap'), color = color(1:(length(morecl)+1)); end

fprintf(1,'Colors:\n');
for i=1:length(color)
    if i <= length(morecl) + 1
        fprintf(1,'Cluster set %3.0f\t%3.2f %3.2f %3.2f\n',i,color{i});
    else
        fprintf(1,'Overlap\t%3.2f %3.2f %3.2f\n',color{i});
    end
end
fprintf(1,'\n')


if strcmp(ovlstr,'nooverlap')
    surfhan = [];
    tor_fig;
     surfhan = make_surface([],cl,morecl,mydistance,color);
%    
     [hh1,hh2,hh3,hl,a1,a2,a3] = make_figure_into_orthviews;
     view(180,-90); [az,el]=view; h = lightangle(az,el);
%     
tor_fig;
 sh = make_surface('brainstem',cl,morecl,mydistance,color);
    surfhan = [surfhan sh];
    scn_export_papersetup(600)
  %  saveas(gcf,'brainstem','fig')
  %  saveas(gcf,'brainstem','png')
    
 tor_fig;
 sh = make_surface('limbic',cl,morecl,mydistance,color);
    surfhan = [surfhan sh];
        scn_export_papersetup(600)
   % saveas(gcf,'left_limbic','fig')
   % saveas(gcf,'left_limbic','png')
    
    sh = make_surface('left',cl,morecl,mydistance,color);
     surfhan = [surfhan sh];
% 
    sh = make_surface('right',cl,morecl,mydistance,color);
     surfhan = [surfhan sh];  
% 
%     sh = make_surface('bg',cl,morecl,mydistance,color);
%     surfhan = [surfhan sh];
    
surf_handles = surfhan;

else

    for i = 1:4
        if length(morecl) < i, morecl{i} = []; end
    end

    % Run surfaces

    surfhan = cluster_surf(cl,cl2,mydistance,color(1:3));
    if(~isempty(cl3) && ~isempty(cl4)), surfhan = cluster_surf(surfhan, cl3,cl4,mydistance,color(4:6));end
    set(gcf,'Color','w'); axis off; camzoom(1.3)

    [hh1,hh2,hh3,hl,a1,a2,a3] = make_figure_into_orthviews;
    view(180,-90); [az,el]=view; h = lightangle(az,el);

    sh2 = cluster_surf(cl,cl2,mydistance,'left',color(1:3));
    if(~isempty(cl3) && ~isempty(cl4)), sh2 = cluster_surf(sh2, cl3,cl4,mydistance,color(4:6));end
    [az,el]=view; h = lightangle(az,el);set(gcf,'Color','w'); axis off;

    sh3 = cluster_surf(cl,cl2,mydistance,'right',color(1:3));
    if(~isempty(cl3) && ~isempty(cl4)), sh3 = cluster_surf(sh3, cl3,cl4,mydistance,color(4:6));end
    [az,el]=view; h = lightangle(az,el);set(gcf,'Color','w'); axis off;

    sh4 = cluster_surf(cl,cl2,mydistance,'bg',color(1:3));
    if(~isempty(cl3) && ~isempty(cl4)), sh4 = cluster_surf(sh4, cl3,cl4,mydistance,color(4:6));end
    lightRestoreSingle(gca);
    [az,el]=view; h = lightangle(az,el);set(gcf,'Color','w'); axis off;

    surf_handles = [surfhan sh2 sh3 sh4];

end


return





 function surfhan = make_surface(typestr,cl,morecl,mydistance,color)
     
     if isempty(typestr), typestr = mydistance; end  % kludgy fix to avoid empty input
     
    surfhan = [];   %cluster_surf(cl,mydistance,color{1});
    for i = length(morecl):-1:1
        if i == length(morecl)
            surfhan = cluster_surf(surfhan,morecl{i},mydistance,color(i+1),typestr);
        else
            cluster_surf(surfhan,morecl{i},mydistance,color(i+1));
        end
        %surfhan = [surfhan sh];
    end
    
    surfhan = cluster_surf(surfhan,cl,mydistance,color(1));
    %surfhan = [surfhan sh];
    
    return
