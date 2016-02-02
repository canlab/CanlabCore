function surf_handles = cluster_surf_batch2(varargin)
% :Usage:
% ::
%
%    surf_handles = cluster_surf_batch(varargin)
%
% :Examples:
% ::
%
%    % Single-map visualization
%    P2 = threshold_imgs('rob_tmap_0001.img', tinv(1-.005, 12), 15, 'pos');
%    cluster_surf_batch(P2);
%    surf_handles = cluster_surf_batch({cl cl2}, {[1 0 0]});
%
%    % Two maps with overlap
%    surf_handles = cluster_surf_batch(cl, {[1 0 0] [0 1 0] [1 1 0]}, cl2);


    OVERLAP_COLOR = [1 1 1];
    DEFAULT_DISTANCE = 5;

    dooverlap = 0;
    overlap_color = OVERLAP_COLOR;
    regions = {[], 'brainstem', 'limbic', 'left', 'right'};

    for i=1:length(varargin)
        if(ischar(varargin{i}))
            switch(varargin{i})
                case 'clusters'
                    cluster_input = varargin{i+1};
                    if(iscell(cluster_input) && isstruct(cluster_input{1}))
                        cl_groups = cluster_input;
                    elseif ischar(cluster_input)
                        cl_groups = mask2clusters(cluster_input);
                    else
                        error('You must input either an image name, or a cell array containing cluster structure arrays in each cell for the clusters.');
                    end
                case 'colors'
                    colors = varargin{i+1};
                case 'overlap'
                    dooverlap = varargin{i+1};
                case {'overlap color', 'overlapcolor'}
                    overlap_color = varargin{i+1};
                case 'regions'
                    regions = varargin{i+1};
            end
        end
    end

    if(~exist('cl_groups', 'var') || isempty(cl_groups))
        disp('Choose thresholded image to get clusters from.')
        imgname = spm_get(1);
        cl_groups = {mask2clusters(imgname)};
    end

    if(~exist('colors', 'var') || isempty(colors))
        colors = {[1 1 0] [1 .5 0] [1 .2 0] [0 0 1] [0 .2 1] [0 .5 1]};
    end
    if(length(colors) < length(cl_groups))
        for i=length(colors):length(cl_groups)
            colors{i} = rand(1,3);
        end
    end

    local_print_colors(colors, cl_groups, dooverlap, overlap_color);

    surf_handles = [];
    for i=1:length(regions)
        surf_handles(end+1) = local_display_cluster_surf(regions{i}, cl_groups, colors, dooverlap, overlap_color, DEFAULT_DISTANCE);

        if(isempty(regions{i}))
            camzoom(1.3);
            make_figure_into_orthviews();
            view(180, -90);
        else
            switch(regions{i})
                case 'bg'
                    lightRestoreSingle(gca);
                case 'brainstem'
                    scn_export_papersetup(600);
                case 'limbic'
                    scn_export_papersetup(600);
            end
        end
        [az, el]=view;
        h = lightangle(az, el);
    end
end




%---------------------------
% Local functions
%---------------------------
function hSurfFig = local_display_cluster_surf(region, cl_groups, colors, dooverlap, overlap_color, mm_depth)
    if(dooverlap)
        hSurfFig = cluster_surf(cl_groups{1}, cl_groups{2}, mm_depth, region, [colors(1:2) overlap_color]);

        if(~isempty(cl_groups{3}) && isempty(cl_groups{4}))
            hSurfFig = cluster_surf(hSurfFig, cl_groups{3}, mm_depth, colors(3));
        elseif(~isempty(cl_groups{3}) && ~isempty(cl_groups{4}))
            hSurfFig = cluster_surf(hSurfFig, cl_groups{3}, cl_groups{4}, mm_depth, [colors(3:4) overlap_color]);
        end

        if(~isempty(cl_groups{5}) && isempty(cl_groups{6}))
            hSurfFig = cluster_surf(hSurfFig, cl_groups{5}, mm_depth, colors(5));
        elseif(~isempty(cl_groups{5}) && ~isempty(cl_groups{6}))
            hSurfFig = cluster_surf(hSurfFig, cl_groups{5}, cl_groups{6}, mm_depth, [colors(5:6) overlap_color]);
        end
    else
        hSurfFig = cluster_surf(cl_groups{1}, mm_depth, region, colors(1));
        for i=2:length(cl_groups)
            cluster_surf(hSurfFig, cl_groups{i}, mm_depth, colors(i));
        end
    end

    set(gcf, 'Color', 'w');
    axis off;
end

function local_print_colors(colors, cl_groups, dooverlap, overlap_color)
    for i=1:length(cl_groups)
        fprintf(1, 'Cluster set %d colors: %3.2f %3.2f %3.2f\n', i, colors{i});
    end
    if(dooverlap)
        fprintf(1, 'Overlap color: %3.2f %3.2f %3.2f\n', overlap_color);
    end
    fprintf(1,'\n')
end
