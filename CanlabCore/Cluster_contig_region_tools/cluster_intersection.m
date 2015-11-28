function cl = cluster_intersection(varargin)
% Computes the intersection of the clusters passed in.
% 
% :Usage:
% ::
%
%     intersect_cl = cluster_intersection(cl1, cl2, cl3, ...)
%
%     % simple intersection
%     cl = cluster_intersection(robust0001_poscl(2), robust0002_poscl(4), robust0003_poscl(17));
%   
%     % intersection between sets of clusters
%     % alternatively, see cluster_set_intersection.m
%     cl = cluster_intersection(clusters2CLU(robust0001_poscl), clusters2CLU(robust0002_poscl));
%
% Note:
%   Only works with single clusters. To compute the intersection between sets of clusters,
%   use cluster_set_intersection()
%
% ..
%    Created by Matthew Davidson, 07/08/24
% ..

    clusters = [];
    parse_params();
    cl = blank_cluster(clusters);

    [cl.voxSize, cl.M] = intersection_cluster_dims(clusters);
    cl.XYZmm = intersection_voxels(clusters)';
    
    cl.XYZ = mm2voxel(cl.XYZmm, cl.M, 1)';
    cl.numVox = size(cl.XYZmm, 2);
    cl.mm_center = mean(cl.XYZmm, 2)';
    cl.center = mean(cl.XYZ, 2)';

    function parse_params()
        for i=1:length(varargin)
            switch(class(varargin{i}))
                case 'struct'
                    clusters = merge_clusters(clusters, varargin{i}); %#ok
                case 'char'
                    %nothing yet - placeholder for future
            end
        end
        if(length(clusters) < 2)
            error('Need 2 or more clusters to compute intersection.');
        end
    end
end



function [voxSize, M] = intersection_cluster_dims(clusters)
    max_cluster_size = prod(clusters(1).voxSize);
    largest_cluster_idx = 1;

    for i=2:length(clusters)
        if(prod(clusters(i).voxSize) > max_cluster_size)
            max_cluster_size = prod(clusters(i).voxSize);
            largest_cluster_idx = i;
        end
    end

    voxSize = clusters(largest_cluster_idx).voxSize;
    M = clusters(largest_cluster_idx).M;
end

% function voxels = intersection_voxels(clusters, voxels)
%     if(~exist('voxels', 'var'))
%         voxels = orient_coord_list(clusters(1).XYZmm);
%     else
%         voxels = intersect(orient_coord_list(clusters(1).XYZmm), voxels, 'rows');
%     end
%     if(length(clusters) > 1)
%         voxels = intersection_voxels(clusters(2:end), voxels);
%     end
% end
function voxels = intersection_voxels(clusters)
    voxels = orient_coord_list(clusters(1).XYZmm);
    for i=2:length(clusters)
        voxels = intersect(orient_coord_list(clusters(i).XYZmm), voxels, 'rows');
    end
end

function coords = orient_coord_list(coords)
    if(size(coords, 1) == 3)
        coords = coords';
    end
end

function cl = blank_cluster(clusters)
    cl.title = 'Intersection cluster';
    cl.name = 'Intersection cluster';
    cl.threshold = 1;
    cl.Z = [];
    cl.numPeaks = [];
    if isfield(clusters, 'P'), cl.P = strvcat(clusters.P); end
    if isfield(clusters, 'imP'),cl.imP = strvcat(clusters.imP); end
end
