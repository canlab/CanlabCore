function subcl = subclusters_from_local_max(cl, dist_thresh)
% Breaks apart a cluster into smaller clusters
%
% :Usage:
% ::
%
%    subcl = subclusters_from_local_max(cl, dist_thresh)

    if isempty(cl), return; end
    subcl = [];

    if ~isfield(cl, 'dim'), cl(1).dim = []; end
    N = fieldnames(cl(1));

    for i = 1:length(cl)
        [xyz, XYZmm, Z, class] = cluster_local_maxima(cl(i), dist_thresh, 0);

        if all(class == 0)  % no peak maxima; just use the existing cluster
            %this_subcl = cl;
            % make same as subclusters by saving only relevant fields

            create_subcl;
            this_subcl = add_fields(N, this_subcl, cl, i);
        else
            this_subcl = cluster2subclusters(cl(i), class);
            this_subcl = add_fields(N, this_subcl, cl, i);
        end

        % keep track of which larger cluster this came from; for cluster_table
        for j = 1:length(this_subcl)
            this_subcl(j).from_cluster = i;
        end

        try
            subcl = [subcl this_subcl];
        catch
            warning('internal error with subcluster consistency. this should work fine, but fix code...');
            subcl = merge_clusters(subcl, this_subcl);
        end
    end

    
    
    function create_subcl()
        
        if ~isfield(cl(i), 'name'), cl(i).name = []; end
        
        % N is field names, cl = clusters, i = which cluster
        this_subcl = struct('title', cl(i).title, 'threshold', cl(i).threshold, 'M', cl(i).M, 'dim', cl(i).dim, 'voxSize', cl(i).voxSize, ...
            'name', cl(i).name, 'numVox', cl(i).numVox, 'XYZ', cl(i).XYZ, 'XYZmm', cl(i).XYZmm, 'Z', cl(i).Z, ...
            'mm_center', cl(i).mm_center, 'from_cluster', i);
    end
end


function this_subcl = add_fields(N, this_subcl, cl, i)
    for j = 1:length(N)
        for k = 1:length(this_subcl)
            if ~isfield(this_subcl(k), N{j}) || isempty(this_subcl(k).(N{j}))
                this_subcl(k).(N{j}) = cl(i).(N{j});
            end
        end
    end
end
