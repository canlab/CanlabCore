function cl = cluster_set_intersection(varargin)
% Computes the intersection of the sets of clusters passed in.
%
% :Usage:
% ::
%
%    intersect_cl = cluster_set_intersection(cls1, cls2, cls3, ...)
%    cl = cluster_intersection(robust0001_poscls, robust0002_poscls, robust0003_poscls);
%
% Note:
%    Designed for sets of clusters. To compute the intersection between individual clusters,
%    use cluster_intersection(). cluster_set_intersection will work, but is not needed.
%
% ..
%    Created by Matthew Davidson, 07/08/24
% ..

    cluster_sets = {};
    other_args = {};
    for i=1:length(varargin)
        switch(class(varargin{i}))
            case 'struct'
                cluster_sets{end+1} = clusters2CLU(varargin{i});
            otherwise
                other_args{end+1} = varargin{i};
        end
    end
    if(length(cluster_sets) < 2)
        error('Need 2 or more cluster sets to compute intersection.');
    end

    cl = cluster_intersection(cluster_sets{:}, other_args{:});
end
