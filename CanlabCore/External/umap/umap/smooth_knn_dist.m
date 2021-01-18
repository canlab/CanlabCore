function [knn_dist, nn_dist] = smooth_knn_dist(distances, k, varargin)
%SMOOTH_KNN_DIST Compute a continuous version of the distance to the kth
% nearest neighbor. That is, this is similar to knn-distance but allows
% continuous k values rather than requiring an integral k. In esscence we
% are simply computing the distance such that the cardinality of fuzzy set
% we generate is k.
% 
% [knn_dist, nn_dist] = SMOOTH_KNN_DIST(distances, k, local_connectivity,
% n_iter, bandwidth)
% 
% Parameters
% ----------
% distances: array of size (n_samples, n_neighbors)
%     Distances to nearest neighbors for each samples. Each row should be a
%     sorted list of distances to a given sample's nearest neighbors.
% 
% k: double
%     The number of nearest neighbors to approximate for.
% 
% n_iter: double (optional, default 64)
%     We need to binary search for the correct distance value. This is the
%     max number of iterations to use in such a search.
% 
% local_connectivity: double (optional, default 1)
%     The local connectivity required -- i.e. the number of nearest
%     neighbors that should be assumed to be connected at a local level.
%     The higher this value the more connected the manifold becomes
%     locally. In practice this should be not more than the local intrinsic
%     dimension of the manifold.
% 
% bandwidth: double (optional, default 1)
%     The target bandwidth of the kernel, larger values will produce
%     larger return values.
% 
% Returns
% -------
% knn_dist: array of size (n_samples, 1)
%     The distance to kth nearest neighbor, as suitably approximated.
% 
% nn_dist: array of size (n_samples, 1)
%     The distance to the 1st nearest neighbor for each point.
%
%   AUTHORSHIP
%   Math Lead & Primary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Secondary Developer: Stephen Meehan <swmeehan@stanford.edu>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Provided by the Herzenberg Lab at Stanford University 
%   License: BSD 3 clause
%
    p=parseArguments();
    parse(p,varargin{:});
    args=p.Results;
    local_connectivity = args.local_connectivity;
    n_iter = args.n_iter;
    bandwidth = args.bandwidth;
    same_set = args.same_set;
    
    SMOOTH_K_TOLERANCE = 1e-5;
    MIN_K_DIST_SCALE = 1e-3;
    
    height = size(distances,1);
    
    target = log2(k)*bandwidth;

    lo = zeros(height,1);
    hi = Inf(height,1);
    mid = ones(height,1);

    zero_dists = sum(distances == 0, 2);
    
    if any(zero_dists==size(distances,2))
        warning('There are at least n_neighbors identical data points in the raw data. Results may be inaccurate.');
    end
        
    index = floor(local_connectivity);
    interpolation = local_connectivity - index;
    aug_dists = [lo distances repmat(max(distances, [], 2), [1 index+1])];
            
    idx = sub2ind(size(aug_dists), (1:height)', zero_dists + index + 1);
    rho = aug_dists(idx) + interpolation*(aug_dists(idx) - aug_dists(idx+height));
    
    if same_set
        d = distances(:,2:end) - rho;
    else
        d = distances(:,1:end) - rho;
    end

    for n = 1:n_iter
        
        summands = exp(-max(0, d./repmat(mid, [1 size(d, 2)])));

        psum = sum(summands, 2);
        
        if all(abs(psum - target) < SMOOTH_K_TOLERANCE)
            break
        end
        
        b = ~(psum > target).*hi;
        b(isnan(b)) = 0;
        
        hi = (psum > target).*mid + b;
        lo = (psum > target).*lo + ~(psum > target).*mid;
        
        c = (psum > target).*((lo + hi)/2);
        c(isnan(c)) = 0;
        
        mid = c + ~(psum > target).*min(2*lo, (lo + hi)/2);

    end

    result = mid;
    
    result = max(result, MIN_K_DIST_SCALE * ((rho > 0).*mean(distances,2) + (rho == 0).*mean(mean(distances))));

    knn_dist = result;
    nn_dist = rho;
end
    
    function p=parseArguments(varargin)
        p = inputParser;
        addParameter(p,'local_connectivity', 1);
        addParameter(p,'n_iter', 64);
        addParameter(p,'bandwidth', 1);
        addParameter(p,'same_set','true');
    end
