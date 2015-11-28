function [xyz, XYZmm, Z, class] = cluster_local_maxima(cl, dthresh, verbose)
% Clusters are chosen so that they must be at least dthresh mm apart
% default is 10 mm
%
% :Usage:
% ::
%
%    [xyz, XYZmm, Z, class] = cluster_local_maxima(cl, [dthresh], [verbose])
%
% verbose output: 1/0, default is 0
%
% additional optional outputs (slower):
%
% class: vector of integers for which subcluster this cluster belongs to
%
% ..
%    tor wager
% ..

    if nargin < 2 || isempty(dthresh), dthresh = 10; end
    if nargin < 3, verbose = 0; end

    [N Z xyz] = spm_max(abs(cl.Z), cl.XYZ);

    XYZmm = voxel2mm(xyz, cl.M);


    ncoord = size(cl.XYZmm, 2);
    class = zeros(1, ncoord);

    if isempty(xyz), return, end   % no peak maximum

    if verbose, fprintf(1, 'Local maxima (initial): %3.0f\n', length(Z)); end

    % find maxima within 10 mm and collapse
    d = pdist(XYZmm');
    nd = length(d);
    tooclose = d < dthresh;

    while any(tooclose)

        [d, nd, tooclose, xyz, XYZmm, Z] = omit_lowz_of_closest_pair(d, nd, tooclose, xyz, XYZmm, Z, dthresh, verbose);

    end

    if verbose, fprintf(1, 'Local maxima (final): %3.0f\n', length(Z)); end

    if nargout > 3
        % assign each voxel a subcluster based on closest local max
        speaks = XYZmm';
        npeaks = length(Z);

        % this is slower even for only 1000 vox
        %d = squareform(pdist([speaks; cl.XYZmm']));
        %d = d(1:npeaks, npeaks+1:end);   % peaks x voxels

        for i = 1:ncoord
            d = distance(cl.XYZmm(:, i)', speaks);
            [mind, whclose] = min(d);
            class(i) = whclose(1);
        end
    end
end



function [d, nd, tooclose, xyz, XYZmm, Z] = omit_lowz_of_closest_pair(d, nd, tooclose, xyz, XYZmm, Z, dthresh, verbose)

    tooclose = tooclose .* d;

    % find minimum distance
    tooclose(tooclose == 0) = Inf;
    [mindist, wh] = min(tooclose);

    % get indices of two local maxima that are too close (in rows)
    closest = zeros(1, nd);
    closest(wh(1)) = mindist;
    closest = squareform(closest);
    [rows, cols] = find(closest);    % need to return cols to get correct behavior

    % find which has the lowest Z-score
    [lowZ, whz] = min(Z(rows));      % rows is sufficient b/c there is only one pair

    wh_to_omit = rows(whz); % which in list to omit

    if verbose
        fprintf(1, 'Omitting index: %3.0f with z-value of %3.2f, which is %3.2f mm from another max whose z is %3.2f\n', ...
            wh_to_omit, lowZ, mindist, max(Z(rows)));
    end

    % eliminate important variables
    xyz(:, wh_to_omit) = [];
    XYZmm(:, wh_to_omit) = [];
    Z(wh_to_omit) = [];

    % find maxima within 10 mm and collapse
    d = pdist(XYZmm');
    nd = length(d);
    tooclose = d < dthresh;

end


function d = distance(u, v)
% d = distance(u, v)
% Euclidean distance between u and v
% u and v should be column vectors or matrices in k-dim space, where k is
% the number of columns of u and v.
%
% if u is a single row, replicates u to size(v,1)

if size(u,1) < size(v,1)
    u = repmat(u,size(v,1),1);
    if size(u,1) < size(v,1), error('Inputs must have same number of rows, or one row for 1st input!'),end
end

delt = u - v;
d = (sum(delt .^ 2,2)) .^.5;

end


