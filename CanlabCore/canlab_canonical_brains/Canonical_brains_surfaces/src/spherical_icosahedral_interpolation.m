function [Cq, coord, vec] = spherical_icosahedral_interpolation(V,F,C,Vq,varargin)
%   My method for selecting the nearest face F to a vertex Vq for
%   barycentric interpolation is also potentially wrong. It doesn't work in
%   all instances, and it should, so I've probably mistaken something along
%   the way, but the results look good enough for government work.

    % V - vertices of source space
    % F - faces of source space
    % C - data in source psace
    % Vq - vertices of target space
    % Fq - faces of target space
    %
    % Optional flag: 'nearest'  - performs nearest neighbor interpolation instead of 
    %      barycentric linear interpolation. Unlike for barycentric interpolation I have
    %      no reservations about the accuracy of my implementation here. It's a trivial
    %      application of built in matlab functions.
    %
    % Cq - data in target space

    dist_fun = @(x1,x2)(sqrt(sum((x1 - x2).^2,2)));
    
    % make sure we're dealing with spheres
    tol = 1e-4;
    r = sqrt(sum(double(V).^2,2));
    rq = sqrt(sum(double(V).^2,2));
    if std(r)/mean(r) > tol
        % if V is not origin centric, center it
        for i = 1:3
            d = dist_fun(V, [0,0,0]);
            furthest = V(d == max(d),:);
    
            r = mean(dist_fun(V, mean(V,1))); % average distance
            V = V - V(d == max(d),:) + mean(r)*furthest/norm(furthest);
        end
    end
    if std(r)/mean(r) > tol
        for i = 1:3
            % if Vq is not origin centric, center it
            d = dist_fun(Vq, [0,0,0]);
            furthest = Vq(d == max(d),:);
    
            rq = mean(dist_fun(Vq, mean(Vq,1))); % average distance
            Vq = Vq - Vq(d == max(d),:) + mean(rq)*furthest/norm(furthest);
        end
    end
    assert(mean(r) - mean(rq) < tol, 'V and Vq do not have the same radius (within tolerance). This function only works for spherical meshes with the same radius.')
    r = mean(r);

    TR = triangulation(double(F),double(V));
    nn = TR.nearestNeighbor(Vq);

    if contains(varargin,'nearest')
        Cq = C(nn);
        coord = ones(length(nn),1);
        vec = nn;

        return
    end

    % our general strategy here involves 3 steps
    % 1) identify nearest neighbor of Vq in V
    % 2) get first-fourth order faces connected to nearest neighbor
    % 3) project Vq to planes containing faces of (2) to define projVq (one for each face)  [edit: commentd out]
    % 4) convert projVq into barycentric coordinates and find face that contains projVq  [edit: just using Vq directly now, not sure which is better]
    % 5) perform barycentric interpolation of projVq based on face form (4)
    Cq = zeros(size(Vq,1),1);
    coord = zeros(size(Vq,1),3);
    vec = zeros(size(Vq,1),3);
    for i = 1:length(Vq)
        first_F_ind = F(cell2mat(TR.vertexAttachments(nn(i))),:);
        second_F_ind = F(unique(cell2mat(TR.vertexAttachments(double(unique(first_F_ind(:))))')),:);
        second_F_ind = F(unique(cell2mat(TR.vertexAttachments(double(unique(second_F_ind(:))))')),:);
        second_F_ind = unique(cell2mat(TR.vertexAttachments(double(unique(second_F_ind(:))))'));



        %{
        %% define planes containing faces
        % normals + radius define the planes of each triangle
        face_normals = TR.faceNormal(double(second_F_ind'));
        % compute offset for each face (~100, but not exactly due to
        % imprecise datatypes
        r = V(F(second_F_ind,:),:)*face_normals';
        r = diag(r'*kron(eye(length(second_F_ind)),1/3*ones(3,1)));

        %% project Vq to planes defined above
        % change basis so Vq(i,:) is origin
        rq = r - face_normals*Vq(i,:)';

        % compute closest point p to Vq(i,:) on each plane
        p = face_normals.*rq;

        % change back to original basis
        projVq = p + Vq(i,:);
        %}
        
        %% identify face containing one of the Vq
        %Bq = cartesianToBarycentric(TR, second_F_ind(:), projVq);
        Bq = cartesianToBarycentric(TR, second_F_ind(:), repmat(Vq(i,:),length(second_F_ind),1));
        enclosing_tri_Bq_ind = find(all(Bq >= 0 & Bq <= 1,2));
        % the triangles won't always inclose points due to numerical
        % interpolation errors, so instead we find nearest triangles
        if isempty(enclosing_tri_Bq_ind)
            enclosing_tri_Bq_ind = find(all(abs(Bq) < 1,2));
        end

        % if there are multiple solutions, find the triangle with all the
        % smallest vertex distance from the original point. This should
        % only happey when a point falls on an edge, in which case the
        % candidates will only differ in their third point, and we use
        % whichever has the closer 3rd point.
        mean_dist = @(x1)(mean(sum((V(F(second_F_ind(x1),:),:) - Vq(i,:)).^2,2)));
        dist = arrayfun(@(x1)mean_dist(x1), enclosing_tri_Bq_ind);
        enclosing_tri_Bq_ind = enclosing_tri_Bq_ind(min(dist) == dist);

        %% get barycentric coordinates of projVq
        neighborhood_V_ind = F(second_F_ind(enclosing_tri_Bq_ind),:);
        Bq = Bq(enclosing_tri_Bq_ind,:);

        %% perform barycentric interpolation
        try
            coord(i,:) = Bq;
            vec(i,:) = neighborhood_V_ind;
            Cq(i) = Bq*C(neighborhood_V_ind);
        catch
            warning('Could not identify face containing vertex %d. Using nearest neighbor Vq.',i);
            Cq(i) = C(nn(i));
            coord(i,:) = [1,0,0];
            vec(i,:) = [nn(i), 1, 1]; % the other two don't matter, we give them zero weight
            %{
                keyboard
                % this code will plot the candidate triangles and your target point to inspect
                % why we're getting an error
                figure;
                cla
                for j = 1:length(second_F_ind)
                    this_F = F(second_F_ind(j),:);
                    these_V = V(this_F,:);
                    these_V = [these_V; these_V(1,:)];
                    plot3(these_V(:,1), these_V(:,2), these_V(:,3),'-');
                    hold on;
                end
                plot3(Vq(i,1),Vq(i,2),Vq(i,3),'+')
            %}
        end
    end
end
