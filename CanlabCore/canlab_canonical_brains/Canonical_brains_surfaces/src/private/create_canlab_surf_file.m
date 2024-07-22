function create_canlab_surf_file(gii_path, basename)
    gii = gifti(gii_path);
    faces = gii.faces;
    vertices = gii.vertices;
    % for grid interpolation to surfaces we need to take the affine
    % translations into account, but for spherical rotations we don't want
    % to do this since we need spheres to be aligned for surface based
    % interpolation between differently tesellated spheres!
    mat = gii.mat;
    %mat = eye(4);
    vertices = mat*[vertices,ones(size(vertices,1),1)]';
    vertices = vertices(1:3,:)';
    save(sprintf('/home/bogdan/.matlab/canlab/CanlabCore/CanlabCore/canlab_canonical_brains/Canonical_brains_surfaces/%s',basename),'faces','vertices','mat');
end