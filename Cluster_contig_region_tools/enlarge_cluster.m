% function cl = enlarge(cl)
%
% enlarge cluster

function cl = enlarge_cluster(cl)

    XYZ = cl.XYZ;
    XYZnew = XYZ;

    for i = 1:3
        Xtmp = XYZ;
        Xtmp(i, :) = Xtmp(i, :) + 1;
        XYZnew = [XYZnew Xtmp];         % add to dim
        Xtmp(i, :) = Xtmp(i, :) - 2;
        XYZnew = [XYZnew Xtmp];      % subtract 1 from dim
    end

    XYZnew = unique(round(XYZnew)', 'rows')';
    cl.XYZ = XYZnew;
    cl.XYZmm = voxel2mm(XYZnew, cl.M);
    cl.Z = ones(1, size(cl.XYZ, 2));

    cl.numVox = size(cl.XYZmm, 2);
end
