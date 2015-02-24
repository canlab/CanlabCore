function [cq mcq Xc]=getmeanquality(X,c,linkagetype)

    %point estimate cluster solution
    Y = pdist1(X);

    if exist('linkage.m', 'file') && exist('cluster.m', 'file')
        Z = linkage(Y,linkagetype);
        Xc = cluster(Z,c);
    else
        Z = linkage_t(Y,linkagetype);
        Xc = cluster_t(Z,c);
    end

    %convert X to X, binary input to doquality
    Xcx=makebinary(Xc)';

    %point estimate quality of each element;
    [cq center]=clustquality(Xcx,X);      %cq is cluster quality of each element

    %mean cq
    if all(cq == 1)
        mcq = 0;
    else
        mcq = mean( cq(cq~=1) );    %don't include single element clusters in your quality estimate
    end

end