function K = triangular(kern,dat1,dat2,ind1,ind2,kerParam),

% k_ij = -||x2_i-x1_j||_2


X2 = get_x(dat2,ind2);
X1 = get_x(dat1,ind1); 


K = [];

[szx1m, szx1n] =size(X1);
[szx2m, szx2n] =size(X2);



if issparse(X1) & issparse(X2)
    % more efficient way to compute the triangular kernel operates on the
    % indices directly
    
    %extract indices and values 
    [ix1,jx1,sx1] = find(X1); 
    [ix2,jx2,sx2] = find(X2); 
    Ix1 = [ix1,jx1,sx1];
    Ix2 = [ix2,jx2,sx2];
    
    % reshape X1 such that each row appears szx2m times consecutively
    [m1,n1] = size(Ix1); [m2,n2] = size(Ix2);
    Ix1 = reshape(repmat(reshape(Ix1,m1*n1,1),1,szx2m)',szx2m*m1,n1);
    Ix1(:,1) = szx2m*Ix1(:,1)  + reshape(repmat(-fliplr([0:(szx2m-1)]),m1,1)',m1*szx2m,1);
    X1 = sparse(Ix1(:,1),Ix1(:,2),Ix1(:,3),szx1m*szx2m,szx1n);
    
    % reproduce X2 szx1m times
    Ix2 = repmat(Ix2,szx1m,1);
    
    Ix2(:,1) = Ix2(:,1) +reshape(repmat(szx2m*[0:(szx1m-1)],m2,1), m2*szx1m,1);%; reshape(repmat(szx1m*[0:(m2-1)],szx1m,1), m2*szx1m,1);
    X2 = sparse(Ix2(:,1),Ix2(:,2),Ix2(:,3),szx1m*szx2m,szx2n);
    
    
    % calc and reshape kernel
    K = reshape(-sqrt(sum((X1-X2).^2,2)),szx2m,szx1m);
else    
    K = -reshape(sqrt(sum((reshape(repmat(reshape(X1,szx1m*szx1n,1),1,szx2m)',szx2m*szx1m,szx1n)-repmat(X2,szx1m,1)).^2,2)),szx2m,szx1m);
end

