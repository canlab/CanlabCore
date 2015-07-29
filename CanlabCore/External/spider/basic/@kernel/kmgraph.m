function K = kmgraph(kern,d1,d2,ind1,ind2,kerparam),
% Marginalized Graph Kernel by Koji Tsuda
% The actual implementation uses the \delta version of the kernel (see publication)
% and a tuned C implementation by Alex Zien.
% Despite the paper we are using a normed K(X1,X2)/sqrt(K(X1,X1)*K(X2,X2))
% version.
% To use this kernel ensure that the function graphkerneldelta.cpp
% is compiled in the 'functions' folder.
% Furthermore note that d1.X must contain !cells! of the form :
% Let x,y,z be a graph, then 
% struct(x)
%     adjacency: [NxN double]
%     vertices: [Nx1 double]
% and d=data({x,y,z}') is a valid data object
%
% gb 03-Mar-2004


X1=get_x(d1,ind1);
X2=get_x(d2,ind2);

K=zeros(length(X2),length(X1));
eps=1e-5;

lambda = kerparam;

if( size(ind1)==size(ind2))
    if( ind1==ind2)
        for i=ind1
            for j=ind2
                Ks(j,i)=graphkerneldelta(X1{i}.adjacency,X2{j}.adjacency,X1{i}.vertices,X2{j}.vertices,lambda,eps);

            end
        end
        
        for i=ind1
            for j=ind2
                K(j,i)=Ks(j,i)/sqrt(Ks(j,j)*Ks(i,i));
            end
        end
    end
else
    
    for i=ind1
        for j=ind2
            k12=graphkerneldelta(X1{i}.adjacency,X2{j}.adjacency,X1{i}.vertices,X2{j}.vertices,lambda,eps);
            k11=graphkerneldelta(X1{i}.adjacency,X1{i}.adjacency,X1{i}.vertices,X1{i}.vertices,lambda,eps);
            k22=graphkerneldelta(X2{j}.adjacency,X2{j}.adjacency,X2{j}.vertices,X2{j}.vertices,lambda,eps);
            K(j,i)=k12/sqrt(1e-5+k11*k22);            
        end
    end
    
    
    
end
