function K = kmgraph_classical(kern,d1,d2,ind1,ind2,kerparam),
% Marginalized Graph Kernel by Koji Tsuda
% The actual implementation uses the \delta version of the kernel (see publication)
% and a tuned C implementation by Alex Zien.
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

verb=kern.algorithm.verbosity>=2;
if( size(ind1)==size(ind2))
    if( ind1==ind2)
        for i=ind1
            if(verb==1)
                fprintf('%d  ',i)
                if(mod(i,30)==0)
                    fprintf('\n');
                end
                
            end
            for j=ind2
                if(i>=j)
                    K(j,i)=graphkerneldelta(X1{i}.adjacency,X2{j}.adjacency,X1{i}.vertices,X2{j}.vertices,lambda,eps);
                    K(i,j)=K(j,i);
                end
            end
        end
    end
    
else
    for i=ind1
        for j=ind2
            
            K(j,i)=graphkerneldelta(X1{i}.adjacency,X2{j}.adjacency,X1{i}.vertices,X2{j}.vertices,lambda,eps);
        end
    end
end