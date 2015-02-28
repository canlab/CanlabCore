function [H,A,c] = svm_multi_init(targets,K)  
% Initialization of the parameters for the optimisation
% problem of multiclass support vector machines.
%
% USAGE:      [H,A,c] = svm_multi_init(targets,K)
% 
% PARAMETERS: targets -> m vector of m Training targets in {1,..,Q}
%             K       -> Gram matrix K_{ij} = x_i.x_j where . is the dot product
%
% SUBROUTINES: none
%
% DESCRIPTION:
%                This procedure computes the matrices H and A that correspond to the dual
%                problem of the multiclass svm optimization problem. H is the hessian of the
%                dual and A is the equality constraint matrix of the dual. It is used only
%                by svm_multi_pred.m.
%                
% Andre Elisseeff, Sep. 2000
% aelissee@eric.univ-lyon2.fr             
%
% NOTES: Les calculs sont exposes dans ma these p.56-57 avec les memes notations (a peu pres).
  
% Description of variables:
%
%        m  (integer) -> number of input and target points
%        Q  (integer) -> number of classes
%        H  (Q*m,Q*m) matrix -> hessian of the dual
%        A  (Q,Q*m) matrix   -> equality constraint matrix
%        c  Q*m vector       -> linear part of the dual
  
  
% init
  m = size(targets,1);
  Q = max(targets);
  H = zeros(m*Q,m*Q);
  A = zeros(Q,m*Q);
  c = ones(m*Q,1);
% Computation of H and c
  for i=1:Q,
    for j=i+1:Q,
      B = zeros(m,Q*m);
      for p=1:m,
        indexe=targets(p);
        % commputes the value for c
        c((p-1)*Q + indexe)=0;
        if (i==indexe)
         B(p,(p-1)*Q+1:(p-1)*Q+Q)=B(p,(p-1)*Q+1:(p-1)*Q+Q)+ones(1,Q);
        end;
        B(p,(p-1)*Q+j)=B(p,(p-1)*Q+j)+1;
        B(p,(p-1)*Q+i)=B(p,(p-1)*Q+i)-1;
        if (j==indexe)
         B(p,(p-1)*Q+1:(p-1)*Q+Q) = B(p,(p-1)*Q+1:(p-1)*Q+Q)-ones(1,Q);
        end;
      end; %% p
      H = H + 2*B'*K*B;
    end; %% j
  end; %% i
  
% Computation of the constraints A
  for i=1:Q,
    for p=1:m,
      if (i==targets(p)) 
        for j=1:Q,
        A(i,Q*(p-1)+j) = 1;
        end; %% j
      end;
       A(i,Q*(p-1)+i) = A(i,Q*(p-1)+i)-1;
    end; %% p
  end; %% i
