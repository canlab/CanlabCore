function [beta, bo] = svm_multi_predK(X,Y,C,K)
% SVM_MULTI_PREDK
%
% Support Vector Multi Classification
%
% USAGE: [beta, bo] = svm_multi_pred(X,Y,C,K)
%
% PARAMETERS:  X      - (m,d) matrix of m Training inputs in R^d
%
%              Y      - m vector of m Training targets in {1,..,Q}
%
%              C      - Trade-off between regularization and empirical error 
%
%              K      - (m,m) Gram matrix
%
%              beta   - (m,Q) Coefficients matrix of the expansion of
%                       the vectors w_i over the x_p's
%                       w_i = \sum_p beta(p,i)x_p
%
%              bo     - (Q,1) bias terms
%
% SUBROUTINES: svm_multi_init.m -> initialization of the optimization problem
%             squadsolve.m          -> optimization 
%             compute_kernel.m -> computation of the gram matrix
%
% DESCRIPTION:
%
%              This procedure is the implementation of a multiclass SVM corresponding to 
%              the following problem
%
%              \min_{w_i,w_j} \sum_{i\neq j} \|w_i - w_j\|^2 + C\sum_{p,j} \xi_{pj}
%
%              subject to :
%               for all p=1,..,m   ( w_{c(p)} - w_j ).x_p + b_{c(p)} - b_j \geq 1 - \xi_{pj}
%                                  \xi_{pj} \geq 0
%
%      where:  m is the number of training points
%              Q is the number of classes
%              i and j are indices between 1 and Q
%              p is an index between 1 and m
%              w.x is the dot product between w and x
%
%              To solve this problem, a direct approach is used. All the gram matrix K (K_{i,j} =x_i.x_j)
%              is computed before beginning the optimization. The memory cost is thus proportionnal to Q^2*m^2.
%              The solution is not approximated but is provided by a quadratic programming method implemented in
%              squadsolve.m. Instead of minimizing the primal objective function, we maximize the dual. Its 
%              formulation is not given here.
%
%              
% ERRORS AND BUGS:
%              It is recommended not to use this code for large problems: unless your computer has a huge memory,
%              the procedure will stop with the message: 'memory exhausted'. Use the code msvm_fw.m or chunk_msvm.m
%              instead. Sometimes, the bias may not be computable because the term 'bias_eps' (see below) is too large.
%              In that case, decrease 'bias_eps' and it should work. In all cases, the program will output the bias and
%              try to give the best solution. 
%              
% EXAMPLE OF USE:
%
%      > [beta,bo] = svm_multi_pred(inputs,targets,10,K);
%
%       -> solve the multiclass optimization problem for a learning set S = \cup_i {(inputs(i,:),targets(i))}
%          with C=10, and with the gram matrix K.
%
% Description of variables:
%
%      kernel = 0,1,2 (integer)    -> type of kernel
%      m              (integer)    -> number of input points
%      Q              (integer)    -> number of classes
%      bias_eps       (real)       -> precision for the computation of the bias
%      reg_eps        (real)       -> to avoid ill conditionning, a diagonal matrix
%                                     is added to the hessian of the minimization problem,
%                                     reg_eps scales this diagonal matrix. (usual < 10^(-6))
%      A              ((Q,Q*m) matrix) -> equality constraint matrix 
%      H            ((Q*m,Q*m) matrix) -> hessian of the objective function
%      c               (Q*m vector)    -> linear part of the objectif function
%      alpha           (Q*m vector)    -> variables of the optimization problem
%      constraints     (Q*m vector)    -> store 1-(w_{c(p)} - w_j).x_p
%      beta            ((m,Q) matrix)  -> vector w_i is computed as w_i = sum_p beta_{pi} x_p
%      b               ((Q,Q) matrix)  -> store the differences between bias: b(i,j) = b_i - b_j
% Avoid warning messages
warning off;
disp('Optimizing...');
% Test if nargin is correct
if (nargin < 2) | (nargin > 6),
help svm_multi_pred;
else
  
% Init
 empty_list_elements_ok = 1;
 m = size(X,1);
 Q = max(Y);
 if (nargin<3),
     C=Inf;
 end;  
 if (C==Inf)
    tol = 1e-5;
  else
    tol = C*1e-6;
 end;
 bias_eps=10^(-8)/C;
 reg_eps = 10^(-6);
 sv_eps=10^(-6);    
% Message
 % disp(sprintf('\nMulti Support Vector Classification\n'));
 % disp(sprintf('-----------------------------------\n'));
   
% Set up the parameters for the Optimisation problem
%disp(sprintf('Initialization...\n'));
[H,A,c]=svm_multi_init(Y,K);
% Add small amount of zero order regularization to
% avoid problems when Hessian is badly conditioned.
H = 1/(2*Q^2)*H+reg_eps*eye(size(H));% the problem is to minimize (1/4Q^2)x'Hx - c.x
c = -c;
% Some variables are irrelevant. Consider only relevant variables
mul = Q*(0:m-1);
indx = mul' + Y;
s1 = (1:Q*m);
s2 = (indx);
indx = setdiff(s1,s2)';
clear s1;clear s2;clear mul;
% Solve the optimization problem
% the solution is stored in alpha
%disp(sprintf('Begin the optimization ...\n'));
bo = zeros(size(A,1),1);
[x1,y] = quadsolve(H(indx,indx),c(indx),A(:,indx),bo,C);
x1=x1(1:length(indx));
obj = 0.5*x1'*(H(indx,indx)-reg_eps*eye(length(indx)))*x1 + c(indx)'*x1;
clear H;clear c;
alpha = zeros(Q*m,1);
alpha(indx)=x1;
clear x1;
% Output the objective value
%disp(sprintf('Final objective function: %f \n',-obj));
% if the alpha's are all greater than C then, C can be considered as infinite
% (for the computation of the bias)
if (max(abs(alpha)) < C*0.95)
    C==Inf;
end;
% Compute the number of support vector
% Can be removed if not desired
M = max(alpha);
nsv=0;
for i=1:m,
  if ~(isempty(find(alpha(Q*(i-1)+1:Q*i)>100*eps)))
    nsv = nsv+1;
  end;
end;
%disp(sprintf('Support vectors : %d (%3.1f)\n',nsv,100*nsv/(m)));
% Compute the coefficients beta_ip of the vector w_i = sum_p beta_ip x_p
% clean alpha
beta=zeros(m,Q);
notvoid = find(alpha>=100*eps);
alpha_tmp=zeros(Q*m,1);
alpha_tmp(notvoid) = alpha(notvoid);
alpha=alpha_tmp;
clear notvoid;clear alpha_tmp;
% compute beta
for i =1:Q,
 for p=1:m,
  tmp=0;
  if (i==Y(p))
   tmp = tmp + sum(alpha((p-1)*Q + 1: (p-1)*Q + Q))/Q;
  end;
  beta(p,i) = tmp - alpha((p-1)*Q + i)/Q;
 end;
end;
% Computation of the bias
%disp(sprintf('\n Computation of the bias........\n'));
% Computation of the outputs 1-(w_{c(i)} - w_j).x_i
errorcache = zeros(Q*m,1);
for p=1:m,
  for i=1:Q,
    temp = ((beta(:,Y(p))-beta(:,i))'*K(:,p));
    errorcache((p-1)*Q+i)=temp;
  end;
end;     
 
  % Computation of the output 1-(w_Y(p)-w_i).x_p
  constraints=1-errorcache;
   if C==Inf,
     b=-Inf*ones(Q,Q);
     for i=1:m,
         for k=1:Q,
            b(Y(i),k) = max([b(Y(i),k),(constraints(Q*(i-1)+k) - 0.05)]);             
         end;
     end;
     for k=1:Q,
         b(k,k)=0;
         b(k,find(b(k,:)==-Inf))=0;
     end;
   end; %if C>max(alpha),
%%% For finite C
if C<Inf,
% computation of alpha_bias
  % here: computing of the bias that minimize the \sum_{ip} \xi_{pi}
  % subject to linear constraints:
  % (w_{c(p)} - w_i).x_p + b_{c(p)} - b_i >= 1 - \xi_{pi}
  % \xi_{pi}
  % This is a linear program whose dual is easy to compute. It can be
  % solved with the slinearsolve method.
  alpha_bias=zeros(Q*m,1);
  c_un = ones(Q*m,1);
  [alpha_bias(indx)] = slinearsolve(errorcache(indx)-c_un(indx),A(1:(Q-1),indx),zeros(Q-1,1),1);
 
  
  
  % clean alpha_bias
  tmp = 0:m-1;
  tmp = tmp*Q;
  tmp = tmp + Y';
  alpha_bias(tmp)=0;
  
  
  % if problem, consider the alpha instead
  if (flag~=1)
    alpha_bias=alpha;
  end;
  
  
  
  b=zeros(Q,Q); % matrix of the difference b_i-b_j  
  % Computation of the output 1-(w_Y(p)-w_i).x_p
  constraints=1-errorcache;
  
  % consider the alpha_ip s.t. 0 < alpha_bias_ip < 1
  % here, i consider alpha_new, because i compute the bias
  % that minimize the l_1 error of the current solution
  % and not of the optimal.
  if C == Inf,
    Cm =1;
  else
    Cm=C;
  end;
  if (flag~=1)
   indice_fix = find(alpha_bias < C*(1-bias_eps) & alpha_bias > Cm*bias_eps);
  else
   indice_fix = find(alpha_bias < (1-bias_eps) & alpha_bias > bias_eps);
  end;
  
  % if there are no such alpha -> problem
  if isempty(indice_fix)
    disp('WARNING : Unable to compute the bias, use another method\n');      
  end;
  % the bias are computed by averaging 
  % the matrix 'compteur' counts when  b_i-b_j is computed
  compteur=zeros(Q,Q);
  for ind=indice_fix',
    % find the index Q and p corresponding to ind
    indQ = rem(ind,Q);
    if (indQ==0), indQ=Q;end;
    indp = (ind-indQ)/Q;
    compteur(Y(indp+1),indQ)=compteur(Y(indp+1),indQ)+1;
    b(Y(indp+1),indQ) = b(Y(indp+1),indQ) + constraints(ind);
  end;
  % average by dividing by compteur
  for i=1:Q,
    b(i,i) = 0;
    for j=1:Q,
      if (compteur(i,j))
       b(i,j)=b(i,j)/compteur(i,j);
      end;
    end;
  end;
  end; %if C...,
  clear indice_fix;clear constraints;
  
  
  % transform b in order to be symetric
  for i=1:Q,
    for j=i+1:Q,
      if (b(i,j)==0)
        b(i,j)=-b(j,i);
      end;
      if (b(j,i)==0)
        b(j,i)=-b(i,j);
      end;                       
    end;
  end;
  
  % average between both symetric elements
  if C==Inf,
   for i=1:Q,
    for j=i+1:Q,
        b(i,j)=max(b(i,j),-b(j,i));
        b(j,i)=-b(i,j);
    end;
  end; 
  else,
  for i=1:Q,
    for j=i+1:Q,
        b(i,j)=(b(i,j)-b(j,i))/2;
        b(j,i)=-b(i,j);
    end;
  end;
  end;
  
  
  % compute the values of b_i s.t. sum_i b_i = 0
  ok=0;% ok > 0, if computation is possible otherwise ok=0
  b_guess = zeros(Q,1); % if no computation is possible still give a guess
  bo = zeros(Q,1); % the bias vector
  for i=1:Q,
    z = find(b(i,:)~=0);
    if (length(z)==Q-1)
      ok=ok+1;
      temp = sum(b(i,:))/Q;
      for j=1:Q,
        bo(j) =bo(j)+temp-b(i,j);
      end;
    else
      %disp(sprintf('Warning : bias %d -> problem... %d\n',i,length(z)));
      if length(z),
        b_guess(i) = sum(b(i,:))/length(z);
      end;
      for j=1:length(z),
        b_guess(z(j)) =b_guess(i) -b(i,z(j));
      end;
    end;                
  end;
  
  if ~ok
    % then one should compute the bias with another method
    disp('WARNING : Problem with the bias, use another method');
    bo=b_guess;
  else
    bo =bo/ok;
  end;
  
end; %if nargin...
