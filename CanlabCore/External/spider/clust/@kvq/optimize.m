function [x,fval] = optimize(a,f,A,b,l,u)

  if nargin < 5
    % set lower and upper bound to -inf resp +inf
    l = -inf*ones(size(A,2),1);  u = -l;
  end
			
  % optimizing with builtin matlab linprog 
  if strcmp(a.optimizer,'linprog')
     [x,fval,exitflag]=linprog(f,A,b,[],[],l,u);
  end
  
  % optimizing with Cplex
  if strcmp(a.optimizer,'cplex')
    [x,y_upd_,how_upd_,fval]= lp_solve(cplex_license, f, sparse(A), b, l, u, 0, 0);
  end
