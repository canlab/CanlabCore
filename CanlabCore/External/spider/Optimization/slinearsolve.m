function [x,y,flag,z,s,w] = slinearsolve(c,A,b,C), 
% Simple Linear Optimization for Franke and Wolfe's multiclass SVM method
% 
% USAGE  :   [x,y,flag,z,s,w] = slinearsolve(c,A,b,C)
%
% PARAMETERS:
%          c = (n,1) vector
%          A = (m,n) matrix 
%          b = (m,1) vector
%          C = scalar
%
%          x = primal variables
%          y = lagrangian coeff of equality constraints
%          z = dual variables of x
%          s = primal auxiliary variable (only if C < Inf)
%          w = dual variable of s
%          flag = integer value, set to 1 if no problem (-1 if unbounded and 0 otherwise)
%            
% DESCRIPTION:
%                minimize c'*x + 0.5*trdoff*(x-alpha)'*(x-alpha)
%
%                subject to A*x=b
%                           0<= x(i) <= C(i)
%                
%            The optimization method is a predictor-corrector method with a logarithmic barrier.
%            For separable problems such as this, this method is quite fast. It tends to minimize
%            the dot product c'*x. 
%            The term trdoff controls the trade-off between both parts. Here, trdoff=0 and alpha=0. 
%
%
% ERRORS AND BUGS:
%            As for the method squadsolve.m, this procedure is iterative and the max. number of iteration
%            is stored in max_iter. If it is reached (which is unlikely to happen), the output of the procedure
%            may give a suboptimal solution.
%
% NOTES:     
%            The heuristics used to set the barrier parameter and the init value are taken from LOQO's method. 
%            The other tricks (such as putting all parameters away from zero) are taken from HOPDM 
%            of Altman and Gondzio.
%
%            This code should be read with the tech. report:
%                "Regularized Symmetric Indefinite Systems in Interior Point
%                 Methods for Linear and Quadratic Optimization", 
%                 A. Altman and J. Gondzio, Logilab Tech. Report 1998.6
%
%            The quadratic term is not used in this implementation. However, if one wants to add it, just 
%            change the value of trdoff and of alpha.
%
% Andre Elisseeff, May 2001
% aelissee@barnhilltechnologies.com
  
  
% init values
  verbose = 0;
  flag=0; 
  n = size(A,2);
  m = size(A,1);
  mu = 1;
  maxiter = 50; 
  itref = 1;
  cont = 1;  
  clean_eps = 10^(-8);
  trdoff = 0;
  alpha = zeros(n,1);
  
% init from LOQO
  x=100*ones(n,1);
  z=100*ones(n,1);
  y=100*ones(m,1);
  if C < Inf,
    s=100*ones(n,1);
    w=100*ones(n,1);
  else
    s=[];w=[];
  end;
 
% Values of the original HOPDM of Gondzio and Altman  
  dinf = 10^(-14);
  smallz = 10^(-14);
  smallt = 2.3*10^(-16);  
  opttol = 10^(-6);
  
% Description of variables:
%
% theta (real) (n,1) vector  -> diagonal matrix (thus stored as vector) to
%                               be inverted to compute the predictor and the
%                               corrector step. See tech. report of Altman and
%                               Gondzio for further details.
% mu         -> barrier parameter
% u          -> vector of upper bounds on x (if C < Inf)
% pobj       -> primal objective function
% pobjo      -> old primal objective function
% dobj       -> dual objective function
% dobjo      -> old dual objective function
% dlgap      -> dual gap
% oldgap     -> old dual gap
% objQ       -> 0.5*x'*x
% n          ->   number of variables in the initial pb (size of x)
% m          ->   number of constraints in A
% Q          ->   number of classes
% x,s  (real) (n,1) vector  ->   primal variables
% z,w  (real) (n,1) vector  ->   dual variables
% dinf       ->   smallest acceptable values for x,y,s etc...
% smallz     ->   smallest acceptable values for z in the inversion
%                 of theta
% smallt     ->   smallest acceptable values for t in the computation
%                 of the matrix theta
% clean_eps  ->   precision to clean x (at the end of the procedure)
% opttol     ->   acceptable tolerance for optimality conditions
% itref      ->   iteration counter
% maxiter    ->   maximum number of iterations
% cont       ->   boolean value to continue in the loop
% flag       ->   set to 1 => computation ok, set to 0 => pb in the computation
%                 set to -1 => solution unbounded.
%                 example: maxiter has been reached (flag=0).
% trdoff     ->   trade-off between the linear and the quadratic part of the objective
%                 function
% ind_tmp_Y  ->   indices corresponding to alpha_{c(p)p}
  
%## AE : scale the problem, has been tested and seems to work
%## better : more stable... 
  if C < Inf,
    u = ones(n,1);
    b = b/C;
  end;
    
    
 % init values of primal and dual objective functions
   objQ = 0.5*trdoff*x'*x;
  
 % add the alpha term in the linear part
   c = c - trdoff*alpha;
  
  
  pobjo = abs(c'*x+objQ) + 1;
  if C < inf,
    dobjo = abs(-objQ+b'*y - u'*w) +1;
  else
    dobjo = abs(b'*y-objQ) +1;
  end;
%###############
%## MAIN LOOP
%###############
  while (cont)&(itref<=maxiter),
      
  % Compute the primal objective function
    objQ = 0.5*trdoff*x'*x;
    pobj = c'*x + objQ;
    if (C<Inf)
      dobj = b'*y - u'*w -objQ;
    else
      dobj = b'*y -objQ;
    end;
    dlgap = pobj - dobj;
    dp = abs(pobj)/(abs(pobjo)+1);
    dd = abs(dobj)/(abs(dobjo)+1);
    dobjo = dobj;
    pobjo = pobj;
    if verbose,
        if C<Inf,
            disp(sprintf('%d - pobj : %f - dobj : %f\n',itref,pobj,dobj));
        else   
            disp(sprintf('%d - pobj : %f - dobj : %f\n',itref,pobj,dobj));
        end;
    end;    
  % Check if the solution are bounded    
    if (dp > 10^6) 
      disp(sprintf('Solution not bounded in the primal. Exit.\n'));
      flag=-1;
      return;
    end;
    if (dd > 10^6) 
      disp(sprintf('Solution not bounded in the dual. Exit.\n'));
      flag=-1;
      return;
    end;
    
  % test if optimality    
    oldgap = dlgap;
    dp = abs(dobj) + 1;
    if (abs(dlgap)/dp) <= opttol & itref > 1
      if verbose,
       disp(sprintf('Optimal solution found. Exit.\n'));
      end;
      cont = 0;
      break;
    end;
    
    dp = dp + abs(pobj);
    T = abs(dlgap)/dp;
    
  % put the variables away from zero (from QHOPDM)
    if (itref <= 3)
      ax = 2*10^(-3);
      az = 10^(-3);
    elseif (T >= 0.8)
      ax = 2*10^(-4);
      az = 10^(-4);
    elseif (T >= 0.1)
      ax = 2*10^(-5);
      az = 10^(-5);
    elseif (T >=0.01)
      ax = 2*10^(-6);
      az = 10^(-6);
    elseif (T>=0.001)
      ax = 2*10^(-7);
      az = 10^(-7);
    elseif (T>=0.0001)
      ax = 2*10^(-8);
      az = 10^(-7);
    elseif (T>=0.00001)
      ax = 2*10^(-9);
      az = 10^(-9);
    else
      ax = T*10^(-5);
      az = ax;
    end;
    
  % consider only variables that can be changed
    x = x + ax;
    z = z + az;
    if C < Inf
      s = s + ax;
      w = w + az;
    end;
  % Compute the values of xi_b, xi_c and xi_u    
  % These variables are not described here. See tech.
  % report of Altman and Gondzio for further details.
  
    xi_b = -A*x + b;    
    xi_c = c - A'*y - z + trdoff*x;
    xi_z =  - x.*z;
    if C < Inf,
      xi_c=xi_c + w;
      xi_u = u - x - s;
      xi_w = - s.*w;
    end;
    
  % Compute theta = xs/(zs + wx)
    
  % for bounded variables
    if C<Inf,
      dp = x;
      if (max(abs(dp))<= smallz)
        disp(sprintf('Conditioning problem to invert theta. Abort.\n'));
        return;
      end;
      dpp= s;
      if (max(abs(dpp))<= smallz)
        disp(sprintf('Conditioning problem to invert theta. Abort.\n'));
        return;
      end;
      theta=z./dp + w./dpp;
    end;
    
  % for unbounded variables
    if C==Inf,
      dp = x;
      if (max(abs(dp))<= smallz) 
        disp(sprintf('Conditioning problem to invert theta. Abort.\n'));
        return;
      end;
      theta = z./dp;
    end;
    tmp = ones(size(theta));
    theta = tmp./(trdoff*(tmp)+theta);
     
  % neglect small elements of theta array    
    neglect = find(theta < smallt);
    if ~isempty(neglect)
      theta(neglect)=zeros(size(neglect));
    end;
    
  % and control large elements of theta    
    neglect = find(theta >= 10^8);
    if ~isempty(neglect)
      theta(neglect)=(10^4)*sqrt(theta(neglect));
    end;
  % solve the system for mu = 0        
    if C < Inf,
      f = xi_c-xi_z./x+(xi_w - xi_u.*w)./s ;
      h = xi_b;
    else
      f = xi_c - xi_z./x;
      h = xi_b;
    end;
    AA = zeros(size(A));
    for i=1:size(A,1),
      AA(i,:) = A(i,:).*theta';
    end;
    to_inv = AA*f + h;
    AAA = AA*A';
    dy = AAA\to_inv;
    dx = AA'*dy - f.*theta;
    dz = (xi_z-z.*dx)./x;
    if C<Inf,
      ds = xi_u - dx;
      dw = (xi_w-w.*ds)./s;
    end;          
  % end of solving the system
    
  % determine the maximum step size apk (primal) and
  % adk (dual) to stay in feasible region
  % (x,s,z,w must be positive and greater than dinf)
    indz = find(dz<0);    
    indx = find(dx<0);
    inds=[];mins=1;
    indw=[];minw=1;
    if C < Inf,
      inds = find(ds<0);
      indw = find(dw<0);
      if ~isempty(inds)
        mins = min(-(s(inds)-dinf)./ds(inds));
      else
    mins = 1;
      end;
      if ~isempty(indw)
        minw = min(-(w(indw)-dinf)./dw(indw));
      else
    minw = 1;
      end;      
    end;
    if ~isempty(indx)
      minx = min(-(x(indx)-dinf)./dx(indx));
    else
      minx = 1;
    end;
    apk = min([minx,mins,1]);
    if ~isempty(indz),
      minz = min(-(z(indz)-dinf)./dz(indz));
    else
      minz = 1;
    end;
    adk = min([minw,minz,1]);
    
    ax = sum(x.*z);
    as = sum((x+apk*dx).*(z+adk*dz));
    az = sum(dx.^2+dz.^2);
    if C < Inf,
      ax = ax + sum(s.*w);
      as = as + sum((s+apk*ds).*(w+adk*dw));
      az = az + sum(ds.^2+dw.^2);
    end;
    
  % check if complementary gap is less than opttol    
    if (as <= opttol)
      if verbose,
        disp(sprintf('Complementary gap is less than %f\n',opttol));
      end;
      cont = 0;
      x = x + apk*dx;
      y = y + adk*dy;
      z = z + adk*dz;
      if C < Inf,
       s = s + apk*ds;
       w = w + adk*dw;
      end;
    % get out of the while loop
      break;
    end;
    
    
  % Set the barrier parameter : LOQO's heuristic
  %  ap = min(apk,adk);
  % mu = (ax/(2*n))*(0.95*(1/ap) -1)^2/(0.95*(1/ap)+10)^2;
  % HOPDM heuristic:
    mu = (as/ax)^2*as/n;
    
  % Compute the new direction (algo. of Mehrotra) of order 1
    xi_z =  - x.*z + mu*ones(size(x)) - dx.*dz;
    f = xi_c - xi_z./x;
    if C < Inf
      xi_w = - s.*w + mu*ones(size(s)) - ds.*dw;
      f=f+ xi_w./s - (w.*xi_u)./s;
    end;
    h = xi_b;
    
  % solve the system for mu different from zero      
    
    AA = zeros(size(A));
    for i=1:size(A,1),
      AA(i,:) = A(i,:).*theta';
    end;
    to_inv = AA*f + h;
    AAA = AA*A';
    dy = AAA\to_inv;
    dx = AA'*dy - f.*theta;
    dz = (xi_z-z.*dx)./x;
    if C<Inf,
      ds = xi_u - dx;
      dw = (xi_w-w.*ds)./s;
    end;      
    
  % end for solving
    
  % determine the maximum step size apk (primal) and
  % adk (dual) to stay in feasible region
  % (x,s,z,w must be positive)
    
    indz = find(dz<0);    
    indx = find(dx<0);
    if C < Inf,
      inds = find(ds<0);
      indw = find(dw<0);
      if ~isempty(inds)
        mins = min(-s(inds)./ds(inds));
      else
        mins=1;
      end;
      if ~isempty(indw)
        minw = min(-w(indw)./dw(indw));
      else
        minw=1;
      end;
    end;
    if ~isempty(indx)
      minx = min(-x(indx)./dx(indx));
    else
      minx=1;
    end;
    alpha_p = min([minx,mins,1]);
    if ~isempty(indz)
      minz = min(-z(indz)./dz(indz));
    else
      minz=1;
    end;
    alpha_d = min([minw,minz,1]);    
    spd = min(alpha_p,alpha_d);
  % Compute step factors     
  %      fp =0.95*spd;
  %      fd =0.95*spd; (from QHOPDM)
    fp = 0.95*spd;
    fd = 0.95*spd;      
    x = x +fp*dx;      
    y = y +fd*dy;
    z= z + fd*dz;
    if C < Inf,
      w = w + fd*dw;
      s = s +fp*ds;      
    end;
    itref=itref+1;
end;
%#####################
%### End Main Loop
%#####################
  if (cont==0)
    flag=1;
  end;
  
% clean x
  tmp = find(x < clean_eps/n | x < 100*eps);
  if ~isempty(tmp)
    x(tmp) = zeros(size(tmp));
  end;
  
% rescale x and y
  if C<Inf,
    x= C*x;
%    y = C*y;
  end;
  %end of function
