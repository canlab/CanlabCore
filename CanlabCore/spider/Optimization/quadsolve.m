function [x,y,z,s,w,flag] = quadsolve(Q,c,A,b,C),  
%SQUADSOLVE
% 
% USAGE:   [x,y,z,s,w,flag] = squadsolve(Q,c,A,b,C)
%
% PARAMETERS:  Q -> (n,n) symetric matrix (definite positive)
%              c -> (n,1) vector
%              A -> (m,n) matrix 
%              b -> (m,1) vector
%              C -> scalar
%
%            x -> primal variables
%            y -> lagrangian coeff of equality constraints
%            z -> dual variables of x
%            s -> primal auxiliary variable (only if C < Inf)
%            w -> dual variable of s
%            flag -> set to 0 => no problem, set to 1 => problem
%
% DESCRIPTION: Primal-dual method for quadratic programming
%                
%            minimize c'*x + 0.5*x'*Q*x
%
%            subject to  A*x=b
%                        0<= x <= C
%            The method used here is a primal dual method with a predictor-corrector
%            approach and a logarithmic barrier. I used the heuristic from two 
%            existing methods LOQO and HOPDM. The method is an iterative method. The 
%            maximal number of iteration is stored in the variable 'max_iter'.
%
% ERRORS AND BUGS: there is no test about the conditionning of the matrix Q. If the iteration
%                  50 has been reached, then the optimization may not be finished and the output
%                  may be wrong.
%
% NOTES: 50 iterations have always been sufficient to solve all problems.
%        This code should be read with the tech. report:
%                "Regularized Symmetric Indefinite Systems in Interior Point
%                 Methods for Linear and Quadratic Optimization", 
%                 A. Altman and J. Gondzio, Logilab Tech. Report 1998.6
%
% Andre Elisseeff, Dec. 1999
% aelissee@eric.univ-lyon2.fr
  
% init    
  verbose = 0;
  n = size(Q,1);
  m = size(A,1);
  H = zeros(n+m,n+m);
  flag=1;
  
% Values of the original HOPDM of Gondzio and Altmann
  dinf = 10^(-14);
  smallz = 10^(-14);
  smallt = 2.3*10^(-16);  
  opttol = 10^(-6);
  itref = 1;
  mu = 1;
  maxiter = 50; 
  
% init values of the primal and dual variables
  x=ones(n,1);
  z=ones(n,1);
  y=ones(m,1);
  if C < Inf,
    s=ones(n,1);
    w=ones(n,1);
  else
    s=[];w=[];
  end;
% Description of variables:
%
%    x,s     -> primal variables
%    z,w     -> dual variables
%    n       -> number of variables in the initial pb (size of x)
%    m       -> number of constraints in A
%
%    dinf     -> smallest value for all variables   
%    smallz   -> smallest value of z
%    smallt   -> smallest value for t in the computation of the matrix theta
%    opttol   -> acceptable tolerance for optimality conditions
%    itref    -> iteration counter
%    maxiter  -> maximum number of iteration
%    DEBUG    -> 0 => no debug, 1 => debug
 
% Analyze the constraints...
  if verbose,
   disp(sprintf('Analyzing the equality constraints...\n'));
  end;
  [QQ,RR]=qr(A',0);
  [mm,nb] = size(QQ); %% number of eq constraints
  ind = 1:1:nb;
  for i=1:nb,
    if abs(RR(i,i)) < 10^(-10),
      disp(sprintf('Constraints %d removed because of dependence\n',i));
      ind(i)=0;
    end;
  end; 
  indice = find(ind >0);
  
  if (isempty(indice))
    disp(sprintf('No equality constraints... \n'));
    A=[];
    m=0;
  else
    A = A(indice,:); %% new independent eq constraints 
    b = b(indice);
    y = y(indice);
    m=length(indice);
  end;
  clear QQ;clear RR;
%%%% AE : scale the problem, has been tested and seems to work
%%%% better : more stable... 
  if C < Inf,
    u = ones(n,1);
    Q = Q*C;
    b = b/C;
    c = c;
    opttol = opttol/C;
  end;
% init values before looping
  cont = 1;
  objQ = 0.5*x'*Q*x;
% init values of primal and dual objective functions
  pobjo = abs(c'*x+objQ) + 1;
  if C < inf,
    dobjo = abs(b'*y - u'*w - objQ);
  else
    dobjo = abs(b'*y - objQ);
  end;
%%%%%%%%%%%%%%%
%% MAIN LOOP
%%%%%%%%%%%%%%%
  while (cont)&(itref<=maxiter)
  % Compute the primal objective function
    objQ = 0.5*x'*Q*x;
    pobj = c'*x+objQ;
    if (C<Inf)
        if ~isempty(A),
            dobj = b'*y - u'*w - objQ;
        else
            dobj = - u'*w - objQ;
        end;
    else
        if ~isempty(A),
            dobj = b'*y - objQ;
        else
            dobj = - objQ;
        end;      
    end;
    dlgap = pobj - dobj;
    dp = abs(pobj)/(abs(pobjo)+1);
    dd = abs(dobj)/(abs(dobjo)+1);
    dobjo = dobj;
    pobjo = pobj;
    if verbose,
        disp(sprintf('%d - pobj : %f - dobj : %f\n',itref,pobj,dobj));
    end;
  % Check if the solution are bounded
    if (dp > 10^6) 
      disp(sprintf('Solution not bounded in the primal. Exit.\n'));
      return;
    end;
    if (dd > 10^6) 
      disp(sprintf('Solution not bounded in the dual. Exit.\n'));
      return;
    end;
    
  % test if optimality
    oldgap = dlgap;
    dp = abs(dobj) + 1;
    if ((abs(dlgap)/dp) <= opttol) & (itref>1)
      if verbose,
       disp(sprintf('Optimal solution found. Exit.\n'));
      end;
      cont = 0;
      break;
    end;
    
    dp = dp + abs(pobj);
    T = abs(dlgap)/dp;
    
  % put the variables away from zero (from HOPDM)
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
    if ~isempty(A),
     xi_b = -A*x + b;    
     xi_c = c - A'*y - z + Q*x;
    else
      xi_b = [];    
      xi_c = c - z + Q*x;
     end;
     
    xi_z =   - x.*z;
    if C < Inf,
      xi_c=xi_c + w;
      xi_u = u - x - s;
      xi_w =  - s.*w;
    end;
  % Compute theta = (z/x + w/s)
    
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
  % factorize H = [-Q-theta^(-1)   A^T]
  %               [ A               0 ]
    
  
    H = zeros(n+m,n+m);
    H(1:n,1:n) = -Q-diag(theta);
    if ~isempty(A)
     H(n+1:n+m,1:n) = A;
     H(1:n,n+1:n+m) = A';
    end
 
  % Compute the predictor step
    if C < Inf,
      f = xi_c-xi_z./x+(xi_w - xi_u.*w)./s ;
      h = xi_b;
    else
      f = xi_c - xi_z./x;
      h = xi_b;
    end;
    delta=H\[f;h];
    dx = delta(1:n);
    dy = delta(n+1:n+m);
    dz = (xi_z-z.*dx)./x;
    if C<Inf,
      ds = xi_u - dx;
      dw = (xi_w-w.*ds)./s;
    end;      
      
 % determine the maximum step size alpha_p (primal) and
 % alpha_d (dual) to stay in feasible region
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
        break;
      end;
      
    % Set the barrier parameter : LOQO's heuristic
      ap = min(apk,adk);
      mu = (ax/(2*n))*(0.95*(1/ap) -1)^2/(0.95*(1/ap)+10)^2;
      
    % Compute the new direction (algo. of Mehrotra) of order 1 (corrector step)
      xi_z =  - x.*z + mu*ones(size(x)) - dx.*dz;
      f = xi_c - xi_z./x;
      if C < Inf
        xi_w = - s.*w + mu*ones(size(s)) - ds.*dw;
        f=f+ xi_w./s - (w.*xi_u)./s;
      end;
      h = xi_b;
      
      delta=H\[f;h];
      dx = delta(1:n);
      dy = delta(n+1:n+m);
      dz = (xi_z-z.*dx)./x;  
      if C<Inf,
        ds = xi_u - dx;  
        dw = (xi_w-w.*ds)./s;
      end;      
      
      
    % determine the maximum step size alpha_p (primal) and
    % alpha_d (dual) to stay in feasible region
    % (x,s,z,w must be positiv)
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
      
    % Compute step factors     
      fp =0.9*min(alpha_p,alpha_d);
      fd =0.9*min(alpha_d,alpha_p);    
      x = x +fp*dx;      
      y = y +fd*dy;
      z= z + fd*dz;
      if C < Inf,
        w = w + fd*dw;
        s = s +fp*ds;      
      end;
      itref=itref+1;
  end;
%%%%%%%%%%%%%%%%%%%%%
%%% End of main loop
%%%%%%%%%%%%%%%%%%%%%
  
    if (cont==0)
        if verbose,
          disp(sprintf('Optimal Solution found after %d iteration.\n',itref));
          disp(sprintf('Value of the objective : %f.\n',pobj));
        end;
      flag=0;
    end;
  % rescale x and y
    if C<Inf,
      x= C*x;
      %y = C*y;
    end;
   
