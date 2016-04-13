function B = Bspline(t,k,u,v,ForceSup)
% :Usage:
% ::
%
%     Bspline(t,k,u[,v])
%
% Create B-Spline basis of order k, with knots u, evaluated at t.
% If control verticies v are specified then then B is the spline
% function instead of the basis.
%
% u must be at least length(k)+1 
%
% ..
%    $Id: Bspline.m,v 1.3 1998/11/09 00:55:06 nicholst Exp $
% ..

global iB

if (nargin<3)
  help Bspline
  return
end

if (k+1>length(u))
  error('u must be at least length k+1')
end

% Silent flag for forcing the support to be defined between u(k)
% and u(end-k+1)
if (nargin<5)
  ForceSup = 1;
end

if ((nargin>3) & (length(v)+k ~= length(u)) & (ForceSup))
    error(sprintf('%d knots requires %d control verticies', ...
	length(u), length(u)-k))
end

% columnize
t=t(:);
u=u(:);

nBasis = length(u)-k;
B = zeros(length(t),nBasis);
iB = zeros(1,nBasis);

for i=1:nBasis
   B(:,i) = recu(t,i,k,u);
   iB(i) = (u(i+k)-u(i))/k;
end

% zero outside of valid range, if there's enough
if (length(u)>=2*k)
  if (ForceSup)
    B(t<u(k) | u(end-k+1)<t,:) = 0;
  end
else
  warning('Insufficient knots to be a proper spline basis')
end

if (nargin>3)
  B = B*v(:);
end


return


function B = recu(t,i,k,u)

if (k==1)

  if (u(i)==u(i+1))
    B = zeros(size(t));
  else
    B = (u(i)<=t) & (t<u(i+1));
  end

else
  
  B = w(t,i,k,u).*recu(t,i,k-1,u)+(1-w(t,i+1,k,u)).*recu(t,i+1,k-1,u);
  
end





function wt = w(t,i,k,u)

if (u(i)~=u(i+k-1))
  wt = (t-u(i))./(u(i+k-1)-u(i));
else
  wt = 0;
end

