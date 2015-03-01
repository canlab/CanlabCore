function Ainv = incremental_invert(Pinv,Q,S)  
%
%   Computes the new inverse of a symmetric matrix from an old submatrix inverse and the new parts. 
%
%   So if   A = [P,Q;Q',S] incremental_invert gets inv(P) , Q and S and
%   returns inv(A)
%                           
%

if isempty(S) && isempty(Q)
    Ainv = Pinv; return
end

Qt = Q';
M = S-Qt*Pinv*Q;
Minv = pinv(M);

Ptmp = Pinv + Pinv*Q*Minv*Qt*Pinv;
Qtmp = -Pinv*Q*Minv;

Ainv = [Ptmp,Qtmp;Qtmp',Minv];