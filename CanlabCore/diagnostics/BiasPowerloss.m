function [b bias pl Pc Pe] = BiasPowerloss(tc, X, c, beta, df, z, pval)
% Calculate the approximate bias and power loss due to mis-modeling
% This works with the Mismodeling Toolbox described by Loh et al. 2008
%
% :Usage:
% ::
%
%     [b bias pl Pc Pe] = BiasPowerloss(tc, X, c, beta, df, z, pval)
%
% :Inputs:
%
%   **tc:**
%        fMRI time course
%
%   **X:**
%        design matrix for multiple regression
%
%   **c:**
%        contrast of interest
%
%   **beta:**
%        (mismodeled) beta value
%
%   **df:**
%        degrees of freedom
%
%   **z:**
%        p-value calculated from ResidScan
%
%   **pval:**
%        cut-off p-value
%
% :Outputs:
%
%
%   **b:**
%        updated (correct) beta value
%
%   **bias:**
%        bias
%
%   **pl:**
%        power loss
%
% :References:
%   Loh, J. M., Lindquist, M. A., Wager, T. D. (2008). Residual Analysis for Detecting Mis-modeling in fMRI. Statistica Sinica, 18, 1421-1448.
%
% ..
%    By Martin Lindquist & Ji-Meng Loh, July 2007
% ..
%



% Update design matrix using correct model
if (z< pval),
    [Gamma b] = EditBasis(X, tc, z, pval);
else
    b = beta;
end;

bias = beta-b;      % Calculate bias

% Power loss
tstar = tinv(1-pval,df);       
Pc = 1- nctcdf(tstar,df,c*b);       % Power for correct model

% Approximate power using incorrect model
Pe = 1- nctcdf(tstar,df,c*(beta));       
        
pl = max(0, Pc-Pe); % Calculate power loss

end

% END MAIN FUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Subfunctions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Gamma b] = EditBasis(X, tc, z, pval)

% use residuals to update design matrix

b = pinv(X)*tc;
XX = X;

maxiter = 5;
niter = 0;
while (z < pval || niter < maxiter),

    e = tc - XX*b;
    [z sres v] = ResidScan(e, 4);

    v(abs(sres)<1.645) =0;
    XX(:,2) = XX(:,2) + v;
    XX(:,2) = XX(:,2)./max(XX(:,2));
    b = pinv(XX)*tc;
    niter = niter +1;
    
end

Gamma = XX-X;

end
