function [z, r] = fisherz(r)
% Fisher's r to z' transform, and the inverse
%
% :Outputs:
%
%   **z:**
%        z = z', treating input r as correlation
%
%  **r:**
%        treating input r as a z' score

    z = .5 .* log((1+r) ./ (1-r));     % Fisher's r-to-z transform

    if nargout > 0
        r = (exp(2.*r) - 1) ./ (exp(2.*r) + 1);    % inverse
    end
end
