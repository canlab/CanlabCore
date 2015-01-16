% Fisher's r to z' transform, and the inverse
% [z, r] outputs:
% z = z', treating input r as correlation
% r = r, treating input r as a z' score
function [z, r] = fisherz(r)
    z = .5 .* log((1+r) ./ (1-r));     % Fisher's r-to-z transform

    if nargout > 0
        r = (exp(2.*r) - 1) ./ (exp(2.*r) + 1);    % inverse
    end
end