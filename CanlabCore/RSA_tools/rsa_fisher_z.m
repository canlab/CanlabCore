function z = rsa_fisher_z(r)
% rsa_fisher_z  Safe Fisher r-to-z transform for averaging/testing correlations.
%
% z = atanh(r), with r clamped to (-1, 1) so that trivial self-correlations
% (r = 1) and tiny round-off overshoots do not produce +/-Inf or complex values.
%
% USE THIS ONLY ON CORRELATION-METRIC similarities. Dot-product and cosine
% similarities are NOT bounded in [-1, 1] (cosine can round slightly past 1;
% dot products are unbounded), so atanh of those is meaningless or complex.
% Compute statistics in z-space; report descriptive means back in r-space with
% tanh(mean_z) if desired.
%
% :Usage:
% ::
%   z = rsa_fisher_z(r)
%
% :Inputs:
%   **r:** scalar / vector / matrix of correlation coefficients.
%
% :Outputs:
%   **z:** Fisher-z transformed values, same size as r. NaNs are preserved.
%
% :Examples:
% ::
%   z = rsa_fisher_z([0.041 0.085 1.0]);   % last value clamped, finite
%   r_back = tanh(mean(z, 'omitnan'));     % average in z-space, report in r-space
%
% ..
%    2026 Michael Sun. GPL v3.
% ..

r = min(max(r, -0.999999), 0.999999);   % clamp away from +/-1
z = atanh(r);

end
