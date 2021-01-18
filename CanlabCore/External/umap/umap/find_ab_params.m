function [a,b] = find_ab_params(spread, min_dist)
%FIND_AB_PARAMS Fit a and b parameters for the differentiable curve used in
% lower dimensional fuzzy simplicial complex construction. We want the
% smooth curve (from a pre-defined family with simple gradient) that best
% matches an offset exponential decay.
%
% [a,b] = FIND_AB_PARAMS(spread, min_dist)
%
% Parameters
% ----------
% spread: double
%     The effective scale of embedded points.
%
% min_dist: double
%     The effective minimum distance between embedded points. Smaller values
%     will result in a more clustered/clumped embedding where nearby points
%     on the manifold are drawn closer together, while larger values will
%     result on a more even dispersal of points. The value should be set
%     relative to the "spread" value, which determines the scale at which
%     embedded points will be spread out.  
% 
% Returns
% -------
% a: double
%     Parameter of differentiable approximation of right adjoint functor.
% 
% b: double
%     Parameter of differentiable approximation of right adjoint functor.
%
%   AUTHORSHIP
%   Math Lead & Primary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Secondary Developer: Stephen Meehan <swmeehan@stanford.edu>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Provided by the Herzenberg Lab at Stanford University 
%   License: BSD 3 clause
%

    curve = @(a,b,x) (1./(1 + a*x.^(2*b)));

    xv = linspace(0, 3*spread, 300);
    yv = (xv < min_dist) + ~(xv < min_dist).*exp(-(xv - min_dist) / spread);

    f = fit(xv', yv', curve);
    params = coeffvalues(f);
    a = params(1);
    b = params(2);
end