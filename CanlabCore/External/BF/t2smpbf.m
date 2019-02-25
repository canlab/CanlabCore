function bf10 = t2smpbf(t,nx,ny,r)
%
% bf10 = t2smpbf(t,n,[r=0.707])
%
% Calculates JZS Bayes Factor for a two-sample t-test given t and sample sizes nx and ny.
% The optional input r is the scale factor which defaults to 0.707.
% This quantifies the evidence in favour of the alternative hypothesis. 
% See Rouder et al, 2009, Psychon Bull Rev for details.
%

% Default scale factor
if nargin < 4
    r = 0.707;
end

% Function to be integrated
F = @(g,t,nx,ny,r) (1+(nx*ny/(nx+ny)).*g.*r.^2).^(-1./2) .* (1 + t.^2./((1+(nx*ny/(nx+ny)).*g.*r.^2).*(nx+ny-2))).^(-(nx+ny-1)./2) .* (2.*pi).^(-1./2) .* g.^(-3./2) .* exp(-1./(2.*g));

% Bayes factor calculation
bf01 = (1 + t^2/(nx+ny-2))^(-(nx+ny-1)/2) / integral(@(g) F(g,t,nx,ny,r),0,Inf);

% Invert Bayes Factor
bf10 = 1 / bf01;
