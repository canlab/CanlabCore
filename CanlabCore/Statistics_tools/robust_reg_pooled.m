function [betas,w] = robust_reg_pooled(X,Y)
% :Usage:
% ::
%
%     [betas,w] = robust_reg_pooled(X,Y)
%     [betas,stats] = weighted_reg(X,Y,'w',w,'uni');

nancols = find(any(isnan(Y) | Y==0,1)); % Get rid of missing values
Y(:,nancols) = [];

[m,n] = size(Y);
w = ones(m,1);

% pool weights across all voxels
[betas,invxwx,bform,fits] = get_betas_singleweight(X,Y,w);

% --------------------------------------
% * Weights and computational setup
% --------------------------------------

% this is done only once in robust regression (not iterated)
lev = diag(X * bform);  % leverage, diagonals of hat matrix
w = 1./lev;                         % starting weights, 1 / leverages

h = min(.9999, sum(lev.*lev,2));

adjfactor = 1 ./ sqrt(1-h);
adjfactor = repmat(adjfactor,1,n);
xrank = rank(X);                    % should recompute (?), but this is rough...

%    tiny_s = 1e-6 * mean(std(Y),2);
% if tiny_s==0
%     tiny_s = 1;
% end

% --------------------------------------
% * Set up iterations
% --------------------------------------

D = sqrt(eps(class(X)));
tol = max(D*max(abs(betas)));
ntol = n; %round(n .* .95);         % number of regressions that should have converged before we proceed

iter = 0;
nconverged = 0;
iterlim = 50;
str = sprintf('Target: %3.0f%% converged.  Iterating %03d',100*ntol./n,0); disp(str);

% --------------------------------------
% * Iterate to convergence
% --------------------------------------

while (iter==0) || (nconverged < ntol)
    iter = iter+1;
    fprintf(1,'\b\b\b\b%03d',iter);
    str2 = sprintf('%3.0f%% converged',round(100*nconverged./n)); disp(str2);
    
    if (iter>iterlim)
        warning('Iteration limit reached.');
        break;
    end

    % Compute residuals from previous fit, then compute weights
    r = Y - fits;
    radjust = r .* adjfactor;
    w = bisquare_weight(r,radjust,xrank);
    w = mean(w,2);      % pool weights over voxels

    oldbetas = betas;
    % re-fit betas based on new weights
    [betas,invxwx,bform,fits] = get_betas_singleweight(X,Y,w);

    nconverged = sum(1 - any(abs(betas-oldbetas) > tol));   % number of regressions converged
    erase_string(str2);
end

erase_string(str);

return





% --------------------------------------
%
% * Sub-functions
%
% --------------------------------------



function [betas,invxwx,bform,fits] = get_betas_singleweight(X,Y,w)

W = diag(w);                    % Weight matrix

%X = repmat(1,m,1);              % Design matrix - 1 column of all ones to calculate average
% and, separately, use bcon if that's entered

invxwx = inv(X'*W*X);
bform = invxwx * X'* W;         % beta-forming matrix.  hat = X * bform

% rows are columns of design (X), cols are Y variables
betas = bform*Y;

if nargout > 3
    fits = X * betas;
end

return



function w = bisquare_weight(r,radjust,xrank)
% r is residuals
% radjust is adjustment factor: DuMouchel & O'Brien
% xrank is rank of weighted X matrix (design)
% w is weights from bisquare function
% n is number of Y variables to replicate weights over
tuneconst = 4.685;

r = r .* radjust;
s = mad_sigma_pooled(r,xrank);
r = r ./ (s*tuneconst);
w = (abs(r)<1) .* (1 - r.^2).^2;
return


function s = mad_sigma_pooled(r,xrank)
%    Compute std estimate using MAD of residuals from 0
rsort = sort(abs(r));
rsort = rsort(xrank:end,:); % eliminate smallest; like reducing df
s = median(rsort(:)) / 0.6745;
return

function str = display_string(str)
str = sprintf(str); fprintf(1,'%s',str);
return


function erase_string(str)

len = length(str);
str2 = repmat('\b',1,len);

fprintf(1,str2);

return
