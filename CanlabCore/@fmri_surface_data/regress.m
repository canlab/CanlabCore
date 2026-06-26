function bobj = regress(obj, varargin)
% regress Grayordinate-wise OLS regression of the maps onto a design matrix.
%
% :Usage:
% ::
%     bobj = regress(obj)              % uses obj.X
%     bobj = regress(obj, X)           % explicit design matrix
%
% Fits an ordinary least-squares regression at each grayordinate: for each
% grayordinate, the values across maps (columns of obj.dat) are regressed onto
% the design matrix X ([nMaps x p], one row per map/observation). Returns an
% fmri_surface_data whose maps are the p regression coefficients, with the
% t-statistics, p-values, and standard errors in .additional_info.statistic.
%
% This is a lightweight native implementation (no external toolbox). For the full
% GLM machinery (contrasts, diagnostics), export the volumetric part with
% to_fmri_data and use fmri_data.regress / glm_map.
%
% :Inputs:
%   **obj:** fmri_surface_data with multiple maps.
%   **X:**   (optional) [nMaps x p] design matrix. Default: obj.X. An intercept
%            is NOT added automatically -- include a column of ones if desired.
%
% :Outputs:
%   **bobj:** fmri_surface_data; .dat is [nGrayordinates x p] coefficients.
%             .additional_info.statistic has .t, .p, .se, .dfe (each [nGray x p]
%             or scalar dfe), .betas.
%
% :Examples:
% ::
%     obj.X = [ones(n,1) age(:)];      % intercept + age
%     b = regress(obj);
%     surface(get_wh_image(b, 2));     % render the age effect (2nd coefficient)
%
% :See also: ttest, predict, to_fmri_data, fmri_data.regress

if nargin >= 2 && ~isempty(varargin{1}) && isnumeric(varargin{1})
    X = varargin{1};
else
    X = obj.X;
end
if isempty(X)
    error('fmri_surface_data:regress:noX', ...
        'No design matrix. Set obj.X or pass X as the second argument.');
end

nMaps = size(obj.dat, 2);
if size(X, 1) ~= nMaps
    error('fmri_surface_data:regress:size', ...
        'Design matrix has %d rows but the object has %d maps.', size(X,1), nMaps);
end

p = size(X, 2);
dfe = nMaps - p;
if dfe < 1
    error('fmri_surface_data:regress:df', 'Not enough maps (%d) for %d regressors.', nMaps, p);
end

Y = double(obj.dat).';            % [nMaps x nGray]
B = X \ Y;                        % [p x nGray]
resid = Y - X * B;
sigma2 = sum(resid.^2, 1) / dfe;  % [1 x nGray]
XtXinv = inv(X' * X);             %#ok<MINV> small p x p
se = sqrt(diag(XtXinv) * sigma2); % [p x nGray]
tvals = B ./ se;                  % [p x nGray]
pvals = 2 * tcdf(-abs(tvals), dfe);

bobj = rebuild_like(obj, B.');    % [nGray x p] coefficients
bobj.intent = 'dscalar';
bobj.image_names = arrayfun(@(k) sprintf('beta%d', k), 1:p, 'UniformOutput', false)';

stat = struct('betas', B.', 't', tvals.', 'p', pvals.', 'se', se.', 'dfe', dfe);
bobj.additional_info.statistic = stat;
bobj.history{end+1} = sprintf('regress: OLS on %d-column design, %d maps (betas in .dat; t/p in additional_info.statistic)', p, nMaps);
end
