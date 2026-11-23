function varargout = rsa_compare_models(dat, varargin)
% rsa_compare_models  Thin alias for rsa_compare_lme_models.
%
% See rsa_compare_lme_models.m for full documentation.
[varargout{1:nargout}] = rsa_compare_lme_models(dat, varargin{:});
end
