function obj = replace_basis_set(obj, condition_num, xBF_hires)
% Replace a basis set in an fmri_model object with another one of your
% choosing.
%
% This allows one to use a custom basis set, and also to use different
% basis sets for different trial types.
%
% Each condition across all sessions must be modeled with the same basis
% set. That is, there can be only one basis set per condition, e.g., one
% for anticipation (used in each session) and one for pain.
%
% :Usage:
% ::
%
%     obj = replace_basis_set(obj, condition_num, xBF_hires)
%
% :Examples:
% ::
%
%    % generate a custom spline basis set and use that for Condition 1,
%    % and the standard one for Condition 2:
%
%    [xBF_hires, xBF] = fmri_spline_basis(2, 'length', 12, 'nbasis', 3, 'order', 3, 'plot');
%
%    %save this to get info that is not typically in basis set until after
%    %model is built.

oldBF = obj.xBF(1); 

% SPM adds these things here, so we will too, for consistency
xBF_hires.T = oldBF.T;
xBF_hires.T0 = oldBF.T0;
xBF_hires.UNITS = oldBF.UNITS;
xBF_hires.Volterra = oldBF.Volterra;

obj.xBF(condition_num) = xBF_hires;

end

