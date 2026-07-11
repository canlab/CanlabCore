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

% Expand the basis-set array to one entry per condition (filled with the
% current default), so that replacing one condition's basis leaves the others
% unchanged regardless of whether the model has been built yet. (build's
% check_model only replicates a length-1 xBF, which would otherwise overwrite
% all conditions with the replacement.)
nconds = condition_num;
if ~isempty(obj.Sess) && isfield(obj.Sess, 'U') && ~isempty(obj.Sess(1).U)
    nconds = max(condition_num, length(obj.Sess(1).U));
end
for c = length(obj.xBF) + 1 : nconds
    obj.xBF(c) = oldBF;
end

% SPM adds these things here, so we will too, for consistency. Inherit them
% from the existing basis when the replacement does not provide them.
for f = {'T', 'T0', 'UNITS', 'Volterra'}
    if ~isfield(xBF_hires, f{1}) && isfield(oldBF, f{1})
        xBF_hires.(f{1}) = oldBF.(f{1});
    end
end

% Assigning into the struct array obj.xBF(condition_num) requires identical
% field sets and order. Different basis-set builders (spm_get_bf,
% fmri_spline_basis, ...) return different fields, so harmonize first: add any
% missing fields (filled with []) to both sides, then match the ordering.
allflds = union(fieldnames(obj.xBF), fieldnames(xBF_hires), 'stable');
obj.xBF    = local_ensure_fields(obj.xBF, allflds);
xBF_hires  = local_ensure_fields(xBF_hires, allflds);
xBF_hires  = orderfields(xBF_hires, obj.xBF(1));

obj.xBF(condition_num) = xBF_hires;

end


% =========================================================================
function s = local_ensure_fields(s, flds)
% Add any missing fields (default []) to every element of struct array s,
% then reorder fields to match flds.
for k = 1:numel(flds)
    if ~isfield(s, flds{k})
        [s.(flds{k})] = deal([]);
    end
end
s = orderfields(s, flds);
end

