function pm = fit_predictive_model(obj, pm, varargin)
% fit_predictive_model  Fit a @predictive_model on this fmri_data.
%
% Bridge from the neuroimaging side (@fmri_data) to the numeric
% predictive_model API. Extracts X = obj.dat' (n_images x n_voxels),
% Y = obj.Y, and optionally a grouping vector; then forwards to
% pm.fit() or pm.crossval() depending on what the caller asked for.
%
% :Usage:
% ::
%     pm = predictive_model('algorithm','svm','task','classification');
%     pm = fit_predictive_model(my_fmri_data, pm);
%     pm = fit_predictive_model(my_fmri_data, pm, 'crossval', true);
%     pm = fit_predictive_model(my_fmri_data, pm, ...
%                                'crossval', true, ...
%                                'cv', cv_splitter.group_kfold(10), ...
%                                'groups', subject_id);
%
% :Optional Inputs (name/value):
%   'crossval'  default false. If true, calls pm.crossval(); otherwise
%               pm.fit() (in-sample).
%   'cv'        cv_splitter object; forwarded if crossval=true.
%   'groups'    grouping vector; if a scalar string 'id' and
%               obj.metadata_table has a column of that name, that
%               column is used.
%   'scoring'   scorer name or cv_scorer object.

    p = inputParser; p.KeepUnmatched = true;
    addParameter(p, 'crossval', false);
    addParameter(p, 'cv',       []);
    addParameter(p, 'groups',   []);
    addParameter(p, 'scoring',  []);
    parse(p, varargin{:});

    docv     = p.Results.crossval;
    cv       = p.Results.cv;
    groups   = p.Results.groups;
    scoring  = p.Results.scoring;

    % --- Extract X, Y ---
    if isempty(obj.dat)
        error('fmri_data:fit_predictive_model:NoData', '.dat is empty.');
    end
    if isempty(obj.Y)
        error('fmri_data:fit_predictive_model:NoY', '.Y is empty.');
    end
    X = double(obj.dat)';     % images x voxels
    Y = double(obj.Y);

    % --- Resolve groups ---
    if ischar(groups) || isstring(groups)
        col = char(groups);
        if isprop(obj, 'metadata_table') && istable(obj.metadata_table) ...
                && any(strcmpi(obj.metadata_table.Properties.VariableNames, col))
            groups = obj.metadata_table.(col);
            if iscell(groups), groups = grp2idx(groups); end
        else
            error('fmri_data:fit_predictive_model:NoGroupCol', ...
                'No metadata_table column named ''%s''.', col);
        end
    end

    % --- Dispatch ---
    if docv
        args = {'groups', groups};
        if ~isempty(cv),      args = [args, {'cv', cv}]; end
        if ~isempty(scoring), args = [args, {'scoring', scoring}]; end
        pm = crossval(pm, X, Y, args{:});
    else
        pm = fit(pm, X, Y, 'id', groups);
    end
end
