function obj = build_single_trial(obj, inputhrf)
% Build a single-trial design matrix (xx) for an fmri_model object
%
% We assume that the same conditions are modeled for each session
% We assume that you have one basis set per condition (this is different
% from SPM, which only allows a single basis set across all conditions)
%
% This is used in single_trial_estimates, which assumes that you have
% estimated an initial model and saved image data.
%
% The idea behind this is somewhat different from other canlab single-trial
% analyses, in that it takes in a single, custom HRF for each condition,
% rather than using a basis set.  In single_trial_estimates, custom HRFs
% are created for each voxel by using the condition- and voxel-specific hrf
% estimates stored during model fitting.
%
% The sequence would be:
%   1. robustfit(my_model), to fit average model and get HRF est for each
%      voxel
%   2. single_trial_estimates(my_model), to use this function to build
%      single-trial design matrices and fit them.
%
% :Usage:
% ::
%
%     obj = build_single_trial(fmri_model_obj, inputhrf)
%
% :Inputs:
%
%   **inputhrf:**
%        should be a cell array of length nconds (number of conditions).


obj = check_model(obj); % Check assumptions and basis set
% We assume that the same conditions are modeled for each session

% ----------------------------------------------
% Define sessions and number of conditions
% ----------------------------------------------

nsess = length(obj.Sess);
nconds = length(obj.Sess(1).U);
TR = obj.xY.RT;

[sess_delta,  sess_conditions,  sess_X, sess_C, sess_B] = deal(cell(1, nsess)); 


for s = 1:nsess

[ons, name] = deal(cell(1, nconds));

[ons{:}] = deal(obj.Sess(s).U(:).ons);
%[name{:}] = deal(obj.Sess(s).U(:).name);

% make sure onsets are in TRs, not secs
switch obj.xBF(1).UNITS
    case 'secs'
        % onsets are in sec, convert
        for i = 1:nconds
            ons{i} = ons{i} ./ TR;
        end
        
    case {'tr', 'TR', 'trs'}
        % ok, do nothing
end

% ----------------------------------------------
% Build predictors for each session
% ----------------------------------------------

delta = onsets2delta(ons, obj.nscan(s));

ns = size(delta, 1);

[trialdelta, cond_assignment] = deal(cell(1, nconds));

% ----------------------------------------------
% For each condition, parse into separate columns
% for each onset (single-trial)
% ----------------------------------------------
for j = 1:nconds
    
    wh = find(delta(:, j));
    
    trialdelta{j} = false(ns, length(wh));
    
    for k = 1:length(wh)
        trialdelta{j}(wh(k), k) = true;
    end
    
    cond_assignment{j} = ones(length(wh), 1);
    
    bf = obj.xBF(j);
    xvals = 0:TR:bf.length - TR;
    hrf = interp1(xvals, inputhrf{j}(1:length(xvals)), 1:bf.dt:bf.length, 'spline', 'extrap');  
    
    condX{j} = getPredictors(trialdelta{j}, hrf, 16);
    
end % conditions

sess_delta{s} = cat(2, trialdelta{:}); % not needed?

sess_conditions{s} = blkdiag(cond_assignment{:});

% design matrix - of interest
sess_X{s} = cat(2, condX{:});

% design matrix - nuisance
sess_C{s} = obj.Sess(s).C.C;

sess_B{s} = ones(obj.nscan(s), 1);

end  % Session

% ----------------------------------------------
% Put the pieces together
% ----------------------------------------------

obj.xX.X = [blkdiag(sess_X{:}) blkdiag(sess_C{:}) blkdiag(sess_B{:})];

obj.xX.cond_assignments = cat(1, sess_conditions{:});

obj.xX.iH = find(any(obj.xX.cond_assignments, 2));
obj.xX.iC = find(~any(obj.xX.cond_assignments, 2));
obj.xX.iB = size(obj.xX.X, 2) - nsess : size(obj.xX.X, 2);


end % function



% ----------------------------------------------
% ----------------------------------------------

% SUB-functions

% ----------------------------------------------
% ----------------------------------------------


function myname = replaceblanks(myname)

myname(myname == ' ') = '-';
tmp = diff(double(myname)) == 0; % fixed 2010a vs. b compat.
myname([false tmp] & myname == '-') = [];

end


function obj = check_model(obj)

nsess = length(obj.Sess);


nconds = length(obj.Sess(1).U);
if nconds == 0
    disp('You must assign conditions and onsets using fmri_model or add (method)')
    disp('before building the design matrix.');
    error('Exiting');
end
for i = 2:nsess
    if any(length(obj.Sess(1).U) ~= nconds)
        disp('The number of conditions must be the same for all sessions.')
        error('Exiting');
    end
end

% replicate xBF if necessary
if length(obj.xBF) == 1
    obj.xBF(2:nconds) = obj.xBF(1);
elseif length(obj.xBF) < nconds
    disp('The number of basis set structures is > 1 but not equal to the number of conditions.')
    error('Exiting');
end

for i = 1:nconds
    obj.xBF(i).name = sprintf('%s for Condition %3.0f', obj.xBF(i).name, i);
end

end



function obj = concat_names(obj)

% names
% ---------------------------------------
[obj.xX(1).name{:}] = deal(obj.Sess(1).U(:).name);

% add PM names
s = 1;
for i = 1:length(obj.Sess(1).U)
    is_pm = isfield(obj.Sess(s).U(i), 'P') && ~isempty(obj.Sess(s).U(i).P) && isfield(obj.Sess(s).U(i).P, 'P') && ~isempty(obj.Sess(s).U(i).P.P);
    if is_pm
        obj.xX(1).name{end + 1} = [obj.Sess(s).U(i).name '-' obj.Sess(s).U(i).P.name];
    end
end

% Remove session labels from names
for i = 1:length(obj.xX(1).name)
    myname = replaceblanks(obj.xX(1).name{i});
    wh = findstr(myname, 'Sess-');
    if ~isempty(wh)
        to_remove = [];
        for j = 1:length(wh)
            to_remove = [to_remove wh(j) : wh(j) + 6];
        end
        myname(to_remove) = [];
    end
    
    obj.xX(1).name{i} = myname;
end

end

