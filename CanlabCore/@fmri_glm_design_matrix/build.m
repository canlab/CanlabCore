function obj = build(obj)
% Build the design matrix (xx) for an fmri_model object
%
% We assume that the same conditions are modeled for each session.
% We assume that you have one basis set per condition (this is different
% from SPM, which only allows a single basis set across all conditions)
%
% :Usage:
% ::
%
%     obj = build(fmri_model_obj)
%


nsess = length(obj.Sess); % Define sessions and number of conditions

% add predictors for each session
[X, C, B, names] = deal(cell(1, nsess));

% Check assumptions and basis set
% We assume that the same conditions are modeled for each session
obj = check_model(obj);

% Do most of the work here; get regressors, and modulators, and names
for s = 1:nsess
    [X{s}, obj.Sess(s).delta, C{s}, B{s}, names{s}]  = get_session_X(obj, s);
end


switch obj.build_method
    case 'Separate sessions'
        % of-interest part of design
        iH = blkdiag(X{:});
        
        obj.xX(1).name = cat(2, names{:});
        
    case 'Concatenated sessions'
        iH = cat(1, (X{:}));
        
        nconds = length(obj.Sess(1).U);
        obj.xX(1).name = cell(1, nconds);
        
        obj = concat_names(obj);
        
        
    otherwise
        error('Unknown build method for fmri_model object. See help for valid strings.');
        
end

% of-interest part of design
%iH = blkdiag(X{:});
iC = blkdiag(C{:});
iG = [];
iB = blkdiag(B{:});

obj.xX(1).X = [iH iC iG iB];

% indices for each partition
niH = size(iH);
niC = size(iC);
niG = size(iG);
niB = size(iB);

sz = [niH; niC; niG; niB];
if any(sz(:, 1) ~= sz(1) & sz(:, 1) ~= 0)
    disp('Number of rows in interest, covs, global, and baseline partitions do not match.');
    disp('Debug me!')
    keyboard
end

wh = [ones(1, sz(1, 2)) 2*ones(1, sz(2, 2)) 3*ones(1, sz(3, 2)) 4*ones(1, sz(4, 2))];
obj.xX(1).iH = find(wh == 1);
obj.xX(1).iC = find(wh == 2);
obj.xX(1).iG = find(wh == 3);
obj.xX(1).iB = find(wh == 4);

% Condition assignments matrix
obj.xX(1).cond_assignments = get_condition_assignments(obj);


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

