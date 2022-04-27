function [Xs, delta, C, B, names] = get_session_X(obj, s)
% Get design matrix (predictors) for one session of fmri_model object, using
% basis functions defined in the object and onsets for one session (s).
%
% :Usage:
% ::
%
%     [Xs, delta, C, B, names] = get_session_X(obj, session number)
%

% ..
%    Define sessions and number of conditions
% ..

nsess = length(obj.Sess);

if s > nsess, error('Session %3.0f does not exist', s); end
    

TR = obj.xY.RT;

nconds = length(obj.Sess(s).U);

[ons, name] = deal(cell(1, nconds));

[ons{:}] = deal(obj.Sess(s).U(:).ons);
[name{:}] = deal(obj.Sess(s).U(:).name);

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
% Predictors
% ----------------------------------------------

delta = onsets2delta(ons, obj.nscan(s));

% Of-interest part of design matrix
% time res is defined as TR / 16, so build and downsample by 16 to TR
% Allow for different basis sets for each condition.

Xs = cell(1, nconds);

for i = 1:nconds
    
    bf = obj.xBF(i).bf;
    
    Xs{i} = getPredictors(delta(:, i), bf, 16);
end

Xs = cat(2, Xs{:});

% covariate part of design matrix
C = [];
if ~isempty(obj.Sess(s).C) && isfield(obj.Sess(s).C, 'C') 
    C = obj.Sess(s).C.C;
end

% baseline
% ----------------------------------------------

B = ones(obj.nscan(s), 1);

% ----------------------------------------------
% modulators
% ----------------------------------------------

Xs_pm = cell(1, nconds);

for i = 1:nconds
    is_pm = isfield(obj.Sess(s).U(i), 'P') && ~isempty(obj.Sess(s).U(i).P) && isfield(obj.Sess(s).U(i).P, 'P') && ~isempty(obj.Sess(s).U(i).P.P);
    
    if is_pm
        
        pm_vals = {obj.Sess(s).U(i).P.P};
        
        model = onsets2parametric_mod_X(ons(i), pm_vals, obj.nscan(s), obj.xBF(i).bf, 16);
        
        Xs_pm{i} = model;

    end

end

Xs_pm = cat(2, Xs_pm{:});
       
% add to Xs
Xs = [Xs Xs_pm];

% ----------------------------------------------
% get names for each BF
% Allow variable number of basis fcns for each condition
% ----------------------------------------------

names = {};

for i = 1:nconds
    
    nbf = size(obj.xBF(i).bf, 2);
    
    for j = 1:nbf
        myname = [name{i} ' BF' num2str(j)];
        myname = replaceblanks(myname);

        names{end+1} = myname;
    end
end

all_pmnames = {};

for i = 1:nconds
    is_pm = isfield(obj.Sess(s).U(i), 'P') && ~isempty(obj.Sess(s).U(i).P) && isfield(obj.Sess(s).U(i).P, 'P') && ~isempty(obj.Sess(s).U(i).P.P);
    
    if is_pm
        % Get param mod names
        % ------------------------------------------------
        pmnames = cell(1, nbf);
        
        for j = 1:nbf
            myname = [obj.Sess(s).U(i).P.name ' BF' num2str(j)];
            myname = replaceblanks(myname);
            
            pmnames{1, j} = myname;
        end
        
        all_pmnames = cat(2, all_pmnames, pmnames{:});
        
    end
end

names = names(:)';

names = [names all_pmnames];

end

% ----------------------------------------------
% ----------------------------------------------


% ----------------------------------------------
% ----------------------------------------------


function myname = replaceblanks(myname)

myname(myname == ' ') = '-';
% this works in 2010b but not a...
%tmp = diff(num2str(myname)) == 0;

tmp = diff(double(myname)) == 0;
myname([false tmp] & myname == '-') = [];

end
        
