function [paths, stats] = mediation(D, xvarname, yvarname, mvarname, varargin)
% Run single or multilevel mediation analysis on a canlab_dataset object.
% Calls mediation.m (see mediation.m) in mediation toolbox
%
% :Usage:
% ::
%
%    [paths, stats] = mediation(D, xvarname, yvarname, mvarname, [optional inputs])
%
% ..
%     Author and copyright information:
%
%     Copyright (C) 2013 Tor Wager
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
% ..
%
% :Inputs:
%
%   **D:**
%        is a canlab_dataset object
%
%   **xvarname:**
%        X, the initial variable (valid variable name in the dataset)
%
%   **yvarname:**
%        Y, the outcome variable (valid variable name in the dataset)
%
%   **mvarname:**
%        M, the potential mediator ((valid variable name in the dataset)
%
%
% :Optional Inputs:
%
%   Takes any optional inputs to mediation.m
%   e.g., 'noverbose', 'dosave', 'names', 'M', 'L2M', 'covs', others
%
%   **wh_keep:**
%        followed by 1/0 vector of subjects to keep.
%           - must be same length as subjects
%           - subjects with value 0 will be excluded
%
%   **rankdata:**
%        ranks all data before mediation; "Nonparametric"
%
%
% :Outputs:
%
%   **paths:**
%        see mediation.m from mediation toolbox
%
%   **stats:**
%        see mediation.m from mediation toolbox
%
% :Examples:
% ::
%
%    [paths, stats] = mediation(D, 'Group', 'DeltaDon', 'DeltaDist', 'M2', 'DeltaTend', 'wh_keep', wh_keep);
%



covstr = 'nocovs';
covvarnames = {};
c = [];

m2str = 'noM2';
m2 = [];
m2names = {};
m2varsdescrip = {};

nboot = 10000;
singlelevel = 0;  % force single-level; update later...

%wh_keep = true(size(D.Subj_Level.data(:, 1)));

%
% % Get levels of main variables
% % ----------------------------------------------
% varlevel = [0 0 0];
% varlevel(1) = get_varlevel(D, xvarname);
% varlevel(2) = get_varlevel(D, mvarname);
% varlevel(3) = get_varlevel(D, yvarname);

[x, xcell, varlevel(1), xvardescrip] = get_var(D, xvarname);
[m, mcell, varlevel(2), mvardescrip] = get_var(D, mvarname);
[y, ycell, varlevel(3), yvardescrip] = get_var(D, yvarname);





% Select level and data
% -------------------------------------------------
if any(varlevel == 1)
    % enforce single-level
    singlelevel = 1;
    x = nanmean(x, 2);
    m = nanmean(m, 2);
    y = nanmean(y, 2);
else
    % Multi-level
    x = xcell;
    m = mcell;
    y = ycell;
end

wh = strcmp(varargin, 'wh_keep');
if any(wh)
    wh_keep = varargin{find(wh)+1};
    
    if length(wh_keep) ~= length(x) || ~islogical(wh_keep);
        error('wh_keep must be logical vector of which subjects to keep.  wrong size or datatype.');
    end
    
    x = x(wh_keep);
    y = y(wh_keep);
    m = m(wh_keep);
    
end

% Get covariates and 2nd mediators
% ----------------------------------------------

wh = find(strcmp(varargin, 'covs'));
if any(wh)
    covstr = 'covs';
    for i = 1:length(wh)
        % WILL HANDLE ONLY ONE COVARIATE NOW...
        covvarnames{i} = varargin{wh+1};
        [c, ccell, cvarlevel(i)] = get_var(D, covvarnames{i});
        
        if singlelevel
            c = nanmean(c, 2); % <- works for both Subject/Event level var
        else
            c = ccell;
        end
        
        if exist('wh_keep', 'var')
            c = c(wh_keep);
        end
        
    end
end

wh = find(strcmp(varargin, 'M2'));
if any(wh)
    m2str = 'M2';
    for i = 1:length(wh)
        m2names{i} = varargin{wh+1};
        [m2, m2cell, m2varlevel(i), m2varsdescrip] = get_var(D, m2names{i});
        
        if singlelevel
            m2 = nanmean(m2, 2); % <- works for both Subject/Event level var
        else
            m2 = m2cell;
        end
        
        if exist('wh_keep', 'var')
            m2 = m2(wh_keep);
        end
    end
end

% Rank data, if requested
% ----------------------------------------------
if any(strcmp(varargin, 'rankdata'))
    
    if singlelevel
        
        x = rankdata(x);
        y = rankdata(y);
        m = rankdata(m);
        c = rankdata(c);
        m2 = rankdata(m2);
        
    else
        cfun = @(x) rankdata(x);
        x = cellfun(cfun, x, 'UniformOutput', 0);
        y = cellfun(cfun, y, 'UniformOutput', 0);
        m = cellfun(cfun, m, 'UniformOutput', 0);
        
        if ~isempty(c)
            c = cellfun(cfun, c, 'UniformOutput', 0);
        end
        
        if ~isempty(m2)
            m2 = cellfun(cfun, m2, 'UniformOutput', 0);
        end
        
        
    end
    
end



% Names
% ----------------------------------------------
names = [{xvardescrip} {yvardescrip} {mvardescrip} {m2varsdescrip{:}}];

[paths, stats] = mediation(x, y, m, 'boot', 'plots', 'verbose', 'bootsamples', nboot, 'names', names, 'covs', c, 'M', m2, varargin{:});


end % function


%
% function varlevel = get_varlevel(D, varname)
%
% [varlevel, varlevel2] = deal(0);
%
% varlevel = any(strcmp(D.Subj_Level.names, varname));
%
% varlevel2 = 2*any(strcmp(D.Event_Level.names, varname));
%
% if any(varlevel & varlevel2)
%     error('Some variables are in both Subject and Event levels.');
% end
%
% varlevel = varlevel + varlevel2;
%
% if any(~varlevel)
%     error(['Variable ' varname ' does not exist in this dataset. Check var input names.'])
% end
%
% end
