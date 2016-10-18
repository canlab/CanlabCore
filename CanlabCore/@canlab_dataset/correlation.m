function scatterplot(D, v1, v2, varargin)
% Correlation of two variables in dataset
%   - can be either event-level or subject-level
%   - two types of correlations in event-level data are returned: across
%   raw datapoints (both within- and  between-subject) and on centered data
%   (within-subject only)
%   - both variables must be valid names (case-sensitive)
%
% :Usage:
% ::
%
%    fig_han = correlation(D, varname1, varname2, [optional inputs])
%
% ..
%     Author and copyright information:
%
%     Copyright (C) 2016 Tor Wager
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
%
% :Inputs:
%
%   **D:**
%        a canlab_dataset object
%
%   **v1:**
%        x variable name
%
%   **v2:**
%        y variable name
%
%
% :Optional Inputs:
%
%   **nofig:**
%        suppress creation of new figure
%
%   **subjtype:**
%        group by the following variable name
%
%   **wh_keep:**
%        followed by logical
%
%   **colors:**
%        followed by colors.
%
%   **dorobust:**
%        do robust corr.  if enabled, colors will not work and subjtype grouping will not work well until
%        the function plot_correlation_samefig is updated, at some point in the future.
%
%
% :Outputs:
%
%   **fig_han:**
%        figure handle
%
% :Examples:
% ::
%
%    scatterplot(D, 'Anxiety', 'Frustration');
%    fig_han = scatterplot(D, D.Subj_Level.names{1}, D.Subj_Level.names{2});
%    scatterplot(D, D.Event_Level.names{1}, D.Event_Level.names{2});
%

% STANDARD CODE FOR CANLAB_DATASET METHODS
% -------------------------------------------------------------------------
% fig_han = [];
% dofig = 1;
%grouping_var_name=[];
wh_keep = true(size(D.Subj_Level.id)); %everyone, by default
colors{1}='k';
%dorobust=0;

for i=1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            %             case 'subjtype'
            %                 grouping_var_name = varargin{i+1};
            case 'wh_keep'
                wh_keep = varargin{i+1};
                %             case 'nofig'
                %                 dofig=0;
                %             case {'robust', 'dorobust'}
                %                 dorobust=1;
        end
    end
end

[dat1, dcell1, whlevel1] = get_var(D, v1, wh_keep, varargin{:});
[dat2, dcell2, whlevel2] = get_var(D, v2, wh_keep, varargin{:});
dat1_level{1}=dat1; %to support grouping
dat2_level{1}=dat2;

if whlevel1 ~= whlevel2
    disp('No plot: Variables are not at same level of analysis.');
    return
end

for i = 1:length(dat1_level)    
    
    switch whlevel1
    case 1  
        x=dat1_level{i}; y= dat2_level{i};
            
        shortstr = correlation_subfunction(whlevel1, x, y, v1, v2);


    case 2
        
        shortstr = correlation_subfunction(whlevel1, dcell1, dcell2, v1, v2);

    otherwise
        error('Illegal level variable returned by get_var(D)');
    end

end

end % main function



function shortstr = correlation_subfunction(whlevel1, dat1, dat2, v1, v2)

shortstr = [];

switch whlevel1
    case 1
        
    [wasnan, dat1, dat2] = nanremove(dat1, dat2);

    [rtotal, ptotal] = corr(dat1, dat2); 

        shortstr = sprintf('r = %3.2f\n', rtotal);
        str = sprintf('(%s, %s): r = %3.2f, p = %3.6f\n', v1, v2, rtotal, ptotal);
        disp(str)
        
    case 2
        Xc = cat(1, dat1{:});
        Yc = cat(1, dat2{:});
        rtotal = corr(Xc, Yc);

        for i = 1:length(dat1)
            x1{i} = scale(dat1{i}, 1); % mean-center
            x2{i} = scale(dat2{i}, 1);
        end
        
        rwithin = corr(cat(1, x1{:}), cat(1, x2{:}));
        
        % sqrt var increase when not removing between-subject random intercepts
        rbetween = sqrt(rwithin .^ 2 - rtotal .^ 2);
        
        % get P-values
        % within var: T-test on slopes for each subject
        for i = 1:length(dat1)
            b(i,:) = glmfit(dat1{i}, dat2{i});      % total
            b2(i,:) = glmfit(x1{i}, x2{i});         % within
        end
        
        [~, ptotal] = ttest(b(:, 2));
        [~, pwithin] = ttest(b2(:, 2));

        shortstr = sprintf('r = %3.2f', rtotal);
        
        str = sprintf('(%s, %s): \nr across all data: %3.2f, p = %3.6f\nr within subjects: %3.2f, p = %3.6f\nr between subjects: %3.2f\n', v1, v2, rtotal, ptotal, rwithin, pwithin, rbetween);
        disp(str)
        
    otherwise ('Illegal value!');
        
end

end
