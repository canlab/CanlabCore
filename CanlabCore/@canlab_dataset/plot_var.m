function [meandat, stedat] = plot_var(D, varname, varargin)
% Plot the mean and standard error of a variable across events.
%
% :Usage:
% ::
%
%    [meandat, stedat] = plot_var(D, varname, [optional inputs])
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
%
% :Inputs:
%
%   **D:**
%        a canlab_dataset object
%
%   **varname:**
%        the name of a valid variable to get from dataset
%           -Looks for var name at either level, returns Event level if exists at both levels
%
% :Optional inputs:
%
%   **subjtype:**
%        followed by name of grouping variable
%           - must be categorical subject-level variable
%           - if entered, plot lines or bars based on these categories
%           - 'eventmeans' will plot bars; without, it will plot line
%             plots across events with standard error shading
%           - the grouping variable's description, if it exists, will
%             be split along commas, and those values will be used as
%             column lables
%
%   **eventmeans:**
%        calculate and plot subject means across event-level variables
%           - if entered, will plot bar plots of means by condition
%
%   **wh_keep:**
%        followed by 1/0 vector of subjects to keep.
%           - must be same length as subjects
%           - subjects with value 0 will be excluded
%
%   **color:**
%        followed by one color for all bars, or cell array with names of colors cell for each line/bar
%
%   **nofig:**
%        don't make a new figure
%
%   **other:**
%        other varargin are passed directly to barplot_columns.  So
%        for example, '95CI' will make 95% confidence interals, instead
%        of SE bars.
%
%
% :Outputs:
%
%   **meandat:**
%        mean values
%
%   **stedat:**
%        standard error values
%
% :Examples:
% ::
%
%    plot_var(D, 'Frustration')
%    plot_var(D, 'RT')
%    plot_var(D, 'RT', 'eventmeans');
%    plot_var(D, 'RT', 'subjtype', 'Placebo');
%    plot_var(D, 'RT', 'eventmeans', 'subjtype', 'Placebo');
%    plot_var(D, 'RT', 'eventmeans', 'subjtype', 'Placebo', 'color', {'r' 'b'});
%


grouping_var_name = ''; % plot variable as a function of event number
event_means = 0;
wh_keep = true(size(D.Subj_Level.id)); %everyone
% colors <- defined below
myfontsize = 22;
nofig=0;

for i=1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            
            case 'subjtype'
                grouping_var_name = varargin{i+1};
                
            case 'eventmeans'
                event_means = 1;
                
            case 'wh_keep'
                wh_keep = varargin{i+1};
            
            case 'nofig'
                nofig=1;
        end
    end
end

if ~nofig, create_figure(varname), end

set(gca, 'FontSize', myfontsize);

[dat, dcell, wh_level, descripVar] = get_var(D, varname, wh_keep);

if any(any(isnan(dat(:,:)))), warning('Some NaNs!'); end
meandat = nanmean(dat);
stedat = ste(dat);

if ~isempty(grouping_var_name)
    % We have a grouping variable
    
    %get a wh_keep for all the levels
    [grouping_var, dum, dum, descripGrp] = D.get_var(grouping_var_name, wh_keep);
    
    levels = unique(grouping_var);
    
    % colors
    colors = scn_standard_colors(length(levels));
    wh = strcmp(varargin, 'colors');
    if any(wh), colors = varargin{find(wh)+1}; end
    
    
    for i=1:length(levels)
        wh_keep_lev{i} = (D.get_var(grouping_var_name,wh_keep)==levels(i));
        
        %plot each level
        dat_level{i} = dat(wh_keep_lev{i},:);
        
        if event_means || wh_level==1
            % prep bars
            ev_means{i} = nanmean(dat_level{i}, wh_level);
            
        else
            % prep plot of events - standard error fills
            h = fill_around_line(nanmean(dat_level{i}), ste(dat_level{i}), colors{i});
            lineh(i) = plot(nanmean(dat_level{i}), 'o-', 'Color', colors{i}, 'LineWidth', 3);
        end
    end
    
    
    groupnames = regexp(descripGrp, ',', 'split'); 
    
    if event_means || wh_level==1
        % Bar plot of groups
        if wh_level==2
            barplot_columns(ev_means, prep_name(descripVar), [],'nofig', varargin{:});
        else
            barplot_columns(dat_level, prep_name(descripVar), [],'nofig', varargin{:});
        end

       
        set(gca, 'FontSize', myfontsize);
        ylabel(prep_name(descripVar));
        xlabel([]);
        
        set(gca, 'XTick', 1:length(levels), 'XTickLabel', groupnames);
        
    else
        set(gca, 'FontSize', myfontsize);

        legend(lineh, groupnames)

        xlabel('Event number');
        ylabel(prep_name(descripVar));
        
    end
    
else  % don't split by subj type
    
    if wh_level==1 %subject level variable
        barplot_columns(dat, prep_name(descripVar), [],'nofig', varargin{:});
        ylabel(descripVar); xlabel('all subjects'); set(gca, 'XTickLabel', []);

    else %event level variable

        % Average across subjects, with std. err across subjects
        h = fill_around_line(meandat, stedat, 'k');

        set(gca, 'FontSize', myfontsize);

        plot(meandat, 'k', 'LineWidth', 3);

        xlabel('Event number');

        ylabel(descripVar);
    end
end

set(gca, 'FontSize', myfontsize);


end % function


function varname = prep_name(varname)

if iscell(varname), varname = varname{1}; end

wh = varname == '_';

varname(find(wh)) = ' ';


end
