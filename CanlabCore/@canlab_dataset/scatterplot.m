function fig_han = scatterplot(D, v1, v2, varargin)
%
% fig_han = scatterplot(D, varname1, varname2, varargin)
%
% Scatterplot of two variables in dataset
% - can be either event-level or subject-level
% - event-level data is plotted as multi-line plot, one line per subject
% - both variables must be valid names (case-sensitive)
%
% Optional inputs:
%  - 'nofig': suppress creation of new figure
%  - 'subjtype': group by the following variable name
%  - 'wh_keep' : followed by logical
%  - 'colors' : followed by colors. 
%  - 'dorobust': do robust corr.  if enabled, colors will not work and subjtype grouping will not work well until
%  the function plot_correlation_samefig is updated, at some point in the future.
%  
%
% Example:
%
% scatterplot(D, 'Anxiety', 'Frustration');
% fig_han = scatterplot(D, D.Subj_Level.names{1}, D.Subj_Level.names{2});
% scatterplot(D, D.Event_Level.names{1}, D.Event_Level.names{2});
%
% Copyright Tor Wager, 2013

fig_han = [];
dofig = 1;
grouping_var_name=[];
wh_keep = true(size(D.Subj_Level.id)); %everyone
colors{1}='k';
dorobust=0;

for i=1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}            
            case 'subjtype'
                grouping_var_name = varargin{i+1};
            case 'wh_keep'
                wh_keep = varargin{i+1};
            case 'nofig'
                dofig=0;
            case {'robust', 'dorobust'}
                dorobust=1;
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

if isempty(dat1) || isempty(dat2)
    % skip
    disp('No plot: Missing variables');
    return
end

if dofig
    fig_han = create_figure([v1 '_vs_' v2]);
else
    fig_han = gcf;
end


if ~isempty(grouping_var_name)  % We have a grouping variable
    
    %get a wh_keep for all the levels
    [grouping_var, dum, dum, descripGrp] = D.get_var(grouping_var_name, wh_keep);
    
    levels = unique(grouping_var);
    
    % colors
    colors = {'r' 'b' 'g' 'k' 'y'};
    colors = colors(1:length(levels));
    wh = strcmp(varargin, 'colors');
    if any(wh), colors = varargin{find(wh)+1}; end

        
    for i=1:length(levels)
        wh_keep_lev{i} = (D.get_var(grouping_var_name,wh_keep)==levels(i));
        
        dat1_level{i} = dat1(wh_keep_lev{i},:);
        dat2_level{i} = dat2(wh_keep_lev{i},:);     
    end
end  
    
  
for i=1:length(dat1_level)    
    
    switch whlevel1
    case 1  
        x=dat1_level{i}; y= dat2_level{i};
            
        if dorobust
            plot_correlation_samefig(x,y,[],[],[],1)
            grid off
        else
            scatter(x,y,65,  'MarkerFaceColor', colors{i}, 'MarkerEdgeColor', colors{i});%, 'within')
            inds = isnan(x) | isnan(y);
            h=refline(polyfit(x(~inds),y(~inds),1))
            set(h, 'Color', colors{i}, 'LineWidth', 2)
        end
        
    case 2
        
        han = line_plot_multisubject(dcell1, dcell2, varargin{:});
        
    otherwise
        error('Illegal level variable returned by get_var(D)');
end

set(gca, 'FontSize', 24)

xlabel(strrep(v1, '_', ' '));
ylabel(strrep(v2, '_', ' '));

rtotal = corr(dat1(:), dat2(:));

switch whlevel1
    case 1
        str = sprintf('r = %3.2f\n', rtotal);
        disp(str)
        
    case 2
        for i = 1:length(dcell1)
            x1{i} = scale(dcell1{i}, 1); % mean-center
            x2{i} = scale(dcell2{i}, 1);
        end
        rwithin = corr(cat(1, x1{:}), cat(1, x2{:}));
        
        str = sprintf('r across all data: %3.2f\nr within subjects: %3.2f', rtotal, rwithin);
        disp(str)
        
    otherwise ('Illegal value!');
        
end


xloc = mean(dat1(:)) + std(dat1(:));
yloc = mean(dat2(:)) + std(dat2(:));

text(xloc, yloc, str, 'FontSize', 24);


end % function

