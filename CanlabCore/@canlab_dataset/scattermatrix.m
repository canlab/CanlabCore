function fig_han = scattermatrix(D, wh_level, wh_vars)
% Scatterplot matrix of pairwise event-level variables
%
% fig_han = scattermatrix(D, wh_level, wh_vars)
%
% wh_level: 1 (Subject) or 2 (Event)
%
% Examples:
% fig_han = scattermatrix(D);
%
% wh = [5:9];
% fig_han = scattermatrix(D, 2, wh);
%
% f = scattermatrix(D, 2, {'Choice' 'RT' 'Pain' 'SwitchNext' 'Frustration' 'Anxiety' 'Control'});
%
% Copyright Tor Wager, 2013

switch wh_level
    case 1
        names = D.Subj_Level.names;
        
    case 2
        names = D.Event_Level.names;
    otherwise
        error('Enter 1 or 2 for wh_level.');
end

if nargin < 3 || isempty(wh_vars)
    wh_vars = 1:length(names);
end

% convert from cell of names of needed
if iscell(names)
    wh = zeros(size(names));
    for i = 1:length(wh_vars)
        wh = wh + strcmp(wh_vars{i}, names);
    end
    wh_vars = find(wh);
end

names = names(wh_vars);

n = length(names);

fig_han = create_figure('scattermatrix', n, n);

for i = 1:n
    
    subplot(n, n, (i-1)*n + i)
    histogram(D, names{i}, 'nofig');
    
    for j = i+1:n
        
        subplot(n, n, (i-1)*n + j)
        
        scatterplot(D, names{i}, names{j}, 'nofig');
        
        title(' ');
        drawnow
        
    end
    
end

end % function



