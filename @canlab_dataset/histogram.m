function fig_han = scatterplot(D, v1, varargin)
%
% fig_han = scatterplot(D, varname1, varargin)
%
% Histogram of one variable in dataset
% - can be either event-level or subject-level
% - event-level data is plotted as concatenated events across subject-level
% - both variables must be valid names (case-sensitive)
%
% Optional inputs:
%  - 'nofig': suppress creation of new figure
%
% Example:
%
% histogram(D, 'Anxiety');
%
% Copyright Tor Wager, 2013

fig_han = [];
dofig = 1;

if any(strcmp(varargin, 'nofig'))
    dofig = 0;
end

[dat1, dcell1, whlevel1] = get_var(D, v1, varargin{:});


if isempty(dat1) 
    % skip
    disp('No plot: Missing variables');
    return
end

if dofig
    fig_han = create_figure([v1 ' Histogram']);
else
    fig_han = gcf;
end

% Set number of bins
nbins = max(10, length(dat1(:)) ./ 100);
nbins = min(nbins, 200);
nbins = min(nbins, length(unique(dat1(:))));

switch whlevel1
    case 1
        hist(dat1(:), nbins);
        
    case 2
        
        hist(dat1(:), nbins);
        
    otherwise
        error('Illegal level variable returned by get_var(D)');
end

        grid off
han = gca;
set(gca, 'FontSize', 24)

xlabel(v1);


end % function

