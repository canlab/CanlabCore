function h = wani_pie(X, varargin)
% Draw a little better pie chart
%
% :Usage:
% ::
%
%    h = wani_pie(X, varargin)
%
% :Inputs:
%
%   **X:**
%        a vector
%
% :Optional Inputs: Enter keyword followed by variable with values
%
%   **'cols' or 'colors':**
%        colors N x 3 {default: using colors from microsoft office}
%        OR cell array of 3-element colors (no text input for colors!)
%
%
%   **'notext':**
%        no text for percentage {default: false}
%
%   **'fontsize':**
%        font size for percentage {default: 15}
%
%   **'hole':**
%        add a hole in the middle of the pie chart {default: no hole}
%
%   **'hole_size':**
%        specify the size of the middle hole {default: 5000}
%
%   **'outline'**
%
%   **'outlinecol'**
%
%   **'outlinewidth'**
%
%   **'outlinewidth'**
%
%   **'labels':**
%       followed by a cell array with string labels for each slice. Can be
%       empty.
%
% :Output:
%
%   **h:**
%        graphic handles
%
% :Examples:
% ::
%
%    % data
%    X = rand(10,1);
%    h = wani_pie(X, 'notext', 'hole')
%
%    savename = 'example_pie.pdf';
%
%    try
%        pagesetup(gcf);
%        saveas(gcf, savename);
%    catch
%        pagesetup(gcf);
%        saveas(gcf, savename);
%    end
%
% ..
%    Copyright (C) 2015 Wani Woo
% ..


% Programmers notes:
% 8/21/2017 - Stephan
% Added option to specify labels for each slice. Varargin 'labels' followed
% by a cell array with labels for each slice. See matlab's pie.m help for
% more info on labels.


cols = [0.0902 0.2157 0.3686 
    0.2157 0.3765 0.5725
    0.5843 0.2157 0.2078
    0.4667 0.5765 0.2353
    0.3765 0.2902 0.4824
    0.1922 0.5216 0.6118
    0.8941 0.4235 0.0392]; % Defaults

n = numel(X);

while size(cols, 1) < n
    newcol = cols(randperm(7),:);
    cols = [cols; newcol(1:3)];
end

dotext = 1;
fs = 15;
hs = 5000;
dohole = 0;
doout = 0;
outlinecol = [0 0 0];
outlinewidth = 1.2;
dolabels = 0;
labels = {};

% Optional inputs
% -------------------------------------------------------------------------

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            % functional commands
            case {'cols', 'colors'}
                cols = varargin{i+1};
            case {'notext'}
                dotext = 0;
            case {'fontsize'}
                fs = varargin{i+1};
            case {'hole'}
                dohole = 1;
            case {'hole_size'}
                hs = varargin{i+1};
            case {'outline'}
                doout = 1;
            case {'outlinecol'}
                outlinecol = varargin{i+1};
            case {'outlinewidth'}
                outlinewidth = varargin{i+1};
            case {'labels'}
                dolabels = 1;
                labels = varargin{i+1};
        end
    end
end

% note: will break with text color inputs, so check...***

if iscell(cols)
    cols = cat(1, cols{:});
end

% avoid exactly-zero values, which will give error
X(X <= 0) = abs(min(X(X ~= 0) ./ 1000)); 

% Create pie
% -------------------------------------------------------------------------

if dolabels
    h = pie(X,zeros(size(X)),labels);
else
    h = pie(X);
end
%set(gcf, 'color', 'w', 'position', [360 393 389 305]);

for i = 1:n
    set(h(2*i-1), 'facecolor', cols(i,:), 'edgecolor', 'none');
    hold on;
    if dotext
        set(h(2*i), 'fontSize', fs);
    else
        set(h(2*i), 'String', '');
    end
end
if dohole
    scatter(0, 0, hs, 'w', 'filled');
end

if doout
    scatter(0, 0, 43500, outlinecol, 'linewidth', outlinewidth);
end

end % function
