function newcm = colormap_tor(lowcolor, hicolor, varargin)
% Create a new colormap of your choosing.
%
% :Usage:
% ::
%
%    newcolormap = colormap_tor(lowcolor, hicolor, [midcolor], [midcolor2], etc.)
%
% :Inputs:
%
%   **lowcolor:**
%        3-element color vector [r g b] 
%
%   **hicolor:**
%        3-element color vector [r g b] 
%
% :Optional Inputs:
%   **'n':**
%        Followed by number of desired colormap elements
%        Must be divisible by # colors entered - 1
%
%   **[r g b]:**
%        Additional [rgb]3-element color vectors, as many as you want
%        These define intermediate colorsﬂ
%
% :Outputs:
%
%   **cm:**
%        A colormap matrix, [n x 3]
%
% :Examples:
% ::
%
%    colormap_tor([.2 .2 .6], [1 1 0]);  % slate to yellow
%    colormap_tor([.9 .5 .2], [1 1 0]);  % orange to yellow
%    colormap_tor([.8 .1 .1], [1 1 0], [.9 .6 .1]);  %red to orange to yellow
%    colormap_tor([.2 .2 .4], [1 1 0], [.9 .6 .1]);  %slate to orange to yellow
%
%    cm = colormap_tor([.2 .2 .6], [1 1 0], 'n', 25);  % with 25 elements in color map
%
%    % Generate a split colormap with 32 grayscale and 32 colored values:
%    cm = colormap_tor([.2 .2 .6], [1 1 0], [.4 .4 .4], 'n', 32);
%    cm = [cmgray; cm];
%
%    % Set the current figure colormap to cm:
%    colormap(cm)
%
% ..
%    tor wager, sept. 2007
% ..

% -------------------------------------------------------------------------
% DEFAULT ARGUMENT VALUES
% -------------------------------------------------------------------------

n = 64;
newcm = NaN .* zeros(n, 3);  %initialize
additional_colors = {};

% -------------------------------------------------------------------------
% OPTIONAL INPUTS
% -------------------------------------------------------------------------

allowable_inputs = {'n'};

% optional inputs with default values - each keyword entered will create a variable of the same name

for i = 1:length(varargin)
    
    if ischar(varargin{i})
        
        switch varargin{i}
            
            case allowable_inputs
                
                eval([varargin{i} ' = varargin{i+1}; varargin{i+1} = [];']);
                
            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
        
    elseif length(varargin{i}) == 3
        
        % We have additional colors
        additional_colors{end+1} = varargin{i};
        
    end
    
end % Input parsing

% -------------------------------------------------------------------------
% MAIN FUNCTION
% -------------------------------------------------------------------------


if isempty(additional_colors)
    % no mid-color
    newcm = [linspace(lowcolor(1), hicolor(1), n)' linspace(lowcolor(2), hicolor(2), n)' linspace(lowcolor(3), hicolor(3), n)' ];
    
else
    % we have more colors
    
    allcolors = [{lowcolor} additional_colors {hicolor}];
    
    % Redefine n based on number of total colors
    % 3 colors turns into 2 segments, 4 to 3, etc.
    n = n ./ (length(allcolors) - 1);
    
    if (n - round(n))
        error('Length of colormap requested must be divisible by total number of colors entered - 1');
    end
    
    newcm = [];
    
    for i = 1:(length(allcolors)-1)
         
        newlocolor = allcolors{i};
        newhicolor = allcolors{i + 1};
        
        newcolorcm = [linspace(newlocolor(1), newhicolor(1), n)' linspace(newlocolor(2), newhicolor(2), n)' linspace(newlocolor(3), newhicolor(3), n)' ];
        newcm = [newcm; newcolorcm];
        
    end % colors loop 
    
end


%         for i = 1:length(varargin) - 1
%             newcm1 = [linspace(varargin{i}(1), varargin{i+1}(1), n/2)' linspace(varargin{i}(2), varargin{i+1}(2), n/2)' linspace(varargin{i}(3),varargin{i+1}(3), n/2)' ];
%             newcm = [newcm; newcm1];
%         end
%
%         newcm1 = [linspace(varargin{end}(1), hicolor(1), n/2)' linspace(varargin{end}(2), hicolor(2), n/2)' linspace(varargin{end}(3),hicolor(3), n/2)' ];
%         newcm = [newcm; newcm1];

% %
% %         % we have a mid-color
% %         newcm1 = [linspace(lowcolor(1), midcolor(1), n/2)' linspace(lowcolor(2), midcolor(2), n/2)' linspace(lowcolor(3), midcolor(3), n/2)' ];
% %         newcm2 = [linspace(midcolor(1), hicolor(1), n/2)' linspace(midcolor(2), hicolor(2), n/2)' linspace(midcolor(3), hicolor(3), n/2)' ];
% %
% %         newcm = [newcm1; newcm2];
end


% %
% %     % ------------------------------------------------------------
% %     % for heatmap option: define color maps - biscale hot/cool
% %     % -------------------------------------------------------------
% %
% %     % color map - hot
% %     % --------------------------------------------
% %     %h1 = (0:1/99:1)';
% %     h1 = linspace(.7,1,300)';
% %     %h2 = ones(size(h1));
% %     h2 = linspace(0,.8,200)';
% %     h2 = [h2; linspace(.8,1,100)'];
% %
% %     h3 = zeros(size(h1));
% %     h = [h1 h2 h3];
% %     %h = [h1 h3 h3; h2 h1 h3; h2 h2 h1];
% %     %h(1:75,:) = []; % take only red values to start
% %     % in new matlab: h = colormap(hot(300));
% %
% %     % color map - winter
% %     % --------------------------------------------
% %     %h1 = (0:1/249:1)';
% %     h1 = linspace(0,.5,250)';
% %     %h2 = (1:-1/(249*2):.2)';
% %     h2 = linspace(1,.7,250)';
% %     h3 = zeros(size(h1));
% %     hc = [h3 h1 h2];
