function newcm = colormap_tor(lowcolor, hicolor, varargin)
    % newcolormap = colormap_tor(lowcolor, hicolor, [midcolor], [midcolor2], etc.)
    %
    % Create a new colormap of your choosing.
    %
    % colormap_tor([.2 .2 .6], [1 1 0]);  % slate to yellow
    % colormap_tor([.9 .5 .2], [1 1 0]);  % orange to yellow
    % colormap_tor([.8 .1 .1], [1 1 0], [.9 .6 .1]);  %red to orange to yellow
    % colormap_tor([.2 .2 .4], [1 1 0], [.9 .6 .1]);  %slate to orange to yellow
    %
    % tor wager, sept. 2007


    n = 64;
    newcm = NaN .* zeros(n, 3);  %initialize


    if nargin < 3
        % no mid-color
        newcm = [linspace(lowcolor(1), hicolor(1), n)' linspace(lowcolor(2), hicolor(2), n)' linspace(lowcolor(3), hicolor(3), n)' ];
    else
        % we have more colors
        newcm = [linspace(lowcolor(1), varargin{1}(1), n/2)' linspace(lowcolor(2), varargin{1}(2), n/2)' linspace(lowcolor(3),varargin{1}(3), n/2)' ];
        
        for i = 1:length(varargin) - 1
            newcm1 = [linspace(varargin{i}(1), varargin{i+1}(1), n/2)' linspace(varargin{i}(2), varargin{i+1}(2), n/2)' linspace(varargin{i}(3),varargin{i+1}(3), n/2)' ];
            newcm = [newcm; newcm1];
        end
        
        newcm1 = [linspace(varargin{end}(1), hicolor(1), n/2)' linspace(varargin{end}(2), hicolor(2), n/2)' linspace(varargin{end}(3),hicolor(3), n/2)' ];
        newcm = [newcm; newcm1];
        
% %         
% %         % we have a mid-color
% %         newcm1 = [linspace(lowcolor(1), midcolor(1), n/2)' linspace(lowcolor(2), midcolor(2), n/2)' linspace(lowcolor(3), midcolor(3), n/2)' ];
% %         newcm2 = [linspace(midcolor(1), hicolor(1), n/2)' linspace(midcolor(2), hicolor(2), n/2)' linspace(midcolor(3), hicolor(3), n/2)' ];
% % 
% %         newcm = [newcm1; newcm2];
    end

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
