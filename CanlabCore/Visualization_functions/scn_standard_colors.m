function colors = scn_standard_colors(varargin)
% Create a set of unique colors in a standardized order.
%
% :Usage:
% ::
%
%    colors = scn_standard_colors(100)
%
% Optional input: minimum number of colors to generate
%
% Repeats after 36 colors

        % unique colors for each blob
        colors = {[1 0 0] [0 1 0] [0 0 1] [1 1 0] [1 0 1] [0 1 1]};
        
        colors = [colors {[1 .5 0] [.5 1 0] [.5 0 1] [1 0 .5] [0 1 .5] [0 .5 1]}];
        
        % next 6 colors, make pastel
        for i = 1:length(colors), colors = [colors {colors{i} .* .5}]; end
        
        % next 6, add blue
        for i = 7:12, colors = [colors colors(i)];  colors{end}(3) = colors{end}(3) + .5; end
        
        % next 6, add red
        for i = 7:12, colors = [colors colors(i)];  colors{end}(1) = colors{end}(1) + .5; end

        if ~isempty(varargin)
            n = varargin{1};
            
            while length(colors) < n, colors = [colors colors]; end
            
            colors = colors(1:n);
        end
        
        for i = 1:length(colors)
            colors{i}(colors{i} < 0) = 0;
            colors{i}(colors{i} > 1) = 1;
        end
        
end
