function obj = legend(obj, varargin)
% Creates legend for fmridisplay object
% Adds legend axis handles to obj.activation_maps{:}
%
% :Usage:
% ::
%
%     obj = legend(obj, varargin)
%
%     obj = legend(obj, 'figure') % new figure
%
% ..
%    Tor Wager
%    8/17/2016 - pkragel updated to accomodate split colormap
% ..
%
% Notes: scaleanchors is min and max values for pos and neg range (for
% splitmap)
%

donewfig = 0;

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            case {'figure' 'newfig'}, donewfig = 1;
                
            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end

if donewfig
    create_figure('legend'); axis off
    
    mypositions = {[0.100    0.1200    0.80    0.04] ...
        [0.100    0.2800    0.80    0.04] ...
        [0.100    0.46    0.80    0.04] ...
        [0.100    0.64    0.80    0.04] ...
        };
    
    myfontsize = 18;
    
else
    mypositions = {[0.05    0.1200    0.20    0.04] ...
        [0.350    0.1200    0.20    0.04] ...
        [0.650    0.1200    0.20    0.04] ...
        [0.05    0.2800    0.20    0.04] ...
        };
    myfontsize = 14;
    
end

for c = 1:length(obj.activation_maps)
    
    currentmap = obj.activation_maps{c};
    
    scaleanchors = currentmap.cmaprange; %for text labels
    
    if isempty(scaleanchors), continue, end
    
    % adjuts in case there are no values on one end
    scaleanchors(isnan(scaleanchors)) = 0;
    
    if ~diff(scaleanchors)
        disp('No variability in mapped values. Not plotting legend.');
        continue
    end
    
    if c > length(mypositions)
        disp('Maximum number of legends exceeded. Not plotting remaining legends.');
        break
    end
    
    obj.activation_maps{c}.legendhandle = axes('Position', mypositions{c});
    hold on;
    
    % mapd = currentmap.mapdata(:); mapd = mapd(mapd ~= 0 & ~isnan(mapd));
    % scaleanchors = [min(mapd) max(mapd)];
    
    % This info should already be in currentmap
    % cmaprange = [prctile(mapd, 20) prctile(mapd, 80)];
    % currentmap.cmaprange = cmaprange;
    % currentmap.maxcolor = [1 1 0];
    % currentmap.mincolor = [1 .3 0];
    
    % fix for multiple blobs
    if size(currentmap.mincolor, 2) > 3 || size(currentmap.maxcolor, 2) > 3
        fprintf('Warning! Extra colors in map...legend will not display correctly with multiple blobs');
        
        currentmap.mincolor = currentmap.mincolor(:, 1:3);
        currentmap.maxcolor = currentmap.maxcolor(:, 1:3);
    end
    
    nsteps = 100;
    
    wvals = linspace(0, 1, nsteps);
    
    % separate plot for split colormap
    issplitmap = size(currentmap.mincolor, 1) - 1; % zero for one row, 1 for 2+
    
    if issplitmap
        
%         legvals = linspace(scaleanchors(1), scaleanchors(4), nsteps);
%         legvals(legvals==0)=sign(mean([scaleanchors(1), scaleanchors(4)])*1E-6);
         
        [legvals, rgb] = get_split_colormap_values_and_colors(nsteps, scaleanchors, currentmap);
        
    else
        legvals = linspace(scaleanchors(1), scaleanchors(2), nsteps);
        
        
    end
    
    
    
    for i = 2:nsteps
        
        if ~issplitmap % || (issplitmap && legvals(i) > 0)
            
            fcolor = (1 - wvals(i)) * currentmap.mincolor(1, :) + wvals(i) * currentmap.maxcolor(1, :);
            
        elseif issplitmap % && legvals(i) < 0
            
            %fcolor = (1 - wvals(i)) * currentmap.mincolor(2, :) + wvals(i) * currentmap.maxcolor(2, :);
            
            fcolor = rgb(i, :);
            
        else
            error('This should not happen...debug me.')
        end
        
        fcolor(fcolor > 1) = 1;
        fcolor(fcolor < 0) = 0;
        
        try
            fill([legvals(i-1) legvals(i-1) legvals(i) legvals(i)], [0 1 1 0], fcolor, 'EdgeColor', 'none');
            
        catch
            disp('problem with legend fill.')
            keyboard
        end
        
    end
    
    if ~issplitmap
        
        ntick=5;
        myxtick = linspace(scaleanchors(1), scaleanchors(2), ntick);
        
    else
        
        ntick=4;
        myxtick = scaleanchors;  % linspace(scaleanchors(1), scaleanchors(4), ntick);
        
    end
    
    % kludgy fix for sig digits
    for i = 1:ntick
        mylabels{i} = sprintf('%3.2f', myxtick(i));
    end
    
    if sum(strcmp(mylabels, mylabels{1})) > 1
        for i = 1:ntick
            mylabels{i} = sprintf('%3.3f', myxtick(i));
        end
    end
    
    % Relabel if we need more significant digits; we have duplicates
    if sum(strcmp(mylabels, mylabels{1})) > 1
        for i = 1:ntick
            mylabels{i} = sprintf('%3.4f', myxtick(i));
        end
    end
    
    % remove duplicates and make sure sorted ascending. needed if missing
    % pos or neg values
    [myxtick, wh] = unique(myxtick);
    mylabels = mylabels(wh);
    
    set(gca, 'YTickLabel', '', 'YColor', 'w', 'FontSize', myfontsize);
    axis tight
    set(gca, 'XLim', [min(scaleanchors) max(scaleanchors)], 'XTick', myxtick, 'XTickLabel', mylabels);
     
end


end % function






function [legvals, rgb] = get_split_colormap_values_and_colors(nsteps, scaleanchors, currentmap)
% Values corresponding to colors
legvals = zeros(nsteps + 3, 1);
legvals = linspace(scaleanchors(1), scaleanchors(2), nsteps./2)';
legvals = [legvals; 0; 0; scaleanchors(3); linspace(scaleanchors(3), scaleanchors(4), nsteps./2)'];

% Colors
rgb = zeros(nsteps + 3, 3);

for i = 1:3
    rgb1(:, i) = linspace(currentmap.mincolor(2, i), currentmap.maxcolor(2, i), nsteps./2)';
end

rgb1 = [rgb1; ones(3, 3)]; % pad with 1 1 1 white

for i = 1:3
    rgb2(:, i) = linspace(currentmap.mincolor(1, i), currentmap.maxcolor(1, i), nsteps./2)';
end

rgb = [rgb1; rgb2];

end
