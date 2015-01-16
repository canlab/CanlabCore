function obj = legend(obj, varargin)
% Creates legend for fmridisplay object
% Adds legend axis handles to obj.activation_maps{:}
%
% obj = legend(obj, varargin)
% obj = legend(obj, 'figure') % new figure
% 
% Tor Wager

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
    mypositions = {[0.100    0.1200    0.20    0.04] ...
        [0.100    0.2800    0.20    0.04] ...
        [0.400    0.1200    0.20    0.04] ...
        [0.400    0.2800    0.20    0.04] ...
        };
    myfontsize = 14;
    
end

for c = 1:length(obj.activation_maps)
    
    currentmap = obj.activation_maps{c};
    
    scaleanchors = currentmap.cmaprange; %for text labels
    
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
    
    legvals = linspace(currentmap.cmaprange(1), currentmap.cmaprange(2), nsteps);
    wvals = linspace(0, 1, nsteps);
    
    % separate plot for split colormap
    issplitmap = size(currentmap.mincolor, 1) - 1; % zero for one row, 1 for 2+
    
    for i = 2:nsteps
        
        if ~issplitmap || (issplitmap && legvals(i) > 0)
            
            fcolor = (1 - wvals(i)) * currentmap.mincolor(1, :) + wvals(i) * currentmap.maxcolor(1, :);
        
        elseif issplitmap && legvals(i) < 0
            
            fcolor = (1 - wvals(i)) * currentmap.mincolor(2, :) + wvals(i) * currentmap.maxcolor(2, :);
        
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
    
    myxtick = linspace(scaleanchors(1), scaleanchors(2), 5);

    % kludgy fix for sig digits
    for i = 1:5
        mylabels{i} = sprintf('%3.2f', myxtick(i));
    end
    
    if sum(strcmp(mylabels, mylabels{1})) > 1
        for i = 1:5
            mylabels{i} = sprintf('%3.3f', myxtick(i));
        end
    end
    
    
    if sum(strcmp(mylabels, mylabels{1})) > 1
        for i = 1:5
            mylabels{i} = sprintf('%3.4f', myxtick(i));
        end
    end
    

    set(gca, 'YTickLabel', '', 'YColor', 'w', 'FontSize', myfontsize);
    set(gca, 'XLim', scaleanchors, 'XTick', myxtick, 'XTickLabel', mylabels);
    axis tight
    
end


end % function