% function lineh = tor_bar_steplot(avg, ste, plotcolor, [x locations], [shift by], [multx by])
%
% for a GROUPED 2 x 2 error bar plot (or 2 x n):
% tor_bar_steplot(means, se, {'k'}, .35, .5, .2)
%
% tor_line_steplot(d(2, :), s(2, :), {'k'}, 1.6, .27)
% tor_line_steplot(d(1, :), s(1, :), {'k'}, .35, .5, .2)
%
% This function has been largely replaced with a better, more flexible one!
% See lineplot_columns (takes in data matrix, not avg/ste), so this
% function may still be useful sometimes.
% out = lineplot_columns(hotopen, 'color', 'r', 'marker', 'none', 'w', 1, 'shade');

function lineh = tor_line_steplot(avg, ste, plotcolor, varargin)
     wid = .1;
     
     lineh = [];
     
    if length(varargin) > 0
        xlocations = varargin{1}; 
        wid = min(diff(xlocations)) .* .3;
    else xlocations = 1:size(ste, 2); 
    end
    if length(varargin) > 1, shift = varargin{2}; else shift = 0; end
    if length(varargin) > 2, multx = varargin{3}; else multx = 1; end
    if length(varargin) > 3, shift2 = varargin{4}; else shift2 = 0; end

    hold on
    if isempty(plotcolor)
        plotcolor{1} = 'b';
    end
    %wid = get(gca, 'XLim');
    %wid = (wid(2) - wid(1))/50;
   


    %if isfield(Op, 'window'), xlocations = Op.window(1):Op.window(2);, end

    if length(xlocations) ~= size(avg, 2), 
        warning('Window is wrong length, using default'); %#ok
        xlocations = 1:size(avg, 2);
    end

    for i = 1:size(ste, 2)
        if length(plotcolor) < i
            plotcolor{i} = plotcolor{1};
        end

        xloc = xlocations(i);

        if avg(i) < 0
            multy = -1;
        else
            multy = 1;
        end

        try
            if mod(i, 2) == 0
                shft = shift - shift2;
            else
                shft = shift;
            end

            lineh = [lineh plot([xloc*multx+shft xloc*multx+shft], [avg(i) avg(i) + multy*ste(i)], plotcolor{i}(1))];
            lineh = [lineh plot([xloc*multx-wid+shft xloc*multx+wid+shft], [avg(i) + multy*ste(i) avg(i) + multy*ste(i)], plotcolor{i}(1))];
            
        catch
            disp('Can''t plot error bar.')
            i
            avg
            ste
            plotcolor
            xloc
            multx
            multy
            wid 
            shft
        end
    end
end
