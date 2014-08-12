% hh = plot_onsets(d1, [color], [miny], [scalef], [duration])
%
% plots non-zero elements of an indicator vector or a list of onset times as vertical lines
%
% d1 is indicator vector; uses entries as scaling for height
%
% optional inputs:
% color, e.g., 'r'
% miny, shifts bars up/down to miny
% scalef, scales height of bars
%
% e.g., 
% plot_onsets(d1, 'k', -.5, .4)
%
%

function hh = plot_onsets(d1, varargin)
    if ~any(d1(2:end) == 0)    % we have an onset list, get indicator
        n = round(max(d1));
        z = zeros(n, 1);
        z(round(d1)+1) = 1;   % 0 is start of 1st onset, 1 of 1st element
        d1 = z;             % convert to indicator
        %wh = d1;
    end



    % make a row vector
    if size(d1, 1) > size(d1, 2)
        d1 = d1';
    end

    % indicator
    wh = find(d1);  


    color = 'k';
    miny = 0;
    scalef = 1;
    dur = 0;

    if length(varargin) > 0, color = varargin{1}; end
    if length(varargin) > 1, miny = varargin{2}; end
    if length(varargin) > 2, scalef = varargin{3}; end
    if length(varargin) > 3, dur = varargin{4}; end

    % adjust by -1 so first element is time 0
    if dur == 0
        % events
        hh = plot([wh; wh]-1, [zeros(size(wh))+miny; miny+scalef.*d1(wh)], color);
    else
        % epochs
        hh = [];
        for i = 1:length(wh)
            hh(end+1) = drawbox(wh(i)-1, dur, miny, scalef.*d1(wh(i)), color);
        end
    end
end


