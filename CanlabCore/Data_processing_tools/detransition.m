function y = detransition(y,varargin)
% For fMRI timeseries that contains large 'jump' artifacts due to motion correction
% or other problems.
%
% :Usage:
% ::
%
%     y = detransition(y,[doplot])
%
% Removes these large spikes.
% Updated version built into spikecorrect in trimts

doplot = 1;
if length(varargin) > 0, doplot = varargin{1};,end

if doplot, figure;subplot(1,2,1);plot(y,'k'); hold on; cols = {'b' 'g' 'r' 'c' 'm'}; subplot(1,2,2),yg = scale(diff(y));,[h,x]=hist(yg,30);,plot(x,h,'k'),drawnow,end

thr = [3 3.3 3.5];

for iter = 1:3
    
    yg = scale(diff(y));
    
    % robust standard deviation on Windsorized data
    sd  = std(trimts(yg,2,[]));
    yg = (yg-mean(yg)) ./ sd;
    
    yw = find(abs(yg) > thr(iter));

    for i = 1:length(yw)
    
        y(yw(i)+1:end) = y(yw(i)+1:end) - (y(yw(i)+1) - y(yw(i)));
    
    end
    
    if doplot, subplot(1,2,1); plot(y,cols{iter}), 
        subplot(1,2,2); hold on; [h]=hist(yg,x);hh=bar(x,h);set(hh,'FaceColor',cols{iter})
        drawnow, ;pause(1), 
    end
end

    y = detrend(y,'linear');
    
    return
    
