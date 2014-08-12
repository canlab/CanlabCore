function equalize_axes(hvec,varargin)
%equalize_axes(hvec,[y only!])
% hvec is vector of axis handles
% tor wager

if length(varargin) > 0, yonly = 1; else yonly = 0; end

hvec(~strcmp(get(hvec, 'Type'), 'axes')) = [];

if isempty(hvec), disp('No valid axes'); return, end

for i = 1:length(hvec)
    
    xm = get(hvec(i),'XLim');
    xmax(i) = xm(2);
    xmin(i) = xm(1);
    
    ym = get(hvec(i),'YLim');
    ymax(i) = ym(2);    
    ymin(i) = ym(1);
    
end

xmax = max(xmax);
xmin = min(xmin);
ymax = max(ymax);
ymin = min(ymin);


for i = 1:length(hvec)
    
    axes(hvec(i))
    if yonly
        set(gca,'YLim',[ymin ymax])
    else
        axis([xmin xmax ymin ymax])
    end
end

end