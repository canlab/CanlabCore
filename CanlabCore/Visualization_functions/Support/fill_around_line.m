function hh = fill_around_line(dat,err,color,varargin)
% h =fill_around_line(dat,err,color,[x indices])
% fills a region around a vector (dat) with + or - err
%
% used in tor_fill_steplot.m

if length(err) == 1
    err = repmat(err,length(dat),1);
end
    
hh = [];

xn = 1:length(err);
if length(varargin) > 0, xn = varargin{1}; end

for i = 1:length(err)-1
    hold on;
    yplus1 = dat(i)+err(i);
    yplus2 = dat(i+1)+err(i+1);
    yminus1= dat(i)-err(i);
    yminus2 = dat(i+1)-err(i+1);
    
%     if ischar(color)

        hh(end+1) = fill([xn(i) xn(i) xn(i+1) xn(i+1)],[yminus1 yplus1 yplus2 yminus2],color,'Edgecolor','none','FaceAlpha',.3);
% 
%     else
%         hh(end+1) = fill([xn(i) xn(i) xn(i+1) xn(i+1)],[yminus1 yplus1 yplus2 yminus2],'FaceColor', color,'Edgecolor','none','FaceAlpha',.3);
%     end

end

return

