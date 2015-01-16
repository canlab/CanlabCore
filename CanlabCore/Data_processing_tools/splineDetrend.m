function [fv,bp,yy,myfft] = splineDetrend(v,varargin)
% [fv,bp,yy,myfft] = splineDetrend(v,'p' [opt])
% input: v, a vector to be detrended
% output: fv, the detrended vector, bp, the knot points,
%           and myfft, the abs(fft) of the detrended vector.
%
% What it does:
% A spline detrend with knot points every 2 s (hard coded number)
% I made this up too - it's not the FDA approved method.
%
% By Tor Wager, 08/04/01
% tmp = cl(1).raw_data(:,1,6); tmp2 = trimts(tmp,3,[],1); [fv,bp,yy]=splineDetrend(tmp2);
% [fv,bp] = splineDetrend(tmp2,'p');
% 
% 'p' means plot
% any other opt argument sets the knots and is treated as an integer, with detrending every n
% elements.

% -------------------------------------------------
% * setup
% -------------------------------------------------

samprate = 1;   % default to 1 to use every knotrate images
knotrate = 10;   % every 100 images by default

m = mean(v);

plotme = 0;

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
        case 'p', plotme = 1;
        otherwise error('Unrecognized argument: acceptable is ''p'' for plot.')
        end
    else
        knotrate = (varargin{i});
        disp(['Setting knot rate to ' num2str(knotrate)])
    end
end

% -------------------------------------------------
% * detrend
% -------------------------------------------------

% knot point every 2 s at 200 Hz
bp = 1:knotrate.*samprate:length(v);

% figure out mean for each segment - 
% this is the interp knot point.
for i = 1:length(bp)-1
    data = v(bp(i):bp(i+1)-1);
    y(i) = mean(data);
    x(i) = mean(bp(i),bp(i+1));
end

% add endpoint
x = [x length(v)];
y = [y y(end)];

xx = 1:length(v);
yy = spline(x,y,xx);
fv = v - yy';

fv = fv + m;



% -------------------------------------------------
% * plot
% -------------------------------------------------
if plotme

figure; subplot(2,1,1)
title('Unfiltered with spline fit and knot points'); 
hold on; plot(v);plot(yy,'r','LineWidth',2)
ylim = get(gca,'YLim');
for i = 1:length(bp),plot([bp(i) bp(i)],ylim,'k','LineWidth',.5),end
subplot(2,1,2); plot(fv); title('Detrended')

% 200 is the sampling rate!
figure;hold on
x = (1:length(fv)) .* (samprate ./ length(fv));
myfft = abs(fft(fv));
plot(x,myfft)
set(gca,'XLim',[0 20])

ylim = get(gca,'YLim');

% find the max in the power spectrum below the nyquist limit
maxx = find(myfft(2:length(myfft./2)) == max(myfft(2:length(myfft./2))));
% add 1 to adjust for fact that you started at 2:length...
maxx = maxx + 1;
% find the frequency of this.
maxx = x(maxx);

plot([maxx(1) maxx(1)],ylim,'r')

end



return

