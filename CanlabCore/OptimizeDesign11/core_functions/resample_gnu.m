% Copyright (C) 2000 Paul Kienzle
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

% usage: y=resample(x,p,q,d)
%
% Change the sample rate of x by a factor of p/q.  Note that p and q do
% not need to be integers since this routine does not use a polyphase
% rate change algorithm, but instead uses bandlimited interpolation,
% wherein the continous time signal is estimated by summing the sinc
% functions of the nearest neighbouring points up to distance d.
%
% This is discussed in:
%     J. O. Smith and P. Gossett (1984). A flexible sampling-rate
%     conversion method. In ICASSP-84, Volume II, pp. 19.4.1-19.4.2. 
%     New York: IEEE Press.
% See the authors page at: http://www-ccrma.stanford.edu/~jos/resample/
%
% Note that the resampling is not yet very fast or very good, but it is
% very flexible.
%
% Example
%    % Speech example
%    [x, fs] = auload(file_in_loadpath("sample.wav"));
%    sound(resample(x,16000,fs), 16000);  % resample at 16 kHz
%
%    % Example from interp1
%    xf=0:0.05:10.95; yf = sin(2*pi*xf/5);
%    xp=0:10;         yp = sin(2*pi*xp/5);
%    r = resample(yp,xp(2),xf(2));
%    plot(xf,yf,';original;',xf,r,';resample;',xp,yp,'*;;');
%
% Note that resample computes all samples up to but not including time
% n+1. If you are increasing the sample rate, this means that it will
% generate samples beyond the end of the time range of the original
% signal. That is why xf must goes all the way to 10.95 in the example.
 
% TODO: Fix so that audible clicking goes away.
% TODO: Change to a faster algorithm.
% TODO: Test on a chirp signal.
   
function y=resample(x,p,q,order,beta)
  if (nargin < 2 | nargin > 5)
    disp('USAGE: y=resample(x,p,q,order)');
  end

  if (nargin < 3), q=1; end
  if (nargin < 4), order = 5; end

  %  % chain to decimate/interpolate if appropriate
  %  if p==1 && q==fix(q)
  %    y=decimate(x,q); order?
  %    return;
  %  elseif q==1 && p==fix(p)
  %    y=interp(x,q); order?
  %    return;
  %  end

  transpose = size(x,1)==1;
  if transpose, x = x.'; end

  % if rate reduction, apply antialiasing filter first
  r=p/q;
  if (r < 1)                 
    b = fir1(2*order+1, r);
    x = fftfilt(b, x);
  end

  % Determine the new sampling times, and their distance to the old
  % ones.  Note that the new series should be the maximum that can
  % be contained in the old series without going over the time
  % allotted to the old series.  In short, you have to go a little
  % beyond the last sample of the old series if your new sampling
  % rate is higher.
  t=[1:1/r:length(x)+1-1/r]';   % the sampling points of the new series
  idx = fix(t);                 % the nearest old point
  t = t-idx;                    % distance to the nearest old point

  % generate the new series by summing the sinc functions of the
  % nearest neighbour points implicit in the continuous time
  % expansion of the old series.  This new series is truncated
  % to +/- order nearest neighbours.  For convenience, the original
  % series is zero-padded before and after, implicitly setting the
  % neighbours at the start of the signal to zero.
  x = [zeros(order,size(x,2)) ; x ; zeros(order,size(x,2))];
  y = zeros(length(idx),size(x,2));        % the new series
  for i=-order:order
    w = sinc(t-i).*(0.5+0.5*cos(pi*(t-i)/(order+0.5))); % hanning window
    y=y + x(idx+i+order,:).*w(:,ones(size(x,2),1));
  end

  if transpose, y=y.'; end

return


%!demo
%! xf=0:0.05:10.95; yf = sin(2*pi*xf/5);
%! xp=0:10;      yp = sin(2*pi*xp/5);
%! r = resample(yp,xp(2),xf(2));
%! oneplot();
%! title("confirm that the resampled function matches the original");
%! plot(xf,yf,';original;',...
%!    xf,r(1:length(xf)),';resample;',...
%!    xp,yp,'*;;');
%! title("");
%! [x, fs] = auload(file_in_loadpath("sample.wav"));
%! sound(resample(x,16000,fs), 16000);
