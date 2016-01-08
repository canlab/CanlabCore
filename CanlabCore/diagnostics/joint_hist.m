function [z, xbins, ybins] = joint_hist(x,y,varargin)
% Create 2-D joint histogram from vectors x and y
%
% :Usage:
% ::
%
%     [z, xbins, ybins] = joint_hist(x,y,[nbins],[noplot])
%
% :Inputs:
%
%     x and y:**
%        are vectors of paired observations on two variables
%
% :Outputs:
%
%   **z:**
%        is the matrix representing the joint histogram
%        cols of z are bins of x, rows are bins of y
%        in plot, X axis is y, Y axis is x
%
% Optional: number of bins, suppress plotting
% 
% :Examples:
% ::
%
%    z = joint_hist(nnmfscores{i}{j}(:, 1),nnmfscores{i}{j}(:, 2), 50, 'noplot');
%    h = plot_joint_hist_contour(z, [0 0 1]);
%
% ..
%    Tor Wager
% ..

nbins = 10;
if length(varargin) > 0, nbins = varargin{1}; end

[h, xbins] = hist(x, nbins);
[h, ybins] = hist(y, nbins);

z = zeros(length(xbins),length(ybins));

for i = 1:length(x)
    
    xx = find(abs(xbins - x(i)) == min(abs(xbins - x(i))));
    yy = find(abs(ybins - y(i)) == min(abs(ybins - y(i))));
    
    z(xx,yy) = z(xx,yy) + (1 ./ length(xx)*length(yy)); 
    % dividing breaks ties by splitting the value assigned evenly among eligible bins

end 

% Transpose to make compatible with x/y for contour, etc.
z = z';

if ~(length(varargin) > 1)
    bar3(z);
    xlabel(inputname(2));
    ylabel(inputname(1));
    set(gca,'YTick',1:length(xbins))
    set(gca,'XTick',1:length(ybins))
    set(gca,'YTickLabel',xbins)
    set(gca,'XTickLabel',ybins)
end
    


return
