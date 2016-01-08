function M = spm_general_hist(hP,mP,textlab,varargin)
% :Usage:
% ::
%
%     function M = spm_general_hist(hP,mP,textlab,[suppress plot - enter anything])
%
% :Inputs:
%
%   **hP:**
%        list of file names to compute histograms from
%
%   **mP:**
%        list of file names to compute masks from
%
%   **textlab:**
%        text string, e.g. 'ventricles' to label output tiffs
% 
% :Output:
%
%   histograms for all input images (usually contrast images from individual
%   subjects) plotted against a normal curve.
%
% The expected output is that each image will have roughly mean 0, with
% bumps or tails in the distribution of there are real activations in some
% parts of the brain.
%
% Looking at these histograms may be helpful for detecting outliers or 
% subjects with strange contrast values.  These may be caused by
% bad scaling, multicolinearity in the design matrix, acquisition artifacts,
% task-correlated head movement, or ???
%
% Histograms (blue) are overlaid on a Gaussian distribution (red)
% with a mean of 0 and a standard deviation equal to that of the observed data.
%
% ..
%    Tor Wager, 10/5/02
% ..

doplot = 1;
if length(varargin) > 0, doplot = 0;, end

mypwd = pwd;

disp([' Running histograms: Text label is ' textlab])

fprintf(1,'\nLoading volumes.\t')

v = spm_read_vols(spm_vol(hP));

    % ----------------------------------------------------------------------------------
    % * load mask, and mask all subjects' contrast values
    % so that we show only those voxels existing for all subjects.
    % ----------------------------------------------------------------------------------

    vm = spm_read_vols(spm_vol(mP));
    vm = (vm ~= 0); vm(vm == 0) = NaN;
    vm = repmat(vm,[1 1 1 size(v,4)]);
    v = v .* vm;


% ----------------------------------------------------------------------------------
% * figure out number of subplots in x and y dims
% ----------------------------------------------------------------------------------

[xd] = floor(sqrt(size(hP,1))); [yd] = ceil(sqrt(size(hP,1)));
while prod([xd yd]) < size(hP,1)
    if xd > yd - 2, yd = yd + 1; , else, xd = xd + 1;, end
end

% ----------------------------------------------------------------------------------
% * overall histogram to get bins
% ----------------------------------------------------------------------------------

fprintf(1,'Getting bins.\t')

if doplot
    [Finter,Fgraph,CmdLine] = spm('FnUIsetup','Display',1);
    figure(Finter);cla; set(gcf,'Color','w')
end

[x,a] = dohist(v,[],textlab,doplot);

% x(a < length(v(:))/18000) = [];     % eliminate relatively empty bins for indiviudal plots

if doplot
    title('All image data input into RFX analysis.')
    saveas(gcf,['all_hist' textlab],'tif')
end

% ----------------------------------------------------------------------------------
% * make individual image histograms
% ----------------------------------------------------------------------------------
if doplot
    figure(Fgraph)
end

fprintf(1,'Making individual histograms.\t')
M = zeros(1,4);

for i = 1:size(hP,1)
    
    subplot(yd,xd,i); hold on
    tmp = v(:,:,:,i);
    [dum,dum,M(i,:)] = dohist(tmp(:),x,textlab,doplot);
    legend off
    nm = deblank(hP(i,:)); if length(nm) > 25, nm = nm(end-25:end);, end, 
    nm(nm == '_') = ' '; nm = nm(1:end-4);
    if doplot, title(nm), end
    
end

if doplot
    saveas(gcf,['Indiv_hist' textlab],'tif')
end

return




% ----------------------------------------------------------------------------------
%
% * sub-functions
% 
% ----------------------------------------------------------------------------------

function [x,a,M] = dohist(v,x,textlab,doplot)

v = v(:);
v(isnan(v)) = [];

if isempty(x),
    [a,x] = hist(v,1000);
else
    [a] = hist(v(:),x);
end

if doplot
    hold on
    plot([0 0],[0 max(a)+1],'r')
end

m = nanmean(v(:)); sd = nanstd(v(:));

if doplot
    plot([m m],[0 max(a)+1],'b')
    legend({'normal' 'data'})

    plot(x,a,'b')


    np = normpdf(x,0,sd);
    sc = sum(a) ./ sum(np);     % scaling area to 1 necessary for hist overlay
    % this scales the same way as used in histfit.m

    %sc = max(a(2:end-1)) ./ max(np);

    h=plot(x,np .* sc,'r');
end

sk = skewness(v(:));
kur = kurtosis(v(:));
M = [m sd sk kur];

if doplot
    str = sprintf('M  %3.2f\nSd %3.2f\nSk %3.2f\nKt %3.2f\n%s', m,sd,sk,kur,textlab);
    text(max(x) - max(x) ./ 2,max(a) - max(a) ./ 2,str)

    drawnow
end

return
