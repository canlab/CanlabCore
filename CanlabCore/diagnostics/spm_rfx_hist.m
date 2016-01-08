function spm_rfx_hist(cwd)
% :Usage:
% ::
%
%     function spm_rfx_hist(cwd)
%
% :Input:
%
%   a directory name where an spm or SnPM random effects analysis lives
%
% :Outputs:
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

mypwd = pwd;
eval(['cd ' cwd])

fprintf(1,'\nLoading volumes.\t')

d = dir('SPMcfg.mat');
d2 = dir('SnPMcfg.mat');
if length(d) > 0
    % then we have an SPM directory
    
    load SPMcfg
    P = str2mat(VY.fname);
    v = spm_read_vols(VY);

    % ----------------------------------------------------------------------------------
    % * load mask, and mask all subjects' contrast values
    % so that we show only those voxels existing for all subjects.
    % ----------------------------------------------------------------------------------

    vm = spm_read_vols(spm_vol('mask.img'));
    vm(vm == 0) = NaN;
    for i = 1:size(v,4)
        v(:,:,:,i) = v(:,:,:,i) .* vm;
    end

elseif length(d2) > 0
    % then its an SnPM directory
    
    load SnPMcfg
    v = spm_read_vols(spm_vol(P));
    
    vm = spm_read_vols(spm_vol('SPMt.img'));
    vm = (vm ~= 0); vm(vm == 0) = NaN;
    for i = 1:size(v,4)
        v(:,:,:,i) = v(:,:,:,i) .* vm;
    end
    
else
    disp('No SPMcfg.mat or SnPMcfg.mat file found in this directory!')
    eval(['cd ' mypwd])
    return
end

% ----------------------------------------------------------------------------------
% * figure out number of subplots in x and y dims
% ----------------------------------------------------------------------------------

[xd] = floor(sqrt(size(P,1))); [yd] = ceil(sqrt(size(P,1)));
while prod([xd yd]) < size(P,1)
    if xd > yd - 2, yd = yd + 1; , else, xd = xd + 1;, end
end

% ----------------------------------------------------------------------------------
% * overall histogram to get bins
% ----------------------------------------------------------------------------------

fprintf(1,'Getting bins.\t')
[Finter,Fgraph,CmdLine] = spm('FnUIsetup','Display',1);
figure(Finter);cla; set(gcf,'Color','w')

[x,a] = dohist(v,[]);

%x(a < length(v(:))/18000) = [];     % eliminate relatively empty bins for indiviudal plots
%x(a < sum(a) ./ 1000) = [];

title('All image data input into RFX analysis.')
saveas(gcf,'all_hist','tif')

% ----------------------------------------------------------------------------------
% * make individual image histograms
% ----------------------------------------------------------------------------------

figure(Fgraph)
fprintf(1,'Making individual histograms.\t')
for i = 1:size(P,1)
    
    subplot(yd,xd,i); hold on
    tmp = v(:,:,:,i);
    dohist(tmp(:),x);
    legend off
    nm = deblank(P(i,:)); if length(nm) > 25, nm = nm(end-25:end);, end, 
    nm(nm == '_') = ' '; nm = nm(1:end-4);
    title(nm)
    
end

saveas(gcf,'Indiv_hist','tif')

eval(['cd ' mypwd])

return




% ----------------------------------------------------------------------------------
%
% * sub-functions
% 
% ----------------------------------------------------------------------------------

function [x,a] = dohist(v,x)

v = v(:);
v(v == 0) = [];

v(isnan(v)) = [];

if isempty(x),
    [a,x] = hist(v,10000);
else
    [a] = hist(v(:),x);
end

hold on
plot([0 0],[0 max(a)+1],'r')
m = nanmean(v(:)); sd = nanstd(v(:));
plot([m m],[0 max(a)+1],'b')
legend({'normal' 'data'})

plot(x,a,'b')

np = normpdf(x,0,sd);
sc = sum(a) ./ sum(np);     % scaling area to 1 necessary for hist overlay
% this scales the same way as used in histfit.m

%sc = max(a(2:end-1)) ./ max(np);

h=plot(x,np .* sc,'r');

str = sprintf('M  %3.2f\nSd %3.2f\nSk %3.2f\nKt %3.2f', m,sd,skewness(v(:)),kurtosis(v(:)));
text(max(x) - max(x) ./ 2,max(a) - max(a) ./ 2,str)

drawnow
return
