function [out,P] = onsets2efficiency(ons,TR,varargin)
% [out,P] = onsets2efficiency(ons,TR,[contrasts])
%
% A summary plot comparing a design with event onsets (ons) to 
% random designs and m-sequences. Tests first contrast entered.
% ** currently, permuted designs are used, and they are sampled
% differently (low-res) than the test design; causes probs!?
%
%
% onsets in seconds, cell array, one cell per condition
% TR in s
%
% Does not account for:
%   Smoothing / filtering
%   Nonlinearity in signal response
%
% Contrasts is a matrix of contasts
% Rows = contrasts, columns = predictors, values = weights
%
% tor wager
%
% calls:
%   efficiency.m
%   hrf_power.m
%   genMseq.m
%
%
% to-do: fixed set of contrasts right now
%
% Examples:
% [out,P] = onsets2efficiency(sac,2,[1 -1; 1 0; 0 1; 1 1]);
% [out,P] = onsets2efficiency(a,2,[eye(5); .33 .33 .33 0 0; 1/3 1/3 1/3 -.5 -.5]);
%
% To go from condition function to onsets (1.5 is TR):
% ons{1} = find(Run == 1) .* 1.5;
% ons{2} = find(Run == 2) .* 1.5;
%
% with ga2.m:
% [evt2,X2,pv2] = paramvec2onsets(mspec,conditions,[],[1 5]);
% [out,P] = onsets2efficiency(evt2,1.5,eye(3));

fprintf(1,'Setup . ')

% ----------------------------------------------
% SETUP
% ----------------------------------------------
nsess = 1;  % number of sessions, could be added as input
res = 16;   % resolution, in samples per second

len = round(max(cat(1,ons{:})) * res);
cf = zeros(len,1);

hrf = spm_hrf(1/res);   % canonical hrf
hrf = hrf ./ max(hrf);

% autocorrelation stuff
% replaced by canonical_autocorrelation.m
%f = inline('A * exp(-a * x)','A','a','x');  % exponential function
%xa = 0:TR:30;
%xc = f(.9475,.5323,xa);    % params at 3T from vnl experiment, n = 10
%xc = xc ./ max(xc);

% ----------------------------------------------
% BUILD DESIGN
% see also getDesign5.m
% onsets2delta.m
% ----------------------------------------------

for i = 1:length(ons)
    
    cf2(:,i) = cf;
    cf2(round(ons{i}*res)+1,i) = 1;
    cf2 = cf2(1:len,:);
    
    x = conv(hrf,cf2(:,i));
    X(:,i) = x(1:length(cf));
    
end

%X = resample(X,1,TR*16);
X = X(1:TR*res:end,:);
 %%%X = modelSaturation(X,2,[]);  % soft threshold 2, exp = .5
 
X(:,end+1) = 1;	

cf2 = [cf2; zeros(size(X,1)*TR*res - size(cf2,1),size(cf2,2))];
P.delta_hires = cf2;

% for HRF shape estimation
len2 = ceil(max(cat(1,ons{:})) ./ TR);
cf = zeros(len2,1);
cf2 = [];

for i = 1:length(ons)
    
    cf2{i} = cf;
    cf2{i}(round(ons{i}./TR)+1) = 1;
    cf2{i} = cf2{i}(1:len2);
    
end

P.delta = cf2;

% ----------------------------------------------
% TEST
% ----------------------------------------------

fprintf(1,'Test . ')

if length(varargin) > 0
    P.contrasts = varargin{1};
else
    %P.contrasts = [1 -1; 1 0; 0 1; 1 1];
    P.contrasts = eye(length(ons));
end

% add 0 contrast values for intercept, if necessary
if size(P.contrasts,2) < length(ons) + 1 ./ nsess, P.contrasts(:,end+1) = 0;,end

% needed for HRF estimation in efficiency.m (to get delta)
P.ISI = TR; 
P.TR = TR;

% autocorrelation
[xc,P.Vi] = canonical_autocorrelation(TR,size(X,1));

[c,h,tmp,X2,P] = efficiency(X,P);

    out.actual.contrasts = c;
    out.actual.contrast_avg = mean(c);
    out.actual.worst_contrast = min(c);
    out.actual.hrf_avg = mean(h);
    out.actual.worst_hrf = min(h);
    
% power: default true, 12 subjects, default intersubject variability
[ipow,gpow] = xpower(X,P.contrasts,[],12,[],P.Vi);
out.actual.ipow = ipow;
out.actual.gpow = gpow;

%
% power for HRF code removed from here - at end
%


%[ipowh,gpowh] = xpower(X2,P.hrfcontrasts,P.truehrf,12,[],P.Vi);
    PP = hrf_power(TR,18,P.delta,P.contrasts);
    ipowh = PP.ipow;
    gpowh = PP.gpow;
out.actual.ipowhrf = ipowh;
out.actual.gpowhrf = gpowh;

% ----------------------------------------------
% BOOTSTRAP ...a sampling distribution of similar designs
% ----------------------------------------------
fprintf(1,'Random . ')

Ptmp = P;
niter = 10;

H = convmtx(PP.hrf,size(X,1));
delta = cell2mat(P.delta);

for i = 1:niter
    Xtmp = [];
    
    % randomize, build model at TR resolution
    dhr = getRandom(delta);
    Ptmp.delta = dhr;
    
    Xtmp = H * dhr;
    Xtmp = Xtmp(1:size(X,1),:);
    Xtmp(:,end+1) = 1;
     
    [ctmp,htmp,tmp,X2tmp] = efficiency(Xtmp,Ptmp);
    [ipow,gpow] = xpower(Xtmp,P.contrasts,[],12,[],P.Vi);

    PP = hrf_power(TR,18,dhr,P.contrasts);
    ipowh = PP.ipow;
    gpowh = PP.gpow;
    
    out.permuted.ipowhrf(:,i) = ipowh;
    out.permuted.gpowhrf(:,i) = gpowh;

    out.permuted.contrasts(:,i) = ctmp;
    out.permuted.contrast_avg(i) = mean(ctmp);
    out.permuted.worst_contrast(i) = min(ctmp);
    out.permuted.hrf_avg(i) = mean(htmp);
    out.permuted.worst_hrf(i) = min(htmp);
    out.permuted.ipow(:,i) = ipow;
    out.permuted.gpow(:,i) = gpow;
    
end

out.summary.contrast_avg = mean(out.permuted.contrast_avg);
out.summary.worst_contrast = mean(out.permuted.worst_contrast);
out.summary.hrf_avg = mean(out.permuted.hrf_avg);
out.summary.worst_hrf = mean(out.permuted.worst_hrf);
    
% ----------------------------------------------
% TEST m-sequences
% ----------------------------------------------

fprintf(1,'M-seq . ')
out.mseq = genMseq(size(X,1),TR,P,10);  % 10 sequences

% ----------------------------------------------
% PLOT
% ----------------------------------------------
colordef white
figure('Color','w'); subplot(2,1,1); hold on;

plot(X),set(gca,'FontSize',16),title('Actual design'),axis off
str = [];
for i = 1:length(ons)
    str{i} = ['Cond' num2str(i)];
end
legend(str)

hh = axes('Position',[.13 .11 .3439 .3439]); hold on;
str = [];

for i = 1:size(P.contrasts,1)
    
    plot(ones(1,niter)*i,out.permuted.contrasts(i,:).^.5,'b.')
    plot(i,out.actual.contrasts(i).^.5,'kd','MarkerFaceColor','y')
    str{i} = ['[' num2str(P.contrasts(i,1:end-1)) ']'];
end

    plot(ones(1,niter)*i+1,out.permuted.hrf_avg,'b.')
    plot(i+1,out.actual.hrf_avg,'kd','MarkerFaceColor','y')
    str{i+1} = 'HRF^2';
    
    set(gca,'FontSize',16)
    title('Efficiency by contrast')
    set(gca,'XTick',1:i+1,'XTickLabel',str,'XLim',[.5 i+1.5])
    xlabel('Contrast')
    ylabel('sqrt(Efficiency)')
    
    
hh = axes('Position',[.55 .11 .3439 .3439]);hold on;

figure('Color','w'); hold on;
tmp1 = max(out.permuted.gpowhrf);
tmp2 = out.permuted.gpow(1,:);
plot(tmp1,tmp2,'bo','MarkerFaceColor','b')
plot(max(out.actual.gpowhrf),out.actual.gpow(1),'kd','MarkerFaceColor','y')
plot(max(out.actual.gpowhrf),out.actual.gpow(1),'ko','MarkerSize',8)

plot(max(out.mseq.gpowhrf),out.mseq.gpow(1,:),'rs','MarkerFaceColor','r')

xlabel('Peak HRF (Z)')
ylabel('Con(1) (Z)')
set(gca,'FontSize',16); title('80% Power Z-scores (N = 12, rand. fx.)')

return





% Extra code

% power for HRF
% ignores downsampling problem to save time
hrf2 = spm_hrf(TR); hrf2 = hrf2 ./ max(hrf2);
truehrf = []; conhrfmult = [];
tt = zeros(1,P.hrfsamples); 
for i = 1:length(P.delta),
    if P.contrasts(1,i)>0, 
        tt2=(tt+1).*hrf2(1:length(tt))';,
        conhrfm = tt + 1;
    else,
        tt2=tt;,
        conhrfm = tt; 
    end,
    truehrf=[truehrf tt2];,
    conhrfmult = [conhrfmult conhrfm];
end
P.hrfcontrasts = eye(size(X2,2) - 1);
P.hrfcontrasts(find(conhrfmult==0),:) = []; % contrasts only for HRFs of interest
P.truehrf = truehrf;




% full code for random iterations, with commented out stuff (older)
for i = 1:niter
    Xtmp = [];
    
    % randomize, build model at TR resolution
    %dhr = getRandom(P.delta_hires);
    dhr = getRandom(delta);
    
    Xtmp = H * dhr;
    %Xtmp = Xtmp(1:TR*res:end,:);   % downsample
    Xtmp = Xtmp(1:size(X,1),:);
    %%%Xtmp = modelSaturation(Xtmp,2,[]);  % soft threshold 2, exp = .5
    
    %for j = 1:length(P.delta)
    %    dtmp{j} = getRandom(P.delta{j});
    %    x = conv(hrf2,dtmp{j});
    %    Xtmp(:,j) = x(1:length(dtmp{j}));
    %end
    Xtmp(:,end+1) = 1;
    
    %Ptmp.delta = dtmp;
    %Ptmp.delta = downsample_delta(dhr,res*TR);
    %Ptmp.delta = [Ptmp.delta; zeros(size(X,1) - size(Ptmp.delta,1),size(Ptmp.delta,2))];
    Ptmp.delta = dhr;
     
    [ctmp,htmp,tmp,X2tmp] = efficiency(Xtmp,Ptmp);
    
    %b = [1 0 0]';
    %y = Xtmp * b + randn(size(Xtmp,1),1);
    %[b,dev,stats] = glmfit(Xtmp,y);
    %out.permuted.t(:,i) = stats.t;
    
    [ipow,gpow] = xpower(Xtmp,P.contrasts,[],12,[],P.Vi);
    
    %[ipowh,gpowh] = xpower(X2tmp,P.hrfcontrasts,P.truehrf,12,[],P.Vi);
    PP = hrf_power(TR,15,dhr,P.contrasts);
    ipowh = PP.ipow;
    gpowh = PP.gpow;
    
    out.permuted.ipowhrf(:,i) = ipowh;
    out.permuted.gpowhrf(:,i) = gpowh;

    out.permuted.contrasts(:,i) = ctmp;
    out.permuted.contrast_avg(i) = mean(ctmp);
    out.permuted.worst_contrast(i) = min(ctmp);
    out.permuted.hrf_avg(i) = mean(htmp);
    out.permuted.worst_hrf(i) = min(htmp);
    out.permuted.ipow(:,i) = ipow;
    out.permuted.gpow(:,i) = gpow;
    
end

