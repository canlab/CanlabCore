function [out,P] = ga3power(pv,mspec,conditions,varargin)
% [out,P] = ga3power(pv,mspec,conditions,[contrasts])
%
% This function creates a power plot of a design with onsets (ons)
% against random designs, block designs, and m-sequences with
% similar design constraints.  Plots are in HRF estimation power (x)
% vs. contrast detection (y)
%
% onsets in seconds, cell array, one cell per condition
% TR in s
%


fprintf(1,'Setup . ')

% ----------------------------------------------
% SETUP
% ----------------------------------------------
nsess = 1;          % number of sessions, could be added as input
res = 4;           % resolution, in samples per second, for testing contrasts
TR = mspec.TR;      % rep time (sampling rate)
% time points to estimate - recommend 30 s worth
truer = [];         % true resp; empty = default in hrf_power, xpower
nsubs = [];          % number of subjects; empty = default
btwvar = [];        % between-subs var, empty = default

% onsets in model
% ----------------------------------------------
if isempty(pv)
    [ons,X2,pv] = paramvec2onsets(mspec,conditions,[]);
else
    [ons] = paramvec2onsets(mspec,conditions,pv);
end
out.actual.onsets = ons;

% time points to estimate - recommend 30 s worth
% ----------------------------------------------
tpest = [];
for i = 1:length(conditions), 
    for j = 1:length(conditions(i).subcond),
        if isfield(conditions(i).subcond,'hrfest')
            tpest(end+1) = conditions(1).subcond(1).hrfest;
        else
            tpest(end+1) = round(30 ./ TR);
        end
    end
end
        
% timeseries length and HRF
% ----------------------------------------------        
%len = round(max(cat(1,ons{:})) * res);
%len2 = ceil(max(cat(1,ons{:})) ./ TR);
len = mspec.numframes * res;        % length of hi-res timeseries
len2 = mspec.numframes;             % length of lo-res timeseries

cf = zeros(len,1);

hrf = spm_hrf(1/res);   % canonical hrf
hrf = hrf ./ max(hrf);
H = convmtx(hrf,len);

% Contrasts
% ----------------------------------------------
if length(varargin) > 0
    P.contrasts = varargin{1};
else
    %P.contrasts = [1 -1; 1 0; 0 1; 1 1];
    P.contrasts = eye(length(ons));
end

% add 0 contrast values for intercept, if necessary
if size(P.contrasts,2) < length(ons) + 1 ./ nsess, P.contrasts(:,end+1) = 0;,end

% needed for HRF estimation (only if inputting delta function later)
%P.ISI = TR; P.TR = TR;

% autocorrelation
% ----------------------------------------------
[xc,P.Vi] = canonical_autocorrelation(TR,len2);



% ----------------------------------------------
% BUILD DESIGN
% see also getDesign5.m
% onsets2delta.m
% ----------------------------------------------

[X,X2] = get_design(conditions,mspec,pv);

%[P,X] = ons2cf(ons,H,TR,res,len2);
% P has hi res and low res delta functions

% ----------------------------------------------
% TEST
% ----------------------------------------------

fprintf(1,'Test . ')


[c] = efficiency(X,P);    % [c,h] = efficiency(X,P);

    out.actual.contrasts = c;
    out.actual.contrast_avg = mean(c);
    out.actual.worst_contrast = min(c);
    %out.actual.hrf_avg = mean(h);
    %out.actual.worst_hrf = min(h);
    
% power: default true, default subjects, default intersubject variability
[ipow,gpow] = xpower(X,P.contrasts,truer,nsubs,btwvar,P.Vi);
out.actual.ipow = ipow;
out.actual.gpow = gpow;

% power for HRF - first contrast only, no intercept
PP = hrf_power(TR,tpest,[],P.contrasts(1,1:length(ons)),[],X2);

ipowh = PP.ipow;
gpowh = PP.gpow;
out.actual.ipowhrf = ipowh;
out.actual.gpowhrf = gpowh;

% ----------------------------------------------
% BOOTSTRAP ...a sampling distribution of similar designs
% ----------------------------------------------
fprintf(1,'Random . ')

Ptmp = P;
niter = 100;

% alt way: get delta with cf2X or efficiency.m, permute
% doesn't produce designs with all the right constraints, tho.
%delta = cell2mat(P.delta);

for i = 1:niter
    Xtmp = [];
    
    [Xtmp,X2tmp] = get_design(conditions,mspec,[]); % random monte-carlo design
    
    
    [ctmp] = efficiency(Xtmp,Ptmp);        % [ctmp,htmp] = efficiency(Xtmp,Ptmp);
    [ipow,gpow] = xpower(Xtmp,P.contrasts,truer,nsubs,btwvar,P.Vi);


    PP = hrf_power(TR,tpest,[],P.contrasts(1,1:length(ons)),[],X2tmp) ;    % no smoothing, use existing X2
    ipowh = PP.ipow;
    gpowh = PP.gpow;
    
    out.permuted.ipowhrf(:,i) = ipowh;
    out.permuted.gpowhrf(:,i) = gpowh;

    out.permuted.contrasts(:,i) = ctmp;
    out.permuted.contrast_avg(i) = mean(ctmp);
    out.permuted.worst_contrast(i) = min(ctmp);
    %out.permuted.hrf_avg(i) = mean(htmp);
    %out.permuted.worst_hrf(i) = min(htmp);
    out.permuted.ipow(:,i) = ipow;
    out.permuted.gpow(:,i) = gpow;
    
end

out.summary.contrast_avg = mean(out.permuted.contrast_avg);
out.summary.worst_contrast = mean(out.permuted.worst_contrast);
%out.summary.hrf_avg = mean(out.permuted.hrf_avg);
%out.summary.worst_hrf = mean(out.permuted.worst_hrf);
    
% ----------------------------------------------
% TEST m-sequences
% m-seqs must be constrained to follow other
% design constraints (min trial time, etc.)
% ----------------------------------------------

fprintf(1,'M-seq . ')

warnflag = 0;
for i = 1:length(conditions), if length(conditions(i).subcond) > 1, warnflag = 1;,end,end
if warnflag, warning('M-sequence power does not apply to designs with multiple parts!'),end
        
%out.mseq = genMseq(mlen,TR,P,10,length(conditions));  % 10 sequences, this
%tests with these params, but we need to constrain the design and do it in
%a separate step.
% this will create sequences of length mlen, randomly divided into trial
% types (conditions)

[tmp,out.mseqs] = genMseq(len2,[],P,10,length(conditions));  % 10 sequences

for i = 1:length(out.mseqs)
    
    % get paramvec, which constrains length to correct values
    [tmp,pvm] = construct_model(mspec,conditions,out.mseqs{i});
    [Xm,X2m] = get_design(conditions,mspec,pvm);

    %[onsm,X2m,pvm] = paramvec2onsets(mspec,conditions,pvm); %[X23] = plotDesign(onsm,[],1.5);   % test

    [ctmp] = efficiency(Xm,X2m);        % [ctmp,htmp] = efficiency(Xtmp,Ptmp);
    [ipow,gpow] = xpower(Xm,P.contrasts,truer,nsubs,btwvar,P.Vi);
    PP = hrf_power(TR,tpest,[],P.contrasts(1,1:length(ons)),[],X2m) ;    % no smoothing, use existing X2
    
    ipowh = PP.ipow;
    gpowh = PP.gpow;
    
    out.mseq.ipowhrf(:,i) = ipowh;
    out.mseq.gpowhrf(:,i) = gpowh;

    out.mseq.contrasts(:,i) = ctmp;
    out.mseq.contrast_avg(i) = mean(ctmp);
    out.mseq.worst_contrast(i) = min(ctmp);
    %out.mseq.hrf_avg(i) = mean(htmp);
    %out.mseq.worst_hrf(i) = min(htmp);
    out.mseq.ipow(:,i) = ipow;
    out.mseq.gpow(:,i) = gpow;
    
end

out.summary.mseq_contrast_avg = mean(out.mseq.contrast_avg);
out.summary.mseq_worst_contrast = mean(out.mseq.worst_contrast);
%out.summary.hrf_avg = mean(out.mseq.hrf_avg);
%out.summary.worst_hrf = mean(out.mseq.worst_hrf);

% ----------------------------------------------
% PLOT
% ----------------------------------------------


colordef white

figure('Color','w'); hold on;

% permuted designs
tmp1 = max(out.permuted.gpowhrf);
tmp2 = out.permuted.gpow(1,:);
plot(tmp1,tmp2,'bo','MarkerFaceColor','b')

% m-sequences
plot(max(out.mseq.gpowhrf),out.mseq.gpow(1,:),'rs','MarkerFaceColor','r')

% actual design
plot(max(out.actual.gpowhrf),out.actual.gpow(1),'kd','MarkerFaceColor','y')
plot(max(out.actual.gpowhrf),out.actual.gpow(1),'ko','MarkerSize',8)


set(gca,'FontSize',16);
xlabel('Peak HRF (Z)')
ylabel('Con(1) (Z)')
title('80% Power Z-scores (N = 12, rand. fx.)')

return




function [P,X] = ons2cf(ons,H,TR,res,len2)
% P has hi res and low res delta functions
% ons = onsets in s
% H = hi-res conv matrix
% res = samps per second
% len2 is design length in TRs

for i = 1:length(ons)
    
    cf2(:,i) = cf;
    cf2(round(ons{i}*res)+1,i) = 1;
    cf2 = cf2(1:len,:);

end

X = cf2X(cf2,H,TR*res,len2);

%cf2 = [cf2; zeros(size(X,1)*TR*res - size(cf2,1),size(cf2,2))];
P.delta_hires = cf2;

% FIR - for HRF shape estimation

cf = zeros(len2,1);
cf2 = [];

for i = 1:length(ons)
    
    cf2{i} = cf;
    cf2{i}(round(ons{i}./TR)+1) = 1;
    cf2{i} = cf2{i}(1:len2);
    
end

P.delta = cf2;

return




function X = cf2X(cf2,H,sampby,len2);
X = H * cf2;
%X = resample(X,1,TR*16);
X = X(1:sampby:end,:);
%%%X = modelSaturation(X,2,[]);  % soft threshold 2, exp = .5

% make sure length is as it should be
X = [X; zeros(len2 - size(X,1), size(X,2))];

X(:,end+1) = 1;	


return


function [X,X2,pv] = get_design(conditions,mspec,pv)

% canonical HRF
for cc = 1:length(conditions)
    for ss = 1:length(conditions(cc).subcond)
        conditions(cc).subcond(ss).convolve = 1;
    end
end
mspec2 = rmfield(mspec,'contrasts'); % get rid of contrasts cause no dx!
[X] = construct_model(mspec2,conditions,pv,mspec2.toverlap);

% FIR - deconvolution
for cc = 1:length(conditions)
    for ss = 1:length(conditions(cc).subcond)
        conditions(cc).subcond(ss).convolve = 0;
    end
end
[X2] = construct_model(mspec,conditions,pv,mspec.toverlap);

return






% Extra code

%randomize - permutation test
    % randomize, build model at TR resolution
    dhr = getRandom(P.delta_hires);
    Ptmp.delta = dhr;
    
    Xtmp = cf2X(dhr,H,TR*res,len2);
     
    [ctmp,htmp,tmp,X2tmp] = efficiency(Xtmp,Ptmp);
    [ipow,gpow] = xpower(Xtmp,P.contrasts,truer,nsubs,btwvar,P.Vi);

    [tmp,dhr2] = downsample_delta(dhr,TR*res);
    PP = hrf_power(TR,tpest,dhr2,P.contrasts);
    ipowh = PP.ipow;
    gpowh = PP.gpow;
    
    
    
    
    
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
