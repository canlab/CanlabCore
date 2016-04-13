function EXPT = wb_hewma_shell(EXPT,varargin)
% :Usage:
% ::
%
%     EXPT = wb_hewma_shell(EXPT,[start slice or vector of slices],[dools],[mask])

% Run EWMA first: wb_multisubject_ewma.
%
% Fields in EXPT needed:
%   - subdir
%   - im_files
%   - FILES.ewma_z
%   - FILES.ewma_var
%   - [cov]
%   - [mask]
%
% ..
%    Tor Wager, 10/2/05
% ..

global xdim
global ydim
global emptyimage
global V
global domask
global lam
global f
global subjnames
global basepts

fprintf(1,'-------------------------------\n')
fprintf(1,'Hierarchical Ewma - whole brain\n')
fprintf(1,'-------------------------------\n')

subjnames = EXPT.subjects;

if isempty(basepts), disp('Assuming 60 baseline timepoints for mixture modeling!'); basepts = 60;, end

% -----------------------------------------------------
% Get sample image to get dimensions
% -----------------------------------------------------
if exist(EXPT.im_files{1}(1,:)) == 2
    V = spm_vol(deblank(EXPT.im_files{1}(1,:)));
    xdim = V.dim(1);
    ydim = V.dim(2);
    zdim = V.dim(3);
    emptyimage = zeros([xdim ydim zdim]);       % for writing slices of stat values
    disp(['Found sample image: ' V.fname]);
else
    error('Enter an image file in EXPT.im_files{1}(1,:) to get dimensions of data images.');
end
   
%load lambda
try
    load([EXPT.subjects{1} filesep 'lam.mat']);
    disp(['Lambda found: lam = ' num2str(lam)]);
catch
    disp(['Cannot find lambda parameter saved in lam.mat in ' EXPT.subjects{1}])
end

% -----------------------------------------------------
% Get slice .mat file names
% -----------------------------------------------------
if ~isfield(EXPT,'FILES'), EXPT.FILES = [];,end
if ~isfield(EXPT.FILES,'ewma_z'), 
    disp('Cannot find EXPT.FILES.ewma_z.  Looking for z_slice*.mat');
    EXPT = getfunctnames2(EXPT,'z_slice*.mat','FILES.ewma_z');  
else
    disp('Found EXPT.FILES.ewma_z.  Assuming z_slice*.mat images are there and using those.');
end

if ~isfield(EXPT.FILES,'ewma_var'), 
    disp('Cannot find EXPT.FILES.ewma_z.  Looking for var_*.mat');
    EXPT = getfunctnames2(EXPT,'var_*.mat','FILES.ewma_var');
else
    disp('Found EXPT.FILES.ewma_var.  Assuming var_*.mat images are there and using those.');
end

if zdim ~= size(EXPT.FILES.ewma_z{1},1)
    disp('Error: Number of slices in EXPT.FILES.ewma_z{1} does not match that of sample image file.');
end


% -----------------------------------------------------
% center and scale predictors
% -----------------------------------------------------
if ~isfield(EXPT,'cov'), EXPT.cov = [];, end
covt = EXPT.cov;

if ~isempty(covt)
    disp('Found behavioral covariate. Regressing against change points using robust IRLS')
    if length(covt) > size(covt,1), covt = covt';, end
    
    covt = scale(covt);
    
    %covt = covt - repmat(mean(covt),size(covt,1),1);
    %for i = 1:size(covt,2)
    %    covt(:,i) = covt(:,i) ./ var(covt(:,i));
    %end
else
    disp('<no behavioral covariates>')
end

% intercept: automatically added if doing robustfit.m
% covt = [ones(1,size(covt,1)) covt];
%covt(:,end+1) = 1;

EXPT.cov = covt;


% ----------------------------------------------------
% * process optional input args
% ----------------------------------------------------

if length(varargin)>0,wh=varargin{1};,else,wh=1:length(EXPT.FILES.ewma_z);,end
if length(varargin)>1,dools = varargin{2};, else, dools = 1;, end
if length(varargin)>2,domask = varargin{3};, else, domask = 0;, end

if length(wh) < 2, wh = wh:size(EXPT.FILES.ewma_z{1},1);, end  % it's a start slice, fill in rest

if domask,disp(['Using mask: ' domask]), else, disp('<no mask image specified.>'),end

% ----------------------------------------------------
% * set up the gui figure
% ----------------------------------------------------
f = [];
%f = figure('Color','w');
%tmp = get(gcf,'Position') .* [1 1 .5 .1];
%set(gcf,'Position',tmp)
%set(gcf,'MenuBar','none','NumberTitle','off')
%figure(f), set(gca,'Xlim',[0 100])



% ----------------------------------------------------
% * create subdirectory
% ----------------------------------------------------
dirindex = 1;
if dirindex < 10, myz = '000';, else, myz = '00';, end
mydir = ['hewma' myz num2str(dirindex)];
eval(['mkdir ' mydir])
cd(mydir)


% ----------------------------------------------------
% * set up vars to save and save setup
% ----------------------------------------------------

% variables containing timeseries values
global grpmean
global grpt
global grpste
global grpp
global grph

% save setup stuff
try
    SETUP.ewma_stat_files = EXPT.FILES.ewma_z;
    SETUP.var_files = EXPT.FILES.ewma_var;
    SETUP.covariates = covt;
    SETUP.V = V;
    save SETUP SETUP
catch
    warning('Error creating SETUP file');
end


% ----------------------------------------------------
% * start slice -- load data if necessary
% ----------------------------------------------------
if wh(1) ~= 1 
    disp('Continuing where we left off...loading hewma_timeseries');
    load hewma_timeseries grpmean grpt grpste lam xdim ydim
end



% ----------------------------------------------------
% * run hewma for selected slices
% ----------------------------------------------------
    %if dools, disp('Running OLS and IRLS comparison (slower) - to turn this off, use 0 as 3rd input argument.'),end
    disp('____________________________________________________________')
    

for i = wh
    % * run robust HLM for each slice 
    % ----------------------------------------------------
    warning off
    fprintf(1,'Slice %3.0f. > ',i);

    if dools
        rob_fit(EXPT.FILES.ewma_z,EXPT.FILES.ewma_var,covt,i,dools,domask);
    else
        rob_fit(EXPT.FILES.ewma_z,EXPT.FILES.ewma_var,covt,i,dools,domask);
    end
    warning on
    fprintf(1,'\n');
end


disp('Saving contiguous clusters in cl');
try
    cl = hewma_save_timeseries('hewma_sig.img',10,grpmean,grpste,xdim,ydim);
    save hewma_cl cl;
catch
end


cd ..



return







% slice-analysis function




function [newP,newP2] = rob_fit(zimgs,varimgs,covt,index,dools,domask)

global xdim
global ydim
global emptyimage
global V
global domask

global grpmean
global grpt
global grpste
global grpp
global grph

global lam
global f
global subjnames

global basepts


% --------------------------------------------
% load data
% --------------------------------------------
t1 = clock; fprintf(1,'Loading data. ');

L = size(zimgs{1},1);   % number of slices
if index > L, error('Slice index > number of slices.');,end
N = length(zimgs);      % number of subjects

for i = 1:N
    tmp = load(deblank(zimgs{i}(index,:))); % the slice ewma stat data
    
    % make sure sizes match
    if i == 1, dims = size(tmp.zdat);, end
    if any(size(tmp.zdat) - dims), disp(['Subject ' num2str(i) ': different slice dims!']);,end
    
    z(:,:,i) = full(tmp.zdat(1:dims(1),1:dims(2)));
    
    tmp = load(deblank(varimgs{i}(index,:))); % the slice variance data
    v(:,:,i) = full(tmp.vardat(1:dims(1),1:dims(2)));
end
    
tp = size(z,2); % time points

fprintf(1,'%3.0f s. ',etime(clock,t1));
    
% --------------------------------------------
% find the in-analysis voxels
% --------------------------------------------

% mask, if spc
if domask
    Vm = spm_vol(domask);, vm = spm_read_vols(Vm);,

    % reshape mask for this slice
    vm = squeeze(vm(:,:,index));
    vm = reshape(vm,prod(size(vm)),1);
    vm = repmat(vm,[1 tp N]);
    
    z = z .* vm;
    v = v .* vm;
end

wh=sum(z,2);
wh = squeeze(wh);
wh = ~any(wh == 0 | isnan(wh),2);
wh = find(wh);      % which voxels are non-zero, non-nan for all subjects

tmp = length(wh); if tmp == 0, disp('No voxels in analysis!'), end


% --------------------------------------------
% set up output arrays
% --------------------------------------------

emptymask =  NaN .* zeros(size(z,1),1);
% to be written in .img files
tvals = emptymask;
pvals = emptymask;
sbmean = emptymask;
hvals = emptymask;
zvals = emptymask;
tthresh = emptymask;
cp = emptymask;
runlen = emptymask;

tvalsdiff = emptymask;
pvalsdiff = emptymask;

hvalsdiff = emptymask;
zvalsdiff = emptymask;
tthreshdiff = emptymask;
cpdiff = emptymask;
runlendiff = emptymask;
        
grpmean{index} = zeros(size(z(:,:,1))) .* NaN;     % save betas
grpt{index} = zeros(size(z(:,:,1))) .* NaN;    % save t values
grpp{index} = zeros(size(z(:,:,1))) .* NaN;    % save p values
grpste{index} = zeros(size(z(:,:,1))) .* NaN;    % save ste values
grph{index}  = zeros(size(z(:,:,1))) .* NaN;    % save hyp test indicator

% mixture modeling
gausscp = emptymask;    % Gaussian-mixture model change point est
gaussnruns = emptymask; % number of runs of OOC points
gausstot = emptymask; 
gausslongest = emptymask;
gaussfirstlen = emptymask;
gausstotaldur = emptymask;

stats = []; stats.cp_ind = []; stats.pdiff = [];

% individual change points
    for subj = 1:length(subjnames)
        cpind{subj} = emptymask;
    end

 cpind3d = NaN .* zeros(xdim,ydim,length(subjnames));

 
 if ~isempty(covt)
     % these output arrays are voxels x predictors
     npreds = size(covt,2);
     covtb = NaN .* zeros(size(z,1),npreds);
     covtt = NaN .* zeros(size(z,1),npreds);
     covtp = NaN .* zeros(size(z,1),npreds);
 end
        
    
    
    
fprintf(1,'. %6.0f voxels. Done:    ',length(wh))


% --------------------------------------------
% perform regression
% --------------------------------------------
et = clock;

for i = 1:length(wh)
    vindx = wh(i);                  % which voxel, x * y (rows), will be reshaped later
    
    dat = squeeze(z(vindx,:,:))';
    vdat = squeeze(v(vindx,:,:))';

    % hewma2, main fitting function
    % with linear detrending; adds between-groups analysis if covt is not
    % empty
    [p,tm,Zcor,sbm,Zpop,ttime,sb,stats] = hewma2(dat,vdat, lam, 0, 1, covt);
    
    % timeseries things to save
    grpmean{index}(vindx,:) = Zpop;     % group timeseries
    grpt{index}(vindx,:) = ttime;       % t-values for each timepoint
    grpste{index}(vindx,:) = sb;        % between-subjects standard dev (sigma)

    
    % summary overall statistics
    tvals(vindx,1) = tm;
    pvals(vindx,1) = p;
    sbmean(vindx,1) = mean(sbm);           % between-subjects st. dev (sigma), pooled over time (sterr)
    hvals(vindx,1) = p < .05;
    zvals(vindx,1) = Zcor;
    tthresh(vindx,1) = stats.tthresh;
    cp(vindx,1) = stats.cp;
    runlen(vindx,1) = stats.maxrunlength;
    
    if ~isempty(covt)
        tvalsdiff(vindx,1) = stats.tmdiff;
        pvalsdiff(vindx,1) = stats.pdiff;

        hvalsdiff(vindx,1) = stats.pdiff < .05;
        zvalsdiff(vindx,1) = stats.Zcordiff;
        tthreshdiff(vindx,1) = stats.tthreshdiff;
        cpdiff(vindx,1) = stats.cpdiff;
        runlendiff(vindx,1) = stats.maxrunlengthdiff;
    end
         
    % individual change points
    for subj = 1:length(stats.cp_ind)
        cpind{subj}(vindx,1) = stats.cp_ind(subj);
    end
    
    % -------------------------------------------------------  
    % do this stuff on significant voxels only
    % -------------------------------------------------------
    
    if p < .05 % | stats.pdiff < .05
        
        % correlation between behavioral covt and change points
        if ~isempty(covt)
            [b,rstat] = robustfit(covt,scale(stats.cp_ind'));
            b = b(2:end); t = rstat.t(2:end); p = rstat.p(2:end);

            % these output arrays are voxels x predictors
            covtb(vindx,:) = b';
            covtt(vindx,:) = t';
            covtp(vindx,:) = p';
        end

        % mixture modeling
        [ind,ind2,gstats] = Gaussian_mix(mean(dat),50,basepts,0,0); 
        
        gausscp(vindx,1) = gstats.cp;    % Gaussian-mixture model change point est
        gaussnruns(vindx,1) = gstats.cnt; % number of runs of OOC points
        gausstot(vindx,1) = gstats.tot;   
        gausslongest(vindx,1) = gstats.longest;   
        gaussfirstlen(vindx,1) = gstats.firstlen;   
        gausstotaldur(vindx,1) = gstats.totaldur;   

    end
    

    if rem(i,10) == 0
        %figure(f), try,barh(100*i / length(x)),catch,end
        %set(gca,'Xlim',[0 100]),set(gca,'YTickLabel',i),drawnow
        %text(5,.5,['Voxel ' num2str(i) ' of ' num2str(length(x))],'Color','r')
        fprintf(1,'\b\b\b%02d%%',round(100*i / length(wh)));
    end
   
end
fprintf(1,'\tDone in %3.0f s\n',etime(clock,et))

save hewma_timeseries grpmean grpt grpste lam xdim ydim



% -------------------------------------------------------
%reshape
% -------------------------------------------------------

% main effects and zero-crossings
tvals = reshape(tvals,xdim,ydim);
pvals = reshape(pvals,xdim,ydim);
sbmean = reshape(sbmean,xdim,ydim);
hvals = reshape(hvals,xdim,ydim);
zvals = reshape(zvals,xdim,ydim);
tthresh = reshape(tthresh,xdim,ydim);
cp = reshape(cp,xdim,ydim);
runlen = reshape(runlen,xdim,ydim);
    
% group differences
if ~isempty(covt)
    tvalsdiff= reshape(tvalsdiff,xdim,ydim);
    pvalsdiff= reshape(pvalsdiff,xdim,ydim);
    hvalsdiff= reshape(hvalsdiff,xdim,ydim);
    zvalsdiff= reshape(zvalsdiff,xdim,ydim);
    tthreshdiff= reshape(tthreshdiff,xdim,ydim);
    cpdiff= reshape(cpdiff,xdim,ydim);
    runlendiff= reshape(runlendiff,xdim,ydim);
end

% mixture modeling
gausscp = reshape(gausscp,xdim,ydim);
gaussnruns = reshape(gaussnruns,xdim,ydim);
gausstot = reshape(gausstot,xdim,ydim);
gausslongest = reshape(gausslongest,xdim,ydim);
gaussfirstlen = reshape(gaussfirstlen,xdim,ydim);
gausstotaldur = reshape(gausstotaldur,xdim,ydim);


% individual change points
for subj = 1:length(stats.cp_ind)
    cpind3d(:,:,subj) = reshape(cpind{subj},xdim,ydim); 
end
  
if ~isempty(covt)
    % these output arrays are voxels x predictors
    covtb = reshape(covtb,xdim,ydim,npreds);  
    covtt = reshape(covtt,xdim,ydim,npreds);  
    covtp = reshape(covtp,xdim,ydim,npreds);  
end       
   

% -------------------------------------------------------
%write image slices
% -------------------------------------------------------

Pt = write_beta_slice(index,V,tvals,emptyimage,'hewma_t');
Pp = write_beta_slice(index,V,pvals,emptyimage,'hewma_p');
Ps = write_beta_slice(index,V,sbmean,emptyimage,'hewma_s');
Ph = write_beta_slice(index,V,hvals,emptyimage,'hewma_sig');
Pz = write_beta_slice(index,V,zvals,emptyimage,'hewma_z');
Pthr = write_beta_slice(index,V,tthresh,emptyimage,'hewma_thresh');
Pc = write_beta_slice(index,V,cp,emptyimage,'hewma_cp');
Pc = write_beta_slice(index,V,runlen,emptyimage,'hewma_runlen');

if ~isempty(covt)
    Pt = write_beta_slice(index,V,tvalsdiff,emptyimage,'hewma_tdiff');
    Pp = write_beta_slice(index,V,pvalsdiff,emptyimage,'hewma_pdiff');
    Ph = write_beta_slice(index,V,hvalsdiff,emptyimage,'hewma_sigdiff');
    Pz = write_beta_slice(index,V,zvalsdiff,emptyimage,'hewma_zdiff');
    Pthr = write_beta_slice(index,V,tthreshdiff,emptyimage,'hewma_threshdiff');
    Pc = write_beta_slice(index,V,cpdiff,emptyimage,'hewma_cpdiff');
    Pc = write_beta_slice(index,V,runlendiff,emptyimage,'hewma_runlendiff');
end

Pt = write_beta_slice(index,V,gausscp,emptyimage,'gausscp');
Pt = write_beta_slice(index,V,gaussnruns,emptyimage,'gaussnruns');
Pt = write_beta_slice(index,V,gausstot,emptyimage,'gausstot');
Pt = write_beta_slice(index,V,gausslongest,emptyimage,'gausslongest');
Pt = write_beta_slice(index,V,gaussfirstlen,emptyimage,'gaussfirstlen');
Pt = write_beta_slice(index,V,gausstotaldur,emptyimage,'gausstotaldur');

Pcp = write_cp_slice(index,V,cpind3d,emptyimage,subjnames);

if ~isempty(covt)
    % define names
    for pred = 1:npreds
        bnms{pred} = ['Covt' num2str(pred) '_cp_correl_b'];
        tnms{pred} = ['Covt' num2str(pred) '_cp_correl_t'];
        pnms{pred} = ['Covt' num2str(pred) '_cp_correl_p'];
    end
    
    Ppredb = write_cp_slice(index,V,covtb,emptyimage,bnms);
    Ppredt = write_cp_slice(index,V,covtt,emptyimage,tnms);
    Ppredp = write_cp_slice(index,V,covtp,emptyimage,pnms);
end


return







function Pw = write_beta_slice(slicei,V,betas,emptyimg,varargin)
% Pw = write_beta_slice(sliceindex,V,data,empty,prefix)
% Slice-a-metric version
warning off % due to NaN to int16 zero conversions
V.dim(4) = 16; % set to float to preserve decimals

prefix = 'image_';
if length(varargin) > 0, prefix = varargin{1};,end

for voli = 1:size(betas,3) % for each image/beta series point
    %if voli < 10, myz = '000';, elseif voli < 100, myz = '00';, else myz = '000';,end
    %V.fname = [prefix myz num2str(voli) '.img'];
    V.fname = [prefix '.img'];
    V.descrip = ['Hewma output image ' num2str(voli)];

    % create volume, if necessary
    if ~(exist(V.fname) == 2), spm_write_vol(V,emptyimg);,end
        
    spm_write_plane(V,betas(:,:,voli),slicei);
end

Pw = which(V.fname);

warning on
return



function Pw = write_cp_slice(slicei,V,betas,emptyimg,varargin)
% Pw = write_con_slice(sliceindex,V,data,empty,prefix)
% Slice-o-matic version
% prefix is a cell array of names

warning off % due to NaN to int16 zero conversions
V.dim(4) = 16; % set to float to preserve decimals

% default names
for coni = 1:size(betas,3)
    prefix{coni} = ['Subject_' num2str(coni)];
end

if length(varargin) > 0, prefix = varargin{1};,end

for voli = 1:size(betas,3)  % for each image/beta series point
 
    V.fname = [prefix{voli} '.img'];
    V.descrip = ['Hewma cp image.'];

    % create volume, if necessary
    if ~(exist(V.fname) == 2), spm_write_vol(V,emptyimg);,end
        
    spm_write_plane(V,betas(:,:,voli),slicei);
    
    Pw{voli} = which(V.fname);
end



warning on
return

