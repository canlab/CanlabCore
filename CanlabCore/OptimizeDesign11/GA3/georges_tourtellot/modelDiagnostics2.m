function M = modelDiagnostics2(finput,varargin)
% M = modelDiagnostics2(finput,optional arguments)
%
% M = modelDiagnostics2(stimList,'stimList','all','ISI',1.2,'TR',2,'contrasts',[1 1 -1 -1;1 -1 1 -1],'screen')
%
%M = 
%        stimlist: your original list of conditions
%           delta: delta matrix sampled at .01 seconds
%           model: model sampled at .01 seconds
%       modelatTR: model sampled at TR
%        smoothed: temporally smoothed model at TR
%        conmodel: predictors for contrasts
%       consmooth: temporally smoothed predictors for contrasts
%           power: power (abs(fft)) matrix for smoothed model
%        conPower: power matrix for smoothed contrasts
%          energy: sum of squared deviations from mean for each predictor
%       conEnergy: sum of squared deviations from mean for each contrast
%             eff: efficiency of design.  1/trace(xtxi)
%  resampleeffinf: efficiency loss due to resampling compared to TR = 1 s.  eff(TR=1) / eff(TR)
%    smootheffinf: loss of efficiency due to smoothing; eff(nosmooth) / eff(smooth)
%         varbeta: variance of each individual beta
%          coneff: efficiency of contrast set
%      convarbeta: efficiency of each contrast beta
%             vif: variance inflation factors for predictors
%          conVif: variance inflation factors for contrasts
%            cond: condition number of model; larger is worse, > 30 is bad
%         conCond: condition number for contrast set
%           colin: correlation coefficients between individual predictors
%        conColin: correlation coefficients between contrasts
%       cbalfreqs: counterbalancing frequencies up to 3rd order
%            freq: frequencies with which each stimulus type occurs
%
% options: inptType is 'stimList', 'sampledModel','plot' (plot only)
%          outptType is 'all', 'short'
%
%
% More examples:
% M = modelDiagnostics2(M.stimlist,'stimList','all','ISI',.85,'TR',2,'contrasts',[1 1 -1 -1],'screen')
% M2 = modelDiagnostics2(M,'plot')
%
% Tor Wager, 11/17/01

%------------------------------------------------------------------------------
% Set up arguments and default values
%------------------------------------------------------------------------------
print2screen = 1;
outptType = 'all';
isfilt = 'n';
doverb = 0;
xc = [];

for i = 1:length(varargin)
    if isstr(varargin{i})
        switch varargin{i}
        case 'ISI', ISI = varargin{i+1};
        case 'TR', TR = varargin{i+1};
        case 'HRF', HRF = varargin{i+1};
        case 'HPlength', HPlength = varargin{i+1};
        case 'LPsmooth',LPsmooth = varargin{i+1};
        case 'contrasts',contrasts = varargin{i+1};
		case 'saturate',satval = varargin{i+1};
		case 'contrastweights',contrastweights = varargin{i+1};
		case 'maxOrder',maxOrder = varargin{i+1};
		case 'xc',xc = varargin{i+1};
		case 'screen',print2screen = 1;
		case 'noplot',print2screen = 0;
		case 'brief',outptType = 'brief';
		case 'plot',outptType = 'plot';
		case 'verbose',doverb = 1;		
        end
    end
end


if doverb,disp('Computing model descriptives.'),end

% -----------------------------------------------------------------------------------
% * determine output type and just plot if only plot is specified
% -----------------------------------------------------------------------------------

switch outptType
case 'all'
    dopower = 1;doenergy = 1;doeff = 1;dovif = 1;docondnum = 1; dohrf = 1;,
    docolin = 1; docbal = 1;dofreq = 1;dotimes = 0;
case 'brief'
    if doverb,disp('Brief output selected.'),end
    dopower = 0; doenergy = 0; doeff = 0; dovif = 0; docondnum = 0; docolin = 1;   
    docbal = 0; dofreq = 1; dotimes = 0; dohrf = 0;
case 'plot'
	M = finput;,if doverb,disp('Using already-computed M and plotting.'),end
	
    M = setup_blockcontrast_fix(M);
    
	if ~isfield(M.ga,'TR'),M.ga.TR = input('M.ga.TR not found: Enter TR in s: ');,end
	if ~isfield(M.ga,'HPlength'),M.ga.HPlength = input('M.ga.HPlength not found: Enter HPlength in s: ');,end
	
	% plot predictors
    plotvarious(M,M.ga.TR,'Predictors')

	% make Q structure with contrast information
    Q.model = M.model * M.contrasts';
    try,Q.smoothed = M.consmooth;,catch,disp('Error: contrasts should be empty but they''re not?'),contrasts,end
	Q.modelatTR = resample(Q.model,1,M.ga.TR*10);

	if isempty(M.ga.HPlength), specstring = 'none';,else specstring = 'specify';,end
	[S,KL,KH] = use_spm_filter(M.ga.TR,size(Q.modelatTR,1),'hrf',specstring,M.ga.HPlength); 
	if isempty(KH),KH = eye(size(Q.modelatTR));,end
    Q.filtered = Q.modelatTR - KH*(KH'*Q.modelatTR);
    Q.delta = M.delta * M.contrasts(:,1:end-1)';

	% plot contrasts
	plotvarious(Q,M.ga.TR,'Contrasts')

    return
	
otherwise eval(['help modelDiagnostics']),error('unknown output type')
end         % end switch


% 1. get every finput type to the stage of model at TR

% -----------------------------------------------------------------------------------
% * determine input type and set up.
% -----------------------------------------------------------------------------------
if 		isstruct(finput),inptType = 'Mstruct';
elseif 	size(finput,2) == 1, inptType = 'stimList';
elseif 	size(finput,2) > 1, inptType = 'model'
else	error('Illegal input.')	
end

switch inptType
case 'Mstruct'
	if isfield(finput,'ga'),opt = finput.ga;
	else opt = finput;
	end
	
	M = finput;

    
    if isfield(M,'contrasts'), contrasts = M.contrasts;, end
    
    N = fieldnames(opt);
    for i = 1:length(N)
        eval([N{i} ' = opt.' N{i} ';'])
    end
    
	if ~isfield(M,'stimlist'), error('M must have stimlist field.'), end
	%if isfield(opt,'ISI'), ISI = opt.ISI;, end
	%if isfield(opt,'TR'), TR = opt.TR;, end
	%if isfield(opt,'nonlinthreshold'), satval = opt.nonlinthreshold;, end
	%if isfield(opt,'HPlength'), HPlength = opt.HPlength;, end
	%if isfield(opt,'LPsmooth'), LPsmooth = opt.LPsmooth;, end
	%if isfield(opt,'contrastweights'), contrastweights = opt.contrastweights;, end
	%if isfield(opt,'maxOrder'), maxOrder = opt.maxOrder;, end
	%if isfield(opt,'xc'), xc = opt.xc;, end
	%if isfield(opt,'contrasts'), contrasts = opt.contrasts;, end
	
    M.stimlist = double(M.stimlist);
	numsamps = ceil(size(M.stimlist,1) * ISI ./ TR);
	
	
case 'model'
	M.model = finput;
	M.modelatTR = M.model;
	if ~(sum(M.model(:,size(M.model,2))) == size(M.model,1)) 
		disp('No intercept detected in last column:  Adding intercept.')   
		M.model(:,size(M.model,2)+1) = 1;
		M.modelatTR(:,size(M.modelatTR,2)+1) = 1;
        isfilt = input('Is the model filtered? (y/n) ','s');
        if strcmp(isfilt,'y')
            M.smoothed = M.modelatTR;
        end    
	end
    numsamps = size(M.modelatTR,1);
    dohrf = 0;      % can''t do hrf without delta function input
    dofreq = 0;
    docbal = 0;

case 'stimList'
	
	numsamps = ceil(size(finput,1) * ISI ./ TR);
	M.stimlist = finput;
    M.stimlist = double(M.stimlist);

end 	% end switch



% -----------------------------------------------------------------------------------
% * build X model matrix, if necessary
% -----------------------------------------------------------------------------------
if ~strcmp(inptType,'model')
	if ~exist('ISI') == 1, ISI = input('Enter ISI in s: ');,end
	if ~exist('TR') == 1, TR = input('Enter TR in s: ');,end

%gst changed line below, because if satval was empty, it should
%remain empty!=> no non-linear thresholding
    if ~exist('satval') == 1,  satval=[]; ,end %satval = 2; ,end    % satval = input('Enter saturation threshold in s: ');,end

    if ~exist('HRF') == 1, 
    	HRF = spm_hrf(.1);
    	HRF = HRF / max(HRF);
	end

	[M.modelatTR,M.model,M.delta] = rna2model(M.stimlist,ISI,HRF,TR,numsamps,satval,[]);
end


% -----------------------------------------------------------------------------------
% Add-on to add block-level (sustained) regressor, if specified
% tor, may 2008
% -----------------------------------------------------------------------------------
M = setup_blockcontrast_fix(M);
if strcmp(inptType, 'Mstruct'), 
    if isfield(M,'contrasts'), contrasts = M.contrasts; end
    contrastweights = M.ga.contrastweights;
end

% 2. set up contrast and smoothing matrices


% -----------------------------------------------------------------------------------
% * set up filtering/smoothing and contrasts, do filtering
% -----------------------------------------------------------------------------------
if ~exist('TR') == 1, TR = input('Enter TR in s: ');,end
if ~exist('HPlength') == 1, HPlength = input('Enter HPlength in s: ');,end
if ~exist('LPsmooth') == 1, LPsmooth = input('Enter LPsmooth value (1,0): ');,end
if ~exist('contrasts') == 1, contrasts = input('Enter contrast matrix: ');,end
if ~exist('contrastweights') == 1, contrasts, contrastweights = input('Enter contrast weight vector: ');,end
	
if doverb,fprintf(1,'setting up filtering...'),end
    
% get smoothing matrix
[S,Vi,svi] = getSmoothing(HPlength,LPsmooth,TR,numsamps,xc);

% set up contrast weights
if isempty(contrastweights),contrastweights = ones(1,size(contrasts,1));,end

% set up contrasts for testing individual predictors (contrasts empty)
if isempty(contrasts), 
	contrasts = eye(size(M.model,2) - 1);
	contrasts(:,end+1) = 0;
end

% make sure contrasts are the right length
if size(M.model,2) > size(contrasts,2)
		% warning('Model is larger than contrasts - extra intercept?  Adding to contrasts...')
		contrasts = [contrasts zeros(size(contrasts,1),size(M.model,2) - size(contrasts,2))]; 
end

if size(M.model,2) < size(contrasts,2)
    warning('Model is smaller than contrasts - extra zero in contrasts?  Truncating...') 
	contrasts = contrasts(:,1:size(M.model,2));
end

if isempty(S), S = 1;,end
if ~(strcmp(inptType,'model')) | strcmp(isfilt,'n')
    M.smoothed = S * M.modelatTR;
end

% -----------------------------------------------------------------------------------
% * assign other output values to M structure
% -----------------------------------------------------------------------------------	

M.contrasts = contrasts;
M.ga.contrasts = contrasts;
M.ga.TR = TR;
if exist('ISI') == 1, M.ga.ISI = ISI;, end
M.ga.HPlength = HPlength;
M.ga.LPsmooth = LPsmooth;
M.ga.contrastweights = contrastweights;
M.ga.xc = xc;




% -----------------------------------------------------------------------------------
% * make contrast-space predictors
% -----------------------------------------------------------------------------------	
if isempty(contrasts)
    contrasts = eye(size(M.model,2));
end

try
	M.conmodel = M.modelatTR * contrasts';
catch
		whos contrasts
		M
		warning('Contrasts are different size from modelatTR.  This should not happen.  Adding intercept 0 to contrast and trying again...')
		contrasts = [contrasts zeros(size(contrasts,1),1)];
		M.contrasts = contrasts;
end
	
M.conmodel = M.modelatTR * contrasts';
M.consmooth = M.smoothed * contrasts';


% -----------------------------------------------------------------------------------
% * calculate individual field values
% -----------------------------------------------------------------------------------

if doenergy
    M.energy = var(M.smoothed) * size(M.smoothed,1); 
    M.conEnergy = var(M.consmooth) * size(M.consmooth,1); 
end

if doeff                                                        		% calculate efficiency
    if doverb,fprintf(1,'eff...'),end	
	xtxitx = pinv(M.smoothed);                                       	% inv(X'S'SX)*(SX)'; pseudoinv of (S*X)
	[M.eff_fitness,M.eff] = calcEfficiency(M.ga.contrastweights,contrasts,xtxitx,svi);
    M.eff = M.eff';
	M.se_contrasts = 1 ./ M.eff;
    
    % ========= calc inefficiency of sampling and smoothing ===============     
	myeff = 1/trace(inv(M.smoothed'*M.smoothed));                 % brief form of model efficiency, not considering autocorrelation
    noresampeff = 1/trace(inv(M.model'*M.model));                 % what would have happened at TR = .1?
    nosmootheff = 1/trace(inv(M.modelatTR'*M.modelatTR));
    M.resample_eff_loss = noresampeff / nosmootheff;
    M.smoothing_eff_loss = nosmootheff / myeff;
	
end

if dohrf
    if doverb,fprintf(1,'DX...'),end
    delta = [];
    for i = 1:max(M.stimlist(:,1))
        delta(:,i) = (M.stimlist == i);
    end   
    %[M.DX] = tor_make_deconv_mtx2(M.delta,round(30 / M.ga.TR),10 * M.ga.TR);
    [M.DX] = tor_make_deconv_mtx2(delta,round(30 / M.ga.TR),M.ga.TR / M.ga.ISI);
    if ~isempty(S), M.DX = S * M.DX; end
    
    % Fix for adding block-level predictor
    if M.ga.trans2block && isfield(M.ga, 'blockcontrast') && M.ga.blockcontrast
        M.DX = [M.blkmodel M.DX];
    end

    xtxitx = pinv(M.DX);                                       		% inv(X'S'SX)*(SX)'; pseudoinv of (S*X)
    [M.hrf_eff_fitness,M.hrf_eff] = calcEfficiency([],[],xtxitx,svi);
    M.hrf_eff_avg = mean(M.hrf_eff(1:end-1));
end
   
if dovif
    if doverb,fprintf(1,'VIFs...'),end
	warning off
    M.vif = getvif(M.smoothed);   
    if ~isempty(M.consmooth)
        M.conVif = getvif(M.consmooth);
    end
	warning on
end
 
if docondnum                                                    % calculate condition numbers
    M.cond = cond(M.smoothed);
    if ~isempty(M.consmooth)
        M.conCond = cond(M.consmooth);
    end
end
warning off
if docolin                                                      %   calculate predictor colinearity
    if doverb,fprintf(1,'colin...'),end
    % sqrt of average squared correlation 
      colin = corrcoef(M.smoothed);
      M.colin = colin;
      if ~isempty(M.consmooth)
        colin = corrcoef(M.consmooth);
        M.conColin = colin;
      end
end
warning on
% calculate counterbalancing
if docbal
    if doverb,fprintf(1,'cbal...'),end
	if ~exist('maxOrder') == 1, maxOrder = input('Enter max counterbalancing order: ');,end

	conditions = 1:max(M.stimlist);
	freqConditions = ones(1,size(conditions,2)) / size(conditions,2);
    [M.cBal M.cbalfreqs] = getCounterBal(M.stimlist, maxOrder,conditions,freqConditions);
    % M.cbalfreqs = round(M.cbalfreqs * 100);
end

if dofreq
    for i = 1:size(conditions,2)
         M.freq(i) = sum(M.stimlist == conditions(i)) / size(M.stimlist,1);
    end
end

% calculate average time between stimulus repetitions in each trial type
if isfield(M,'stimlist')
    for i = 1:max(M.stimlist),M.timeBtwn(i) = mean(diff(find(M.stimlist==i))) .* M.ga.ISI;, end
end

if print2screen
    format short g
    format compact

    plotvarious(M,TR,'Predictors')
    if ~isempty(M.consmooth)
        Q.model = M.model * contrasts';
		Q.modelatTR  = M.modelatTR * contrasts';
        Q.smoothed = M.consmooth;
        Q.delta = M.delta * contrasts(:,1:end-1)';
        plotvarious(Q,TR,'Contrasts')
    end
end


if doverb
    disp(' ')
    disp('Correlated columns:')
    [x,y] = find(abs(corrcoef(xX.X)) > .5 & abs(corrcoef(xX.X)) < 1); [x y]
end
return


% subfunctions
% ================================================================== %
function plotvarious(M,TR,textlbl)


    if strcmp(textlbl,'Contrasts'), contit = 1;, else, contit = 0;, end
    
    % image the model
    % --------------------------------------


   figure; set(gcf,'Color','w')
    imagesc(M.modelatTR);colormap(copper)
    numsamp = size(M.smoothed,1);

    if contit
        title('Optimized Contrasts Heatmap, after filtering','FontSize',14)
    else
        title('Optimized Model Heatmap, after filtering','FontSize',14)
    end
    
    if isfield(M,'delta')
    	numregs = size(M.delta,2);
	% subtract one b/c last one is assumed to be intercept? sometimes...not always.
	
    else
	numregs = size(M.smoothed,2);	% subtract one b/c last one is assumed to be intercept? sometimes...not always.
    end

    nummodel = size(M.model,1);
    numsecs = numsamp * TR;
    deltascale = .1;
    x = 1:deltascale:numsecs+1;
    x = x(1:nummodel);
    x2 = 1:TR:numsecs+1;
    x2 = x2(1:numsamp);
    ymin = min(min(M.smoothed));
    ymax = max(max(M.smoothed));
    

    % Plot all columns in model after resampling and smoothing
    % --------------------------------------

    figure;set(gcf,'Color','w')
    for i = 1:numregs
  	  subplot(numregs,1,i)
      hold on
      plot(x,M.model(:,i),'k','LineWidth',1.5)
      plot(x2,M.modelatTR(:,i),'b','LineWidth',1.5)
      plot(x2,M.smoothed(:,i),'r','LineWidth',1.5)
      if isfield(M,'delta')
        for j = 1:size(M.delta,1)
            if M.delta(j,i) == 1
                plot([j*deltascale j*deltascale],[ymin 0],'k')
            end
        end
      end
      set(gca,'Ylim',[ymin ymax]);
      %axis([0 100 ymin ymax])
      set(gca,'XColor',[0 0 1])
      set(gca,'YColor',[0 0 1])
    end
    subplot(numregs,1,1);
    if contit
        title('Optimized Model Contrasts and Onsets','FontSize',14)
    else
        title('Optimized Model Regressors and Onsets','FontSize',14)
    end
        
    legend({'high-resolution' 'sampled at TR' 'filtered'})
    %title(textlbl,'FontSize',14)
    %title([textlbl ' Green: high-resolution   Blue: resampled at TR   Red: sampled and HP/LP filtered (smoothed)'])   
    drawnow
    
    % Plot first column at various stages
    % --------------------------------------

    figure;set(gcf,'Color','w')
        subplot(3,1,1)
        ymin = min(min(M.model));
        ymax = max(max(M.model));
        plot(x,M.model(:,1),'k','LineWidth',2)
        title([textlbl '(#1): Original model - high samp rate'])
        set(gca,'Ylim',[ymin ymax]);
        if contit
            title('Optimized Model - first Contrast','FontSize',14)
        else
            title('Optimized Model - first Regressor','FontSize',14)
        end
        
        %axis([0 100 -1 2])
        subplot(3,1,2)
        plot(x2,M.modelatTR(:,1),'k','LineWidth',2)
        title('Resampled at TR')
        set(gca,'Ylim',[ymin ymax]);
        %axis([0 100 -1 2])
        subplot(3,1,3)
        plot(x2,M.smoothed(:,1),'k','LineWidth',2)
        title('Filtered')
        set(gca,'Ylim',[ymin ymax]); 
        % axis([0 100 -1 2])
        
    drawnow
    
return



% -----------------------------------------------------------------------------------
% Add-on to add block-level (sustained) regressor, if specified
% tor, may 2008
% -----------------------------------------------------------------------------------
function M = setup_blockcontrast_fix(M)

if M.ga.trans2block && isfield(M.ga, 'blockcontrast') && M.ga.blockcontrast

    disp('Adding block predictor.');
    numsamps = size(M.modelatTR, 1);
    
    a = [ones(M.ga.restevery, 1); zeros(M.ga.restlength, 1)]; ab = [a; -a];
    blkstimlist = repmat(ab, ( ceil(2 * size(M.stimlist, 1)) ./ length(ab) ), 1);
    blkmodel = sampleInSeconds(blkstimlist,M.ga.ISI);
    
    hrf = spm_hrf(1 ./ (M.ga.TR * 10)); hrf = hrf ./ max(hrf);
    
    blkmodel = conv(blkmodel, hrf);              % convolve with canonical
    
    % add to hi-res
    if isfield(M, 'model')
        M.model = [blkmodel(1:size(M.model, 1)) M.model];
    end
    
    if isfield(M, 'delta')
        tmp = zeros(size(M.delta, 1), 1);
        tmp(1:length(a):length(tmp)) = 1;
        M.delta = [tmp M.delta];
    end
    
    blkmodel = resample(blkmodel,1,M.ga.TR*10); % downsample to TR
    blkmodel = blkmodel(1:numsamps, :);
    
    M.modelatTR = [blkmodel M.modelatTR];
    
    % filter
    if isempty(M.ga.HPlength), specstring = 'none';,else specstring = 'specify';,end
    [S,KL,KH] = use_spm_filter(M.ga.TR,size(M.modelatTR,1),'hrf', specstring, M.ga.HPlength); 
    M.smoothed = M.modelatTR - KH*(KH'*M.modelatTR);
    
    % contrasts
    M.conmodel = [blkmodel M.conmodel];
    M.consmooth = [blkmodel M.consmooth];
    
    % put in contrasts
    if size(M.contrasts, 2) < size(M.model, 2)
        M.contrasts = blkdiag(1, M.contrasts);
    end
    
    if size(M.ga.contrastweights, 2) < size(M.contrasts, 1)
        M.ga.contrastweights = [1 M.ga.contrastweights];
    end
    
    M.DX = [blkmodel M.DX]; % not always used; may never be used?
    M.blkmodel = blkmodel;
    
    disp('Remember to put in extra element in contrastweights for the block predictor.');
    
end

return

