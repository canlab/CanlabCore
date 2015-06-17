function [fitness,convar,xtxi,models] = testlist(stimLists,GA,HRF,varargin)
% [fitness,convar,xtxi,models] = testlist(stimLists,GA,HRF,[opt args])
%
% stimLists can have multiple lists (columns)
%
% 9/12/01 Tor Wager

clockStart = cputime;			  	% keeps track of starting time

for i = 1:length(varargin)
    if isstr(varargin{i})
        switch varargin{i}
        case 'rest', restMatrix = varargin{i+1};
        case 'svi', svi = varargin{i+1};
        case 'S', S = varargin{i+1};
        end
    end
end

% ----------------------------------------------------------------
% * vars setup
% ---------------------------------------------------------------- 	

trans2switch = GA.trans2switch;
trans2block = GA.trans2block;
maxOrder = GA.maxOrder;
conditions = GA.conditions;
freqConditions = GA.freqConditions;
cbalColinPowerWeights = GA.cbalColinPowerWeights;
contrastweights = GA.contrastweights;
NumStimthresh = GA.NumStimthresh;
maxCbalDevthresh = GA.maxCbalDevthresh;
maxFreqDevthresh = GA.maxFreqDevthresh;
ISI = GA.ISI;
TR = GA.TR;
nonlinthreshold = GA.nonlinthreshold;
scanLength = GA.scanLength;
restevery = GA.restevery;
restlength = GA.restlength;
HPlength = GA.HPlength;
xc = GA.xc;
contrasts = GA.contrasts;

% ----------------------------------------------------------------
% * contrast setup
% ----------------------------------------------------------------
  
% work out whether to use contrasts.
if isempty(contrasts)
   nconds = sum(conditions > 0);
   %contrasts = zeros(nconds,nconds);
    %for i = 1:nconds
    %    contrasts(i,i) = 1;
    %end
    %disp('  ...no contrasts entered; using original model')
    docontrasts = 0;
 else 
    %disp(' ...contrasts entered')   
    docontrasts = 1;
    %contrasts
end
    
%if ~(size(contrasts,1) == sum(conditions > 0)), 
    %disp('	...fyi: Num contrasts does not equal num conditions')
%end


% ----------------------------------------------------------------
% * setup freqConditions and rests
% ---------------------------------------------------------------- 

% normalize freqConditions
if ~sum(freqConditions) == 1,
	freqConditions(1:end-1) = (1 - freqConditions(end)) / (size(freqConditions,2)-1);
	disp(['	...freqConditions does not sum to 1: normalizing to ' num2str(freqConditions)])
end

% flag for rest length.
	if ~isempty(restevery) &  ~isempty(restlength)
		disp(['	...Using rests of length ' num2str(restlength) ' and resting every ' num2str(restevery)])
		dorests = 1;
	else disp('	...No rests specified.'),dorests = 0;
	end

if ~(mod(scanLength / ISI,1) == 0),disp('Warning: Scan length in s is not an even multiple of ISI!'),end


% ----------------------------------------------------------------
% * initial computation of list lengths and output
% ----------------------------------------------------------------  

numStim = ceil(scanLength / (ISI));
if dorests,
	numRestStim = (ceil(numStim/(mean(restevery)+restlength)) - 1) * restlength;
else
	numRestStim = 0;
end
numStimEachCond = ceil((numStim-numRestStim) * freqConditions);			            % row vector of stim in each cond
numsamps = ceil(numStim*ISI/TR);


% ----------------------------------------------------------------
% * get smoothing matrix and autocorrelation matrix
% ----------------------------------------------------------------
if isempty(HPlength),HPlength = 'none';,end
if ~isfield(GA,'LPsmooth'),GA.LPsmooth = 1;,end

dofilter = 1;
if strcmp(HPlength,'none') & ~GA.LPsmooth
    dofilter = 0;
end

%disp(['	...setting smoothing, hrf, and autocorrelation matrices, HPlength = ' num2str(HPlength) ', Smoothing = ' num2str(GA.LPsmooth)])

if ~(exist('S') == 1)
%use spm_filter
if strcmp(HPlength,'none') & GA.LPsmooth									% LP only
	%disp('		...LP only')
	[S,KL] = use_spm_filter(TR,numsamps,'hrf','none',[]);
elseif GA.LPsmooth															% HP and LP
	[S,KL,KH] = use_spm_filter(TR,numsamps,'hrf','specify',HPlength);
	%disp('		...HP and LP')
elseif strcmp(HPlength,'none') & ~GA.LPsmooth								% neither
	%disp('		...no filtering')
	S = eye(numsamps);
	dofilter = 0;
else [S,KL,KH] = use_spm_filter(TR,numsamps,'none','specify',HPlength);		% HP only
	%disp('		HP only')
end

end % if ~exist S

if ~(exist('svi') == 1)
    
if isempty(xc),Vi = eye(numsamps);,disp('Using white noise autocorrelation')
else Vi = getv('make',xc,numsamps);
end

svi = S * Vi * S';

end % if exist svi

clear KL, clear KH, clear Vi, clear myscannerV



% ----------------------------------------------------------------
% * criterion measures setup
% ---------------------------------------------------------------- 
if ~isempty(NumStimthresh) | ~isempty(maxCbalDevthresh) | ~isempty(maxFreqDevthresh)
	docriterion = 1;
	if isempty(NumStimthresh),NumStimthresh = 10000;,end
 	if isempty(maxCbalDevthresh),maxCbalDevthresh = 10000;,end
	if isempty(maxFreqDevthresh),maxFreqDevthresh = 10000;,end
else
	docriterion = 0;
end

maxrestthresh = 2;					% max number of rests in a row.

% ----------------------------------------------------------------
% * saturation setup
% ----------------------------------------------------------------
if isempty(nonlinthreshold) | strcmp(nonlinthreshold,'none'),dosaturation = 0;,nonlinthreshold = 'none';,else dosaturation = 1;,end


% add intercept
% ------------------------------
if docontrasts, contrasts(:,end+1) = 0;, end





% ----------------------------------------------------------------
% * do it.
% ---------------------------------------------------------------- 

for z = 1:size(stimLists,2)

	stimList = stimLists(:,z);

% The Switch Hack: transform list of stimuli into list of switches/no switches for each stimtype (R or L)
      if trans2switch, stimList = transform2switches(stimList);,end
      if trans2block,
	     restMatrix = varargin{1};
	     restlength = GA.restlength; 

             restList = restMatrix(:,z);
             stimList = transform2block(stimList,restList,restlength);
      end

      
      % ===== counterbalancing ============================================  
      %if (cbalColinPowerWeights(1) > 0)    removed 3/30 to implement criterion rather than weights.
   		[cBal,dummy,maxDev] = getCounterBal(stimList, maxOrder,conditions,freqConditions);
		 % ...if optimizing this measure, add it to fitness scores.
         if (cbalColinPowerWeights(1) > 0),fitnessMatrix(1,z) = cBal;,end
	  %end
   
	  
	  % ===== frequency discrepancy ========================================
	  % calculate max discrepancy between actual frequencies and requested frequencies 
	  %if size(cbalColinPowerWeights,2) > 3                  
      	%if (cbalColinPowerWeights(4) > 0)    	 
            for i = 1:size(freqConditions,2) 
         		freqMat(i) = sum(stimList == conditions(i)) / size(stimList,1);
      	    end
            freqMat = freqMat ./ sum(freqMat);   % adjust for probe ''0''s
			maxFreqDev = max(abs(freqConditions - freqMat));
			% ...if optimizing this measure, add it to fitness scores.
      		if (cbalColinPowerWeights(4) > 0), dummy = 1 - maxFreqDev;,end 
		%end
	  %end
	  
	  
	  % criterion measures
	  % ====================================================================
	    go = 1;
	    if docriterion
	  	[maxNumStim,maxrest] = getMaxInARow(stimList,dojitter); 
	  	if maxNumStim > NumStimthresh | maxrest > maxrestthresh | maxDev > maxCbalDevthresh | maxFreqDev > maxFreqDevthresh
		  go = 0;
		  error('Criteria not met.')
          	end
	   end
	
	 
		% =====power calculation =============================================
    if (cbalColinPowerWeights(2) > 0 | cbalColinPowerWeights(3) > 0) & go   % build predictor set of vectors and convolve with HRF
        model = sampleInSeconds(stimList,ISI);
        model = getPredictors(model,HRF);
        model = resample(model,1,TR*10);
	if dosaturation,model = modelSaturation(model,nonlinthreshold);,end						% saturation (nonlinear responses)
        if size(model,1) > numsamps, model = model(1:numsamps,:);,end

	% add intercept
	% ------------------------------
	model(:,end+1) = 1;
   
	

        if dofilter
		try
            model = S * model;                                                  % temporal smoothing and HP filter
        catch
            whos model
            whos S
            error('model smoothing: wRoNg Sizes!')
        end
	end
        xtxitx = pinv(model);                                       		% inv(X'S'SX)*(SX)'; pseudoinv of (S*X)

	
	if docontrasts,
		try
            con1 =  model * contrasts(1,:)';
        catch
            whos model
            whos contrasts
            error('contrast multiplication: wRoNg Sizes!')
        end
	end
	%figure; subplot(2,1,1)
	%plot(model(:,1))
	%title(['Model ' num2str(z) ' Predictor 1'])
	%subplot(2,1,2)
	%plot(con1)
	%title(['Model ' num2str(z) ' Contrast 1'])
	%drawnow; pause(5)

	xtxi{z} = inv(model'*model);
	models{z} = model;

        if docontrasts											% compute the variance of the contrasts
           try
              fitness(z) = 1./(contrastweights * diag(contrasts*xtxitx*svi*xtxitx'*contrasts'));
	      convar(z) = (contrasts(1,:)*xtxitx*svi*xtxitx'*contrasts(1,:)');
              % it's 1/ because we want to minimize the variance, but we maximize fitnessMatrix value.
           catch
             disp('using contrasts: wrong sizes!!!')
             contrasts	
             whos xtxitx
             whos svi
             diag(xtxitx*svi*xtxitx')
             error('Exiting...')
           end
           
        else	
           try
               fitness(z) = 1./(contrastweights * diag(xtxitx*svi*xtxitx')); % variance of chosen regressors.
		myconvar = diag(xtxitx*svi*xtxitx');
	       convar(z) = myconvar(1);
           catch
               disp(['no contrasts version: wrong sizes!!!'])
               whos contrastweights
               whos xtxitx
               whos svi
               diag(xtxitx*svi*xtxitx')
               error('Exiting...')
           end
	end
    end

end % loop through designs

%figure;for i = 1:length(models),subplot(1,length(models),i),imagesc(models{i}),
%	title(['Model ' num2str(i) ' v = ' num2str(convar(i)) ' f = ' num2str(fitness(i))])
%end

disp(['testlist done: ' num2str(cputime-clockStart) ' s']) 


return
