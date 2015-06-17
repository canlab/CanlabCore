function RESULTS = runsim(SIM)
% function RESULTS = runsim(SIM)
%
% 	SIM is a structure with parameters
%	e.g., SIM.params is passed in from testMvsOther.m
%
%	RESULTS is a structure with results
%	e.g., stored as SIM(simindex).results in testMvsOther.m
%
%	Made by Tor Wager, 05/24/01

format compact


% ----------------------------------------------------------------
% * get all variable names out of the structure (for compatibility)
% ----------------------------------------------------------------  

HRF_estlen=SIM.HRF_estlen;  %gst (FEATURE)
NULL_flag=SIM.NULL_flag;  %gst (FEATURE)

HRF = SIM.HRF;
powerTvoxels = SIM.powerTvoxels;
powerTdf = SIM.powerTdf;
conditions = SIM.conditions; 
ISI = SIM.ISI;
noise_var = SIM.noise_var;
c = SIM.contrasts;  
niterations = SIM.niterations;  
freqConditions = SIM.freqConditions;
scanLength = SIM.scanLength;
TR = SIM.TR;
cbalColinPowerWeights = SIM.cbalColinPowerWeights;
numGenerations = SIM.numGenerations;
sizeGenerations = SIM.sizeGenerations;
lowerLimit = SIM.lowerLimit;
maxOrder = SIM.maxOrder;
plotFlag = SIM.plotFlag;  
beta = SIM.beta;
HPlength = SIM.HPlength;
xc = SIM.xc;
restlength = SIM.restlength;
restevery = SIM.restevery;
trans2switch = SIM.trans2switch;
trans2block = SIM.trans2block;  
NumStimthresh = SIM.NumStimthresh;  
maxCbalDevthresh = SIM.maxCbalDevthresh;
maxFreqDevthresh = SIM.maxFreqDevthresh; 
nonlinthreshold = SIM.nonlinthreshold;
maxrestthresh = SIM.maxrestthresh;
if ~isfield(SIM,'LPsmooth')
    SIM.LPsmooth = 1;
else 
    LPsmooth = 1; %gst
%    LPsmooth = SIM.LPsmooth;  %gst
end
if ~isfield(SIM,'dofirst')
	dofirst = 0;
else 
	dofirst = SIM.dofirst;
end


% ----------------------------------------------------------------
% * setup freqConditions and rests
% ---------------------------------------------------------------- 

numStim = ceil(scanLength / (ISI));

% normalize freqConditions
if ~sum(freqConditions) == 1,
	freqConditions(1:end-1) = (1 - freqConditions(end)) / (size(freqConditions,2)-1);
	disp(['	...freqConditions does not sum to 1: normalizing to ' num2str(freqConditions)])
end

% flag for rest length.
	if ~isempty(restevery) &  ~isempty(restlength)
		disp(['	...Using rests of length ' num2str(restlength) ' and resting every ' num2str(restevery)])
		numRestStim = (ceil(numStim/(mean(restevery)+restlength)) - 1) * restlength;
		dorests = 1;
	else 
		disp('	...No rests specified.'),dorests = 0; numRestStim = 0;
	end

if ~(mod(scanLength / ISI,1) == 0),disp('Warning: Scan length in s is not an even multiple of ISI!'),end


% ----------------------------------------------------------------
% * general parameters
% ----------------------------------------------------------------  

fitnessMatrix = zeros(4,sizeGenerations);
numStimEachCond = ceil((numStim-numRestStim) * freqConditions);			            % row vector of stim in each cond
numsamps = ceil(numStim*ISI/TR);


% ----------------------------------------------------------------
% * saturation setup
% ----------------------------------------------------------------
if isempty(nonlinthreshold) | strcmp(nonlinthreshold,'none'),dosaturation = 0;,nonlinthreshold = 'none';,else dosaturation = 1;,end



% ----------------------------------------------------------------
% * get filtering matrix, if any
% ---------------------------------------------------------------- 
if isempty(HPlength),HPlength = 'none';,end

dofilter = 1;
disp(['	...filtering: HPlength = ' num2str(HPlength) ', LP Smoothing = ' num2str(LPsmooth)])

%use spm_filter
if strcmp(HPlength,'none') & LPsmooth										% LP only
	[S,KL] = use_spm_filter(TR,numsamps,'hrf','none',[]);
elseif LPsmooth																% HP and LP
	[S,KL,KH] = use_spm_filter(TR,numsamps,'hrf','specify',HPlength);
elseif strcmp(HPlength,'none') & ~LPsmooth									% neither
	S = 0;												% for 'no filter' in designsim - only thing diff from optimizeGA.m
	dofilter = 0;
else [S,KL,KH] = use_spm_filter(TR,numsamps,'none','specify',HPlength);		% HP only
end


switch SIM.type

case 'ga'
% ****************************************************************
% * for genetic algorithm simulations
% ****************************************************************


% ----------------------------------------------------------------
% * simulation - setup ga simulation with different noise vectors
% ----------------------------------------------------------------	
stimList = SIM.stimlist;

%gst commented out line below
%if ISI > TR,stimList = stimList(1:numsamps+1,1);,end   		% cut off to save time if it's too long 
if TR > ISI,stimList = stimList(1:numsamps+1,1);,end   		% gst

stimList = sampleInSeconds(stimList,ISI,.1);				% hi-res sample at .1 s, matches HRF
dmodel = getPredictors(stimList, HRF);
dmodel = resample(dmodel,1,TR*10);
if dosaturation,
	dmodel = modelSaturation(dmodel,nonlinthreshold);
end

%%%%%%%%%%%% gst added lines below: to find Efficiency %%%%%%%%%%%%%%
        HRF = spm_hrf(.1);						% SPM HRF sampled at .1 s
        HRF = HRF/ max(HRF);
        nonlinthreshold=[]; S=[];
                
%        dmodel = designvector2model(stimList,ISI,HRF,TR,numsamps,nonlinthreshold,S);
        dmodel = rna2model(stimList,ISI,HRF,TR,numsamps,nonlinthreshold,S);
        
        xtxitx=pinv(dmodel);
        contrastweights = 1; %ones(1,size(xtxitx,1)-1);
        svi = eye(size(xtxitx,2));

        [RESULTS.eff_fitness, RESULTS.eff]=calcEfficiency(contrastweights,c,xtxitx,svi);

        dmodel=dmodel(:,1:size(dmodel,2)-1);  %gst a hack: rna2model (called above) adds a column of 1's, 
                                            %but designsim (called below) 
                                            %expects that no such column has be added yet.

%%% hrf eff below 
    delta = [];
    for jj = 1:max(stimList(:,1))
        delta(:,jj) = (stimList == jj);
    end
    [model] = tor_make_deconv_mtx2(delta,round(HRF_estlen/ TR), TR / ISI); %from modelDiagnostics2
     xtxitx = pinv(model);   % a-optimality   % inv(X'S'SX)*(SX)'; pseudoinv of (S*X)
     svi = eye(size(xtxitx,2));
     [RESULTS.hrf_eff_fitness,RESULTS.hrf_eff] = calcEfficiency([],[],xtxitx,svi,0);                                            
                                            
%%%%%%%%%%%%

if size(dmodel,1) > numsamps
	disp('	...truncating dmodel to number of samples')
	dmodel = dmodel(1:numsamps,:);
elseif size(dmodel,1) < numsamps
	warning('	...stimlist for GA is too short for this ISI!  Results may not be accurate for this ISI!')
	dmodel = [dmodel; zeros(numsamps - size(dmodel,1),size(dmodel,2))];
end

% ----------------------------------------------------------------
% * simulation - run ga simulation with different noise vectors
% ----------------------------------------------------------------

disp(['	*> simulating GA results with noisy  at variance ' num2str(noise_var)])     

[dummy,RESULTS.se,RESULTS.t] = designsim('myxc',dmodel,HRF,ISI,TR,noise_var,c,beta,S,xc,niterations, NULL_flag); %gst
%[dummy,RESULTS.se,RESULTS.t] = designsim('myxc',dmodel,HRF,ISI,TR,noise_var,c,beta,S,xc,niterations);

	% Note:
	% no filter:[dummy,SIM.ga.HPse,SIM.ga.HPt] = designsim('myxc',dmodel,HRF,ISI,TR,noise_var,c,beta,0,xc,niterations);	
	% (S = 0 if no filter)






case 'rnd'
% ****************************************************************
% * for random design simulations
% ****************************************************************


disp(['	*> simulating random designs with ' num2str(niterations) ' iterations...']);

% ----------------------------------------------------------------
% * jitter setup
% ----------------------------------------------------------------  

dojitter = 0;
if size(freqConditions,2) == size(conditions,2) + 1
    conditions = [conditions max(conditions)+1];   % insert 0 as last condition if jitter
    dojitter = 1;
elseif size(freqConditions,2) == size(conditions,2)
else error('Frequency and condition vectors must be same length, or freq must be 1 longer for jitter')
end


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


% ----------------------------------------------------------------
% * create unsorted list
% ----------------------------------------------------------------

% this is overwritten if trans2switch is used.
unsortedList = [];
for i = 1:size(conditions,2)									% unsorted list of all stims in proper frequencies
   startIndex = size(unsortedList,1)+1;
%   unsortedList(startIndex:(startIndex + 2 * numStimEachCond(i)-1),1) =   conditions(i);	%gst
   unsortedList(startIndex:(startIndex + numStimEachCond(i)-1),1) = conditions(i);  %gst
end

% ----------------------------------------------------------------
% * simulation - setup random simulation - trans2switch
% ----------------------------------------------------------------  

if trans2switch
    if size(conditions,2) > 4
	NEWfreqConditions = [.5 .5 freqConditions(5)]; NEWconditions = [1 2 5];
    else NEWfreqConditions = [.5 .5];,NEWconditions = [1 2];
    end
	NEWnumStimEachCond = ceil((numStim-numRestStim) * NEWfreqConditions);			            % row vector of stim in each cond

	disp(['		...USING TRANS2SWITCH > num stim each cond = ' num2str(NEWnumStimEachCond)])

	for i = 1:size(NEWfreqConditions,2)									% unsorted list of all stims in proper frequencies
   		startIndex = size(unsortedList,1)+1;
   		unsortedList(startIndex:(startIndex + 2 * NEWnumStimEachCond(i)-1),1) = NEWconditions(i);
	end
end


% ----------------------------------------------------------------
% * simulation - setup random simulation - trans2block
% ---------------------------------------------------------------- 

if trans2block
	disp(['		...USING TRANS2BLOCK >'])
        freqConditions = [freqConditions./2 freqConditions./2];
        disp(['     	freqConditions = ' num2str(freqConditions)])
   
        if dojitter
           maxcond = max(conditions)-1;
           conditions = [conditions(1:maxcond) maxcond+1:2*length(conditions)-1];
           disp(['     CONDITION ' num2str(conditions(end)) ' is blank jitter.'])
	else
           maxcond = max(conditions);
           conditions = [conditions maxcond+1:2*length(conditions)];
        end
        disp(['     	conditions = ' num2str(conditions)])
	
	% ----------------------------------------------------------------
	% * setup special first trial
	% ----------------------------------------------------------------
	if isempty(dofirst) 
		dofirst = 0;
	else 
		disp(['		...block design: adding special predictor for first trial in block.'])
		try
			% this is really kludgy, but should give a rough estimate of how many first trials to expect.
			numfirst = numStim ./ mean(restevery);
			proportionfirst = numfirst./numStim;
		catch
			error('restevery not found...must use rest periods with transform2block!')
		end
		freqConditions = freqConditions - proportionfirst./size(freqConditions,2);
		freqConditions = [freqConditions proportionfirst];
		disp('		...first trial special condition - adjusted expected frequencies:')
	end
end



% ----------------------------------------------------------------
% * create random listMatrix and restMatrix
% ----------------------------------------------------------------

clear restMatrix, clear listMatrix
if dorests,
	% make rest matrix
	atleast = min(restevery);
	restPeriods = numStim / atleast + 5;	% make this many rows in restMatrix
	unsortedrestList = [];
	for i = 1:size(restevery,2):restPeriods
		unsortedrestList = [unsortedrestList;restevery'];
	end
else
	restMatrix = [];
end

for i = 1:niterations
	if dorests, restMatrix(:,i) = getRandom(unsortedrestList);,end
	listMatrix(:,i) = getRandom(unsortedList);
end

if dorests,
	listMatrix = insert_rests(listMatrix,restMatrix,restlength,numStim);
end		
		
% insert rests before transforming
for i = 1:niterations
	if trans2switch, 
		listMatrix(:,i) = transform2switches(listMatrix(:,i));
	end
	if trans2block, 	
		if ~dorests,error('trans2block is not (yet) compatible with using no rests!'),end 
             	listMatrix(:,i) = transform2block(listMatrix(:,i),restMatrix(:,i),restlength,dofirst);
      	end
end



		

% ----------------------------------------------------------------
% * simulation - setup random simulation - criterion measures
% ----------------------------------------------------------------

              if docriterion
                disp('   ...Making sure all lists meet criteria...')

                % criterion measures
	        % ====================================================================
	        for i = 1:niterations
                	stimList = listMatrix(:,i);
			restList = restMatrix(:,i);

                	go = 0;tryindex = 1;

                	while ~go 
                  		[maxNumStim,maxrest] = getMaxInARow(stimList);
 
                  		% ====== Frequency Deviation =======================================
                  		for j = 1:size(conditions,2) 
         				freqMat(j) = sum(stimList == conditions(j)) / size(stimList,1);
      		  		end
                  		freqMat = freqMat ./ sum(freqMat);   % adjust for probe ''0''s
		  		maxFreqDev = max(abs(freqConditions(1:size(conditions,2)) - freqMat));

                  		% ====== Counterbalancing Deviation =======================================
                  		[cBal,dummy,maxDev] = getCounterBal(stimList, maxOrder,conditions,freqConditions);

                  		% ====== Test list for criteria satisfaction ==============================
		  		if maxNumStim < NumStimthresh & maxrest < maxrestthresh & maxDev < maxCbalDevthresh & maxFreqDev < maxFreqDevthresh
		     		go = 1;
		  		else
                     		restList        = getRandom(unsortedrestList);
			 		stimList		 = getRandom(unsortedList);
			 		stimList        = insert_rests(stimList,restList,restlength,numStim);
                     		if trans2switch,
                        		stimList = transform2switches(stimList);
                     		end
		     		if trans2block, 
             					stimList = transform2block(stimList,restList,restlength,dofirst);
      			 		end

                     		tryindex = tryindex+1;
                  		end

		     	if mod(tryindex,100) == 0,disp(['Tried 100 times on model ' num2str(i)])
                          disp(['Avg. max stim in a row = ' num2str(maxNumStim)])
                          disp(['Avg. deviation from freqs = ' num2str(maxFreqDev)])
                          disp(['Avg. deviation from counterbal = ' num2str(maxDev)])   
                          disp(['Avg. frequencies are = ' num2str(freqMat)])
                          disp(['Max restcond is ' num2str(maxrest)])
					conditions
					freqConditions
					freqMat
                     	end           
				if mod(tryindex,500) == 0
					disp('Cannot meet criteria after 500 tries.')
					conditions
					freqConditions
					freqMat
					cBal
					error('Exiting now.')
				end
                	end  % end while

                	restMatrix(:,i) = restList;
			listMatrix(:,i) = stimList;
                	if mod(i,50) == 0, fprintf('%d ',i),end
                end  % end for ... iterations
	     end  % if docriterion

% ----------------------------------------------------------------
% * simulation - setup random simulation - make models
% ----------------------------------------------------------------


for i = 1:niterations

%%%%%%%%%% gst added lines below %%%%%%%%
	if i==1  %%%% calculate hrf_eff  (before do sampling)

                stimList = listMatrix(:,i);

                delta = [];
                for jj = 1:max(stimList(:,1))
                    delta(:,jj) = (stimList == jj);
                end
	%             
	%             % how much HRF estimation we do (up to 30 s) matters
	%             % we want at least as long as the expected actual response
    [model] = tor_make_deconv_mtx2(delta,round(HRF_estlen / TR), TR / ISI); %from modelDiagnostics2
%    [model] = tor_make_deconv_mtx2(delta,round(30 / TR),TR / ISI);  % hard-coded for 30 s of HRF estimation
	
	%             if ~isempty(S), model = S * model;,end
	%             
	%             if dflag
	%                 xtxitx = model'*svi*model;  % d-optimality
	%                 fitnessMatrix(3,z) = calcEfficiency([],[],xtxitx,[],dflag);
	%             else
	
                  xtxitx = pinv(model);   % a-optimality   % inv(X'S'SX)*(SX)'; pseudoinv of (S*X)
                  svi = eye(size(xtxitx,2));
                 [RESULTS.hrf_eff_fitness,RESULTS.hrf_eff] = calcEfficiency([],[],xtxitx,svi,0);
        
	end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Test Random models
%gst commented out lines below

% 	  	% Prepare Stimulus List
%       	        stimList = listMatrix(:,i);
% 		SIM.stimList = stimList;
% 		stimList = sampleInSeconds(stimList,ISI,.1);
% 		dmodel = getPredictors(stimList, HRF);
% 		dmodel = resample(dmodel,1,TR*10);
% 		if dosaturation, dmodel = modelSaturation(dmodel,nonlinthreshold);,end
% 
% 		if size(dmodel,1) > numsamps
% 			dmodel = dmodel(1:numsamps,:);
% 		elseif size(dmodel,1) < numsamps
% 			dmodel = [dmodel; zeros(numsamps - size(dmodel,1),size(dmodel,2))];
% 		end
		
		
%%%%%%%%%%%%%gst added lines below: %%%%%%%%%%%%%%
        HRF = spm_hrf(.1);						% SPM HRF sampled at .1 s
        HRF = HRF/ max(HRF);
        nonlinthreshold=[]; S=[];

      	        stimList = listMatrix(:,i); %gst
                
%        dmodel = designvector2model(stimList,ISI,HRF,TR,numsamps,nonlinthreshold,S);
        dmodel = rna2model(stimList,ISI,HRF,TR,numsamps,nonlinthreshold,S);
        
        xtxitx=pinv(dmodel);
        contrastweights = 1; %ones(1,size(xtxitx,1)-1);
        svi = eye(size(xtxitx,2));

        [RESULTS.eff_fitness(i,:), RESULTS.eff(i,:)]=calcEfficiency(contrastweights,c,xtxitx,svi);

        dmodel=dmodel(:,1:size(dmodel,2)-1);  %gst a hack: rna2model (called above) adds a column of 1's, 
                                            %but designsim (called below) 
                                            %expects that no such column has be added yet.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ----------------------------------------------------------------
% * simulation - run random simulation
% ----------------------------------------------------------------

	if i == 1   
		[dummy,RESULTS.se(:,i),RESULTS.t(:,i),GLM] = ... 		%added 4th output to plot
		designsim('myxc',dmodel,HRF,ISI,TR,noise_var,c,beta,S,xc, 1, NULL_flag);				
%        designsim('myxc',dmodel,HRF,ISI,TR,noise_var,c,beta,S,xc);
        subplot(4,1,4); plot(dmodel(:,1));title('First predictor of model')
		drawnow
	else
		[dummy,RESULTS.se(:,i),RESULTS.t(:,i)] = ...
		designsim('myxc',dmodel,HRF,ISI,TR,noise_var,c,beta,S,xc, 1, NULL_flag);
%		designsim('myxc',dmodel,HRF,ISI,TR,noise_var,c,beta,S,xc);
    end

	RESULTS.listMatrix = listMatrix;


    
end	% loop through iterations

end	% end switch simulation type

RESULTS.meant = mean(RESULTS.t);
RESULTS.meanse = mean(RESULTS.se);

% ----------------------------------------------------------------
% * power calculation for all sims
% ----------------------------------------------------------------

RESULTS.powerTcrit = tinv(1-.05/powerTvoxels,powerTdf);
RESULTS.power = sum(abs(RESULTS.t) > RESULTS.powerTcrit) / length(RESULTS.t);


return

