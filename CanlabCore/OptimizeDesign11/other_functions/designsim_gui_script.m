% designsim_gui_script
% Last modified 5/20/01 by Tor Wager		add optional trans2switch from GUI.  still need block.
% also added power calculation.

% sets defaults, reads GUI values from optimizeGUI, and runs simulation

disp('Designsim_gui_script')
disp('===============================')
disp('  ...Setting up variables from GUI.')
% returns vector of t values as t in workspace.

if ~exist('M') | ~isstruct(M),disp('	...No GA results (M struct) detected - using only random designs.'),end

clear SIM
clear H;

% ----------------------------------------------------------------
% * setup defauls
% ----------------------------------------------------------------
% fixed parameters
% --------------------
powerTvoxels = 30000;
powerTdf = 19;
format compact

% parameters from GUI
% --------------------
nonlinthreshold = 4;
noise_var = 0.3;
beta = [.5 .5 0 0;1 1 0 0];	% 20
c = [1 1 -1 -1];					% 13
niterations = 500;
conditions = [1 2 3 4];
ISI = 1.2;
restevery = 9;						% 26
restlength = 2;					% 25
HPlength = 60;						% 23
xc = [1 0];							% 24
trans2switch = 0;

% ----------------------------------------------------------------
% * setup input from GUI
% ----------------------------------------------------------------
% read values from GUI for GA
H = findobj('Tag', 'E1');                                                                                                             if ~isempty(H)              
conditions = str2num(get(H(1), 'String'));   

H = findobj('Tag', 'E4');   
ISI = str2num(get(H(1), 'String')); 

H = findobj('Tag', 'E12');  
noise_var = str2num(get(H(1), 'String')); 

H = findobj('Tag', 'E13');
c = str2num(get(H(1), 'String'));  

H = findobj('Tag', 'E14');  
niterations = str2num(get(H(1), 'String')); 

H = findobj('Tag', 'E15'); 
respMean = str2num(get(H(1), 'String')); 

H = findobj('Tag', 'E16');
respSD = str2num(get(H(1), 'String')); 

H = findobj('Tag', 'C1');  
addResponseVariability = get(H(1), 'Value'); 

% additional values required for random sim.

H = findobj('Tag', 'E2');  
freqConditions = str2num(get(H(1), 'String')); 

H = findobj('Tag', 'E3');
scanLength = str2num(get(H(1), 'String'));   

H = findobj('Tag', 'E5'); 
TR = str2num(get(H(1), 'String')); 

H = findobj('Tag', 'E6');
cbalColinPowerWeights = str2num(get(H(1), 'String'));

H = findobj('Tag', 'E7'); 
numGenerations = str2num(get(H(1), 'String')); 

H = findobj('Tag', 'E8');
sizeGenerations = str2num(get(H(1), 'String'));  

H = findobj('Tag', 'E9'); 
lowerLimit = str2num(get(H(1), 'String'));  

H = findobj('Tag', 'E10');  
maxOrder = str2num(get(H(1), 'String')); 

H = findobj('Tag', 'E11');  
plotFlag = str2num(get(H(1), 'String'));  

H = findobj('Tag', 'E20');  
beta = str2num(get(H(1), 'String'));  

H = findobj('Tag', 'E23');  
HPlength = str2num(get(H(1), 'String'));  

H = findobj('Tag', 'E24');  
eval(['load ' get(H(1), 'String')]);  

H = findobj('Tag', 'E25');  
restlength = str2num(get(H(1), 'String')); 

H = findobj('Tag', 'E26');  
restevery = str2num(get(H(1), 'String')); 

H = findobj('Tag', 'E27');  
trans2switch = str2num(get(H(1), 'String')); 

H = findobj('Tag', 'E28');  
trans2block = str2num(get(H(1), 'String')); 

H = findobj('Tag', 'E30');  
NumStimthresh = str2num(get(H(1), 'String')); 
 
H = findobj('Tag', 'E31');  
maxCbalDevthresh = str2num(get(H(1), 'String')); 	% counterbalancing deviation

H = findobj('Tag', 'E32');  
maxFreqDevthresh = str2num(get(H(1), 'String')); 

H = findobj('Tag', 'E33');  
nonlinthreshold = str2num(get(H(1), 'String')); 

H = findobj('Tag', 'E34');  
dofirst = str2num(get(H(1), 'String')); 

end % ~isempty(H) : means found the GUI


% ----------------------------------------------------------------
% * setup lists, HRF, and other parameters
% ----------------------------------------------------------------

HRF = spm_hrf(.1);
HRF = HRF/ max(HRF);

SIM.beta = beta;
SIM.contrasts = c;
SIM.noise_var = noise_var;
SIM.xc = xc;
SIM.ISI = ISI;
SIM.HPlength = HPlength;
SIM.restevery = restevery;
SIM.restlen = restlength;
SIM.niterations = niterations;
SIM.noise_var = noise_var;
SIM.freqConditions = freqConditions;
 
SIM.NumStimthresh = NumStimthresh;   
SIM.maxCbalDevthresh = maxCbalDevthresh; 
SIM.maxFreqDevthresh = maxFreqDevthresh;
SIM.nonlinthreshold = nonlinthreshold;

disp('  ...Generating design list.')
% for building RANDOM designs to test against GA design
fitnessMatrix = zeros(4,sizeGenerations);
numStim = ceil(scanLength / (ISI));
numRestStim = (ceil(numStim/(mean(restevery)+restlength)) - 1) * restlength;
numStimEachCond = ceil((numStim-numRestStim) * freqConditions);			            % row vector of stim in each cond
numsamps = ceil(numStim*ISI/TR);

SIM.numsamps = numsamps;
if exist('M') & isstruct(M), SIM.numsamps = size(M.modelatTR,1);,end
SIM.numStimEachCond = numStimEachCond;

% this is overwritten if trans2switch is used.
unsortedList = [];
for i = 1:size(conditions,2)									% unsorted list of all stims in proper frequencies
   startIndex = size(unsortedList,1)+1;
   unsortedList(startIndex:(startIndex + 2 * numStimEachCond(i)-1),1) = conditions(i);	
end

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

maxrestthresh = 2;					% max number of rests in a row.


% ----------------------------------------------------------------
% * stochastic variability in response times - now not used
% ---------------------------------------------------------------- 
if addResponseVariability == 1

else

% ----------------------------------------------------------------
% * simulation - get smoothing matrix
% ---------------------------------------------------------------- 
	figure;
    disp(['  ...loading smoothing matrix with hrf LP and HP filter length of ' num2str(HPlength)])
    [hpS] = use_spm_filter(TR,SIM.numsamps,'none','specify',HPlength); 
	[fullS] = use_spm_filter(TR,SIM.numsamps,'hrf','specify',HPlength);  


% ----------------------------------------------------------------
% * simulation - setup ga simulation with different noise vectors
% ----------------------------------------------------------------	
	if exist('M') & isstruct(M),
	    if ISI > TR,stimList = M.stimlist(1:numsamps+1,1);,end   	% cut off to save time if it's too long
                %stimList = insert_rests(stimList,M.restlist,restlength,numStim);not necessary;done in GA
                stimList = sampleInSeconds(M.stimlist,ISI,.1);				% hi-res sample at .1 s, matches HRF
		dmodel = getPredictors(stimList, HRF);
		dmodel = resample(dmodel,1,TR*10);
		dmodel = modelSaturation(dmodel,nonlinthreshold);

		if size(dmodel,1) > SIM.numsamps
			disp('	...truncating dmodel to number of samples')
			dmodel = dmodel(1:SIM.numsamps,:);
		elseif size(dmodel,1) < SIM.numsamps
		        warning('	...stimlist for GA is too short for this ISI!  Results may not be accurate for this ISI!')
	                dmodel = [dmodel; zeros(SIM.numsamps - size(dmodel,1),size(dmodel,2))];
	        end

% ----------------------------------------------------------------
% * simulation - run ga simulation with different noise vectors
% ----------------------------------------------------------------

  		disp(['  ...simulating GA results with noisy data at variance ' num2str(noise_var)])     
	        % Test idealized GA - no response variability
		fprintf('		...original	model')
    	[dummy,SIM.ga.se,SIM.ga.t] = designsim('myxc',dmodel,HRF,ISI,TR,noise_var,c,beta,0,xc,niterations);
 		fprintf('...high-pass filtered')
    	[dummy,SIM.ga.HPse,SIM.ga.HPt] = designsim('myxc',dmodel,HRF,ISI,TR,noise_var,c,beta,hpS,xc,niterations);	
		fprintf('...smoothed\n')   
		[dummy,SIM.ga.HLse,SIM.ga.HLt] = designsim('myxc',dmodel,HRF,ISI,TR,noise_var,c,beta,fullS,xc,niterations);
	end

    disp(['  ...simulating random designs with ' num2str(niterations) ' iterations...']);
 
SIM

% ----------------------------------------------------------------
% * simulation - setup random simulation - trans2switch
% ----------------------------------------------------------------  

% This is part of the switching hack. uses 2 conditions, then figures out 4 columns from that.
% 5 is the rest condition.
if trans2switch
    if size(conditions,2) > 4
	NEWfreqConditions = [.5 .5 freqConditions(5)]; NEWconditions = [1 2 5];
    else NEWfreqConditions = [.5 .5];,NEWconditions = [1 2];
    end
	NEWnumStimEachCond = ceil((numStim-numRestStim) * NEWfreqConditions);			            % row vector of stim in each cond

	disp(['	...USING TRANS2SWITCH > num stim each cond = ' num2str(NEWnumStimEachCond)])

	for i = 1:size(NEWfreqConditions,2)									% unsorted list of all stims in proper frequencies
   		startIndex = size(unsortedList,1)+1;
   		unsortedList(startIndex:(startIndex + 2 * NEWnumStimEachCond(i)-1),1) = NEWconditions(i);
	end
end


% ----------------------------------------------------------------
% * simulation - setup random simulation - trans2block
% ---------------------------------------------------------------- 
if trans2block
	disp(['	...USING TRANS2BLOCK >'])
        freqConditions = [freqConditions./2 freqConditions./2];
        disp(['     freqConditions = ' num2str(freqConditions)])
   
        if dojitter
           maxcond = max(conditions)-1;
           conditions = [conditions(1:maxcond) maxcond+1:2*length(conditions)-1];
           disp(['     CONDITION ' num2str(conditions(end)) ' is blank jitter.'])
	else
           maxcond = max(conditions);
           conditions = [conditions maxcond+1:2*length(conditions)];
        end
        disp(['     conditions = ' num2str(conditions)])
	
	% ----------------------------------------------------------------
	% * setup special first trial
	% ----------------------------------------------------------------
	if isempty(dofirst) 
		dofirst = 0;
	else 
		disp(['	...block design: adding special predictor for first trial in block.'])
		try
			% this is really kludgy, but should give a rough estimate of how many first trials to expect.
			numfirst = numStim ./ mean(restevery);
			proportionfirst = numfirst./numStim;
		catch
			error('restevery not found...must use rest periods with transform2block!')
		end
		freqConditions = freqConditions - proportionfirst./size(freqConditions,2);
		freqConditions = [freqConditions proportionfirst];
		disp('	...first trial special condition - adjusted expected frequencies:')
		freqConditions
	end
end


% ----------------------------------------------------------------
% * simulation - setup random simulation - list and rest matrix
% ----------------------------------------------------------------
		disp('	...Setting up lists.')
                clear restMatrix, clear listMatrix
		% make rest matrix
		atleast = min(restevery);
		restPeriods = numStim / atleast + 5;	% make this many rows in restMatrix
		unsortedrestList = [];
		for i = 1:size(restevery,2):restPeriods
			unsortedrestList = [unsortedrestList;restevery'];
		end

		for i = 1:niterations
			restMatrix(:,i) = getRandom(unsortedrestList);
			listMatrix(:,i) = getRandom(unsortedList);
		end

		listMatrix = insert_rests(listMatrix,restMatrix,restlength,numStim);
		
		% insert rests before transforming
		for i = 1:niterations
			if trans2switch, 
				listMatrix(:,i) = transform2switches(listMatrix(:,i));
		    end
			if trans2block, 
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
      	% Test Random models

	  	% Prepare Stimulus List
      	        stimList = listMatrix(:,i);
		SIM.stimList = stimList;
		stimList = sampleInSeconds(stimList,ISI,.1);
		dmodel = getPredictors(stimList, HRF);
		dmodel = resample(dmodel,1,TR*10);
		dmodel = modelSaturation(dmodel,nonlinthreshold);

		if size(dmodel,1) > SIM.numsamps
			dmodel = dmodel(1:SIM.numsamps,:);
		elseif size(dmodel,1) < SIM.numsamps
			dmodel = [dmodel; zeros(SIM.numsamps - size(dmodel,1),size(dmodel,2))];
		end
		
		
% ----------------------------------------------------------------
% * simulation - run random simulation - no smoothing
% ----------------------------------------------------------------
% this gets the average over [numnoise] simulations with different noise vectors
   		if i == 1   
	  		[dummy,SIM.rnd.se(:,i),SIM.rnd.t(:,i),GLM] = ... 
			designsim('myxc',dmodel,HRF,ISI,TR,noise_var,c,beta,0,xc);	%added 4th output to plot
			set(gcf,'Position',[88         514        1185         420])	;drawnow
			subplot(4,1,4); plot(dmodel(:,1));title('First predictor of model')
			drawnow
		else
			[dummy,SIM.rnd.se(:,i),SIM.rnd.t(:,i)] = ...
			designsim('myxc',dmodel,HRF,ISI,TR,noise_var,c,beta,0,xc);
		end


% ----------------------------------------------------------------
% * simulation - run random simulation - HP filter
% ----------------------------------------------------------------
	% high-pass only
	[dummy, SIM.rnd.HPse(:,i),SIM.rnd.HPt(:,i)] = ...
	designsim('myxc',dmodel,HRF,ISI,TR,noise_var,c,beta,hpS,xc);


% ----------------------------------------------------------------
% * simulation - run random simulation - high and low filter
% ----------------------------------------------------------------

      	% now with smoothing
      	[dummy, SIM.rnd.HLse(:,i),SIM.rnd.HLt(:,i)] = designsim('myxc',dmodel,HRF,ISI,TR,noise_var,c,beta,fullS,xc);
      	fprintf('.')
	  if mod(i,30) == 0, fprintf(' %d ',i),end
	  if mod(i,90) == 0, fprintf('\n'),end
	
    end
	% for niterations
  
end
		% simulation with fixed response times.

% ----------------------------------------------------------------
% * power calculation for all sims
% ----------------------------------------------------------------

powerTcrit = tinv(1-.05/powerTvoxels,powerTdf);
for i = {'SIM.ga.t' 'SIM.ga.HPt' 'SIM.ga.HLt' 'SIM.rnd.t' 'SIM.rnd.HPt' 'SIM.rnd.HLt'}
	eval([i{1} 'power = sum(abs(' i{1} ') > powerTcrit) / length(' i{1} ');'])
end

	try
		plotSim(SIM)
	catch
		disp('Error plotting. command: plotSim(SIM)')
	end


save latestSIM SIM




