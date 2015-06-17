% designsim_gui_script
% Last modified 3/21/01 by Tor Wager		fix HRF bug / various
name = input('Enter output file name: ');

% sets defaults, reads GUI values from optimizeGUI, and runs simulation

disp('Designsim_gui_script')
disp('===============================')
disp('  ...Setting up variables from GUI.')
% returns vector of t values as t in workspace.

if ~exist('M') | ~isstruct(M),disp('No GA results (M struct) detected - using only random designs.'),end

clear SIM
clear H;

format compact

% defaults
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



% read values from GUI for GA
H = findobj('Tag', 'E1');                                                                                                                                 
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

HRF = spm_hrf(.1);
HRF = HRF/ max(HRF);

for myISI = 1:16
ISI = myISI;
disp(['Starting simulation with ISI ' num2str(myISI)]);

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

disp('  ...Generating design list.')
% for building RANDOM designs to test against GA design
fitnessMatrix = zeros(4,sizeGenerations);
numStim = ceil(scanLength / (ISI));
numRestStim = (ceil(numStim/restevery) - 1) * restlength;
numStimEachCond = ceil((numStim-numRestStim) * freqConditions);			            % row vector of stim in each cond
numsamps = ceil(numStim*ISI/TR);

SIM.numsamps = numsamps;
if exist('M') & isstruct(M), SIM.numsamps = size(M.modelatTR,1);,end
SIM.numStimEachCond = numStimEachCond;

unsortedList = [];
for i = 1:size(conditions,2)									% unsorted list of all stims in proper frequencies
   startIndex = size(unsortedList,1)+1;
   unsortedList(startIndex:(startIndex + numStimEachCond(i)-1),1) = conditions(i);	
end



if addResponseVariability == 1
    % THIS ISN'T FINISHED, AND I'M NOT SURE I SEE THE POINT, REALLY.

else
    disp(['  ...loading smoothing matrix with hrf LP and HP filter length of ' num2str(HPlength)])
    [hpS] = use_spm_filter(TR,SIM.numsamps,'none','specify',HPlength); 
	[fullS] = use_spm_filter(TR,SIM.numsamps,'hrf','specify',HPlength);  
	
	if exist('M') & isstruct(M),
	        stimList = M.stimlist(1:numsamps+1,1);   % cut off to make more efficient
                stimList = sampleInSeconds(stimList,ISI,.1);
		dmodel = getPredictors(stimList, HRF);
		dmodel = resample(dmodel,1,TR*10);
		% dmodel = dmodel / max(max(dmodel));                     % scale so that max = 1;
		if size(dmodel,1) > SIM.numsamps
			dmodel = dmodel(1:SIM.numsamps,:);
		elseif size(dmodel,1) < SIM.numsamps
		        warning('stimlist for GA is too short for this ISI!  Results may not be accurate for this ISI!')
	                dmodel = [dmodel; zeros(SIM.numsamps - size(dmodel,1),size(dmodel,2))];
	        end

  		disp(['  ...simulating GA results with noisy data at variance ' num2str(noise_var)])
                SIM        
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
whos unsortedList

	% to do:
% cut off model at numscans  done, i think
% rebuild ga at different isis  done
% cut off stimlist beforehand cause it's too long.   done
% warning if stimlist in GA is too short for selected isi - then pad with zeros.  done


	% computing random SIM model
	% ============================================================================================
	
	for i = 1:niterations
      % Test Random models

	  	% Prepare Stimulus List
      	stimList = getRandom(unsortedList);
      	% ====== insert rests, if needed, and trim list ===================
        stimList = insert_rests(stimList,restevery,restlength);
        if size(stimList,1) < numStim, stimList(end+1:numStim,:) = 0;
        elseif size(stimList,1) > numStim,stimList = stimList(1:numStim,:);
        end
		stimList = sampleInSeconds(stimList,ISI,.1);
		dmodel = getPredictors(stimList, HRF);
		dmodel = resample(dmodel,1,TR*10);
		% dmodel = dmodel / max(max(dmodel));                     % scale so that max = 1;
		if size(dmodel,1) > SIM.numsamps
			dmodel = dmodel(1:SIM.numsamps,:);
		elseif size(dmodel,1) < SIM.numsamps
			dmodel = [dmodel; zeros(SIM.numsamps - size(dmodel,1),size(dmodel,2))];
		end
		
		% this gets the average over [numnoise] simulations with different noise vectors
	  % =====================================================================================================
      % no smoothing or filtering
      [dummy,SIM.rnd.se(:,i),SIM.rnd.t(:,i)] = designsim('myxc',dmodel,HRF,ISI,TR,noise_var,c,beta,0,xc);
	  % high-pass only
	  [dummy, SIM.rnd.HPse(:,i),SIM.rnd.HPt(:,i)] = designsim('myxc',dmodel,HRF,ISI,TR,noise_var,c,beta,hpS,xc);
      % now with smoothing
      [dummy, SIM.rnd.HLse(:,i),SIM.rnd.HLt(:,i)] = designsim('myxc',dmodel,HRF,ISI,TR,noise_var,c,beta,fullS,xc);
      fprintf('.')
	  if mod(i,30) == 0, fprintf(' %d ',i),end
	  if mod(i,90) == 0, fprintf('\n'),end
	  
	  % this gets the average over [numnoise] simulations with different noise vectors
	  % =====================================================================================================
      % no smoothing or filtering
      %[SIM.rnd.t(:,i), SIM.rnd.se(:,i)] = designsim('myxc',dmodel,HRF,ISI,TR,noise_var,c,beta,0,xc,100);
	  % high-pass only
	  %[SIM.rnd.HPt(:,i), SIM.rnd.HPse(:,i)] = designsim('myxc',dmodel,HRF,ISI,TR,noise_var,c,beta,hpS,xc,100);
      % now with smoothing
      %[SIM.rnd.HLt(:,i), SIM.rnd.HLse(:,i)] = designsim('myxc',dmodel,HRF,ISI,TR,noise_var,c,beta,fullS,xc,100);
      %fprintf('.')
      %if mod(i,30) == 0, fprintf(' %d\n',i),end
      %if mod(i,50) == 0, disp(['done ' num2str(i)]),end;

	 
    end
      
end


      MyISIMean.ga.t(myISI) = mean(SIM.ga.t);
      MyISIMean.rnd.t(myISI) = mean(SIM.rnd.t);
      MyISIMean.ga.HPt(myISI) = mean(SIM.ga.HPt);
      MyISIMean.rnd.HPt(myISI) = mean(SIM.rnd.HPt);
      MyISIMean.ga.HLt(myISI) = mean(SIM.ga.HLt);
      MyISIMean.rnd.HLt(myISI) = mean(SIM.rnd.HLt);
      MyISIMean.ISI(myISI) = ISI;
      MyISIMean.ga.se(myISI) = std(SIM.ga.t)/sqrt(size(SIM.ga.t,2));
      MyISIMean.rnd.se(myISI) = std(SIM.rnd.t)/sqrt(size(SIM.rnd.t,2));
      MyISIMean.ga.HPse(myISI) = std(SIM.ga.HPt)/sqrt(size(SIM.ga.HPt,2));
      MyISIMean.rnd.HPse(myISI) = std(SIM.rnd.HPt)/sqrt(size(SIM.rnd.HPt,2));
      MyISIMean.ga.HLse(myISI) = std(SIM.ga.HLt)/sqrt(size(SIM.ga.HLt,2));
      MyISIMean.rnd.HLse(myISI) = std(SIM.rnd.HLt)/sqrt(size(SIM.rnd.HLt,2));

end % end loop through simulations
      eval(['save ' name ' MyISIMean'])






