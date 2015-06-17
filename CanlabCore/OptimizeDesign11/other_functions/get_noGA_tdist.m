% defaults

conditions = [1 2 3 4];				% vector of condition #s (1 for each trial type)

freqConditions = [.35 .35 .15 .15];	% relative freq of conditions - should sum to 1

scanLength = 360;					% length in s of the scan,must be multiple of ISI! 

numCondNoInt = 0;					% 1 = insert a condition of no interest for jitter

ISI = 1200;							% inter-stimulus interval

TR = 2;								% TR of your experiment

cbalColinPowerWeights = [0 5 5]; 	% rel. importance of three factors, scale from 0 to 1+

respjitter = 284;					% st. dev of msec value to jitter estimated response + or - by

numGenerations = 10;				% number of generations for genetic algorithm

sizeGenerations = 26;				% number of organisms per generation, must be even number!

lowerLimit = 0.05;					% lower limit of acceptable power frequency, in Hz

maxOrder = 1;						% maximum order of trial sequence counterbalancing to consider

plotFlag = 0;						% flag to subfunctions to plot output; turn on after design is chosen.



H = findobj('Tag', 'E1');                                                                                                                                 

conditions = str2num(get(H, 'String'));     

H = findobj('Tag', 'E2');                                                                                                                                 

freqConditions = str2num(get(H, 'String')); 

H = findobj('Tag', 'E3');                                                                                                                                 

scanLength = str2num(get(H, 'String'));  

H = findobj('Tag', 'E4');                                                                                                                                 

ISI = str2num(get(H, 'String'));  

H = findobj('Tag', 'E5');                                                                                                                                 

TR = str2num(get(H, 'String')); 

H = findobj('Tag', 'E6');                                                                                                                                 

cbalColinPowerWeights = str2num(get(H, 'String'));

H = findobj('Tag', 'E7');                                                                                                                                 

numGenerations = str2num(get(H, 'String')); 

H = findobj('Tag', 'E8');                                                                                                                                 

sizeGenerations = str2num(get(H, 'String'));  

H = findobj('Tag', 'E9');                                                                                                                                 

lowerLimit = str2num(get(H, 'String'));  

H = findobj('Tag', 'E10');                                                                                                                                 

maxOrder = str2num(get(H, 'String')); 

H = findobj('Tag', 'E11');                                                                                                                                 

plotFlag = str2num(get(H, 'String'));  

H = findobj('Tag', 'E12');                                                                                                                                 

noise_var = str2num(get(H, 'String')); 

H = findobj('Tag', 'E13');                                                                                                                                 

c = str2num(get(H, 'String'));  

H = findobj('Tag', 'E14');                                                                                                                                 

niterations = str2num(get(H, 'String')); 

H = findobj('Tag', 'E15');                                                                                                                                 

respMean = str2num(get(H, 'String')); 

H = findobj('Tag', 'E16');                                                                                                                                 

respSD = str2num(get(H, 'String')); 

H = findobj('Tag', 'C1');                                                                                                                                 

addResponseVariability = get(H, 'Value'); 



HRF = spm_hrf(ISI/1000);

HRF = HRF/ max(HRF);

clear nullt;



% generate an unsorted list of conditions to use for all random generations

numStim = scanLength / (ISI/1000);

numStimEachCond = numStim * freqConditions;  			% row vector of stim in each cond

unsortedList = [];

for i = 1:size(conditions,2)									% unsorted list of all stims in proper frequencies

   startIndex = size(unsortedList,1)+1;

   unsortedList(startIndex:(startIndex + numStimEachCond(i)-1),1) = conditions(i);	

end



disp(['Starting design simulation of ' num2str(niterations) ' iterations...']);

if addResponseVariability == 1

	for i = 1:niterations

		stimList = getRandom(unsortedList); 

 	   MODEL = getRandRespPredictors(stimList,conditions, ISI, respMean, respSD);

       [nullt(i), nullr(i)] = designsim(MODEL,ISI,noise_var,c,1);

	   if mod(i,50) == 0, disp(['done ' num2str(i)]),end;

    end

else

   for i = 1:niterations

      stimList = getRandom(unsortedList); 

 	   MODEL = getPredictors(stimList,conditions, HRF);

      [nullt(i), nullr(i)] = designsim(MODEL,ISI,noise_var,c,1);

	  if mod(i,50) == 0, disp(['done ' num2str(i)]),end;

   end

end



figure;

hist(nullt,niterations/25); 

title(['Distribution of t-values with random designs using your params']);

figure;

hist(nullr,round(niterations/25)); 

title(['Distribution of predictor r-values for contrast [' num2str(c) '] using your design']);   

