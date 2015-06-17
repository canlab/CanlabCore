% User-Specified Parameters
% ===============================================================

% Conditions and stimuli
% ---------------------------------------------------------------
GA.conditions = [1];
GA.freqConditions = [.5];
GA.scanLength = 480;  
GA.ISI = 2;  
GA.TR = 2; 

% Genetic algorithm parameters
% ---------------------------------------------------------------
nmodels = 1; 
GA.cbalColinPowerWeights = [1 1 1 0];	% 1 = cbal, 2 = eff, 3 = hrf shape, 4 = freq
GA.numGenerations = 2; 
GA.sizeGenerations = 20;  
GA.maxTime = 600;						% max time to run in s, or Inf
GA.alph = 2.1;  
GA.plotFlag = 0; 

% Filtering, counterbalancing, and design tolerance
% ---------------------------------------------------------------
GA.lowerLimit = [];  
GA.HPlength = []; 
GA.LPsmooth = []; 
GA.maxOrder = 1; 
GA.NumStimthresh = []; 
GA.maxCbalDevthresh = []; 	% counterbalancing deviation
GA.maxFreqDevthresh = []; 

% Contrast setup
% ---------------------------------------------------------------
GA.contrasts = [1]; 
GA.contrastweights = [1];	% or predictor weights, if no contrasts


% Autocorrelation and special options
% ---------------------------------------------------------------
AutocorrelationFileName = 'myscannerxc';
GA.restlength = []; 
GA.restevery = []; 
GA.trans2switch = 0; 
GA.trans2block = 0; 
GA.dofirst = 0; 
GA.nonlinthreshold = []; 

eval(['load ' AutocorrelationFileName]);  
GA.xc = myscannerxc;



% =============================================
% * vary by parameter - set this in this script
% * if running multiple models, nmodels > 1
% ---------------------------------------------

varyparam = 0;
fieldname = 'alph';		% or 'freqConditions(end), etc.
incrementby = 0.3;
incrementevery = 5;		% every n models
contrastsofinterest = 1;% contrasts or predictors, if no cons specified

% =============================================


if varyparam
	eval(['paramvalue = GA.' fieldname ';'])
	disp(' ');disp('* *********************************************************************')
	disp('*  You have selected to vary a parameter across models in ga_gui_script')
	disp(['*  GA.' fieldname ' starting at ' num2str(paramvalue)]) 
	disp('* *********************************************************************')
end


for nm = 1:nmodels

   eval(['diary model' num2str(nm) '.log'])
   disp(['Starting Model ' num2str(nm)])
    % freqConditions

   M = optimize_rand_search(GA);


   Models(:,nm) = M.stimlist;
   MM{nm} = M;
   save GAworkspace
   
   if varyparam
   	if isfield(M,'consebeta'),
		Fit(1,nm) = mean(M.consebeta(contrastsofinterest));
	else 
		Fit(1,nm) = mean(M.sebeta(contrastsofinterest));
	end
   	eval(['Fit(2,nm) = GA.' fieldname ';'])
   	FC(nm,:) = GA.freqConditions;
   	ALPH(nm,:) = GA.alph;
	str = ['if mod(nm,' num2str(incrementevery) ') == 0,GA.' fieldname ' = GA.' fieldname '+' num2str(incrementby) ';,end']
   	eval(str);
	eval(['paramvalue = GA.' fieldname ';'])
	disp('* *********************************************************************')
	disp('*  Varying parameter option on.')
	disp(['*  GA.' fieldname ' is now ' num2str(paramvalue)]) 
	disp('* *********************************************************************')
   end

   eval(['save model' num2str(nm) '_' M.date ' M'])
   save GAworkspace
   diary off
   fprintf('\n')
end



