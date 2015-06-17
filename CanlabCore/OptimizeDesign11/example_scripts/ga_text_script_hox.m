clear Models, clear MM

% User-Specified Parameters
% ===============================================================

% Conditions and stimuli
% ---------------------------------------------------------------
GA.conditions = [1 2];
GA.freqConditions = [.5 .5];
GA.scanLength = 480;  
GA.ISI = 8;  
GA.TR = 2; 
GA.doepochs = 0;        % build epochs instead of events

% Hox optimization parameters
% A hox gene is a master gene.  These numbers control the stimulus
% parameters, which can be optimized within the GA rather than
% pre-specified.
% Hox elements in this GA go in the stimlist as the first several elements.
% hox sequence currently codes, in number order: 
% ISI TR cuelen cuerest stimlen stimrest resplen resprest (in s)
% TR does not work right now
% use zeros to use default rather than variable hox parameters
% ---------------------------------------------------------------
GA.numhox = 1;
GA.hoxrange = [1 16];
    % rows are each hox gene's allowable range, cols index hox elements

% Genetic algorithm parameters
% ---------------------------------------------------------------
nmodels = 1; 
GA.cbalColinPowerWeights = [0 1 1 1];	% 1 = cbal, 2 = eff, 3 = hrf shape, 4 = freq
GA.numGenerations = 100000; 
GA.sizeGenerations = 40;  
GA.maxTime = 20;						% max time to run in s, or Inf
GA.alph = 2.1;  
GA.plotFlag = 0; 

% Filtering, counterbalancing, and design tolerance
% ---------------------------------------------------------------
GA.lowerLimit = [];  
GA.HPlength = [120]; 
GA.LPsmooth = []; 
GA.maxOrder = 1; 
GA.NumStimthresh = []; 
GA.maxCbalDevthresh = []; 	% counterbalancing deviation
GA.maxFreqDevthresh = []; 

% Contrast setup
% ---------------------------------------------------------------
GA.contrasts = [1 -1];                            % three elements per condition for epochs (3 epochs)!
GA.contrastweights = [1];	                      % or predictor weights, if no contrasts


% Autocorrelation and special options
% ---------------------------------------------------------------
AutocorrelationFileName = 'myscannerxc';
GA.restlength = []; 
GA.restevery = []; 
GA.trans2switch = 0; 
GA.trans2block = 0; 
GA.dofirst = 0; 
GA.nonlinthreshold = [1.2]; 

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

   M = optimizeGA_epochs(GA);


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



