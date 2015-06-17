clear Models, clear MM

% User-Specified Parameters
% ===============================================================

% Conditions and stimuli
% ---------------------------------------------------------------
GA.conditions = [1 2 3 4 5];
	% a vector of integers 1:# of event-related conditions (e.g., [1 2 3 4])
	% to add "rest" intervals, add an extra condition and set its contrast weight
	% to zero in all contrasts.
	% e.g., condition 3 in this design could be passive rest.

GA.freqConditions = [0.20 0.20 0.20 0.20 0.20];
	% vector of frequencies of each trial type; should sum to 1

GA.scanLength = 544;
	% how long your run is, in seconds.  
        % This and the ISI determine how many stimuli are in each design vector

GA.ISI = 2;  
	% how long between stimulus presentations? (you can also include "rest" presentations)
	% also the time resolution of stimulus condition function (list of stimuli)
	% designs will be constructed in time units of this resolution
    
GA.TR = 2; 
    % the TR (sampling resolution) of your experiment; time for volume acquisition

GA.doepochs = 0;        % build epochs instead of events
	% use this if you have epochs instead of events; alpha version.  use with caution.
	% if you use this, use the GA function optimizeGA_epochs.m
	% also, the contrast matrix will be different, with 3 periods per trial explicitly modeled
	% this wasn't intended to be a general solution, but a specific optimization for a certain design.

% Hox optimization parameters
% A hox gene is a master gene.  These numbers control the stimulus
% parameters, which can be optimized within the GA rather than
% pre-specified.
% Hox elements in this GA go in the stimlist as the first several elements.
% hox sequence currently codes, in number order: 
% ISI TR cuelen cuerest stimlen stimrest resplen resprest (in s)
% the last 6 parameters are for use with doepochs.  
% They code epoch lengths for 6 periods within each trial.
% TR does not work right now
% use zeros to use default rather than variable hox parameters
% ---------------------------------------------------------------
GA.numhox = 0;
GA.hoxrange = [];
    % rows are each hox gene's allowable range, cols index hox elements
	% use this if you want to determine the ISI, etc. using the GA
	% also may require significant user input/program modification at this stage.	

% Genetic algorithm parameters
% ---------------------------------------------------------------
nmodels = 1; 
	% how many runs do you want to optimize?
	% GA runs one separate optimization for each model.

GA.cbalColinPowerWeights = [0 1 1 .5];	% 1 = cbal, 2 = eff, 3 = hrf shape, 4 = freq
	% first element: counterbalancing of stimuli
	% second: contrast detection efficiency
	% third: hrf shape estimation efficiency
	% maintenance of input frequencies for each trial type
	% fitness scores for each measure are multiplied by these values
	% before collapsing to a single, final fitness measure for each design.
	% does not have to sum to 1.


GA.numGenerations = 100000; 
	% how many iterations of the GA to run.
GA.sizeGenerations = 20;  
	% how many designs to test  for each generation?  Population size.
GA.maxTime = Inf;						
	% max time to run in s, or Inf for infinite time
	% The GA stops when either numGenerations or maxTime is reached.

GA.alph = 2.1; 
	% "selection pressure": higher is more extreme selection; always pick the best 50% of designs
	%1 is no selection pressure at all, or random selection of designs for recombination
	% selection pressure is like in evolution; refers to whether the best designs are selected to
	% continue on to the next generation.
	% an intermediate value is best, to keep the population heterogeneity high.

GA.plotFlag = 0; 
	% plot results after each run.  
	% Recommended to leave this off if nmodels > 1, and then plot your final results later.

% Filtering, counterbalancing, and design tolerance
% ---------------------------------------------------------------
% Empty brackets indicate the option is not to be used.

GA.lowerLimit = [];  
GA.HPlength = [120]; 
	% high-pass filter length, in s, for analysis; [] for no HP filter
	% Used to filter out noise at lower frequencies than design during data analysis
	% The cutoff you should use depends on how much power in your design is below the cutoff -
	% you want all the power in your design to be above the cutoff.

GA.LPsmooth = ['hrf']; 
	% low-pass filter type for analysis; [] for no LP filter
	% choices are 'hrf' or []

GA.maxOrder = 1;   
    % order of counterbalancing to use
    % 1 = each trial type follows each other with equal probabilities, adjusted for base frequencies
    % 2 = one-back and two-back counterbalancing.  3 = 1+2+3 back, etc.


% Hard constraints
% ---------------------------------------------------------------
% Here you can specify tolerances for designs
% If designs do not meet these criteria, they will receive -infinity fitness scores
% Empty brackets indicate the option is not to be used.

GA.NumStimthresh = [];      % maximum number of repeats of any one event type in a row
GA.maxCbalDevthresh = []; 	% maximum acceptable counterbalancing deviation (actual freq - expected freq)
GA.maxFreqDevthresh = [];   % maximum acceptable deviation from input frequencies (GA.freqConditions)

% Contrast setup
% ---------------------------------------------------------------
% Here you specify contrasts across conditions, for use with the contrast estimation fitness measure
% Each contrasts should be a row in the matrix GA.contrasts
% 
GA.contrasts = [1 0 0 0 0  -1  0  0  0  0; ...
                0 1 0 0 0   0 -1  0  0  0; ...
                0 0 1 0 0   0  0 -1  0  0; ...
                0 0 0 1 0   0  0  0 -1  0];


%GA.contrasts = [0 1 0 0 0 0 -1 0 0 0; 0 0 1 0 0 0 0 -1 0 0];
%[1 0 0 0 -1 0 0 0;  0 1 0 0 0 -1 0 0;  0 0 1 0 0 0 -1 0; 0 0 1 0 0 0 -1 0];

    % there should be one column per condition in your design, not including the intercept.
    %if trans2switch or trans2block = 1, double the number of columns.
    % when using epoch design, three elements per condition for epochs (3 epochs)!
GA.contrastweights = [1 1 1 1];
    % Weighting function for contrasts (rows of GA.contrasts)
    % Contrast efficiencies will be multiplied by these weights before computing overall design fitness
    % If no contrasts are specified, this vector can specify weights for predictors (columns)


% Autocorrelation and special options
% ---------------------------------------------------------------
AutocorrelationFileName = 'myscannerxc';
    % This is the name of a mat file without the .mat extension
    % In the mat file, there should be a variable called myscannerxc
    % which is a row vector containing the autocorrelation fucntion you wish to use.
    % This function is used as the intrinsic noise autocorrelation estimate.
    % Included with the GA scripts are these functions:
    % myscannerxc   (U of M 3T, 7 subjects)
    % hiautocorr    (1/f)
    % hiautocorr3   (another 1/f)
    % noautocorr    (identity; no autocorrelation)
    
GA.restlength = [1]; 
    % if inserting probe or instruction periods, modeled with periodic rests in the simulation,
    % this is the length of the rest/probe/instr periods to insert, in units of the ISI
GA.restevery = [15]; 
    % This is how many regular events you want to have (in ISIs) between rests.
    % Leave restlength and restevery blank to avoid using these options.
GA.trans2switch = 0; 
    % 1 or 0.  This option works if the conditions field is 2 elements long.
    % This option creates new conditions 3 and 4 that occur whenever 1 is followed by 2 or 2 by 1
    % So this option creates trial history-dependent predictors.  1 and 2 are repeats, 3 and 4 are switches
    % Useful when studying habituation, switching between items, etc.
GA.trans2block = 1; 
    % 1 or 0.  This option doubles the input length of conditions - so [1 2] is transformed to [1 2 3 4]
    % This option creates alternating ABAB blocks of your input event types (e.g., 1 2) with the new
    % event types (e.g., 3 4).  This is useful when simulating mixed block/event related designs.
    % Trial types 1-2, for example, can be two events within block A, and 3-4 can be the same (or different)
    % events within block B.  The alternation frequency (in ISIs) is specified by GA.restevery.
    % If you use this option, the contrasts field must contain twice as many columns as you have input conditions.
    % for example, if conditions = [1 2 3 4], contrasts = [1 1 1 1 -1 -1 -1 -1] specifies block A vs block B.
    % contrasts = [1 1 1 1 -1 -1 -1 -1; 1 1 -1 -1 1 1 -1 -1] specifies A vs B and events [1 2 5 6] - [3 4 7 8].
GA.dofirst = 0; 
    % 1 or 0.  This creates a separate trial type for the first trial following rests in GA.restlength.
    % It only works if you're using trans2switch and trans2block together; unreliable under other circumstances.
GA.nonlinthreshold = [2]; 
    % If the value here is x, predictor heights are clipped (thresholded) at x times the unit HRF height.
    % If empty, no thresholding is performed.
    % This is useful for modeling nonlinearities in the BOLD signal, which is important if your design uses
    % short ISIs (generally, 2 s or less).  The clipping acts as a 'saturation' factor.
    % I like to use 2, because data in our lab suggests the the nonlinearities in the HRF can be
    % roughly approximated by a clipping function with this parameter value.
 
% ---------------------------------------------------------------
%
% * E N D   U S E R   I N P U T
%
% ---------------------------------------------------------------
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

   M = optimizeGA(GA);


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



