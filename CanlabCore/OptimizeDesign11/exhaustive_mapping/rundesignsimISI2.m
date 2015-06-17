function out = param_map(GA)

% rundesignsimISI2
% Tor Wager, 11/24/01
% if your input var name is already a variable in the wkspace, uses that...
% otherwise creates one with zeros.

warning off
index = 1;

% -------------------------------------------------
% load GA field names
% -------------------------------------------------
fN = fieldnames(GA);
for i = 1:size(fN,1)
    eval([fN{i} ' = GA.' fN{i}])
end

% -------------------------------------------------
% common user variables
% -------------------------------------------------
tickRes = [.1 .1];
ISIrange = [.1 16]
FCrange = [0 0];
GA.sizegenerations = 100;

myISI = ISIrange(1):tickRes(1):ISIrange(2);	
myFC = FCrange(1):tickRes(2):FCrange(2);

ISIsize = length(myISI);
FCsize = length(myFC);

ISItick = 1:max(round(ISIsize./25),1):ISIsize;
ISIticklab = myISI(ISItick);
FCtick = 1:2:FCsize;
FCticklab = myFC(FCtick);

disp(['Starting ' num2str(prod([ISIsize FCsize])) ' surface points.'])

out.ga = GA;
out.myISI = myISI;
out.myFC = myFC;

% ----------------------------------------------------------------
% * HRF GAMMA FUNCTION - spm 99
% ----------------------------------------------------------------

   HRF = spm_hrf(.1);						% SPM HRF sampled at .1 s
   HRF = HRF/ max(HRF);

% ----------------------------------------------------------------
% * make the map
% ----------------------------------------------------------------

for RR = 1:size(myISI,2)
	
    % set ISI
    GA.ISI = myISI(RR);
    numsamps = ceil(numStim*GA.ISI/TR);
    
    out.numsamps(find(myISI == GA.ISI)) = numsamps;
    
    % ----------------------------------------------------------------
    % * get smoothing matrix and autocorrelation matrix
    % ----------------------------------------------------------------
    if ~isfield(GA,'LPsmooth'),GA.LPsmooth = 1;,end
    [S,Vi,svi] = getSmoothing(HPlength,GA.LPsmooth,TR,numsamps,xc);
    clear Vi


	for CC = 1:size(myFC,2)
		disp(['Starting model ' num2str(index)])		
	
		% set rest frequency
		GA.freqConditions(end+1) = myFC(CC);
		GA.freqConditions(1:end-1) = (1 - myFC(CC)) / (size(GA.freqConditions,2)-1);

		disp(['FREQC and ISI are:	' num2str(GA.freqConditions) ' and ' num2str(GA.ISI)])
	
        % ----------------------------------------------------------------
        % * generate an unsorted list of conditions
        % ----------------------------------------------------------------   
        numStim = ceil(scanLength / (GA.ISI));
        if ~isempty(restevery) & ~isempty(restlength),
            numRestStim = (ceil(numStim/(mean(restevery)+restlength)) - 1) * restlength;
        else
	        numRestStim = 0;
        end
        numStimEachCond = ceil((numStim-numRestStim) * freqConditions);			            % row vector of stim in each cond
        numBlanks = numStim - numRestStim - sum(numStimEachCond);
           
        unsortedList = [];
        for i = 1:size(conditions,2)	
            unsortedList = [unsortedList; conditions(i) * ones(numStimEachCond(i),1)];
        end
        unsortedList = [unsortedList; zeros(numBlanks,1)];
            
        % ----------------------------------------------------------------
        % * randomize order and test each list
        % ----------------------------------------------------------------

	    clear listMatrix
        for z = 1:GA.sizeGenerations
            
  	        stimList = getRandom(unsortedList);

            % -------------------------------------------------------------------------------------------------
		    % * counterbalancing
		    % -------------------------------------------------------------------------------------------------
	
            cbal(z) = getCounterBal(stimList, maxOrder,conditions,freqConditions);
            
            % -------------------------------------------------------------------------------------------------
		    % * efficiency
		    % -------------------------------------------------------------------------------------------------
	
			model = designvector2model(stimList,ISI,HRF,TR,numsamps,nonlinthreshold,S);
        	xtxitx = pinv(model);                                       		% inv(X'S'SX)*(SX)'; pseudoinv of (S*X)
			eff(z) = calcEfficiency(contrastweights,contrasts,xtxitx,svi);
            
            % -------------------------------------------------------------------------------------------------
		    % * HRF shape estimation efficiency
		    % -------------------------------------------------------------------------------------------------
            
            delta = [];
            for i = 1:max(stimList(:,1))
                delta(:,i) = (stimList == i);
            end
            
			[model] = tor_make_deconv_mtx2(delta,5,TR / ISI);
            if ~isempty(S), model = S * model;,end
            
        	xtxitx = pinv(model);                                       		% inv(X'S'SX)*(SX)'; pseudoinv of (S*X)
			hrf_eff(z) = calcEfficiency([],[],xtxitx,svi);
            

            % -------------------------------------------------------------------------------------------------
		    % calculate average time between stimulus repetitions in each trial type
		    % -------------------------------------------------------------------------------------------------
            
            for i = 1:max(stimList),mytimeBtwn(i) = mean(diff(find(stimList==i))) .* GA.ISI;, end
            timeBtwn(z) = mean(mytimeBtwn);
            
            % old method:
            %[fitness,convar,xtxi,models] = testlist(listMatrix,GA,HRF,'svi',svi,'S',S);
            %convar = 1./convar;
	        %eval([name '(RR,CC) = mean(convar,2);']);
            
        end % testing at this parameter combination
        
        out.cbal(RR,CC) = mean(cbal);
        out.eff(RR,CC) = mean(eff);
        out.hrf_eff(RR,CC) = mean(hrf_eff);
		
        out.cbal_se(RR,CC) = std(cbal) / sqrt(sizeGenerations);
        out.eff_se(RR,CC) = std(eff) / sqrt(sizeGenerations);
        out.hrf_eff_se(RR,CC) = std(hrf_eff) / sqrt(sizeGenerations);
        out.timeBtwn(RR,CC) = mean(timeBtwn);
        out.timeBtwn_se(RR,CC) = std(timeBtwn) / sqrt(sizeGenerations);
        
		%eval(['save ' name num2str(index) ' ' name ' SIM'])
		index = index + 1;
	end
end


