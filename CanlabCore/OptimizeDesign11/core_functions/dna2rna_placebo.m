function [rna,allmodels,alldelta] = dna2model(listMatrix,restMatrix,restlength,numStim,dorests,trans2switch,trans2block,dofirst,HRF,nonlinthreshold,numsamps,S)
% function [rna,allmodels,alldelta] = dna2model(listMatrix,restMatrix,restlength,numStim,dorests,trans2switch,trans2block,dofirst,HRF,nonlinthreshold,numsamps,S)
%
% transforms 'genetic' code used in crossbreeding (dna)
% into 'expressed' form (rna) used to build the design matrix 
% in high-frequency sampling.
%
% Tor Wager, 12/29/01


% before start of generation - ops on overall list
% ---------------------------------------------------------------------------------------------
if dorests,
		listMatrix = insert_rests(listMatrix,restMatrix,restlength,numStim);
end


% within generation - ops on each design vector
% ---------------------------------------------------------------------------------------------

for z = 1:size(listMatrix,2) 		% do this for each organism in the generation
	stimList = listMatrix(:,z);

    % The Switch Hack: transform list of stimuli into list of switches/no switches for each stimtype (R or L)
    if trans2switch, stimList = transform2switches(stimList);,end
    if trans2block, 
		  if isempty(restMatrix),error('trans2block will not work without some rest intervals specified.'),end
          stimList = transform2block(stimList,restMatrix(:,z),restlength,dofirst);
    end
	
    
    % stretch out trials to number of actual events
    
    % hox: a list of parameters added to listMatrix
    % numhox: number of 'hox' parameters added
    % hox are added to front end of stimlist
    % hox codes, in number order: ISI TR cuelen cuerest stimlen stimrest resplen resprest
    hox = stimList(1:numhox);
    
    % get ISI and TR
    ISI = hox(1);
    TR = hox(2);
    
    % figure out trial length in seconds
    tl = sum(hox(3:end));       % sum lengths of each part of the trial - cue, stim, response
    
    % get number of trials in scan length
    numtrials = ceil(scanLength / tl); 
    
    % check to make sure stimList is long enough
    if length(stimList) < numtrials, error(['stimlist is too short: ' num2str(length(stimList))]), end
    
    % get trial length in samples (high sampling freq)
    % sampres = the resolution of sampling (.1 s default)
    sampres = .1;
    tl = ceil(tl / sampres);
    
    % build trials
    newstimList = [];
    for i = 1:numtrials
        mytrial = ones(tl,1);
        hoxindex = 2;
        for myhox = 3:numhox
            mytrial(hox(myhox)./sampres+1:hox(myhox+1)./sampres) = hoxindex;
            hoxindex = hoxindex + 1;
        end
        mytrial= mytrial + max(mytrial) * (stimList(i)-1);  % shift numbers up to reflect condition number
        newstimList = [newstimList; mytrial];
    end
    
	rna(:,z) = newstimList;

    % this part takes the place of designvector2model

    [model,delta] = getPredictors(newstimList,HRF);
    model = resample(model,1,TR*10);
    if ~isempty(nonlinthreshold),model = modelSaturation(model,nonlinthreshold);,end		% saturation (nonlinear responses)      
    if size(model,1) > numsamps, model = model(1:numsamps,:);,end    
    if ~isempty(S)
    		model = S * model;                                                  			% temporal smoothing and HP/LP filter
    end
    model(:,end+1) = 1;																		% add intercept

    allmodels{i} = model;
    alldelta{i} = delta;

end

return