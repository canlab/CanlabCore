function [rna,allmodels,alldelta,all_stimlist,numtrials] = dna2model_epochs(listMatrix,restMatrix,hox,TR,restlength,scanLength,numStim,dorests,trans2switch,trans2block,dofirst,HRF,nonlinthreshold,numsamps,S,varargin)
% function [rna,allmodels,alldelta,all_stimlist,numtrials] = dna2model(listMatrix,restMatrix,hox,TR,restlength,scanLength,numStim,dorests,trans2switch,trans2block,dofirst,HRF,nonlinthreshold,numsamps,S,[opt] verbose)
%
% transforms 'genetic' code used in crossbreeding (dna)
% into 'expressed' form (rna) used to build the design matrix 
% in high-frequency sampling.
%
% works with one or matrix of lists, returns models in cell array
%
% Tor Wager, 12/29/01


verbose = 0;
if ~isempty(varargin) % exist('varargin') == 1
    verbose = varargin{1};
end

% before start of generation - ops on overall list - don't do, bec no rests allowed in epochs
% ---------------------------------------------------------------------------------------------
%if dorests,
%		listMatrix = insert_rests(listMatrix,restMatrix,restlength,numStim);
%end


% within generation - ops on each design vector
% ---------------------------------------------------------------------------------------------

for z = 1:size(listMatrix,2) 		% do this for each organism in the generation
	stimList = double(listMatrix(:,z));

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
    numhox = length(hox);
    
    % get ISI and TR
    ISI = hox(1);
    % TR = hox(2);  % TR cannot vary because S depends on it.  fixed.
    
    % figure out trial length in seconds
    tl = sum(hox(3:end));       % sum lengths of each part of the trial - cue, stim, response
    
    % get number of trials in scan length
    numtrials = ceil(scanLength / tl); 
    
    % check to make sure stimList is long enough
    if size(stimList,1) < numtrials, 
        warning(['stimlist is too short: ' num2str(size(stimList,1)) ' stimuli, ' num2str(numtrials) ' trials']), 
        stimList = [stimList; stimList];    
    end
    
    % truncate stim list and save in truncated listMatrix of stim lists
    stimList = stimList(1:numtrials);
    all_stimlist(1:numtrials,z) = stimList;
    
    % get trial length in samples (high sampling freq)
    % sampres = the resolution of sampling (.1 s default)
    sampres = .1;
    tl = ceil(tl / sampres);
    
    % build trials
    newstimList = [];
    for i = 1:numtrials
        mytrial = [];
        hoxindex = 1;
        for myhox = 3:2:numhox                                              % must be even! set up for 3 periods.
            mytrial = [mytrial; ones(round(hox(myhox)./sampres),1) * hoxindex];   % add condition number to trial
            mytrial = [mytrial; zeros(round(hox(myhox+1)./sampres),1)];            % add rest length (every other hox)
            hoxindex = hoxindex + 1;
        end
        mytrial(mytrial > 0) = mytrial(mytrial>0) + max(mytrial) * (stimList(i)-1);  % shift numbers up to reflect condition number
        if any(mytrial < 0), error('Stimlist should not contain negative numbers or zeros'),end
        newstimList = [newstimList; mytrial];
    end
    
    % pad hi-resolution stim list with zeros if necessary
    num2add = scanLength ./ sampres - length(newstimList);
    if num2add > 0, newstimList = [newstimList; zeros(num2add,1)];, end
    
    newstimList = newstimList(1:scanLength ./ sampres);
    
	rna(:,z) = newstimList;
    
    % this part takes the place of designvector2model

    [model,delta] = getPredictors(newstimList,HRF);
    model = resample(model,1,round(TR*10));
    if ~isempty(nonlinthreshold),model = modelSaturation(model,nonlinthreshold);,end		% saturation (nonlinear responses)      
    if size(model,1) > numsamps, model = model(1:numsamps,:);,end    
    if ~isempty(S)              
    		model = S * model;                                                              % temporal smoothing and HP/LP filter
    end
    model(:,end+1) = 1;																		% add intercept

    allmodels{z} = model;
    alldelta{z} = delta;

end

if verbose
    disp(['Epoch model parameters'])
    disp(['------------------------------------------'])
    disp(['Num hox genes = ' num2str(numhox)])
    disp(['Num trials for this model = ' num2str(numtrials)])
end

all_stimlist = uint8(all_stimlist);

return