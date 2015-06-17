function [X,paramvec,conditions] = construct_model(mspec,conditions,paramvec,varargin)
%[X,paramvec,conditions] = construct_model(mspec,conditions,paramvec,[dooverlap])
%
% paramvec is a cell array that completely specifies the choices of random variables
% to re-create a design, use the original conditions,mspec, and paramvec
%
% if you specify variable trial onsets in condition(i).onsets, then any subsequent
% sub-trial parts you specify as fixed onset may be randomly chosen within the range
% of the first two sub-part onsets.  Why?  The algorithm creates a varying number of
% trials, so must have as many sub-trials of each type.
%
% dooverlap = 1 or 0, last (optional) argument; constrains events to not overlap with
% other trials
% If a 2nd element is entered, this is the min time in TRs between onsets of 
% different trials - this performs a selective assigment of all trial times to conditions
% before calculating subparts, so may not work well with multiple trial types,
% but works well with multiple conditions.
%
% Tor Wager

paramindex = 1;
X = []; alltons = [];
if length(varargin) > 0, nooverlap = varargin{1};,else, nooverlap = 0;,end

if isfield(mspec,'toverlap'), nooverlap = mspec.toverlap;,end

% 2nd element of varargin is minimum abs trial distance between all trials, in TRs
if length(nooverlap) > 1, mintrialtime = nooverlap(2);, nooverlap = nooverlap(1);, 
else mintrialtime = []; 
end


    
% RANDOMIZE ALLTONS INTO CONDITIONS
if nooverlap & ~isempty(mintrialtime) & isempty(paramvec)
    alltons = get_all_trialonsets(conditions(1),mspec,paramvec,1,mintrialtime);
    numc = length(conditions);
    paramvec = assign_random(alltons,numc);
    
elseif nooverlap & ~isempty(paramvec) & ~isempty(mintrialtime)
    
        paramvec = constrain_toverlap(paramvec,mintrialtime,mspec.numframes);
        nooverlap = 0; % taken care of
        
elseif nooverlap
        
    alltons = get_all_trialonsets(conditions,mspec,paramvec,1,mintrialtime);
    
end



% REGULAR (PREVIOUS) WAY
    
for i = 1:length(conditions)
    
    if isempty(conditions(i).parts), conditions(i).parts = 0;,end
    
    if conditions(i).parts == 1
        % make fixed regressor or set of regressors
        if nooverlap
            [x,paramvec,paramindex,ons] = cond2reg(conditions(i),mspec,paramvec,paramindex,alltons);
        else
            [x,paramvec,paramindex,ons] = cond2reg(conditions(i),mspec,paramvec,paramindex);
        end
            if isempty(X), X = x; else, X = [X x];,end
        
    else
        
        % get trial onsets
        trialons = conditions(i).onsets;
        if size(conditions(i).onsets,2) == 2   % variable range of trial onsets
            if length(paramvec) < paramindex
                % build new ones if they're not input
                paramvec{paramindex} = onsetbuilder(conditions(i).onsets,mspec.numframes);
            end
            trialons = paramvec{paramindex};
            paramindex = paramindex + 1;
        end
        
        % get subpart onsets    
        for j = 1:length(conditions(i).subcond)
            if conditions(i).subcond(j).rel_to == 1
                conditions(i).subcond(j).baseons = trialons;
            elseif conditions(i).subcond(j).rel_to == 2
                conditions(i).subcond(j).baseons = trialons;
                conditions(i).subcond(j).add2base = ons;
            else
                error('rel_to field must be 1 or 2')
            end
            
            % store onset parameters in paramvec{paramindex}
            % return ons so that onsets for prev parts are stored in add2base and added to later event onset times
            if nooverlap
                [x,paramvec,paramindex,ons] = cond2reg(conditions(i).subcond(j),mspec,paramvec,paramindex,alltons);
                
                
                % tmp special for Emily partial trial design
               %if j == 2
               %    tmppv = paramvec; 
               %    tmppv{i} = tmppv{i}(3:end); % does not matter

                %   tmpcond = conditions(i).subcond(j);
                %   whtmp = getRandom((1:length(tmpcond.baseons))'); whtmp = whtmp(1:10);
                %    tmpcond.baseons(whtmp)=[];tmpcond.add2base(whtmp) =[];
                %    [x,paramvec,paramindex,ons] = cond2reg(tmpcond,mspec,paramvec,paramindex);
               %end
                
                
            else
                [x,paramvec,paramindex,ons] = cond2reg(conditions(i).subcond(j),mspec,paramvec,paramindex);
            end
            if isempty(X), X = x; else, X = [X x];,end  
        end
        
    end
    
end

% add intercept
X(:,end+1) = 1;


% contrasts across dx columns, if entered
if isfield(mspec,'contrasts')
    numframes = [];,    
    for i = 1:length(conditions),
        for j=1:length(conditions(i).subcond), numframes(end+1) = conditions(i).subcond(j).hrfest;,
        end
    end
    X = apply_dx_contrast(X,numframes,mspec.contrasts);
end


return







% -----------------------------------------------------------------------------------
% * sub-functions
% -----------------------------------------------------------------------------------

function [x,paramvec,paramindex,ons] = cond2reg(c,mspec,paramvec,paramindex,varargin)
% given a condition structure, returns regressor
%
% x is reg or set of regs to add to model
% paramvec is cell array of onsets (relative to...) for design reconstruction
% paramindex is index of which cell in paramvec to use for reconstruction of which condition
% ons is vector of trial onsets, for updating conditions, if necessary
%
% varargin is list of all trials, to constrain to no overlap

    if isfield(c,'parts')
        if c.parts ~= 1
            % can't build regressor from trial type with parts - return
            x = [];
           return
       end
    end
   
    % -----------------------------------------------------------------------------------
    % determine onset times
    % -----------------------------------------------------------------------------------
    delta = zeros(mspec.numframes,1);
    if size(c.onsets,1) == 1
        % ------------------------------------------------------------
        % fixed onsets for condition, but sub-trial may be variable
        % ------------------------------------------------------------
        
        if isfield(c,'baseons')     % for sub-parts of trials
            
            % if no range, just make it the same for all trials, no paramvec required
            if ~any(diff(c.onsets)), c.onsets=repmat(c.onsets(1),1,size(c.baseons,2));,end
            
            if size(c.onsets,2) ~= size(c.baseons,2) %== 2    
            % sub-trial, either range is given or length does not match variable trial onsets; get random
            % uses range of c.onsets(1:2), so if spec as fixed, returns trials within this range
                if length(paramvec) < paramindex,
                    paramvec{paramindex} = round(randrange(c.onsets(1:2),length(c.baseons)));
                end
                c.onsets = paramvec{paramindex};
                paramindex = paramindex + 1;
            end
            
            len = min(length(c.baseons),length(c.onsets));
            
            ons = c.onsets(1:len) + c.baseons(1:len);
            
            % replace baseons with onsets for previous trial for relative onsets
            % leave baseons as the trial onsets in case we need to compute trial overlaps
            if isfield(c,'add2base')
                if ~isempty(c.add2base)
                    ons = c.onsets(1:len) + c.add2base(1:len);
                end
            end
            
            % constrain overlap, if specified
            if length(varargin) > 0, 
                tons = varargin{1};, tons = [tons Inf];
                for i = 1:length(c.baseons)
                    wh = abs(c.baseons(i) - tons);
                    wh=find(wh == min(wh));
                    wh = wh(1);
                    ons(i) = min(ons(i),tons(wh+1)-1);
                end
            end
            
        else % no baseons field
            ons = c.onsets;
            
        end
        
        ons = round(ons);
        
        
        ons(ons < 0) = 0;
        
        if any(1+ons == 0), warning('Some onsets are zero!'),ons,keyboard,end
        if any(1+ons < 0), warning('Some onsets are less than zero!'),ons,keyboard,end
        if ~isreal(1+ons), warning('Ons is not real!'),ons,keyboard,end
        %if any(1+ons > length(delta)), warning('Onsets exceed length of delta!'),end
        
        try,delta(1 + ons) = 1;,catch, disp('Unknown error!'),keyboard,end
        
        % add repeated elements to delta, if specified
        if c.stimlength > 1,
            reps = repeatindex(ons,c.stimlength);
            delta(1 + reps) = 1;
        end
        
    else
        % ------------------------------------------------------------
        % variable onsets
        % ------------------------------------------------------------
        % variable (random) onsets
        if length(paramvec) < paramindex
            % build new ones if they're not input
            paramvec{paramindex} = onsetbuilder(c.onsets,mspec.numframes);
        end
        
        if isfield(c,'baseons')     % for sub-parts of trials
            % I don't think this should ever happen.
            ons = paramvec{paramindex} + c.baseons;
            
            % constrain overlap, if specified
            if length(varargin) > 0, 
                tons = varargin{1};, tons = [tons Inf];
                for i = 1:length(c.baseons)
                    wh = abs(c.baseons(i) - tons);
                    wh=find(wh == min(wh));
                    wh = wh(1);
                    ons(i) = min(ons(i),tons(wh+1)-1);
                end
            end
            
        else 
            ons = paramvec{paramindex};
        end
        delta(1 + (paramvec{paramindex})) = 1;   
        paramindex = paramindex + 1;
        
		% if c is the main condition vector with only 1 part, then
		% we don't have a stimlength field, so get it from the subpart
		if ~isfield(c,'stimlength'), c.stimlength = c.subcond(1).stimlength;, end
		
        % add repeated elements to delta, if specified        
		if c.stimlength > 1,
            reps = repeatindex(ons,c.stimlength);
            delta(1 + reps) = 1;
        end
    end

    
    
    
    % -----------------------------------------------------------------------------------
    % convolve or shift
    % -----------------------------------------------------------------------------------
    
    if isfield(c,'parts')
        if c.parts == 1
            doconv = c.subcond(1).convolve;
            hrfest = c.subcond(1).hrfest;
        else
            error('This should never ever happen.')
        end
    else
        % this c IS a subpart
        doconv = c.convolve;
        hrfest = c.hrfest;
    end
   
    if ~isstr(doconv) & doconv == 1,
        % Canonical SPM hrf
         x = conv(mspec.hrf,delta);
     elseif isstr(doconv)
         x = getPredictors(stimList, mspec.hrf) %*************
    elseif ~isempty(hrfest)
         sf{1} = delta;
         [x] = tor_make_deconv_mtx(sf,round(hrfest),1);
    else
         x = delta;
    end
     
    x = x(1:ceil(mspec.numframes),:);
    
return
    



% -----------------------------------------------------------------------------------
% * fill the run with as many trials as possible, given range of random lengths and delays
% -----------------------------------------------------------------------------------

function [onsets] = onsetbuilder(range,numframes)
% builds list of trials with random onset times
% range is range of onset times in condition.onsets, a 2 x 2 matrix

%if length(range) > size(range,1), range = range';, end

onsets = []; tend = 0;  % 0 is first time point in run
dlen = round(randrange(range(1,:),1));
dnaindex = 2;

while tend + dlen < numframes     % while onset of next trial is < end of run
    
    tlen = round(randrange(range(2,:),1));  % trial length
    if tend + dlen + tlen > numframes, break,end  % trial will not fit in run
    
    % add a trial
    onsets(end+1) = tend + dlen;            % specify onset after start delay of dlen
    tend = onsets(end) + tlen;              % specify ending point of trial
    dlen = round(randrange(range(1,:),1));  % re-randomize delay for next trial
    
end




    
% -----------------------------------------------------------------------------------
% * get a random number or vector within a specified range
% -----------------------------------------------------------------------------------    
function rval = randrange(range,num)
    rval = rand(1,num) * (range(1,2) - range(1,1)) + range(1,1);
return



% -----------------------------------------------------------------------------------
% * given onsets and length of stim, find indices of 'on' elements following onsets
% -----------------------------------------------------------------------------------    
function a = repeatindex(ons,stimlength);
            stimlength = round(stimlength);
            a = repmat(ons,stimlength-1,1); 
            b = (1:stimlength-1)'; b = repmat(b,1,size(ons,2)); 
            a = a + b; a = a(:);
return







% -----------------------------------------------------------------------------------
% * given conditions, get concatentated onsets of all trial types (conditions)
% ----------------------------------------------------------------------------------- 

function alltons = get_all_trialonsets(conditions,mspec,paramvec,paramindex,varargin)
    
alltons = [];


for i = 1:length(conditions)
    
    if isempty(conditions(i).parts), conditions(i).parts = 0;,end
    
    if conditions(i).parts == 1
        % make fixed or variable regressor or set of regressors
        %[d1,d2,d3,ons] = cond2reg(conditions(i),mspec,paramvec,paramindex);
        ons = onsetbuilder(conditions(i).onsets,mspec.numframes);
    else
        
        % get trial onsets
        ons = conditions(i).onsets;
        if size(conditions(i).onsets,2) == 2   % variable range of trial onsets
            if length(paramvec) < paramindex
                % build new ones if they're not input
                paramvec{paramindex} = onsetbuilder(conditions(i).onsets,mspec.numframes);
            end
            ons = paramvec{paramindex};
            
            % add one for each subpart not constant
            toadd = 1;
            for j = 1:conditions(i).parts
                if any(diff(conditions(i).subcond(j).onsets))  
                    toadd = toadd + 1;
                end
            end
            paramindex = paramindex + toadd;    % jump to next trial start
        end
    end % fixed or variable onsets
    
    alltons = [alltons ons];
            
end % loop thru conditions

alltons = sort(alltons);
alltons = alltons(alltons < mspec.numframes);

% if varargin, no overlap, so constrain trial onsets to minimum time
if length(varargin) > 0
    if ~isempty(varargin{1})
        mymin = varargin{1};
        tmp = [Inf diff(alltons)];
        wh = find(tmp < mymin);
        while ~isempty(wh) & alltons(wh(1)) < mspec.numframes
            alltons(wh(1):end) = alltons(wh(1):end) + (mymin - tmp(wh(1)));
            tmp = [Inf diff(alltons)];
            wh = find(tmp < mymin);
        end
    end
end
   
alltons = alltons(alltons < mspec.numframes);

return

    


function paramvec = assign_random(alltons,numc)
%
% alltons : onsets, m-sequence 1's
% numc : number of conditions
% returns paramvec, cell array of numc cells with onsets for each condition
% type.

ind = []; 
tmp1 = floor(length(alltons)./numc);    % at least this many stim per condition (even distribution)

    % create ind vectors of length = num onsets of each type, value =
    % condition number
    for i = 1:(numc), tmp = ones(tmp1,1) .* i; ind = [ind; tmp];,end
    
    % do this first, make sure critical stim don't appear at the end
    ind = getRandom(ind);                   
    ind = [ind; round(randrange([1 5],length(alltons)-length(ind)))'];

    for i = 1:(numc),
        paramvec{i} = alltons(ind == i);
    end
    
return

