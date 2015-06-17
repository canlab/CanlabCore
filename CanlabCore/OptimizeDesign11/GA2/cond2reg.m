function [x,paramvec,paramindex,ons] = cond2reg(c,mspec,paramvec,paramindex)
% given a condition structure, returns regressor
%
% x is reg or set of regs to add to model
% paramvec is cell array of onsets (relative to...) for design reconstruction
% paramindex is index of which cell in paramvec to use for reconstruction of which condition
% ons is vector of trial onsets, for updating conditions, if necessary
%
% This function is not used in GA2: it's an extra function that mirrors a subfunction of
% construct_model.m
%
% THIS IS OLD! DOES IT WORK?
%
% Tor Wager

    if isfield(c,'parts')
        if c.parts ~= 1
            % can't build regressor from trial type with parts - return
            x = [];
           return
       end
    end
   
    % determine onset times
    % -----------------------------------------------------------------------------------
    delta = zeros(mspec.numframes,1);
    if size(c.onsets,1) == 1
        % fixed onsets for condition, but sub-trial may be variable
        
        if isfield(c,'baseons')     % for sub-parts of trials
            if size(c.onsets,2) == 2    % sub-trial, range is given; get random
                if length(paramvec) < paramindex,
                    paramvec{paramindex} = round(randrange(c.onsets,length(c.baseons)));
                end
                c.onsets = paramvec{paramindex};
                paramindex = paramindex + 1;
            end
            ons = c.onsets + c.baseons;
            
        else 
            ons = c.onsets;
            
        end
        ons = round(ons);
        delta(1 + (ons)) = 1;
        
    else
        % variable (random) onsets
        if length(paramvec) < paramindex
            % build new ones if they're not input
            paramvec{paramindex} = onsetbuilder(c.onsets,mspec.numframes);
        end
        
        if isfield(c,'baseons')     % for sub-parts of trials
            ons = paramvec{paramindex} + c.baseons;
        else 
            ons = paramvec{paramindex};
        end
        delta(1 + (paramvec{paramindex})) = 1;   
        paramindex = paramindex + 1;
    end

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
   
    if doconv,
         x = conv(mspec.hrf,delta);
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

    
