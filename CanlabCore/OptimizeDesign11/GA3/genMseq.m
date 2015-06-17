function [OUT,delta] = genMseq(samples,TR,P,varargin)
% function [OUT,delta] = genMseq(samples,TR,P,[number of sequence groups],[seqs per group])
%
% Generates a number of truncated m-sequences
% Tests their power with pre-set parameters, unless this option is turned
% off.
%
% samples   : length of sequence
% TR        : repetition time of sampling.  If empty, no testing is done.
% P         : a structure of options to specify testing parameters.
%             See hrf_power.m.  If empty, no testing is done.
% [seqs]    : optional, number of sequences to generate.
% [seqsperg]: optional, divide "on" values into conditions
%             use only if TR and P are empty.  Returns onsets in paramvec
%             cell array in delta variable (for compatibility with ga2).
%
% uses the m-sequence toolbox.
%
% tor wager, last update 9 / 3/ 04.
% 
% examples: 
% generates 10 m-sequences truncated to length 100
%mseqs = genMseq(100,[],[],10);
%
% generates m-sequences based on mspec and conditions in ga2.m
% [mseqs,newvec(1:nmseq)] = genMseq(mspec.numframes,[],[],nmseq,length(conditions));


% -----------------------------------------------------------------
% Set up input arguments and defaults
% -----------------------------------------------------------------

it = 100; nconds = 1; spergroup = [];
    
if length(varargin) > 0, it = varargin{1};,end
if length(varargin) > 1, spergroup = varargin{2};,end

if ~isempty(TR) & ~isempty(P)
    hrf = spm_hrf(TR);
    hrf = hrf ./ max(hrf);
    
    nconds = size(P.contrasts,2) - 1;   % exclude intercept
end

% should be correct, but mseqs are sometimes 1 shorter...
expon = ceil(log(samples) ./ log(2));


% -----------------------------------------------------------------
% Generate m-sequences
% -----------------------------------------------------------------

for i = 1:it
    
    warning off
    
    for j = 1:nconds
        ms = mseq(2,expon,1); 
        if length(ms) < samples, expon=expon+1;,ms = mseq(2,expon,1);,end   % fix for 1-too-short mseqs 
        ms = ms(1:samples);
    
        mb{i}(:,j) = m2bin(ms);
        delta{i}{j} = mb{i}(:,j);
    end    

    % divide, if necessary
    if (isempty(TR) | isempty(P)) & ~isempty(spergroup)
        delta{i} = assign_random(find(mb{i})',spergroup);
    end
    
end

% -----------------------------------------------------------------
% Test fitness of m-sequences
% -----------------------------------------------------------------

if ~isempty(TR) & ~isempty(P)
    
for i = 1:it
    
    X = getPredictors(mb{i},hrf);
    X(:,end+1) = 1;
    
    tmp = sum(diff(mb{i},2));   % test for empty conditions?
    if any(tmp(1:end-1) == 0), 
        for j = 1:size(P.contrasts,2) - 1
            ms = mseq(2,expon,1); ms = ms(1:samples);
    
            mb{i}(:,j) = m2bin(ms);
            delta{i}{j} = mb{i}(:,j);
        end    
        X = getPredictors(mb{i},hrf);
        X(:,end+1) = 1;
        
    end
    
    trueresp = zeros(1,size(P.contrasts,1)); trueresp(1) = 1;
    [ipow(:,i),gpow(:,i),OUT] = xpower(X,P.contrasts,trueresp,12,[],P.Vi);

    PP = hrf_power(TR,18,delta{i},P.contrasts(1,:));
    % save ipow from this output, get distribution, plot (onsets2eff code
    % use)
    ihrfpow(i) = max(PP.ipow);
    ghrfpow(i) = max(PP.gpow);
    
    warning on
end

OUT.ipow = ipow; OUT.gpow = gpow; OUT.mseq = mb;
OUT.ipowhrf = ihrfpow; OUT.gpowhrf = ghrfpow;

else    % no testing
    
    OUT = mb;
end

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
        
% -----------------------------------------------------------------------------------
% * get a random number or vector within a specified range
% -----------------------------------------------------------------------------------    
function rval = randrange(range,num)
    rval = rand(1,num) * (range(1,2) - range(1,1)) + range(1,1);
return


