function [DX,delta] = onsets2dx(onsets, TR, scansperrun, numseconds)
% :Usage:
% ::
%
%     [DX,delta] = onsets2dx(onsets, TR, scansperrun, numseconds)
%
% :Inputs:
%
%   **onsets:**
%        cell array whose length is num runs * num conditions,
%        e.g., {run 1 ant onsets, run 1 stim onsets, run 2 ant onsets,
%        run 2 stim onsets}
%
%   **TR:**
%        TR in seconds
%
%   **scansperrun:**
%        number of volumes in each run
%
%   **numseconds:**
%        number of seconds after event onsets to generate regressors for [default: 32]
%
% :Examples:
% ::
%
%    TR = 2;
%    scansperrun = [192 196 196 184 190 192];
%    numseconds = 30;
%    [DX,delta] = onsets2dx(onsets, TR, scansperrun, numseconds)
%    EXPT.FIR.model{subjectnumber} = DX;

    if ~iscell(onsets)
        for i = 1:size(onsets,2), tmp{i} = onsets(:,i); end 
        onsets = tmp;
    end
    
    if(~exist('numseconds', 'var') || isempty(numseconds))
        numseconds = 32;
    end
    
    evtsperrun = length(onsets) ./ length(scansperrun);

    delta = {};

    for i = 1:length(scansperrun)
        % onsets for this run
        st = (i - 1) .* evtsperrun + 1;
        wh = st:(st + evtsperrun - 1);

        % delta for this run
        [x,d] = onsets2delta(onsets(wh), TR, scansperrun(i)*TR - 1);

        delta = [delta; d];
    end

    DX = tor_make_deconv_mtx3(cell2mat(delta), round(numseconds./TR), 1, 0, 1, 0, scansperrun);
end
