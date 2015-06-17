function [sf,varargout] = downsample_delta(sf,eres)
% [sf_cell,sf_matrix] = downsample_delta(sf,eres)
% 
% Tor Wager, 2/25/04
% Dowsample an indicator matrix (sf) of ones and zeros, taking every eres-th
% element.  
% This function will preserve the number of ones in sf.  We're assuming
% that a one represents a trial onset, and that you want ones in the
% nearest bin of the downsampled matrix.
%
% 1st output argument is a cell array, 2nd is a matrix with the same
% information, which is the downsampled indicator matrix.

if ~iscell(sf)
    for i = 1:size(sf,2), sf2{i} = sf(:,i);, end
    sf = sf2;
end

% -------------------------------------------------------------------
% * downsample sf to number of TRs
% -------------------------------------------------------------------
numtrs = ceil(length(sf{1}) ./ eres);
myzeros = zeros(numtrs,1);
origsf = sf;

for i = 1:length(sf)
    Snumtrs = length(sf{i}) ./ eres;
    if Snumtrs ~= round(Snumtrs), warning(['sf{ ' num2str(i) '}: length not evenly divisible by eres.']),end
    if numtrs ~= Snumtrs, warning(['sf{ ' num2str(i) '}: different length than sf{1}.']),end

    inums = find(sf{i} > 0);
    inums = inums ./ eres;      % convert to TRs
    inums = ceil(inums);        % nearest TR
                                % i'm getting weird effects with round
                                % and stimuli onsets in between TRs -
                                % if the 1 in the matrix occurs before the
                                % onset of the actual event, the est. HRF is inverted -
                                % thus i'm now using ceiling.
    inums(inums == 0) = 1;      % never use 0th element
    sf{i} = myzeros;
    sf{i}(inums) = 1;           % always use 1 for sf
end

if nargout > 1
    varargout{1} = cell2mat(sf);
end

return
