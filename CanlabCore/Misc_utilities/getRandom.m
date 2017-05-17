function stimList = getRandom(stimList)
% Randomizes rows of a matrix, preserving dependencies between columns
%
% :Usage:
% ::
%
%     stimList = getRandom(stimList)
%
% :Input:
%   a col. vector or matrix of stimulus conditions,
%
%   e.g. [1 1 1 1 2 2 2 2 3 3 4 4]'
%
%output: a randomized permutation of this vector or matrix
% all columns are resorted with the same order
%
% ..
%    Tor Wager, created 2 / 01, last modified 2/25/04 to sort whole matrix
%    Modified 8/7/2015 by Tor to speed and simplify code
% ..


n = size(stimList, 1);

wh = randperm(n)';

stimList = stimList(wh, :);

% Old
%[n, k] = size(stimList);
% randvector = rand(n, 1);						% create random vector
% stimList(:, k + 1) = randvector;
% stimList = sortrows(stimList, k);				% sort the rows by random seed
% stimList = stimList(:, 1:end-1);

end % function

