function stimList = getRandom(stimList)
% stimList = getRandom(stimList)
%
% Randomizes rows of a matrix, preserving dependencies between columns
% Tor Wager, created 2 / 01, last modified 2/25/04 to sort whole matrix

%input: a col. vector of stimulus conditions, e.g. 1 1 1 1 2 2 2 2 3 3 4 4

%output: a randomized permutation of this vector

randvector = rand(size(stimList,1),1);						% create random vector

stimList(:,size(stimList,2)+1) = randvector;

stimList = sortrows(stimList,size(stimList,2));				% sort the rows by random seed

stimList = stimList(:,1:end-1);

return

