% does not consider relationships between blocks - won't work.

% k = number of conditions
% n = counterbalancing order
% numB is number of building blocks = permutations to include
% seqLength is the minimum length of the counterbalanced sequence

% add 1 to n, so that n is the total length of trial sequences to be considered.
n = n+1;

numB = k.^n;
seqLength = k.^(n+1);

% * get matrix of all possible permutations over n trials
blk = ones(n,k);

for i = 1:n
	veclen = numB ./ (k * i);

	vec = [];	
	for j = 1:k
		vec = [vec k*ones(1,veclen)];
	end

	whos vec

	vec2 = [];
	while length(vec2) < numB
		vec2 = [vec2 vec];
	end

	whos vec2

	blks(i,:) = vec2;
end

