function [cBal,fM,msd,eN] = getCounterBal(stimList, maxOrder,conditions,freqConditions)
% [cBal,fM,msd,eN] = getCounterBal(stimList, maxOrder,conditions,freqConditions)
%
% this function determines the counterbalancing fitness of a design vector
%
% Inputs
% column design vector
% max. order to consider in counterbalancing
% row vector of integer codes for conditions of interest, e.g. [1 2 3 4]
% row vect of condition frequencies, e.g. [.35 .35 .15 .15]
%
% Outputs
% cBal: 	counterbalancing fitness value
% fM: 		observed frequency matrix of each condition following each other one, nc x nc x time lag
% msd:		the square root of the maximum squared deviation from counterbalancing across all cells.
%
% by Tor Wager.  Last edit 11/17/01.

% this should be a row vector

warning off

nc = size(conditions,2);
fc = freqConditions(1:nc);

% eN is expected frequency matrix
eN = fc' * fc;
eN = eN ./ sum(sum(eN));		% normalize to sum of 1 so it's comparable with observed, excluding rest intervals
for i = 2:maxOrder,eN(:,:,i) = eN(:,:,1);,end

fM = zeros(nc,nc,maxOrder);

for i = 1:nc
	a = find(stimList == i);
		
	% fill in diagonals, fast
	% -----------------------------------------------------------
	for k = 1:maxOrder
		fM(i,i,k) = sum(diff(a) == k);
	end
			
	% off-diagonals
	% -----------------------------------------------------------
	for j = 1:nc
		
		if j ~= i
			b = find(stimList == j);
			
			%if isempty(b),
			%	warning(['No events of type ' num2str(j) ' in stimList'])
			% this is ok, so i'm turning warning off in this function to avoid this.
			%end
	
			for k = 1:maxOrder
					
				mysum = 0;
				for m = 1:length(a)
					mysum = mysum + sum(b-a(m) == k);	
				end
					
				fM(i,j,k) = mysum;
					
			end
		end
		
	end
end
		
% normalize to sum 1 for each time lag
% each element contains proportion of each contingency
% -----------------------------------------------------------

fM = fM / sum(sum(fM(:,:,1)));  	

dM = fM - eN;

% msd is the square root of the mean squared deviation from counterbalancing across all cells.
msd = sqrt(mean(mean(mean(dM.^2))));

cBal = 100*(1 - msd);

warning on

return



