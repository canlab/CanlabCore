function [listMat,restMat] = breedorgs2(ovf,listMat,alph,varargin)
% function [listMat,restMat] = breedorgs2(ovf,listMat,alph,mutations [opt], restMatrix [optional])
%
% Input: ovf = row vector of fitness z-scores for each agent, matrix of condition lists
%			(cols = agents, rows = stim. presentations)
%			gamma = bias for probabilistic breeding - 1 = no bias?  0 = always choose best
% Output: new listMat of best organisms interbred.
%
% Notes: swaps restMatrix (varargin) along with stimList.  in theory, this could be
%			anything else you want to optimize - i.e., another 'chromosome'.
%
% by Tor Wager, last edit 3/30/01 to implement swapping rest intervals (restMatrix).

if nargin > 3, mutations = varargin{1};, end
if nargin > 4, mutations = varargin{1};,restMat = varargin{2};,swaprests = 1;,else swaprests = 0;,end


% keep the best one to ensure hill climbing
bestOne = listMat(:,(ovf == max(ovf)));
bestOne = bestOne(:,1);					    % in case of multiple matches to max
if swaprests, bestRest = restMat(:,(ovf == max(ovf)));, bestRest = bestRest(:,1);, end

numOrgs = size(listMat,2);
numStim = size(listMat,1);
if nargin < 4, mutations = ceil(numStim*numOrgs / 100);, end 

% ====================================================================
% selection of top half: choices are cross,tourney,sigmoid
% breeding: choices are probabilistic,neighbor,pickrandom
% ====================================================================

% Select top half
% ====================================================================
%[listMat,ovf] = tourney(ovf,listMat,numOrgs);
% [listMat,ovf] = besthalf(ovf,listMat,numOrgs);
if swaprests
	[listMat,restMat] = sigmoid(ovf,listMat,numOrgs,alph,restMat);
else
	[listMat] = sigmoid(ovf,listMat,numOrgs,alph);
end



% Crossover top half
% ====================================================================
% listMat = probbreed(ovf,listMat,numOrgs,numStim,gam);
%listMat = neighbor(listMat,numOrgs,numStim);
if swaprests
	[listMat,restMat] = pickrandom(listMat,numOrgs,numStim,restMat);
else
	listMat = pickrandom(listMat,numOrgs,numStim);
end

% Mutation
% ====================================================================
% .001% chance of swap mutation across whole population
for z = 1:mutations
   a = ceil(rand*numStim);
   b = ceil(rand*numStim);
   org = ceil(rand*numOrgs);
   buffer = listMat(a,org);
 	listMat(a,org) = listMat(b,org);
	listMat(b,org) = buffer;
end

% reinsert best one - spare from mutation and crossover, but still allow to breed.
% reserve slot 1 so that multiple copies of best one don't proliferate if ga fails to improve.
listMat(:,1) = bestOne;
if swaprests,restMat(:,1) = bestRest;,end
return



% ================================================================================== %
% SUB-FUNCTIONS
% ================================================================================== %
function [newMat,newovf] = tourney(ovf,listMat,numOrgs)
% ovf,listMat,numOrgs
for n = 1:numOrgs
   % Tournament: spot 1 is reserved for best player from previous generation
   picker = rand(1,numOrgs);
   if ovf(picker == max(picker)) > ovf(picker == min(picker))
      newMat(:,n) = listMat(:,(picker == max(picker)));
      newovf(1,n) = ovf(1,(picker == max(picker)));
   else
      newMat(:,n) = listMat(:,(picker == min(picker)));
      newovf(1,n) = ovf(1,(picker == min(picker)));
   end
end
return

% ================================================================================== %
function [listMat,ovf,newRestMat] = besthalf(ovf,listMat,numOrgs,varargin)
    if nargin > 3, restMat = varargin{1};,swaprests = 1;,else swaprests = 0;,end
    listMat(:,(ovf < median(ovf))) = [];
    ovf(:,(ovf < median(ovf))) = [];
    listMat = [listMat listMat];
    ovf = [ovf ovf];
    listMat = listMat(:,1:numOrgs);
    ovf = ovf(:,1:numOrgs);
    
    if swaprests
		restMat(:,(ovf < median(ovf))) = [];
    	restMat = [restMat restMat];
    	newRestMat = restMat(:,1:numOrgs);
	end
    
return

% ================================================================================== %
function newMat = probbreed(ovf,listMat,numOrgs,numStim,gam)
% probabilistic breeding
% gamma,ovf,listMat
% gamm = .995;        % 1/rank, adjusted by sigmoid filter based on distance from median.
sigmoidadjust = 1 ./(ovf - gam*(abs((ovf-median(ovf)/median(ovf)))));
newMat = zeros(numStim,numOrgs);

for i = 1:2:numOrgs-1
    fitProbs = rand(1,numOrgs) .* sigmoidadjust; 
    a = listMat(:,fitProbs == max(fitProbs));
    fitProbs(1,fitProbs == max(fitProbs)) = 0;
    b = listMat(:,fitProbs == max(fitProbs));

   % 2 children for each set of parents
   crossPoint = floor(rand * numStim);	                                 % set crossover point
   if crossPoint == 0,crossPoint = 1;,end;
   newMat(1:crossPoint,i) = a(1:crossPoint,1);
   newMat(crossPoint+1:end,i) = b(crossPoint+1:end,1);                   % swap chromosome pieces
   newMat(1:crossPoint,i+1) = b(1:crossPoint,1);
   newMat(crossPoint+1:end,i+1) = a(crossPoint+1:end,1);       
end
return

% ================================================================================== %
function [newMat,newRestMat] = sigmoid(ovf,listMat,numOrgs,alph,varargin)
 	
	if nargin > 4, restMat = varargin{1};,swaprests = 1;,else swaprests = 0;,end
	scaleovf = (5/max(ovf)) * (ovf-mean(ovf));                                            % scale ovf to +/- 5, so sigmoid has an effect
	y = rand(1,numOrgs) + 1./(1 + (2.371.^(-alph*scaleovf)));
    
	% testing the choosing, to see plots.
	%figure; subplot(3,1,1);plot(ovf,'r');title('original fitness')   
	%subplot(3,1,2); plot(scaleovf,'b'); title('squashed to +- 5')
	%subplot(3,1,3);hold on;plot(ovf,y,'ro'); plot([min(ovf) max(ovf)],[median(y) median(y)],'k');title('choosing score (y) vs fitness (x)')
	%chosen = (y >= median(y))  
	%corrcoef(ovf,y)
	%subplot(3,1,1); hold on; for i = 1:size(chosen,2),if chosen(i) == 1,plot(i,ovf(i),'yd','MarkerFaceColor','y'),end,end
	
	listMat(:,(y < median(y))) = [];
    listMat = [listMat listMat];
    newMat = listMat(:,1:numOrgs);
	if swaprests
		restMat(:,(y < median(y))) = [];
    	restMat = [restMat restMat];
    	newRestMat = restMat(:,1:numOrgs);
	end
return

% ================================================================================== %
function listMat = neighbor(listMat,numOrgs,numStim)
% pair off and swap upper halves
for z = 1:2:numOrgs - 1										% pair off
   crossPoint = floor(rand * numStim);	% set crossover point
   if crossPoint == 0,crossPoint = 1;,end;
   % disp(['switching cols ' num2str(z) ' and ' num2str(z+1)])
   a = listMat(1:crossPoint,z);
   b = listMat(1:crossPoint,z+1); % swap chromosome pieces
   listMat(1:crossPoint,z) = b;
   listMat(1:crossPoint,z+1) = a;
end
return

% ================================================================================== %
function [newMat,newRestMat] = pickrandom(listMat,numOrgs,numStim,varargin)

if nargin > 3, restMat = varargin{1};,swaprests = 1;,else swaprests = 0;,end

% random without replacement.  
newMat = zeros(numStim,numOrgs);
if swaprests,newRestMat = zeros(size(restMat,1),numOrgs);,end
fitProbs = rand(1,numOrgs); 

for i = 1:2:numOrgs-1
    a = listMat(:,fitProbs == max(fitProbs));
    if swaprests, c = restMat(:,fitProbs == max(fitProbs));,end
	fitProbs(1,fitProbs == max(fitProbs)) = 0;
    
	b = listMat(:,fitProbs == max(fitProbs));
 	if swaprests, d = restMat(:,fitProbs == max(fitProbs));,end	
	fitProbs(1,fitProbs == max(fitProbs)) = 0;   
		
   % 2 children for each set of parents
   crossPoint = floor(rand * numStim);	                                 % set crossover point
   if crossPoint == 0,crossPoint = 1;,end;
   newMat(1:crossPoint,i) = a(1:crossPoint,1);
   newMat(crossPoint+1:end,i) = b(crossPoint+1:end,1);          % swap chromosome pieces
   newMat(1:crossPoint,i+1) = b(1:crossPoint,1);
   newMat(crossPoint+1:end,i+1) = a(crossPoint+1:end,1);       
   
   	if swaprests
		crossPoint = floor(rand * size(restMat,1));	                                 % set crossover point
   		if crossPoint == 0,crossPoint = 1;,end;
   		newRestMat(1:crossPoint,i) = c(1:crossPoint,1);
   		newRestMat(crossPoint+1:end,i) = d(crossPoint+1:end,1);          % swap chromosome pieces
   		newRestMat(1:crossPoint,i+1) = d(1:crossPoint,1);
   		newRestMat(crossPoint+1:end,i+1) = c(crossPoint+1:end,1);       
	end

end
return

% ================================================================================== %
function test_dummy
index = 1;
for gamma = .8:.05:1.15
    for i = 1:1000
        a(i*10-9:i*10,1) = breedorgs([1 2 3 4 5 6 7 8 9 10],[1 2 3 4 5 6 7 8 9 10],gamma)';
    end
    subplot(2,4,index)
    hist(a,1:10)
    title(['Gamma = ' num2str(gamma)]) 
    index = index+1;
end
figure
index = 1;
for gamma = .9:.01:1
    for i = 1:1000
        a(i*10-9:i*10,1) = breedorgs([1 2 3 4 5 6 7 8 9 10],[1 2 3 4 5 6 7 8 9 10],gamma)';
    end
    subplot(3,4,index)
    hist(a,1:10)
    title(['Gamma = ' num2str(gamma)]) 
    index = index+1;
end

% sigmoid function
for a = 0:.1:1
y = 1./(1 + (2.371.^(-a*x))); hold on;plot(x,y)
end
return