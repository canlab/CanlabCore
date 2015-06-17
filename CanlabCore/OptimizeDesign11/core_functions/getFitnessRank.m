function overallFitnessRank = getFitnessRank(fitnessMatrix, sizeGenerations,cbalColinPowerWeights) 

% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ % 

% GETFITNESSRANK: GET RANK ORDER OF INDIVIDUAL FITNESS SCORES

% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ % 

% Input: fitness matrix, with rows as fitness categories and cols as agents, entries = fitness scores

%			size of generations (# of agents), row vect of weights for importance of fitness categories.

% Output: row vector of overall fitness rank within that generation for each agent.

% rank order fitness to decide relative fitness for breeding

   tempFM = fitnessMatrix;				% do this to scale relative importance of values

   for h = 1:3  							% for Cbal, then Colin, then Power

   	for i = sizeGenerations:-1:1

      	myMax = max(tempFM(h,:));

      	for j = 1:sizeGenerations

         	if tempFM(h,j) == myMax & not(tempFM(h,j) == -1) 

               rankMatrix(h,j) = i;

            	tempFM(h,j) = -1;

         	end

      	end

      end

      % scale ranks by relative importance

      rankMatrix(h,:) = rankMatrix(h,:) * cbalColinPowerWeights(h);

   end

   overallFitnessRank = sum(rankMatrix,1); % sum the rows in each col

return

