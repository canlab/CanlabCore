function [overallFitness,mostFit,f] = getOverallFitness(cbalColinPowerWeights,fitnessMatrix)
% function [overallFitness,mostFit,f] = getOverallFitness(cbalColinPowerWeights,fitnessMatrix)
%
%
% Tor Wager, 11/17/01

   	% Determine overall fitness
   	for i = 1:size(cbalColinPowerWeights,2)				% convert to z scores
      if not(sum(fitnessMatrix(i,:)) == 0)              % only for rows with > 0 weight
         if std(fitnessMatrix(i,:)) == 0
            % disp(['Warning: All fitness scores for row ' num2str(i) ' are the same. Sum ' num2str(sum(fitnessMatrix(i,:)))])
            zfitnessMatrix(i,:) = (fitnessMatrix(i,:) - mean(fitnessMatrix(i,:))) / 1;
         else
      		zfitnessMatrix(i,:) = (fitnessMatrix(i,:) - mean(fitnessMatrix(i,:))) / std(fitnessMatrix(i,:));      
         end
      end
   	end
   zfitnessMatrix((cbalColinPowerWeights == 0),:) = 0;
   overallFitness = (cbalColinPowerWeights * zfitnessMatrix) / sum(cbalColinPowerWeights);  % weighted avg sum z-scores 
   
   % determine the most fit organism overall
   mostFit = find(overallFitness == max(overallFitness));
   
   % avoid duplicates
   mostFit = mostFit(1);	
   
   % return fitness scores for most fit
   f = fitnessMatrix(:,mostFit);
   
return
   