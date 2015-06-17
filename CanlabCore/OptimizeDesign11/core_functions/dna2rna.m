function rna = dna2rna(listMatrix,restMatrix,restlength,numStim,dorests,trans2switch,trans2block,dofirst)
% function rna = dna2rna(listMatrix,restMatrix,restlength,numStim,dorests,trans2switch,trans2block,dofirst)
%
% transforms 'genetic' code used in crossbreeding (dna)
% into 'expressed' form (rna) used to build the design matrix 
%
% Tor Wager, 11/17/01

listMatrix = double(listMatrix);
restMatrix = double(restMatrix);

% before start of generation - ops on overall list
% ---------------------------------------------------------------------------------------------
if dorests,
		listMatrix = insert_rests(listMatrix,restMatrix,restlength,numStim);
end


% within generation - ops on each design vector
% ---------------------------------------------------------------------------------------------

for z = 1:size(listMatrix,2) 		% do this for each organism in the generation
	stimList = listMatrix(:,z);

    % The Switch Hack: transform list of stimuli into list of switches/no switches for each stimtype (R or L)
    if trans2switch, stimList = transform2switches(stimList);,end
    if trans2block, 
		  if isempty(restMatrix),error('trans2block will not work without some rest intervals specified.'),end
          stimList = transform2block(stimList,restMatrix(:,z),restlength,dofirst);
    end
	
	rna(:,z) = stimList;
end


return