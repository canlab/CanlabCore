% rundesignsimMAP2
% Tor Wager, 8/29/01
% if your input var name is already a variable in the wkspace, uses that...
% otherwise creates one with zeros.

warning off
index = 1;
name = input('Enter the name of the output variable and files ','s');
% eval([name ' = zeros(10,10);'])

% -------------------------------------------------
% common user variables
% -------------------------------------------------
tickRes = [.1 .1];
ISIrange = [.1 16]
FCrange = [.01 .81];
GA.sizegenerations = 100;

myISI = ISIrange(1):tickRes(1):ISIrange(2);	
myFC = FCrange(1):tickRes(2):FCrange(2);

ISIsize = length(myISI);
FCsize = length(myFC);

ISItick = 1:max(round(ISIsize./25),1):ISIsize;
ISIticklab = myISI(ISItick);
FCtick = 1:2:FCsize;
FCticklab = myFC(FCtick);

eval([' if ~(exist(''' name ''') == 1), disp(''Creating new variable''),' name ' = zeros(ISIsize,FCsize);, else,disp(''Using existing variable''),end'])

disp(['Starting ' num2str(prod([ISIsize FCsize])) ' surface points.'])

% ----------------------------------------------------------------
% * HRF GAMMA FUNCTION - spm 99
% ----------------------------------------------------------------

   HRF = spm_hrf(.1);						% SPM HRF sampled at .1 s
   HRF = HRF/ max(HRF);
   


% ----------------------------------------------------------------
% * make the map
% ----------------------------------------------------------------

for RR = 1:size(myISI,2)
	% set ISI
    GA.ISI = myISI(RR);
	
    % ----------------------------------------------------------------
    % * get smoothing matrix and autocorrelation matrix
    % ----------------------------------------------------------------
    numStim = ceil(GA.scanLength / (GA.ISI));
    numsamps = ceil(numStim*GA.ISI/GA.TR);
	disp('  ...getting S and svi')
	[S,Vi,svi] = getSmoothing(GA.HPlength,GA.LPsmooth,GA.TR,numsamps,GA.xc);
	if isempty(S), S = 1;, end

   	clear Vi


	for CC = 1:size(myFC,2)
		disp(['Starting model ' num2str(index)])		
	
		% set rest frequency
		GA.freqConditions(end) = myFC(CC);
		GA.freqConditions(1:end-1) = (1 - myFC(CC)) / (size(GA.freqConditions,2)-1);

		disp(['FREQC and ISI are:	' num2str(GA.freqConditions) ' and ' num2str(GA.ISI)])
	
        eval(['mymat = ' name ';'])
        if mymat(RR,CC) == 0        % only do it if it's not done yet.
            
        % ----------------------------------------------------------------
        % * randomize initial set of organisms
        % ----------------------------------------------------------------
        unsortedList = [];
        numStimEachCond = ceil((numStim) * GA.freqConditions);
	    clear listMatrix
    
	    for i = 1:size(GA.freqConditions,2)   									
	    	% unsorted list of all stims in proper frequencies
   		    startIndex = size(unsortedList,1)+1;
		    % multiplied by 1.5 because lists were coming out too short.
   		    unsortedList(startIndex:(startIndex + 2*numStimEachCond(i)-1),1) = GA.conditions(i);
        end


        disp('  ...randomizing organism start state')
        for z = 1:GA.sizeGenerations % generate random ordering of x conditions
  	        stimList = getRandom(unsortedList);
            listMatrix(:,z) = stimList; % a row for each subject
        end


            [fitness,convar,xtxi,models] = testlist(listMatrix,GA,HRF,'svi',svi,'S',S);
            convar = 1./convar;
	        eval([name '(RR,CC) = mean(convar,2);']);
        end % if mymat...
		
		%eval(['save ' name num2str(index) ' ' name ' SIM'])
		index = index + 1;
	end
end

eval(['save ' name])
figure;
set(gcf,'Color','w')
eval(['surf(' name ')']); colormap(copper)
eval(['set(gca,''XTick'',1:size(' name ',2))'])
set(gca,'XTickLabel',myFC)
%eval(['set(gca,''YTick'',1:size(' name ',1))'])
%set(gca,'YTickLabel',myISI)
xlabel('Proportion of blank intervals','FontSize',14)
ylabel('ISI in s','FontSize',14)
title('Surface map of ISI vs. proportion of blank intervals','FontSize',18)
zlabel('Efficiency of contrast [1 -1]','FontSize',14)

set(gca,'YTick',ISItick)
set(gca,'YTickLabel',ISIticklab)

warning on

break
%for line plots
figure;plot(MAP1)
plot(MAP1,'LineWidth',2)
legend({'.01' '.21' '.41' '.61' '.81'})
set(gca,'XTick',1:2:40)
set(gca,'XTickLabels',.1:.4:7.9)
title('Efficiency as a function of ISI and jitter: Linear HRF function','FontSize',14)
grid on

eval(['saveas(gcf,''' name ''',''fig'')'])
eval(['saveas(gcf,''' name ''',''jpg'')'])