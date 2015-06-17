% rundesignsimMAP2
% Tor Wager, 1/2/02
% if your input var name is already a variable in the wkspace, uses that...
% otherwise creates one with zeros.
% run ga_text_script_epochs and break first.

N = fieldnames(GA);
for i = 1:length(N)
    eval([N{i} ' = GA.' N{i} ';'])
end

warning off
index = 1;
name = input('Enter the name of the output variable and files ','s');
% eval([name ' = zeros(10,10);'])
listMatrix = [];

% -------------------------------------------------
% common user variables
% -------------------------------------------------
tickRes = [.2 .2];
ISIrange = [.1 16];     % used to increment hox8
FCrange = [.5 .5];
GA.sizegenerations = 100;
hox = [0 0 8 0 8 0 4 1];  % last one is variable in this script


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
% * Add intercept to contrasts
% ----------------------------------------------------------------   
contrasts(:,end+1) = 0;

% ----------------------------------------------------------------
% * make the map
% ----------------------------------------------------------------

for RR = 1:size(myISI,2)
	% set ISI
    hox(8) = myISI(RR);
	
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

            [rna,allmodels,alldelta,all_stimlist,numtrials] = dna2model_epochs(stimList,[],hox,2,[],480,240,0,0,0,0,HRF,2,240,S);
            model = allmodels{1}; 
            delta = alldelta{1};
            
            % efficiency
            xtxitx = pinv(model);                                       		% inv(X'S'SX)*(SX)'; pseudoinv of (S*X)
            if size(model,2) == size(contrasts,2)
                go = 1;
 			    coneff(z) = calcEfficiency(contrastweights,contrasts,xtxitx,svi);
            else 
                coneff(z) = NaN;
                warning on
                warning('HRF Model is wrong size'); whos model
                keyboard
                go = 0;
            end
            
            % hrf estimation
            [model] = tor_make_deconv_mtx2(delta,round(12 / TR),TR * 10);       % modified for .1 s sampling res delta
            
            if size(model,1) > numsamps, 
                warning on
                warning('Model longer than num samps - this should not happen!'),
                whos model
                numsamps
                model = model(1:numsamps,:);,
            end
            if ~isempty(S), model = S * model;,end
            
        	xtxitx = pinv(model);                                       		% inv(X'S'SX)*(SX)'; pseudoinv of (S*X)
			if go,
                hrfeff(z) = calcEfficiency([],[],xtxitx,svi);
            else
                hrfeff(z) = NaN;
                disp('HRF Model:'); whos model
            end
            
            %[fitness,convar,xtxi,models] = testlist(listMatrix,GA,HRF,'svi',svi,'S',S);
            %convar = 1./convar;
	        eval([name '(RR,CC) = nanmean(coneff);']);
            eval([name '_hrf(RR,CC) = nanmean(hrfeff);']);
            
            if mod(z,10) == 0, fprintf(1,'.'),end
        end % loop through subjects
		
        fprintf(1,'%3.1f\n',z)
        
    end % if mymat...
		%eval(['save ' name num2str(index) ' ' name ' SIM'])
		index = index + 1;
	end
end

eval(['save ' name])
figure;
set(gcf,'Color','w')
try 
    eval(['surf(' name ')']); colormap(copper)
catch
    eval(['plot(' name ')']); 
end
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

figure;
set(gcf,'Color','w')
try 
    eval(['surf(' name '_hrf)']); colormap(copper)
catch
    eval(['plot(' name '_hrf)']); 
end
eval(['set(gca,''XTick'',1:size(' name ',2))'])
set(gca,'XTickLabel',myFC)
%eval(['set(gca,''YTick'',1:size(' name ',1))'])
%set(gca,'YTickLabel',myISI)
xlabel('Proportion of blank intervals','FontSize',14)
ylabel('ISI in s','FontSize',14)
title('Surface map of ISI vs. proportion of blank intervals','FontSize',18)
zlabel('Efficiency of HRF estimation','FontSize',14)

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