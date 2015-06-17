% rundesignsim
clear;
figure;
index = 1;
name = input('Enter the name of the output variable and files ');
eval([name ' = zeros(10,10);'])
	
% set initial frequency of rests
	H = findobj('Tag', 'E2');  
	freqConditions = str2num(get(H(1), 'String'));
	freqConditions = num2str([freqConditions(1:4) .01]);
	set(H,'String',freqConditions);

	% set initial ISI
	H = findobj('Tag', 'E4'); % ISI
	ISI = str2num(get(H(1), 'String')); 
	ISI = num2str(1);
	set(H,'String',ISI);

myISI = 1:.2:16;	
myFC = .01:.05:.9;

for RR = 1:size(myISI,2)
	% set ISI
	H = findobj('Tag', 'E4'); % ISI
	ISI = str2num(get(H(1), 'String')); 
	ISI = num2str(myISI(RR));
	set(H(1),'String',ISI);
	
	for CC = 1:size(myFC,2)
		disp(['Starting model ' num2str(index)])		
	
		% set rest frequency
		H = findobj('Tag', 'E2');  
		freqConditions = str2num(get(H(1), 'String'));
		freqConditions = num2str([freqConditions(1:4) myFC(CC)]);
		set(H(1),'String',freqConditions);
		
		disp(['FREQC and ISI are:	' get(H(1),'String') ' and ' get(findobj('Tag','E4'),'String')])
		designsim_gui_script
	
		H = findobj('Tag','SIMplot'); 
		if ~isempty(H)
			figure(H)
			close
		end
		
		disp(['mean is ' num2str(mean(SIM.rnd.t))])
		eval([name '(RR,CC) = mean(SIM.rnd.t);']);
		eval([name '(1:10,1:10)'])
		eval(['save ' name num2str(index) ' ' name ' SIM'])
		index = index + 1;
	end
end

eval(['save ' name])
eval(['surf(' name ')']); colormap(copper)
eval(['set(gca,''XTick'',1:size(' name ',2))'])
set(gca,'XTickLabel',myFC)
eval(['set(gca,''YTick'',1:size(' name ',1))'])
set(gca,'YTickLabel',myISI)
xlabel('Proportion of blank intervals','FontSize',18)
ylabel('ITI in s','FontSize',18)
title('Surface map of ITI vs. proportion of blank intervals','FontSize',18)