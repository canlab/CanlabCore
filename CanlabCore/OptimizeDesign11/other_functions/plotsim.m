function SIM = plotSim(SIM,varargin)	
% function SIM = plotSim(SIM,'add' [opt],'color',[.4 .3])	
% plotting the SIM model
% ============================================================================================
% use OptimizeDesign scripts to create a SIM struct, then plot it using this.
% 'add' adds to current plot instead of making new figure
% 'color', followed by 4 numbers btwn 0 and 1, modifies color scheme.
%		colors are: 1-2 = [r g] for random.  3-4 = [g b] for ga.
% 'nonoise' : suppress noise text.
% last modified 5/21/01 Tor Wager
	
	addtoplot = 0; 
	noisetext = 1;
	color = [.2 .3 .5 .5];
	for i = 1:nargin-1
		if isstr(varargin{i})
			switch varargin{i}
			case 'add', addtoplot = 1;
			case 'color',color = varargin{i+1};	
			case 'nonoise',noisetext = 0;	
			end
		end
	end
	
	divisor = 12;
	if SIM.niterations < 10, divisor = 1;,warning('Low number of iterations (< 10).'),end
	
   % determine bins
   tempvar = [];
   for i = 1:size(SIM.rnd.t,1), tempvar = [tempvar SIM.rnd.t(i,:)];,end
   for i = 1:size(SIM.rnd.HPt,1), tempvar = [tempvar SIM.rnd.HPt(i,:)];,end
   for i = 1:size(SIM.rnd.HLt,1), tempvar = [tempvar SIM.rnd.HLt(i,:)];,end
   if isfield(SIM,'ga')
     for i = 1:size(SIM.ga.t,1), tempvar = [tempvar SIM.ga.t(i,:)];,end
     for i = 1:size(SIM.ga.HPt,1), tempvar = [tempvar SIM.ga.HPt(i,:)];,end
     for i = 1:size(SIM.ga.HLt,1), tempvar = [tempvar SIM.ga.HLt(i,:)];,end
   end
   [temphist,SIM.x] = hist(tempvar,round(SIM.niterations/divisor));
   
   tempvar = [];
   for i = 1:size(SIM.rnd.se,1), tempvar = [tempvar SIM.rnd.se(i,:)];,end
   for i = 1:size(SIM.rnd.HPse,1), tempvar = [tempvar SIM.rnd.HPse(i,:)];,end
   for i = 1:size(SIM.rnd.HLse,1), tempvar = [tempvar SIM.rnd.HLse(i,:)];,end
   [temphist,SIM.sex] = hist(tempvar,round(SIM.niterations/divisor));
   
   % No smoothing
   % ===============================================================================================
   % ==== plot no smoothing, t =======
   if ~addtoplot,figure; set(gcf,'Tag','SIMplot'),else H = findobj('Tag','SIMplot'); figure(H),end
   subplot(3,2,1);hold on; 
   for i = 1:size(SIM.rnd.t,1)
        SIM.rnd.t_hist(i,:) = hist(SIM.rnd.t(i,:),SIM.x);
        plot(SIM.x,SIM.rnd.t_hist(i,:),'Color',[color(1) color(2) 1/i],'LineWidth',2) % 'LineStyle','--')
        ylim = get(gca,'YLim');
        plot([mean(SIM.rnd.t(i,:)) mean(SIM.rnd.t(i,:))],[0 ylim(2)],'b');
  		if isfield(SIM,'ga')   
                        ylim = get(gca,'YLim');   
			plot([mean(SIM.ga.t(i,:)) mean(SIM.ga.t(i,:))],[0 ylim(2)],'r');
        	[histo,histox] = hist(SIM.ga.t(i,:),SIM.x);
        	plot(histox,histo,'Color',[1/i color(3) color(4)],'LineWidth',2)
		end
   end
   title('t results: no smoothing.')

   % ==== plot no smoothing, r =======
   subplot(3,2,2); hold on;
   %hist(SIM.rnd.se,[min(SIM.rnd.se):(max(SIM.rnd.se)-min(SIM.rnd.se))*10/SIM.niterations:max(SIM.rnd.se)]);
   	if size(SIM.rnd.se,1) > 1
   		for i = 1:size(SIM.rnd.se,1)
   			plot(hist(SIM.rnd.se(i,:),SIM.sex),'Color',[color(1) color(2) 1/i],'LineWidth',2);
		end
	else
		hist(SIM.rnd.se,SIM.sex)
	end
   SIM.rnd.sehist = hist(SIM.rnd.se,SIM.sex);
   	if isfield(SIM,'ga') 
                ylim = get(gca,'YLim');
   		plot([SIM.ga.se SIM.ga.se],[0 ylim(2)],'r'); 
	end

   title('Contrast variance: no smoothing.')

   if noisetext,text(.7,max(SIM.rnd.sehist(:,1))+7,['noise var. = ' num2str(SIM.noise_var)]);,end

   
   % High-pass
   % ===============================================================================================
   % determine bins
   %tempvar = [];
   %for i = 1:size(SIM.rnd.HPt,1), tempvar = [tempvar SIM.rnd.HPt(i,:)];,end
   %[temphist,SIM.x] = hist(tempvar,round(SIM.niterations/10));
   
   % ==== plot high-pass, t =======
   subplot(3,2,3);hold on;
   for i = 1:size(SIM.rnd.t,1)
        SIM.rnd.HPt_hist(i,:) = hist(SIM.rnd.HPt(i,:),SIM.x);
        plot(SIM.x,SIM.rnd.HPt_hist(i,:),'Color',[color(1) color(2) 1/i],'LineWidth',2)
        ylim = get(gca,'YLim');
        plot([mean(SIM.rnd.HPt(i,:)) mean(SIM.rnd.HPt(i,:))],[0 ylim(2)],'b');
  		if isfield(SIM,'ga')         
                        ylim = get(gca,'YLim');
			plot([mean(SIM.ga.HPt(i,:)) mean(SIM.ga.HPt(i,:))],[0 ylim(2)],'r');
        	[histo,histox] = hist(SIM.ga.HPt(i,:),SIM.x);
        	plot(histox,histo,'Color',[1/i color(3) color(4)],'LineWidth',2)
		end		
	end
   title('t results: high-pass filtered.')
   
   % ==== plot high-pass, r =======
   subplot(3,2,4); hold on;
   %hist(SIM.rnd.HPse,[min(SIM.rnd.HPse):(max(SIM.rnd.HPse)-min(SIM.rnd.HPse))*10/SIM.niterations:max(SIM.rnd.HPse)]);
   if size(SIM.rnd.HPse,1) > 1
   		for i = 1:size(SIM.rnd.HPse,1)
   			plot(hist(SIM.rnd.HPse(i,:),SIM.sex),'Color',[color(1) color(2) 1/i],'LineWidth',2);
		end
	else
		hist(SIM.rnd.HPse,SIM.sex)
	end
   SIM.rnd.HPsehist = hist(SIM.rnd.HPse,[min(SIM.rnd.HPse):(max(SIM.rnd.HPse)-min(SIM.rnd.HPse))*10/SIM.niterations:max(SIM.rnd.HPse)]);
   if isfield(SIM,'ga')  
           ylim = get(gca,'YLim');
   	   plot([SIM.ga.HPse SIM.ga.HPse],[0 ylim(2)],'r'); 
   end
   title('Contrast variance: high-pass filtered.')
    if noisetext,text(.7,max(SIM.rnd.sehist(:,1))+7,['noise var. = ' num2str(SIM.noise_var)]);,end
   
   
   
   
   % Low-pass
   % ===============================================================================================
   
   % determine bins
   %tempvar = [];
   %for i = 1:size(SIM.rnd.HLt,1), tempvar = [tempvar SIM.rnd.HLt(i,:)];,end
   %[temphist,SIM.x] = hist(tempvar,round(SIM.niterations/10));
   
   % ==== plot smoothed, t =======
   subplot(3,2,5);hold on;
   for i = 1:size(SIM.rnd.t,1)
        SIM.rnd.HLt_hist(i,:) = hist(SIM.rnd.HLt(i,:),SIM.x);
        plot(SIM.x,SIM.rnd.HLt_hist(i,:),'Color',[color(1) color(2) 1/i],'LineWidth',2)
        ylim = get(gca,'YLim');
        plot([mean(SIM.rnd.HLt(i,:)) mean(SIM.rnd.HLt(i,:))],[0 ylim(2)],'b');
  		if isfield(SIM,'ga')         
                        ylim = get(gca,'YLim');
			plot([mean(SIM.ga.HLt(i,:)) mean(SIM.ga.HLt(i,:))],[0 ylim(2)],'r');
        	[histo,histox] = hist(SIM.ga.HLt(i,:),SIM.x);
        	plot(histox,histo,'Color',[1/i color(3) color(4)],'LineWidth',2)
   		end 
	end
   title('t results: HP and LP smoothed.')
   
   % ==== plot smoothed, r =======
   subplot(3,2,6); hold on;
   %hist(SIM.rnd.HLse,[min(SIM.rnd.HLse):(max(SIM.rnd.HLse)-min(SIM.rnd.HLse))*10/SIM.niterations:max(SIM.rnd.HLse)]);
   if size(SIM.rnd.HLse,1) > 1
		for i = 1:size(SIM.rnd.HLse,1)
   			plot(hist(SIM.rnd.HLse(i,:),SIM.sex),'Color',[color(1) color(2) 1/i],'LineWidth',2);
		end
	else
		hist(SIM.rnd.HLse,SIM.sex)
	end
   %SIM.rnd.HLsehist = hist(SIM.rnd.HLse,[min(SIM.rnd.HLse):(max(SIM.rnd.HLse)-min(SIM.rnd.HLse))*10/SIM.niterations:max(SIM.rnd.HLse)]);
   if isfield(SIM,'ga')   
           ylim = get(gca,'YLim');
   	   plot([SIM.ga.HLse SIM.ga.HLse],[0 ylim(2)],'r');
   end
   title('Contrast variance: HP and LP smoothed')
    if noisetext,text(.7,max(SIM.rnd.sehist(:,1))+7,['noise var. = ' num2str(SIM.noise_var)]);,end
  
 
% ===================================
% * power plot 
% ===================================

plotfields = {'tpower' 'HPtpower' 'HLtpower' 'tpower' 'HPtpower' 'HLtpower'};
myxticklabel = {'  ' 'GA' 'random' 'GA - HP filter' 'random - HP' 'GA - smoothed' 'random - smoothed'};
doplot = 0;maxpower = 0;
for i = plotfields
	if isfield(SIM.rnd,i{1}) | isfield(SIM.ga,i{1}),doplot = 1;,end
end

if doplot
	figure; hold on; index = 1;
	for i = {'SIM.ga.tpower' 'SIM.rnd.tpower' 'SIM.ga.HPtpower' 'SIM.rnd.HPtpower' 'SIM.ga.HLtpower' 'SIM.rnd.HLtpower'}
		if mod(index,2) == 0,pcolor = 'b';, else pcolor = 'r';,end
		eval(['if isfield(SIM.ga,''' plotfields{index} ''') & isfield(SIM.rnd,''' plotfields{index} '''),plotthis = 1;,else plotthis = 0;,end'])
		if plotthis
			eval(['fill([index-.4 index-.4 index+.4 index+.4],[0 ' i{1} ' ' i{1} ' 0],pcolor)'])
			eval(['if ' i{1} ' > maxpower, maxpower = ' i{1} ';,end'])
		end
		index = index + 1;
	end
	set(gca,'YLim',[0 maxpower + .2])
	set(gca,'XTickLabel',myxticklabel)
	set(gcf,'Position',[51   190   878   408])
end
		

return





