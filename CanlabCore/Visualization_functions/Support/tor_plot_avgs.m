function tor_plot_avgs(AVG,STE,varargin)
% function tor_plot_avgs(AVG,STE,varargin)
%
% ------ plot the selective average --------
% Tor Wager, 10/25/01


plotcolor = {'ro-' 'r^--' 'bo-' 'b^--' 'k^--' 'ko-' 'g^--' 'go-' 'm^--' 'mo-'  };
steplot = 1;
for i = 1:length(AVG)
    mylegend{i} = ['Condition ' num2str(i)];
end
mytitle = 'Selective Average';

if nargin > 2
    Op = varargin{1};
    if isfield(Op,'nosteplot'),if Op.nosteplot == 1, steplot = 0;,end, end
    if isfield(Op,'colors'),plotcolor = Op.colors;,end
    if isfield(Op,'legend'),mylegend = Op.legend;,end
    if isfield(Op,'title'),mytitle = Op.title;,end
end

cla
hold on

	for g = 1:length(AVG)

            avg = AVG{g};
            
            % set x axis values
            % ---------------------------------------------------------------
            myx = 1:size(avg,2);
            if isfield(Op,'window'),myx = Op.window(1):Op.window(2);,end
            if length(myx) ~= size(avg,2),
                warning('Window is wrong length, using default')
                myx = 1:size(avg,2);
            end
                
      		try
			plot(myx,avg,plotcolor{g},'LineWidth',2)
		catch
			plotcolor{g} = 'k';
			plot(myx,avg,plotcolor{g},'LineWidth',2)
		end
      		
      		wid = get(gca,'XLim');
      		wid = (wid(2) - wid(1))/50;
    end
    
    
if steplot
    for g = 1:length(AVG)
            avg = AVG{g};
            ste = STE{g};
            
            % set x axis values
            % ---------------------------------------------------------------
            myx = 1:size(avg,2);
            if isfield(Op,'window'),myx = Op.window(1):Op.window(2);,end
            if length(myx) ~= size(avg,2),
                warning('Window is wrong length, using default')
                myx = 1:size(avg,2);
            end
            
      			for i = 1:size(avg,2)
                    j = myx(i);
         			try
                        plot([j j],[avg(i)-ste(i) avg(i)+ste(i)],plotcolor{g}(1))
         			    plot([j-wid j+wid],[avg(i)-ste(i) avg(i)-ste(i)],plotcolor{g}(1))
         			    plot([j-wid j+wid],[avg(i)+ste(i) avg(i)+ste(i)],plotcolor{g}(1))
                    catch
                        disp('Can''t plot error bar.')
                        i
                        avg
                        ste
                        plotcolor
                    end
      			end
    end
end 
 
grid on
title(mytitle,'FontSize',14)
legend(mylegend,0)

return
    