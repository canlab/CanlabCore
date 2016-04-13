function [h,t,w] = mvroi_plot_firs(DATA,r)
% ::
%
%    [h,t,w] = mvroi_plot_firs(DATA,region)
%
% :Examples:
% ::
%
%    DATA.SPEC.firnames = {'Antic (C)'  'Pain (C)' 'Response (C)' 'Antic (P)' 'Pain (P)' 'Response (P)'};
%    DATA.SPEC.firconditions = [1 4]

c = DATA.SPEC.firconditions;
nms = DATA.SPEC.firnames;
regname = DATA.SPEC.names{r};
dat = DATA.DATA.fir(r,c);

colors = {'ro' 'go' 'bo' 'yo' 'co' 'mo' 'ko' 'r^' 'g^' 'b^' 'y^' 'c^' 'm^' 'k^'};

tor_fig(1,2); subplot(1,2,1); 


fprintf(1,'\n\nRegion %3.0f, %s\n----------------------------------\n',r,regname);


% smoothing and baseline subtraction

for i = 1:length(dat)
    for j = 1:size(dat{1},1)
    
        dat{i}(j,:) = smooth_timeseries(dat{i}(j,:),3);
        %dat{i}(j,:)  = dat{i}(j,:) - nanmean(dat{i}(j,1:2));
        
    end
end



% get average standard error for diffs within subs
for i = 2:length(dat)
    stes(i-1,:) = ste(dat{i} - dat{i-1});
    m = mean(dat{i} - dat{i-1});
end
m = mean(m,1);
s = mean(stes,1);
t = m ./ s; disp(['Max t is ' num2str(max(t))])
clear t



for i = 1:length(dat)
    
    m = nanmean(dat{i});
    %s = ste(dat{i});
    
    h(i) = plot(m,[colors{i} '-'],'LineWidth',2);
    
    fill_around_line(m,s,colors{i}(1));
    
end

legend(h,nms(c))
title(regname,'FontSize',24);

% get H, T, W

hconstraint = round(length(dat{i}(1,:)).*.66);

for i = 1:length(dat)
    for j = 1:size(dat{1},1)
    
        y = dat{i}(j,:);
        %figure;plot(dat{i}(j,:),'k','LineWidth',2);
        [h(j,i) t(j,i) w(j,i)] = fir2htw2(y,hconstraint,0);
        
    end
end

subplot(1,2,2);
barplot_columns([diff(h,1,2) diff(t,1,2) diff(w,1,2)],'Param Diffs (1-2)',[],'noind','nofig');

set(gca,'XTickLabel',{'H' 'T' 'W'});

%str = {};
%for i = 1:length(dat), str = [str {['H' num2str(i)]}];,end
%for i = 1:length(dat), str = [str {['T' num2str(i)]}];,end
%for i = 1:length(dat), str = [str {['W' num2str(i)]}];,end
%set(gca,'XTickLabel',str);



% text string for sig
m = mean(diff(h,1,2)); s = ste(diff(h,1,2)); t = m./s;
df = size(h,1) - 1; p = tcdf(t,df); p = 2*min(p,1 - p);

str = sprintf('height diff = %3.2f, t(%3.0f) = %3.2f, p = %3.4f', ...
   m, df, t, p);
    
title(str,'FontSize',14);


ssize = get(0,'screensize');
set(gcf,'Position',ssize ./ [0.0270    0.0029    2.0915    2.9147]);


return

