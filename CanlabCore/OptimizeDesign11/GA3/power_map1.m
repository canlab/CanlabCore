function power_map1

is = 0:.01:.99;
js = 10:10:600;
ks = 5:5:60;

 jind = 1;
 for j = js,
   iind = 1; 
   for i = is
       x = mvnrnd(zeros(j,3),[1 i i; i 1 i; i i 1]); x(:,end+1) = 1;
       [t80_ind(iind,jind),t80_grp(iind,jind),OUT] = xpower(x,[1 0 0],[1 0 0],12); all_n(iind,jind) = OUT.n;
       iind = iind + 1;
   end
   jind = jind + 1;
 end


figure('Color','w');subplot(2,2,1); set(gca,'FontSize',14)
surf(js,is,t80_ind,'EdgeColor','none'); 
ylabel('Correlation between predictors'),xlabel('N per individual'),view(-160,24)
zlabel('Noncentrality parameter')
subplot(2,2,2); hold on; set(gca,'FontSize',14)
plot(js,t80_ind(1,:)),plot(js,t80_ind(40,:),'g'),plot(js,t80_ind(100,:),'r')
xlabel('N per individual'),legend({'r = 0','r = .40','r = .99'})
ylabel('Noncentrality parameter')
mp = [.05 .01 .001 .0001 .00001 .000001];
mt = norminv(1-mp); mt = [mt;mt];
plot(repmat([0 max(js)+1],size(mt,2),1)',mt,'k')
for i = 1:length(mp),text(js(end-1),mt(1,i)+.2,num2str(mp(i))),end
title('Power in individual analysis')

subplot(2,2,3); set(gca,'FontSize',14)
surf(js,is,t80_grp,'EdgeColor','none'); 
ylabel('Correlation between predictors'),xlabel('N per individual'),view(-160,24)
zlabel('Noncentrality parameter')
subplot(2,2,4); hold on; set(gca,'FontSize',14)
plot(js,t80_grp(1,:)),plot(js,t80_grp(40,:),'g'),plot(js,t80_grp(100,:),'r')
xlabel('N per individual'),legend({'r = 0','r = .40','r = .99'})
ylabel('Noncentrality parameter')
mp = [.05 .01 .001 .0001 .00001 .000001];
mt = norminv(1-mp); mt = [mt;mt];
plot(repmat([0 max(js)+1],size(mt,2),1)',mt,'k')
for i = 1:length(mp),text(js(end-1),mt(1,i)+.2,num2str(mp(i))),end
title('Power in group analysis, N = 12')

saveas(gcf,'power_by_correl_and_Nind','fig')

 jind = 1;
 for j = js
   iind = 1; 
   for i = ks
       x = mvnrnd(zeros(j,3),[1 .3 .3; .3 1 .3; .3 .3 1]); x(:,end+1) = 1;
       [t80_ind2(iind,jind),t80_grp2(iind,jind),OUT] = xpower(x,[1 0 0],[1 0 0],i); all_n2(iind,jind) = OUT.n;
       iind = iind + 1;
   end
   jind = jind + 1;
 end
 
 
figure('Color','w');subplot(2,2,1); set(gca,'FontSize',14)
surf(js,ks,t80_ind2,'EdgeColor','none'); 
ylabel('N per group'),xlabel('N per individual'),view(-160,24)
zlabel('Noncentrality parameter')
subplot(2,2,2); hold on; set(gca,'FontSize',14)
plot(js,t80_ind2(1,:)),plot(js,t80_ind2(6,:),'g'),plot(js,t80_ind2(12,:),'r')
xlabel('N per individual'),legend({'r = 5','r = 30','r = 60'})
ylabel('Noncentrality parameter')
mp = [.05 .01 .001 .0001 .00001 .000001];
mt = norminv(1-mp); mt = [mt;mt];
plot(repmat([0 max(js)+1],size(mt,2),1)',mt,'k')
for i = 1:length(mp),text(js(end-1),mt(1,i)+.2,num2str(mp(i))),end
title('Power in individual analysis')

subplot(2,2,3); set(gca,'FontSize',14)
surf(js,ks,t80_grp2,'EdgeColor','none'); 
ylabel('N per group'),xlabel('N per individual'),view(-160,24)
zlabel('Noncentrality parameter')
subplot(2,2,4); hold on; set(gca,'FontSize',14)
plot(js,t80_grp2(1,:)),plot(js,t80_grp2(6,:),'g'),plot(js,t80_grp2(12,:),'r')
xlabel('N per individual'),legend({'r = 5','r = 30','r = 60'})
ylabel('Noncentrality parameter')
mp = [.05 .01 .001 .0001 .00001 .000001];
mt = norminv(1-mp); mt = [mt;mt];
plot(repmat([0 max(js)+1],size(mt,2),1)',mt,'k')
for i = 1:length(mp),text(js(end-1),mt(1,i)+.2,num2str(mp(i))),end
title('Power in group analysis')

saveas(gcf,'Power_by_ind_grp_N','fig')

save power_map1_out


return


 jind = 1;
 for j = 5:30,
   iind = 1; 
   for i = 1:20
       x = mvnrnd(zeros(j,3),[1 0 0; 0 1 0; 0 0 1]); x(:,end+1) = 1;
       [t802(iind,jind),t80_g(iind,jind),OUT] = xpower(x,[1 0 0],[1 0 0],12); alln(iind,jind) = OUT.n;
       
       y = x(:,1) + randn(size(x,1),1);
       [b,dev,stat] = glmfit(x(:,1:3),y);
       bet1(iind,jind) = b(2);
       tstat1(iind,jind) = stat.t(2);
       p1(iind,jind) = log(stat.p(2));
       
       iind = iind + 1;
   end
   jind = jind + 1;
 end

 
 
