clear betas1, clear betas0, clear b t s se p

ind = 1;
for i = 0:.01:.99
    
    for j = 1:100
        
        x = mvnrnd(zeros(30,3),[1 i i; i 1 i; i i 1]);
    
        x(:,end+1) = 1;
    
        y = x * [1 0 0 0]' + randn(30,1); % response to one var
    
        [b(:,j),dev,stat] = glmfit(x(:,1:3),y);
        t(:,j) = stat.t;
        s(j) = stat.s;
        se(:,j) = stat.se;
        p(:,j) = stat.p;
        
    end
    
    betas1(ind,:) = b(2,:); % true effect
    betas0(ind,:) = b(3,:); % null effect
    tval1(ind,:) = t(2,:);
    tval0(ind,:) = t(3,:);
    sex1(ind,:) = se(2,:) ./ s;
    sex0(ind,:) = se(3,:) ./ s;
    sigma(ind,:) = s;
    sig1(ind,:) = p(2,:) < .05;
    sig0(ind,:) = p(3,:) < .05;
    
    clear b t se p s
    
    ind = ind + 1;
    
end
    
    
%figure;hold on; plot(0:.01:.99,std(betas1'),'r','LineWidth',2)

% power and false positives
%figure('Color','w')
%plot(0:.01:.99,sum(sig0') ./ size(sig0,2),'k--')
%hold on;
%plot(0:.01:.99,sum(sig1' & betas1' > 0) ./ size(sig1,2),'k')
%xlabel('Correlation between predictors'),legend({'False positive rate' 'Power'})

% t-values vs. efficiency
figure;hold on
plot(0:.01:.99,mean(sigma'),'g')
plot(0:.01:.99,mean(tval1'),'r')
plot(0:.01:.99,1 ./ mean(sex1'),'b')
plot(0:.01:.99,mean(betas1'),'m--')
xlabel('Correlation between predictors'),legend({'Sigma' 'E(t-value)' 'Efficiency^(1/2)' 'E(b)'})

% power plot
% sigma = 1, beta = 1; 
figure;hold on;set(gca,'FontSize',14)
mp = [.05 .01 .001 .0001 .00001 .000001];
mt = tinv(1-mp,stat.dfe); mt = [mt;mt];
plot(repmat([0 1],size(mt,2),1)',mt,'k')
for i = 1:length(mp),text(.9,mt(1,i)+.2,num2str(mp(i))),end

mint = prctile(tval1',20);
plot(0:.01:.99,mint,'r','LineWidth',2)
title('Individual subject power (80% level)')
xlabel('Correlation between predictors'),ylabel('T-value (20th percentile)')

% power by efficiency
% sigma = 1, beta = 1; 
meff = 1 ./ mean(sex1');
figure('Color','w');hold on;set(gca,'FontSize',14)
mp = [.05 .01 .001 .0001 .00001 .000001];
mt = tinv(1-mp,stat.dfe); mt = [mt;mt];
plot(repmat([0 max(meff)+1],size(mt,2),1)',mt,'k')
for i = 1:length(mp),text(.9,mt(1,i)+.2,num2str(mp(i))),end

mint = prctile(tval1',20);
plot(meff,mint,'ro','LineWidth',2)
title('Individual power (80% level) by efficiency')
xlabel('Sqrt(efficiency)'),ylabel('T-value (20th percentile)')

bp = glmfit(meff,mint);
str = sprintf('Power = %3.2f + %3.2f(Efficiency)^1^/^2',bp(1),bp(2));
legend(str)
refline

% group se is se(g) + se(ind)

    gse = 1;    % group standard deviation (between subjects)
    trueb = 1;
    
    n = 12;
    gsei = b1std ./ (sqrt(n));
    seg = gsei + gse;
    tg = trueb ./ seg;  % expected t-value for group
    
    % group std dev
    % this is not quite right
    for i = 1:size(betas1,1), B(:,i) = bootstrp(1000,'std',betas1(i,:));,end
    colors = {'b' 'g' 'm' 'c'}; ind = 1;
    for n = [10 20 30 40];
        
        for i = 1:size(betas1,1), t20(i) = trueb ./ ((prctile(B(:,i),80) + gse) ./ sqrt(n));,end
        hold on; plot(meff,t20,[colors{ind} 's']);
        ind = ind + 1;
        
    end

% individual diffs - st. error of individual beta
% group - standard error of group activation by n

% variables:  n, b, contrast weights, sigma
% group vs. individual

[h,x]=hist(tval0(:)); figure;hold on; for i = 1:size(tval0,1); h=hist(tval0(i,:),x);plot(x,h);,end
plot(x,h,'r')
title('Distribution of H0 t-values at different correlations with true predictor')

[h,x]=hist(tval1(:)); figure;hold on; for i = 1:size(tval1,1); h=hist(tval1(i,:),x);plot(x,h);,end
plot(x,h,'r')
title('Distribution of H1 t-values at different correlations with true predictor')

[h,x]=hist(betas0(:)); figure;hold on; for i = 1:size(betas0,1); h=hist(betas0(i,:),x);plot(x,h);,end
plot(x,h,'r')
title('Distribution of H0 betas at different correlations with true predictor')

[h,x]=hist(betas1(:)); figure;hold on; for i = 1:size(betas1,1); h=hist(betas1(i,:),x);plot(x,h);,end
plot(x,h,'r')
title('Distribution of H1 betas at different correlations with true predictor')


    
