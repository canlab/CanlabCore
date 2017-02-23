function regression_table(y,X,names,varargin)
%regression_table(y,X,names,[do figures],[do robust])
%
% performs regression using the regress command, and prints an output table
% X should not have an intercept
%
% Tor Wager

dofigs = 1; dorobust = 0;
if length(varargin) > 0, dofigs = varargin{1}; end
if length(varargin) > 1, dorobust = varargin{2}; end

[n, k] = size(X);

if dorobust
    [b,stat]=robustfit(X,y); 
    BINT = NaN * zeros(n + 1,2);
    STATS = [NaN NaN NaN];
    partialr = NaN * zeros(1, k + 1);
    
else
    [b,dev,stat]=glmfit(X,y); 
    [B,BINT,R,RINT,STATS] = regress(y,[ones(n ,1) X]);
    
    % partial correlations
    partialr(1) = NaN;  % intercept
    for i = 1:k, [x,y,partialr(i+1),partialp] = partialcor(X,y,i); end
end


names = [{'Intercept'}, names];

fprintf(1,'Parameter\tb-hat\tt\tp\tConf. interval\t\tPartial Corr\tRob. Partial Corr.\n')

if dofigs 
    create_figure('Regression Table output', k, 2); 

end
    
for i = 1:length(b)
    
    if i == 1, %intercept
            fprintf(1,'%s\t%3.2f\t%3.2f\t%3.4f\t%3.2f\t%3.2f\t\n', ...
        names{i},b(i),stat.t(i),stat.p(i),BINT(i,1),BINT(i,2));
    
    else
        if dofigs
            subplot(k, 2, 2 * (i - 2) + 1); [r,str,sig] = prplot(y,X,i-1,0); xlabel([names{i}])
            subplot(k, 2, 2 * (i - 2) + 2); [r2,str,sig] = prplot(y,X,i-1,1); xlabel([names{i} ': Robust IRLS'])
        else
            r = partialr(i); 
            r2 = NaN;
        end
        
        fprintf(1,'%s\t%3.2f\t%3.2f\t%3.4f\t%3.2f\t%3.2f\t%3.2f\t%3.2f\t\n', ...
        names{i},b(i),stat.t(i),stat.p(i),BINT(i,1),BINT(i,2),r,r2);
    end
end


fprintf(1,'\nOverall R2 = %3.2f, F = %3.2f, p = %3.4f\t\n',STATS(1),STATS(2),STATS(3));


return