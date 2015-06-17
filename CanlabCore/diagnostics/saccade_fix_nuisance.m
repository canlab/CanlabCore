load Nuisance_covariates.mat; load SPMcfg 
 
X = Xn;
px = X * pinv(X);

fprintf(1,'BEFORE\n')
for i = 1:length(xX.iC)

    r = xX.X(:,xX.iC(i)) - px * xX.X(:,xX.iC(i));       % residuals

    pve = 1 - ((r' * r) ./ (xX.X(:,xX.iC(i))' * xX.X(:,xX.iC(i)))); % percentage of variance explained

    fprintf(1,'\n%s\t%3.2f%%',xX.Xnames{xX.iC(i)},100*pve)

end
fprintf(1,'\n')

mX = xX.X(:,xX.iC); mX(:,end+1) = 1; Xn = orthogonalize(mX,X,1);
X = Xn;
px = X * pinv(X);

fprintf(1,'AFTER\n')
for i = 1:length(xX.iC)

    r = xX.X(:,xX.iC(i)) - px * xX.X(:,xX.iC(i));       % residuals

    pve = 1 - ((r' * r) ./ (xX.X(:,xX.iC(i))' * xX.X(:,xX.iC(i)))); % percentage of variance explained

    fprintf(1,'\n%s\t%3.2f%%',xX.Xnames{xX.iC(i)},100*pve)

end
fprintf(1,'\n')


save Nuisance_covariates.mat Xn
reset_SPMcfg
add_nuisance_to_SPMcfg(Xn)


%check one more time

load SPMcfg
mX = xX.X(:,xX.iC); mX(:,end+1) = 1; Xn = xX.X(:,xX.iB);
figure;subplot 121; imagesc(mX); subplot 122; imagesc(Xn)
px = Xn * pinv(Xn);
fprintf(1,'FINAL\n')
for i = 1:length(xX.iC)

    r = xX.X(:,xX.iC(i)) - px * xX.X(:,xX.iC(i));       % residuals

    pve = 1 - ((r' * r) ./ (xX.X(:,xX.iC(i))' * xX.X(:,xX.iC(i)))); % percentage of variance explained

    fprintf(1,'\n%s\t%3.2f%%',xX.Xnames{xX.iC(i)},100*pve)

end
fprintf(1,'\n')
