function [tthr,rthr,tnegthr,rnegthr,tmax,rmax,tmin,rmin,textmax,rextmax] = npm_ttest(x,niter,beh,varargin)
% [tthr,rthr,tnegthr,rnegthr,tmax,rmax,tmin,rmin,textmax,rextmax] = ...
% npm_ttest(x,niter,beh,[p for cluster extent test],[XYZ])
% tor wager
% oct 2003
%
% does cluster extent threshold, given p value for initial cutoff
% and XYZ voxel coordinates (3 x n)
%
% textmax: maximum cluster size for intercept, used for cluster extent
% rextmax: maximum cluster size for regressor, used for cluster extent
%
% to get corrected p-values, e.g.,
% rcorp = 1 - (sum(rmax <= rm) ./ length(rmax));

tthr = [];,rthr=[];tnegthr=[];rnegthr=[];tmax=[];rmax=[];

if length(varargin) > 0, 
    % set up cluster extent threshold
    pthr = varargin{1};, 
    XYZ = varargin{2};, 
    tthr = abs(tinv(pthr,size(x,1) - min(size(beh)) + 1));  % t threshold
else, pthr = [];, 
end

signv = sign(randn(size(x,1),niter));

figure('Color','w');

    
for i = 1:niter
    
    % t-test
    
    xtmp = x .* repmat(signv(:,i),1,size(x,2));
    
    ttmp = mean(xtmp) ./ (std(xtmp)./sqrt(size(xtmp,1)));
    
    tmax(i) = max(ttmp);
    tmin(i) = min(ttmp);
    
    % cluster extent
        
    if ~isempty(pthr)
        tmp = XYZ(:,ttmp > tthr);
        if isempty(tmp), textmax(i) = 0;,
        else, 
            tmp = spm_clusters(tmp);
            textmax(i)=max(hist(tmp,max(tmp)));                % get the max cluster size
        end
    end
        
    
    if mod(i,100)==0, [pltmp,xx] = hist(tmax); plot(xx,pltmp,'r','LineWidth',2),title(['It ' num2str(i)]),drawnow,end
    
    % simple regression
    
    if ~isempty(beh)
        
        btmp = getRandom(beh);
        
        % vectorize simple regression
        
        XX = [ones(size(btmp,1),1) btmp];
        XXpinv = pinv(XX);
        XTXI = diag(inv(XX' * XX));
        y = x;
        %se = sqrt(diag(inv(XX' * XX)) .* (r'*r) / 2)
        %b = pinv([ones(size(X,1),1) X]) * y;
        b = XXpinv * y;
        r = y - XX * b;
        dfe = size(XX,1) - size(XX,2);
        se = sqrt(repmat(XTXI,1,size(y,2)) .* repmat(var(r) .* (size(XX,1)-1) / dfe,size(XX,2),1));
        rtmp = b ./ se;
        
        rmax(i) = max(rtmp(2,:));
        rmin(i) = min(rtmp(2,:));
        
        % cluster extent
        
        if ~isempty(pthr)
            tmp = XYZ(:,rtmp(2,:) > tthr);
            if isempty(tmp), rextmax(i) = 0;,
            else, 
                tmp = spm_clusters(tmp);
                rextmax(i)=max(hist(tmp,max(tmp)));                % get the max cluster size
            end
        end
    
        if mod(i,100)==0, hold on, pltmp = hist(rmax); plot(xx,pltmp,'b','LineWidth',2),drawnow,end
    end
        
end

if any(~isnan(tmax))
    tthr = prctile(tmax,95);
    tnegthr = prctile(tmin,5);
    if ~isempty(beh),
        rthr = prctile(rmax,95);,
        rnegthr = prctile(rmin,5);
    end
else
    return
end

close

return
