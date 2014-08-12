function [r,str,sig,ry,rx,h,rr] = prplot(yy,X,k,varargin)
% [r,str,sig,ry,rx,h,rr] = prplot(y,X,col,[dorobust],[colors])
% Partial residual plot of one column of X against y.
% Uses IRLS estimation to downweight outliers 
% if you enter a 4th argument
%
% Partial residual plot of y ~ X for column k
% tor wager
%
% if y contains multiple columns, different colors
% and symbols will be used, with a separate regression
% for each.
%
% colors: e.g., {'ro' 'bs' 'gd' 'y^' 'cv' 'mx'}

mycols = {'ko'};
robopt = 'IRLS';    %'IRLS' or 'MCD'

if size(yy,2) > 1, 
    
        legstr = cell(1);
        mycols = {'ro' 'bs' 'gd' 'y^' 'cv' 'mx'};
        
end

dorobust = 0;
if length(varargin) > 0, dorobust = varargin{1}; end
if length(varargin) > 1, mycols = varargin{2}; end

% ----------------------------------------
% do this separately for every column of y
% ----------------------------------------

for i = 1:size(yy,2)
    
    y = yy(:,i);
    
    Xm = X; 

    if dorobust && strcmp(robopt,'MCD')
        [res]=fastmcd_noplot([Xm y]);
        % remove n most extreme outliers and recompute correlation
        wh = res.flag==0; nout = sum(res.flag==0);
    
        %figure;plot_correlation_samefig(X(:,2),y,[],'ro');
    
        y(wh,:) = []; Xm(wh,:) = [];
    
        %plot_correlation_samefig(X(:,2),y,[],'kx');
        %nout
    end

    sel = Xm(:,k);
    xmean = mean(sel);  % save mean to add back in later
    
    X2 = Xm; X2(:,k) = [];
    
   
    % problem if cols are duplicated (as in script i'm running now)
    % just a fix - changes nothing if no duplicate columns
    % temporary fix.
    %X2 = unique(X2','rows')';
    %for j = 1:size(X2,2), tmp=corrcoef(X2(:,j),sel); tmpa(j)=tmp(1,2);,end
    %X2(:,tmpa > .99) = [];

    % ----------------------------------------
    % regress out columns X2 from selected X(:,k)
    % leave rx, the final adjusted x predictor
    % ----------------------------------------
        
    % Xm is full design matrix
    % X2 is design matrix without column k
    % sel is column k of design matrix
    
    if dorobust && strcmp(robopt,'IRLS')
        
        % fit model and save overall residuals
        if nargout > 6
            [b,stat] = robustfit(Xm,y);
            rr = y - [ones(size(Xm,1),1) Xm] * b;
            w = stat.w;    % weights
        end
        
        % regress out (control for) columns of no interest (x2) from sel
        [b, stat] = robustfit(X2,sel);
        rx = sel - [ones(size(X2,1),1) X2] * b;              % residual X component of interest
        
    else
        % non-robust method
    
        Xm(:,end+1) = 1;
        
        % save overall residuals
        if nargout > 6
            try rr(:,i) = y - Xm * pinv(Xm) * y; catch end
        end
        
        X2(:,end+1) = 1;    % add intercept
        b = pinv(X2) * sel; 
        rx = sel - X2 * b;
    end
    
    % y done in next lines

   
    
    % ----------------------------------------
    % regress out columns X2 from y
    % leave ry, the final adjusted y data
    % ----------------------------------------
    
    if dorobust && strcmp(robopt,'IRLS')
            b = robustfit(X2,y);    % no intercept above; added by robustfit
            b(end+1) = b(1);        % move intercept to end, for compatibility with OLS
            b = b(2:end);
            X2(:,end+1) = 1;        % add intercept to end, for resid getting later
    else
            b = pinv(X2) * y;       % intercept added above
    end
    
    if size(yy,2) > 1
        
        ry{i} = y - X2 * b;
        ry{i} = ry{i} + b(end); % add intercept back in
        rx = rx + xmean;
        
        [r(i),str{i},sig(i),hh] = plot_correlation_samefig(rx,ry{i},[],mycols{i},0,dorobust);
        h(i) = hh(1);
        legstr{i} = ['DV ' num2str(i)];
        if exist('nout') == 1, legstr{i} = [legstr{i} ' nout=' num2str(nout)]; end
        
    else

        ry = y - X2 * b;
        ry = ry + b(end); % add intercept back in
        rx = rx + xmean;
        
        if length(unique(rx)) > 2
            [r,str,sig,hh] = plot_correlation_samefig(rx,ry,[],mycols{1},0,dorobust);
            
        else
            % t-test: bar plot, adjusting for covariates
            rx = scale(rx, 1); %mediansplit(rx);
            [H,p,ci,stats] = ttest2_printout(ry(rx>0),ry(rx<0),1);
            r = stats.tstat; str = 't-test'; sig = H; hh = [];         
        end
        
            
        h = hh;
    
    end

end % loop through yy

if size(yy,2) > 1
        legend(h,legstr);
end
    
return
