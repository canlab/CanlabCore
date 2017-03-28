function [corr,t,p, fdrsig, fdrthresh] = correlation(meth,x,varargin)
% Multiple types of correlations, including Spearman's rho
% (nonparametric) and phi (dichotomous)
%
% :Usage:
% ::
%
%     [corr,t,p,fdrp, fdrthresh] = correlation(method,x,[y],['matrix'])
%
%     IN PROGRESS : Warning : Use at your own risk.
%     Some methods are not adequately tested yet.
%     Spearman's rho does not correct for ties
%
% :Inputs:
%
%   **Methods:**
%        String indicating the method for computing the correlation
%        coefficient
%        - Pearson's r.          Enter: {'r','pearson',[]}
%        - IRLS                  Enter: {'irls','robust'}
%        - Phi                   Enter: {'phi'}
%        - Spearman's rho        Enter: {'rho','spearman'}
%        - Kendall's Tau (a)     Enter: {'taua','kendalla'}
%        - Tau (b)               Enter: {'tau','kendall','taub','kendallb'}
%        - Gamma                 Enter: {'gamma','kruskal'}
%
%   **x:**
%        Matrix of observations (n instances by p varianbles)
%
% :Optional Inputs:
%
%   **varargin:**
%        To be documented
%
% :Outputs:
%
%   **OUT:**
%        Output stats structure
%
% :Examples:
% ::
%
%    % Corelation between two variables, Pearson's
%    x = rand(10,1); y = rand(10,1);
%    [corr,t,p] = correlation('r',x,y);
%
%    % Correlation matrix of 10 variables, phi correlation:
%    studybyroi = magic(10);
%    [corr,t,p] = correlation('phi',studybyroi);
%
% :See Also: correlation_fast_series.m
%
% ..
%    tor wager, november 2006, Jan 2007
% ..

corr = []; t = []; p = [];


% get y data, if entered
% --------------------------------
for i = 1:length(varargin)
    if ~ischar(varargin{i}) && ~any(size(varargin{i}) - size(x))
        y = varargin{i};
    end
end

% get flag for whether to compute correls among all possible pairs
% --------------------------------
matrixFormFlag = 0;
if any(strcmp(varargin,'matrix')) || ~exist('y','var')
    matrixFormFlag = 1;
end

[n,nvars] = size(x);

% Get indices
% --------------------------------
if matrixFormFlag
    % input is matrix
    [rows,cols,npairs] = corrcoef_indices(nvars);

else
    % input is columns of consecutive pairs to correlate
    npairs = size(x,2); % number of pairs of correlations to compute
end

doverbose = 0;
if npairs > 1000, doverbose = 1; fprintf('%03d%%', 0); end

% Compute correlations
% --------------------------------
for i = 1:npairs

    if doverbose, fprintf('\b\b\b\b%03d%%', round(i*100 ./ npairs)), end
    
    % get data
    if matrixFormFlag
        x1 = x(:,rows(i)); x2 = x(:,cols(i));
    else
        x1 = x(:,i); x2 = y(:,i);
    end

    if nargout > 1
        [c,tt,pp] = compute_single_correl(meth,x1,x2,n);
    else
        c = compute_single_correl(meth,x1,x2,n);
    end

    corr(i,1) = c;

    if nargout > 1
        if ~isempty(tt), t(i,1) = tt; end
        if ~isempty(pp), p(i,1) = pp; end
    end
end

% Reconstruct into matrix, if needed
% --------------------------------
if matrixFormFlag
    corr = reconstruct(corr,nvars,npairs,rows,cols);

    corr = corr + eye(nvars);

    if ~isempty(t), t = reconstruct(t,nvars,npairs,rows,cols); end

    if ~isempty(p), p = reconstruct(p,nvars,npairs,rows,cols); end

    % FDR correction, if requested
    % only works for matrices, otherwise error
    if nargout > 3
        [fdrthresh,fdrsig] = fdr_correct_pvals(p,corr);
    end

elseif nargout > 3
    error('FDR output is only allowed for matrix form output.  Try requesting fewer outputs')
end




return





% sub-functions



function [est,t,p] = compute_single_correl(meth,x,y,n)

est = [];
t = [];
p = [];

[wasnan, x, y] = nanremove(x, y);

switch lower(meth)

    case {'irls','robust'}
        [b,stats]=robustfit(x,y,'bisquare'); %,[],'off');
        est = weighted_corrcoef([x y], stats.w);
        est = est(1,2);

        t = stats.t(2);
        p = stats.p(2);

    case {'r','pearson',[]}

        if islogical(x)
            warning('Logical vectors not appropriate for Pearson''s r.');
            x = double(x);
            y = double(y);
        end

        [est,p] = corrcoef([x y]);
        est = est(1,2);
        p = p(1,2);

    case 'phi'
        % Phi is for dichotomous (2-level) variables
        % formula from Robert Yaffee, NYU, online.


        tab = make_crosstabs(x,y);
        num = det(tab);                 % same as: num = tab(1,1)*tab(2,2) - tab(1,2)*tab(2,1);
        den = ( prod(sum(tab,2)) .* prod(sum(tab,1)) ).^.5;
        est = num ./ den;

        % F = varexp/dfexp / varunexp/dferror; t = sqrt(F)
        if nargout > 1 % for speed
            t = est .* sqrt((n - 2) ./ (1 - est.^2));
            p = 2 .* (1 - tcdf(abs(t),n-2));        % two-tailed p-value
        end
        %wts = [];
        %nonpar = 0;
        % this works,but is unsigned
        %[chi2,df,p,sig] = chi2test([x y],'obs'); %,wts,nonpar);
        %est = (chi2 ./ n) .^ .5;


    case {'rho','spearman'}

        % This method uses midrank method to handle ties; needs checking
        % vs. SPSS
        D = rankdata(x) - rankdata(y);
        est = 1 - (6*(D'*D)) ./ (n * (n^2-1));

        % same as : corrcoef(rankdata(x),rankdata(y));
        % but faster
        % tested against SPSS 11.04
        % problem with large n?

        % % n = 1000;
        % % tic, for i = 1:100, x = rand(n,1); y = x + rand(n,1);
        % % corrcoef(rankdata(x),rankdata(y));
        % % end, toc
        % % tic, for i = 1:100, x = rand(n,1); y = x + rand(n,1);
        % %     n = size(x,1);
        % % D = rankdata(x) - rankdata(y);
        % % 1 - (6*(D'*D)) ./ (n * ((n^2)-1));
        % % end, toc


    case {'taua','kendalla'}

        % checked 1/13/07 against http://www.wessa.net/rwasp_kendall.wasp
        % should probably be checked against SPSS

        [r,i] = sort(x);
        y = y(i);               % get ordered y; don't really need to rank

        % count number of ranks above each successive value of y
        for i = 1:(n - 1)
            P(i) = sum( y(i+1:end) > y(i) );
        end
        est = 4 * sum(P) ./ (n * (n - 1)) - 1;

        if nargout > 1 % for speed
            % actually a z-score
            t = 3 * est * (n * (n - 1)) .^.5 ./ (2 * (2 * n + 5)) .^.5;
            p = 1 - normcdf(abs(t));
            
            if p == 0, p = 10*eps; end
        end

    case {'tau','kendall','taub','kendallb'}

        % checked 1/13/07 against http://www.wessa.net/rwasp_kendall.wasp
        % should probably be checked against SPSS

        [r,i] = sort(x);
        y = y(i);               % get ordered y

        % count number of ranks above each successive value of y
        for i = 1:(n - 1)
            wh = r > r(i);  % find pts with x value greater than r(i); handle x ties
            P(i) = sum( y(wh) > y(i) ); % concordances
            N(i) = sum( y(wh) < y(i) ); % discordances
        end

        %% *** maybe ties have to include all values that are tied"??
        % http://www.statsdirect.com/help/nonparametric_methods/kend.htm
        xties = get_ties(x);
        yties = get_ties(y);

        num = sum(P) - sum(N);

        n1 = n * (n - 1) ./ 2;

        % From http://www.unesco.org/webworld/idams/advguide/Chapt4_2.htm
        % and http://www.statsdirect.com/help/nonparametric_methods/kend.htm
        den = sqrt(  (n1 - sum(xties)) * (n1 - sum(yties))  );

        % Tau b : from http://www.unesco.org/webworld/idams/advguide/Chapt4_2.htm
        %den = sqrt((sum(P) + sum(N) + (n - sum(xties))) * (sum(P) + sum(N) + (n - sum(yties))));

        est = num ./ den;

        if nargout > 1 % for speed
            % actually a z-score
            t = 3 * est * (n * (n - 1)) .^.5 ./ (2 * (2 * n + 5)) .^.5;
            p = 1 - normcdf(abs(t));
        end


    case {'gamma','kruskal'}

        [r,i] = sort(x);
        y = y(i);               % get ranks of y

        % count number of ranks above each successive value of y
        for i = 1:(n - 1)
            wh = r > r(i);  % find pts with x value greater than r(i); handle x ties
            P(i) = sum( y(wh) > y(i) ); % concordances
            N(i) = sum( y(wh) < y(i) ); % discordances
        end

        num = sum(P) - sum(N);
        den = sum(P) + sum(N);
        est = num ./ den;

    otherwise
        error('Unknown correlation option');
end



return






function tab = make_crosstabs(x,y)
% make table (crosstabs)
u1 = unique(x);
u2 = unique(y);
len1 = length(u1);
len2 = length(u2);
tab = zeros(len1,len2);
w = 1;                  % weights

for i = 1:length(u1)                   % rows are 0 then 1 on the first var, "Nos" then "Yesses"
    for j = 1:length(u2)               % for each column
        tab(i,j) = sum( (x == u1(i) & y == u2(j)) .* w );
    end
end

return



function ties = get_ties(x)
r = rankdata(x, 'nomidrank');
for i = 1:length(unique(r))
    ties(i) = sum(r == i);
end
ties = sum(ties .* (ties - 1) ./ 2);
return



function [rows,cols,npairs] = corrcoef_indices(nvars)

    if nvars > 1000
        fprintf('Setting up indices of matrix.');
        % added to deal with large datasets
        tmp = triu(true(nvars));

        for i = 1:nvars
            tmp(i, i) = 0;
        end
        
%         % create logical identity matrix of large size
%         t2 = eye(1000);
%         t2 = logical(t2);
%         eyemtx = t2;
%         while size(eyemtx, 2) < nvars
%             eyemtx = blkdiag(eyemtx, t2);
%         end
%         eyemtx = eyemtx(1:nvars, 1:nvars);
% 
%         tmp = tmp - eyemtx;

        fprintf(1,'Done.\n');

    else
        % upper triangle only
        tmp = triu(ones(nvars));
        tmp = tmp - eye(nvars);
    end

    [rows,cols] = find(tmp);
    npairs = length(rows);
    
    return



function valmat = reconstruct(vals,nvars,npairs,rows,cols)

valmat = zeros(nvars);
for i = 1:npairs
    valmat(rows(i),cols(i)) = vals(i);
end
valmat = valmat + valmat';

return


function [pthr,sig] = fdr_correct_pvals(p,r)

    psq = p; psq(find(eye(size(p,1)))) = 0;
    psq = squareform(psq);
    pthr = FDR(p,.05);
    if isempty(pthr), pthr = 0; end

    sig = sign(r) .* (p < pthr);

    return
    
