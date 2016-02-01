function [b,stats,yadj] = repeated_ancova(X,Y,wicons,btwnnames,winames,ynames,varargin)
% Repeated measures ANCOVA with table and plot
% uses Robust IRLS
%
% :Usage:
% ::
%
%     [b,stats,yadj] = repeated_ancova(X,Y,wicons,btwnnames,winames,ynames,varargin)
%
% :Examples:
% ::
%
%    Y = rand(15,2);
%    X = Y + rand(15,2);
%    cons = [-1 1 -1 1; 1 1 -1 -1];  % placebo vs control, hot vs. warm
%
%    X = R.X(:,1:2);
%    Y = cl(2).CONTRAST.data;
%    cons = [-1 1 0 0; 0 0 -1 1];  % placebo vs control for heat then warm
%    repeated_ancova(X,Y,cons,{'Reported Placebo (C - P)' 'Order'},{'P-C Heat' 'P-C Warm'},{'CH' 'PH' 'CW' 'PW'});
%
% DOES NOT WORK WITH FIXED BTWN-SUBJECTS COVARIATES
%
% X must be a random variable that is observed multiple times for each
% subject, as does Y
%
% ..
%    tor wager, Aug. 06
% ..

dotable = 1;
doplot = 1;

for i = 1:length(varargin)
    if isstr(varargin{i})
        switch varargin{i}
            % reserved keywords
            case 'notable', dotable = 0;
            case 'noplot', doplot = 0;
            
            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end


[n,r] = size(Y);

y = Y(:);

obs = size(y,1); % total num observations

% get predictors for between-subjects effects
% x
% ---------------------------------------------------
if all(size(X) - [n r] == 0)
    % X is obs. x conditions for a single predictor
    nbtwn = 1;
    x = X(:);   %repmat(X,r,1);
    btwnnames = repmat(btwnnames,1,r);
    
elseif size(X,1) == n
    disp('Warning!!!!!! Results are not valid for fixed within-subjects contrasts.')
    % X is a between-Ss effect replicated for each measure
    nbtwn = size(X,2);
    for i = 1:nbtwn
        x(:,i) = repmat(X(:,i),r,1);
    end
end
    
x = scale(x,1);  % center x

% get predictors for repeated measures effects
% c
% ---------------------------------------------------
% effect of repeated measure on y (repeated condition difference in ANCOVA)
% at x = 0
c = expand_contrasts(wicons,n,r);
nwithin = size(c,2);



% get preds for interaction between continuous x and c
% xc
% ---------------------------------------------------
% difference in x-y slopes for different repeated measures in c
if nbtwn
    [xc,xcnames,ninteract] = make_interactions(nbtwn,nwithin,x,c,n,r,btwnnames,winames);
else
    xc = []; xcnames = {}; ninteract = 0;
end



% Subject effects, obs x n-1

%nused = 1 + nbtwn ;  % warning!  not sure if this is right; only needed if btwn effects are replicated for all conditions
% THE ABOVE IS BAD; CANNOT ESTIMATE MODEL WHEN X IS FIXED.

S = repmat([eye(n-1); zeros(1,n-1)], r, 1);
%S = S(:,nused+1:end);

% y is regressed on x, controlling for S, and testing c

% full model
M = [x c xc S]; M(:,end+1) = 1;
p = size(M,2);              % total num params estimated
xtxitx = pinv(M);

res = y - M * xtxitx * y;   % residuals
SSr = res' * res;           % sums of squared residuals
%dfe = obs - p;              % error df
%dfb = p - 1;                % full model df

% reduced model: Subject effects only
Mr = S; Mr(:,end+1) = 1;
res = y - Mr * pinv(Mr) * y;   % residuals
SSr2 = res' * res;              % sums of squared residuals (2)

pr = size(Mr,2);              % total num params estimated, including intercept = n
dfe = obs - pr;              % error df
MSerr = SSr2 ./ dfe; 

% model comparison
SSrdiff = SSr2 - SSr;         % sums of squares explained by model params x c xc
dfdiff = p - pr;              % parameters used to explain SSrdiff
MSmodel = SSrdiff ./ dfdiff;    % mean square for non-subject, non-intercept effects in model


F = MSmodel ./ MSerr;
p = 1 - fcdf(F,dfdiff,dfe);   % all this is very approx now!!! it's late.

[b,stats] = robustfit(M,y,'bisquare',[],'off');

wh_reg = 1;

if nargout > 2
    % get which regs are of interest for fitted data
    z = get_which_of_interest(wh_reg,nbtwn,nwithin,ninteract,S);

    yadj = get_yadj(wh_reg,x,y,M,z,b,n,r);
end

% table
if dotable
    regtable(b,stats,nbtwn,nwithin,ninteract,btwnnames,winames,xcnames);
end


% plot
if doplot
    %regplot(wh_reg,b,x,y,X,Y,M,S,nbtwn,nwithin,ninteract,n,r,c,btwnnames,winames,ynames);
    regplot2(X,Y,r,stats.w,btwnnames,ynames);
    ylabel('Y');
    
    z = get_which_of_interest(wh_reg,nbtwn,nwithin,ninteract,S);
    yadj = get_yadj(wh_reg,x,y,M,z,b,n,r);
    regplot2(X,yadj,r,stats.w,btwnnames,ynames);

end

return






% setup: get within-subjects contrasts

function c = expand_contrasts(cons,n,r)

ncons = size(cons,1);
c = zeros(n*r,ncons);

for i = 1:ncons
    mycon = zeros(n*r,1);
    for j = 1:r
        mycon((j-1)*n+1 : j*n) = cons(i,j);
    end
    c(:,i) = mycon;
end

return



function [xc,xcnames,ninteract] = make_interactions(nbtwn,nwithin,x,c,n,r,btwnnames,winames)
ninteract = nbtwn * nwithin;
xc = zeros(n*r,ninteract);
xccounter = 1;
for i = 1:nbtwn
    for j = 1:nwithin
        xc(:,xccounter) = x(:,i) .* c(:,j);
        xcnames{xccounter} = [btwnnames{i} '*' winames{j}];
        xccounter = xccounter + 1;
    end
end
return



function regtable(b,stats,nbtwn,nwithin,ninteract,btwnnames,winames,xcnames)
fprintf(1,'Parameter\tb-hat\tt\tp\t    \n')

fprintf(1,'Between-subjects effects    \n')
wh = 1:nbtwn;       % which columns for these effects
for i=1:nbtwn
    fprintf(1,'%s\t%3.2f\t%3.2f\t%3.4f\t    \n', ...
        btwnnames{i},b(wh(i)),stats.t(wh(i)),stats.p(wh(i)));
end

fprintf(1,'Within-subjects condition effects    \n')
wh = nbtwn+1:nbtwn+nwithin;
for i=1:nwithin
    fprintf(1,'%s\t%3.2f\t%3.2f\t%3.4f\t    \n', ...
        winames{i},b(wh(i)),stats.t(wh(i)),stats.p(wh(i)));
end

fprintf(1,'Between * Within interactions    \n')
wh = nbtwn+nwithin+1 : nbtwn+nwithin+ninteract;
for i=1:ninteract
    fprintf(1,'%s\t%3.2f\t%3.2f\t%3.4f\t    \n', ...
        xcnames{i},b(wh(i)),stats.t(wh(i)),stats.p(wh(i)));
end

return





function z = get_which_of_interest(wh_reg,nbtwn,nwithin,ninteract,S)
z = zeros(1,nbtwn);
z(wh_reg) = 1;              % between of interest
z = [z ones(1,nwithin)];    % all within are of interest
zxc = zeros(1,ninteract);   % interaction effects
xccounter = 1;
for i = 1:nbtwn
    for j = 1:nwithin
        if i == wh_reg, zxc(xccounter) = 1; end
        xccounter = xccounter + 1;
    end
end
z = [z zxc];
z = [z zeros(1,size(S,2))]; % subject effects
z = [z 1];                  % intercept
z = logical(z);
return


function yadj = get_yadj(wh_reg,x,y,M,z,b,n,r);
% x (of interest) and adjusted y
xdata = reshape(x(:,wh_reg),n,r);
yadj = y - M(:,~z) * b(~z);
yadj = reshape(yadj,n,r);



function regplot2(X,Y,r,w,btwnnames,ynames)

tor_fig; set(gca,'FontSize',24)
colors = {'ro' 'gs' 'cv' 'mh'};
for i = 1:r
    h = plot_correlation(X(:,i),Y(:,i),'colors',colors(i),'robust','noprint','weights',w);
    han(i) = h{1}(1);
end

han2 = findobj(gcf,'Type','text');
delete(han2)

% range over which to get fits for line
xrange = max(X(:,1)) - min(X(:,1));
minx = min(X(:,1)) - .1*xrange;
maxx = max(X(:,1)) + .1*xrange;

set(gca,'XLim',[minx maxx]);

% Y range
yrange = max(Y(:)) - min(Y(:));
miny = min(Y(:)) - .1*yrange;
maxy = max(Y(:)) + .1*yrange;
set(gca,'YLim',[miny maxy]);

xlabel(btwnnames{1});
ylabel('Adjusted data');
%ynames = {'HC' 'HA' 'LC' 'LA'};
legend(han,ynames)

scn_export_papersetup(400);


return



% OLD FUNCTION--but a goodie.  this was for fixed regressors btwn.
function regplot(wh_reg,b,x,y,X,Y,M,S,nbtwn,nwithin,ninteract,n,r,c,btwnnames,winames,ynames)

npar = nbtwn+nwithin+ninteract;

% only significant betas
% sig = stats.p < .10;
% sig = sig(1:npar);
% bsig = b(1:npar) .* sig;
% bsig = [bsig; b(end)];

% all betas for btwn of interest and w/i conditions
% bsig = zeros(1,npar)';
% bsig(wh_reg) = b(wh_reg);  % regressor of interest
% wh = nbtwn+1:nbtwn+nwithin; % condition, within Ss effects
% bsig(wh) = b(wh);
% wh = nbtwn+nwithin+1:npar; % interactions with reg of interest
% wh = wh(1:2);           %%%%%kludgy fix for reg 1 only
% bsig(wh) = b(wh);
% bsig = [bsig; b(end)];

% range over which to get fits for line
xrange = max(X(:,1)) - min(X(:,1));
minx = min(X(:,1)) - .1*xrange;
maxx = max(X(:,1)) + .1*xrange;

% get predicted values for this btwn effect (xf)
x2 = linspace(minx,maxx,n)';
x2 = repmat(x2,1,r);
x2 = x2(:);
xf = x;
xf(:,wh_reg) = x2;

% get interactions
[xcf,xcnamesf,ninteract] = make_interactions(nbtwn,nwithin,xf,c,n,r,btwnnames,winames);

% Model for pred values
Mf = [xf c xcf S]; 
Mf(:,end+1) = 1;

% get which regs are of interest for fitted data
z = get_which_of_interest(wh_reg,nbtwn,nwithin,ninteract,S);

%fit lines
f = Mf(:,z) * b(z);          % fitted data (for lines)
f = reshape(f,n,r);         % fits
xfit = reshape(x2,n,r);   % x values for fits


% x (of interest) and adjusted y
xdata = X;
%xdata = reshape(x(:,wh_reg),n,r);
yadj = y - M(:,~z) * b(~z);
yadj = reshape(yadj,n,r);

% make figure
tor_fig; set(gca,'FontSize',24)
colors = {'ro' 'bs' 'kv' 'gp'};
for i = 1:r
    plot(xdata(:,i),yadj(:,i),colors{i},'MarkerFaceColor',colors{i}(1));  % adjusted data
    han(i) = plot(xfit(:,i),f(:,i),[colors{i}(1) '-']);    % fit line
end

xlabel(btwnnames{wh_reg});
ylabel('Adjusted data');
%ynames = {'HC' 'HA' 'LC' 'LA'};
legend(han,ynames)

set(gca,'XLim',[minx maxx]);

% Y range
yrange = max(yadj(:)) - min(yadj(:));
miny = min(yadj(:)) - .1*yrange;
maxy = max(yadj(:)) + .1*yrange;
set(gca,'YLim',[miny maxy]);

scn_export_papersetup(400);

return
