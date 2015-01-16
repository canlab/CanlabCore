function stats = rmanova2(data,alpha,doplot,ttst)
%
%Repeated-measures two-way ANOVA
%
%USAGE: stats = rmanova2(data,[alpha],[doplot],[ttst]);
%
%INPUTS: 
%    data: can be one of two formats:
%          1. Cell array - each row represents a level of factor 1, and
%          each column represents a level of factor 2. Each cell contains a
%          vector of values of the dependent variable for each subject.
%          2. Matrix - each row represents a trial, with the following
%          columns:
%               column1 - dependent variable
%               column2 - grouping variable for subject
%               column3 - grouping variable for factor 1
%               column4 - grouping variable for factor 2
%    alpha (optional): p-value threshold (default: 0.05)
%    doplot (optional): if 1, will produce a line plot.  
%                       Works only for cell input data (default: 1)
%    ttst (optional): if 1, will perform pairwise t-tests (default: 0)
%
% Aaron Schurger (2005.02.04)
%   Derived from Keppel & Wickens (2004) "Design and Analysis" ch. 18
% Modified by Sam Gershman (2006.11.16)
%

if nargin < 2; alpha = 0.05; doplot = 0; ttst = 0; end
if nargin < 3; doplot = 1;  ttst = 0; end
if nargin < 4; ttst = 0; end

if iscell(data);
    Y = []; S = []; F1 = []; F2 = [];
    [m n] = size(data);
    for i = 1:m;
        for j = 1:n;
            Y = cat(1,Y,data{i,j});
            nsubs = length(data{i,j});
            S = cat(2,S,1:nsubs);
            F1 = cat(1,F1,zeros(nsubs,1)+i);
            F2 = cat(1,F2,zeros(nsubs,1)+j);
        end
    end
    S = S';
else
    Y = data(:,1);
    S = data(:,2);
    F1 = data(:,3); F2 = data(:,4);
    clear data
    su = unique(S); f1u = unique(F1); f2u = unique(F2);
    for a = 1:length(f1u);
        for b = 1:length(f2u);
            ii = intersect(find(F1==f1u(a)),find(F2==f2u(b)));
            d = [Y(ii) S(ii)];
            d = sortrows(d,2);
            data{a,b} = d(:,1);
        end
    end
end

F1_lvls = unique(F1);
F2_lvls = unique(F2);
Subjs = unique(S);

a = length(F1_lvls); % # of levels in factor 1
b = length(F2_lvls); % # of levels in factor 2
n = length(Subjs); % # of subjects

INDS = cell(a,b,n); % this will hold arrays of indices
CELLS = cell(a,b,n); % this will hold the data for each subject X condition
MEANS = zeros(a,b,n); % this will hold the means for each subj X condition

% Calculate means for each subject X condition.
% Keep data in CELLS, because in future we may want to allow options for
% how to compute the means (e.g. leaving out outliers > 3stdev, etc...).
for i=1:a % F1
    for j=1:b % F2
        for k=1:n % Subjs
            INDS{i,j,k} = find(F1==F1_lvls(i) & F2==F2_lvls(j) & S==Subjs(k));
            CELLS{i,j,k} = Y(INDS{i,j,k});
            MEANS(i,j,k) = mean(CELLS{i,j,k});
        end
    end
end

% make tables (see table 18.1, p. 402)
AB = reshape(sum(MEANS,3),a,b); % across subjects
AS = reshape(sum(MEANS,2),a,n); % across factor 2
BS = reshape(sum(MEANS,1),b,n); % across factor 1

A = sum(AB,2); % sum across columns, so result is ax1 column vector
B = sum(AB,1); % sum across rows, so result is 1xb row vector
S = sum(AS,1); % sum across columns, so result is 1xs row vector
T = sum(sum(A)); % could sum either A or B or S, choice is arbitrary

% degrees of freedom
dfA = a-1;
dfB = b-1;
dfAB = (a-1)*(b-1);
dfS = n-1;
dfAS = (a-1)*(n-1);
dfBS = (b-1)*(n-1);
dfABS = (a-1)*(b-1)*(n-1);

% bracket terms (expected value)
expA = sum(A.^2)./(b*n);
expB = sum(B.^2)./(a*n);
expAB = sum(sum(AB.^2))./n;
expS = sum(S.^2)./(a*b);
expAS = sum(sum(AS.^2))./b;
expBS = sum(sum(BS.^2))./a;
expY = sum(Y.^2);
expT = T^2 / (a*b*n);

% sums of squares
ssA = expA - expT;
ssB = expB - expT;
ssAB = expAB - expA - expB + expT;
ssS = expS - expT;
ssAS = expAS - expA - expS + expT;
ssBS = expBS - expB - expS + expT;
ssABS = expY - expAB - expAS - expBS + expA + expB + expS - expT;
ssTot = expY - expT;

% mean squares
msA = ssA / dfA;
msB = ssB / dfB;
msAB = ssAB / dfAB;
msS = ssS / dfS;
msAS = ssAS / dfAS;
msBS = ssBS / dfBS;
msABS = ssABS / dfABS;

% f statistic
fA = msA / msAS;
fB = msB / msBS;
fAB = msAB / msABS;

% p values
pA = 1-fcdf(fA,dfA,dfAS);
pB = 1-fcdf(fB,dfB,dfBS);
pAB = 1-fcdf(fAB,dfAB,dfABS);

% return values
 stats.SS1 = ssA;
 stats.df1 = dfA;
 stats.MS1 = msA;
 stats.F1 = fA;
 stats.P1 = pA;
 stats.SS2 = ssB;
 stats.df2 = dfB;
 stats.MS2 = msB;
 stats.F2 = fB;
 stats.P2 = pB;
 stats.SS12 = ssAB;
 stats.df12 = dfAB;
 stats.MS12 = msAB;
 stats.F12 = fAB;
 stats.P12 = pAB;
 
 stats.alpha = alpha;
 
 %decision rule
if stats.P1 < alpha; stats.significant1 = 'Yes'; else; stats.significant1 = 'No'; end;
if stats.P2 < alpha; stats.significant2 = 'Yes'; else; stats.significant2 = 'No'; end;
if stats.P12 < alpha; stats.significant12 = 'Yes'; else; stats.significant12 = 'No'; end;

%make line plots
if doplot
    anova_line_plot(data);
    
%     [m n] = size(data); plotcounter = 0;
%     figure;
%     for i = 1:m;
%         for j = 1:n;
%             plotcounter = plotcounter + 1;
%             subplot(m,n,plotcounter);
%             boxplot(data{i,j});
%             xlabel([num2str(i),',',num2str(j)]);
%             set(gca,'XTickLabel',[]);
%         end
%     end
end

%pairwise t-tests
if ttst;
    disp('Performing post-hoc t-tests...');
    tcounter = 0; a = 0; b = 0;
    [m n] = size(data);
    for t1 = 1:m;
        for t2 = 1:n;
            a = a + 1; b = 0;
            for t3 = 1:m;
                for t4 = 1:n;
                    b = b + 1;
                    if b > a;
                        tcounter = tcounter + 1;
                        [h,p,ci,stat] = ttest(data{t1,t2},data{t3,t4});
                        stats.ttests(tcounter).P = p;
                        stats.ttests(tcounter).comparison = [num2str(t1),',',num2str(t2),' > ',num2str(t3),',',num2str(t4)];
                        stats.ttests(tcounter).means = [mean(data{t1,t2}) mean(data{t3,t4})];
                        stats.ttests(tcounter).ci = ci';
                    end
                end
            end
        end
    end
end

end
  
  
function anova_line_plot(anova_dat)
    
    if ~iscell(anova_dat)
        error('For line plots, enter cells rather than matrix input!');
    end
  
for i = 1:size(anova_dat, 1)
for j = 1:size(anova_dat, 2)
linedat(i, j) = nanmean(anova_dat{i, j});
lineste(i, j) = ste(anova_dat{i, j});
end
end

% Transpose so F2 is on X-axis, F1 is lines
linedat = linedat';
lineste = lineste';

create_figure('Lines');

plot(linedat, 'o-', 'LineWidth', 3, 'MarkerSize', 12, 'MarkerFaceColor', [.5 .5 .5])

set(gca, 'XLim', [.5 size(anova_dat, 2)+.5], 'XTick', [1:size(anova_dat, 2)]); xlabel('Factor 2');

for j = 1:size(anova_dat, 2), f1names{j} = ['F1: Level ' num2str(j)]; end

legend(f1names)

end

