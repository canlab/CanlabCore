function out = ICC(cse,typ,dat)
    % iccvalue = ICC([1 to 6],['single' or 'k'], data matrix)
    %
%function to work out ICCs according to shrout & fleiss' schema (Shrout PE,
%Fleiss JL. Intraclass correlations: uses in assessing rater reliability.
%Psychol Bull. 1979;86:420-428).
%
% Modified 10/09 by Tor Wager; minor bug fix and changes to documentation
%
% 'dat' is data whose *columns* represent k different raters (judges) & whose
% *rows* represent n different cases or targets being measured. Each target
% is assumed to be a random sample from a population of targets.
%
% 'cse' is either 1,2,3. 'cse' is: 1 if each target is measured by a
% different set of raters from a population of raters, 2 if each target is
% measured by the same raters, but that these raters are sampled from a
% population of raters, 3 if each target is measured by the same raters and
% these raters are the only raters of interest.
%
% 'typ' is either 'single' or 'k' & denotes whether the ICC is based on a
% single measurement or on an average of k measurements, where k = the
% number of ratings/raters.
%
% This has been tested using the example data in the paper by shrout & fleiss.
% 
% Example: out = ICC(3,'k',S_Fdata)
% returns ICC(3,k) of data 'S_Fdata' to double 'out'.
%
% Kevin Brownhill, Imaging Sciences, KCL, London kevin.brownhill@kcl.ac.uk
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Additional documentation
% iccvalue = ICC([1 to 6],['single' or 'k'], data matrix)
%
% Here, columns are 'judges', or more generally, 'measures' that are
% usually ideally intercorrelated.  Rows are items being assessed.
% The ICC assesses the proportion of variance attributed to the items, 
%  shared across measures.
% As the correlation between the measures grows, the icc grows.
% Another way of saying this is that if the rows are consistently different
% across measures, the icc will be high.
% 
% Think of rows as criminals, and columns as judges. The data values are 'guilt scores',
% where higher is more guilty.  If all the judges agree, the most guilty
% cases will be rated as most guilty by all judges, and the icc will be
% high. This is actually consistent with Case 2 or 3 in Shrout and Fleiss.
%
% In Case 1, the columns don't have any real meaning, as there are
% different 'judges' for each row, and variance components due to judge
% cannot be separated from error and the judge x target interaction.
% In Case 2 and 3, they are crossed. Case 2 treats judge as a random
% effect, whereas Case 3 treats judge as a fixed effect.
%
% If the data were an individual differences study of cognitive performance, 
% then the rows would be subjects, and the columns would be tests.
% A high icc would indicate a high correlation across the tests, which
% indicates that subjects are reliably different from one another, i.e.,
% that a large proportion of the total variance is related to subject.
% In such a case, as tests are fixed entities, then Case 3 might be
% appropriate.  
%
% Cronbach's alpha is equal to ICC(3, k) - case 3, k
% This assumes no target x rater interaction
%
% More examples:
%
%dat = mvnrnd([1 1 1], [1 .5 .5; .5 1 .5; .5 .5 1], 50); whos dat
% corrcoef(dat)
% ri = ICC(2, 'k', dat)
% dat = mvnrnd([1 1 1], [1 .9 .9; .9 1 .9; .9 .9 1], 50); whos dat
% corrcoef(dat)
% ri = ICC(2, 'k', dat)
%
% In the example below, judges (measures) have systematically different
% means, and the ICC values are different.  ICC(1, 1) is low because judge
% is not considered as a source of variance.  ICC(2, 1) is higher, but
% intermediate, because judge is considered as a random effect and modeled,
% but we want to generalize to new judges.  ICC(3, 1) is highest, because
% judge is modeled 
% dat = mvnrnd([1 2 3], [1 .5 .5; .5 1 .5; .5 .5 1], 50); whos dat
% ri = ICC(1, 'single', dat)
% ri = ICC(2, 'single', dat)
% ri = ICC(3, 'single', dat)


%number of raters/ratings
k = size(dat,2);
%number of targets
n = size(dat,1);
%mean per target
mpt = mean(dat,2);
%mean per rater/rating
mpr = mean(dat);
%get total mean
tm = mean(mpt);
%within target sum sqrs
WSS = sum(sum(bsxfun(@minus,dat,mpt).^2));
%within target mean sqrs
WMS = WSS / (n * (k - 1));
%between rater sum sqrs
RSS = sum((mpr - tm).^2) * n;
%between rater mean sqrs
RMS = RSS / (k - 1);
% %get total sum sqrs
% TSS = sum(sum((dat - tm).^2));
%between target sum sqrs
BSS = sum((mpt - tm).^2) * k;
%between targets mean squares
BMS = BSS / (n - 1);
%residual sum of squares
ESS = WSS - RSS;
%residual mean sqrs
EMS = ESS / ((k - 1) * (n - 1));
switch cse
    case 1
        switch typ
            case 'single'
                out = (BMS - WMS) / (BMS + (k - 1) * WMS);
            case 'k'
                out = (BMS - WMS) / BMS;
            otherwise
               error('Wrong value for input typ') 
        end
    case 2
        switch typ
            case 'single'
                out = (BMS - EMS) / (BMS + (k - 1) * EMS + k * (RMS - EMS) / n);
            case 'k'
                out = (BMS - EMS) / (BMS + (RMS - EMS) / n);
            otherwise
               error('Wrong value for input typ') 
        end
    case 3
        switch typ
            case 'single'
                out = (BMS - EMS) / (BMS + (k - 1) * EMS);
            case 'k'
                out = (BMS - EMS) / BMS;
            otherwise
               error('Wrong value for input typ') 
        end
    otherwise
        error('Wrong value for input cse')
end