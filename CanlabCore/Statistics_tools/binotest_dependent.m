function [varargout] = binotest_dependent(X, Po)
% This function runs several different types of tests on dependent binomial data.
%
% :Usage:
% ::
%
%     [varargout] = binotest_dependent(X, Po)
%
% Overall, it tests the number of "hits" for each subject (row in X) against a null-hypothesis
% proportion p, across all subjects using a Z-test (two-tailed).
% The second level null hypothesis should be approximated by a normal
% distribution with a mean of p.  This approach assumes that each subject
% has an equal number of independent Bernoulli trials (columns in X) and
% that the number of subjects exceeds n=20 the test will be more accurate as n -> infinity.
% Also, calculates tests for each separate trial (e.g., subject columns),
% and the difference between proportions (two proportion z-test).
%
% :Inputs:
%
%   **X:**
%        X is a matrix of "hits" and "misses", coded as 1s and 0s.
%        where rows = subjects and columns = observations within
%        subject
%
%   **Po:**
%        Po is the null hypothesis proportion of "hits", e.g., often p = 0.5
%
% :Outputs:
%
%   **RES [1:5]:**
%        a structure containing the output of the stats for the
%        z-test, includes the number of subject (N), number of overall 
%        hits (hits), the overall proportion of hits (prop), the standard 
%        deviation (SE), z-statistic (Z), and the two tailed p-value (pval) 
%        trial across all subjects.  Assumes independence
%
%   **RES1:**
%        Independent Single Interval Test for Column 1 (Column 1 only against Po)
%
%   **RES2:**
%        Independent Single Interval Test for Column 2 (Column 2 only against Po)
%
%   **RES3:**
%        Two proportion dependent difference z-test (Column 1 minus Column 2 against 0) 
%
%   **RES4:**
%        Dependent single-interval test (Mean of Column 1 and Column 2 against Po)
%
%   **RES5:**
%        Two proportion dependent addition z-test (Column 1 plus Column 2 against 2 * Po) 
%              (Similar to mean, not sure what this will be used for)
%
% :Examples:
% ::
%
%    [RES1, RES2, RES3, RES4, RES5] = binotest_dependent([1,1,1,1,0; 1,0,1,0,1]',.5)
%
% ..
%     Author and copyright information:
%
%     Copyright (C) 2014  Luke Chang & Tor Wager
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
% ..

% ..
%    Process Data
% ..

X = double(X); % just in case
[N, k] = size(X);  % number of subjects x dependent obs per person
A = X(:, 1);
B = X(:, 2);

% -------------------------------------------------------------------------
% Test if Data is in Correct Format
% -------------------------------------------------------------------------
if k > 2, error('NOT IMPLEMENTED FOR MORE THAN 2 DEP OBSERVATIONS PER UNIT'); end

u = unique(X(:));
if any(u ~= 1 & u ~= 0), error('X MUST HAVE VALUES OF 1 OR 0 ONLY.'); end

if N < 20
    warning(['Running dependent binomial test requires many subjects to approximate normal distribution.  You are only using ' num2str(length(unique(subject_id))) ' subjects.  Interpret results with caution.'])
end

% -------------------------------------------------------------------------
% Calculate Variance and Covariance
% -------------------------------------------------------------------------
% First, calculate variance of the quantity we want for a single trial,
% using Pa and Pb as estimates for the binomial parameters
% (This is not the SD of the SAMPLING distribution, but the SD of the
% distribution of interest (e.g., variance of sum, difference)
% Later, we will use this to get the SD for the SAMPLING distribution,
% which is used to construct confidence intervals and tests.

%Calculate Probabilities of each Subject Trial
Pa = sum(A) ./ N;  % P-hat for A
Pb = sum(B) ./ N;  % P-hat for B

% Calculate Variance for each Subject Trial
Va = Pa * (1 - Pa); % ./ N;
Vb = Pb * (1 - Pb); % ./ N;

% Calculate Covariance between Trials
% C = ( (sum(A .* B) ./ N) * (sum(~A .* ~B) ./ N) ) - ( (sum(A .* ~B) ./ N) * (sum(B .* ~A) ./ N) ); % % Agresti, p. 410 SAME as below
C = ( (sum(A .* B) ./ N) - (Pa * Pb) );   % Covariance - expected prop under independence

% Calculate Variance for Test Statistics adjusting for covariance
Vaplusb = Va + Vb - (2 * C);
Vameanb = ( Va + Vb - (2 * C) ) / 4;  % / 4 because V(aX) = a^2Var(X)
Vaminusb = Va + Vb + (2 * C);

% -------------------------------------------------------------------------
% Define Functions
% -------------------------------------------------------------------------

sdfunction = @(V, N) (V ./ N) .^ .5;  % Calculate Standard Deviation
zfunc = @(Px, Po, SD) (Px - Po) ./ SD; % Calculate Z-Value
pfunction_gauss = @(z) 2 .* min(normcdf(z, 0, 1), normcdf(z, 0, 1, 'upper')); %Calculate p-value for z-distribution

% pfunction_bino = @(N, Po, hits) 2 * min(binocdf(hits, N, Po), (1 - binocdf(hits - 1, N, Po))); %Calculate p-value for binomial distribution

% -------------------------------------------------------------------------
% Define Various Confidence Intervals - TO DO - See Agresti Book
% -------------------------------------------------------------------------

% switch  cimethod
%     case 'Wald'
% 
%         cifunc = @(z)
%         
%     case 'scoring'  % default
%         
%     otherwise
%         error('Unknown CI method. Check inputs.');
% end

% Now, we will use this to get the SD for the SAMPLING distribution,
% which is used to construct confidence intervals and tests.
% e.g., SDa = (Va ./ N) .^ .5;

% -------------------------------------------------------------------------
% Calculate test statistics
% -------------------------------------------------------------------------
% Use z-test (large-sample approximation; ok for N = 20+)

names = {'P1' 'P2' 'P1minusP2' 'PavgP1P2' 'P1plusP2' };
est = [Pa Pb Pa-Pb mean([Pa Pb]) Pa+Pb];
Vs = [Va Vb Vaminusb Vameanb Vaplusb];
Pnull = [Po Po 0 Po 2*Po];

for i = 1:length(names)
    SDs(i) = sdfunction(Vs(i), N);
    Z(i) = zfunc(est(i), Pnull(i), SDs(i));
    pval(i) = pfunction_gauss(Z(i));
    %ci(i) = add confidence intervals

    RES{i} = struct('n', N, 'hits', est(i).*N, 'prop', est(i), 'SE', SDs(i), 'Z', Z(i), 'p_val', pval(i));

    varargout{i} = RES{i};
end

end % function
