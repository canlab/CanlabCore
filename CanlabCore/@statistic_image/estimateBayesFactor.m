function BF = estimateBayesFactor(stat_image,varargin)
% Compute voxel-wise Bayes Factors from a statistic image object.
%
% This code is a wrapper function for code from Sam Schwarzkopf at UCL
% (see https://figshare.com/articles/Bayes_Factors_Matlab_functions/1357917),
% it computes Bayes Factors for t-tests (Rouder et al., 2009), simple 
% correlations (Wetzels et al., 2012), and proportions (http://pcl.missouri.edu/bf-binomial)
% See below for more information. 
%
% Usage:
% ::
%
%    BF = estimateBayesFactor(stat_image,varargin)
%
% ..
%     Author and copyright information:
%
%     Copyright (C) 2018 Philip Kragel
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
%
% :Background:
%
% Bayesian t-tests
% -----------------------------------------------------------------
% We use the BF Bayes' Factor Matlab toolbox, based on:
% 
%     Rouder et al., 2009     % one-sample t-test, t1smpbf.m
%     Rouder et al., 2009     % binomial test, binombf.m
%     Wetzels et al., 2012    % Pearson's r, corrbf.m
%     Boekel et al., 2014     % test for replicating Pearson's r in same direction
% 
% Implemented by Sam Schwarzkopf, UCL
%
% Rouder (2009) derived a formula to calculate Bayes Factors for a
% one-sample t-test, a common test statistic in neuroimaging, particularly 
% when testing contrast or 2nd-level (across-participant) summary statistics 
% in a two-level hierarchical model. They also provided a web application:
% http://pcl.missouri.edu/bf-one-sample
%
% Gönen et al. (2005) provided the corresponding equation
% for the unit-information Bayes factor. Liang et al. (2008)
% provided the corresponding JZS Bayes factors for testing
% slopes in a regression model.
% 
% As Rouder et al. point out: "Researchers need only provide the sample size N and
% the observed t value. There is no need to input raw data.
% The integration is over a single dimension and is computationally
% straightforward."
%
% This function iterates over voxels to calculate bayes factors for a map.
%
% Assumptions about prior distributions and parameter choices
% -----------------------------------------------------------------
% With Bayesian analysis, one must specify the distributional form of the
% prior belief about the effect size, which is integrated with evidence
% from the data to estimate a posterior probabilities of both null and
% alternative hypotheses. The Bayes Factor (BF) is a ratio of these (see below).
% 
% Different choices of prior distribution and effect size will thus yield
% different results for BFs, but there are some standard, reasonable
% choices. In addition, BFs are often not very sensitive to reasonable variation
% in priors, so it is reasonable to use a default choice for many applications.
% 
% The tests in the BF toolbox use default scaling values for prior distributions, and the
% Jeffrey-Zellner-Siow Prior (JZS, Cauchy distribution on effect size).
% This is standard, widely used prior. The JZS prior has heavier tails than
% the normal distribution, so does not penalize very large effect sizes as
% much as the Normal prior (large effects can also be unlikely given an
% assumption of a particular alternative distribution with moderate effect
% sizes). The JZS prior is there more noninformative than Normal prior.
%
% One additional choice is the choice of prior effect size under the
% alternative hypothesis, which is determined by the scale factor (r in the
% BF toolbox). This is a noncentrality parameter that governs the expected
% effect. If the observed effect is much *larger* than the belief about the
% alternative, evidence for the alternative will actually go down!
% However, in this case, the BFs will likely still very strongly favor the
% alternative, so this is not a problem with the JZS prior.
% In Rouder et al. 2009, the default r was 1.0, but it was changed in 2015
% to be 0.707, which is a reasonable choice. We use this default here.
%
% Interpreting Bayes Factors:
% -----------------------------------------------------------------
% Bayes Factors > 1 provide evidence in favor of the Alternative (an
% effect), and those < 1 provide evidence in favor of the Null (no effect).
%
% For example, bf = t1smpbf(3, 50); yields bf = 7.92, or about 8:1 in favor of the alternative.
% bf = bf = t1smpbf(2, 50); yields bf = 0.96, or about 1.04:1 in favor of the null.
%
% The BF toolbox returns BF values in their original scaling. the
% fmri_data.estimateBayesFactor method scales the BFs by 2*ln(BF),
% so that values of 0 indicate equal support for Null and Alternative,
% positive values support the Alternative, and negative values support the
% Null. (See Kass and Raftery 1995)
%
% These are returned in a statistic_image object BF, whose .dat field
% contains 2*ln(BF) values for each voxel.  A value of about 4.6 indicates
% a BF of 10, or 10:1 evidence in favor of the Alternative, which is a typical cutoff. 
% A value of about 6 indicates 20:1 evidence in favor of the Alternative.
% 
% 
% :Inputs:
%
%   **stat_image:**
%        A statistic image that either contains 1) a t-map, 2) a
%        correlation map, or 3) a map of counts for testing proportions
%
%   **input_type:**
%        A string specifying the type of map provided for now, 
%        options are 't', 'r', and 'prop'
%
%
% :Outputs:
%
%   **BF:**
%        fmri_data object with voxel-wise bayes factors (scaled as 2 ln BF)
%        these can be thresholded according to Kass and Raftery 1995
%
% :Examples:
% ::
%
%  % Example 1: Emotion Regulation data and perform 1-sample t-test
%
%        dat=load_image_set('emotionreg');
%        t=ttest(dat);
%        BF_tstat=estimateBF(t,'t');
%        BF_tstat_th = threshold(BF_tstat, [-6 6], 'raw-outside');
%        orthviews(BF_tstat_th);
%
%
%
%  % Example 1: Emotion Regulation data and test of proportions
%
%        dat=load_image_set('emotionreg');
%        t=ttest(dat);
%        prop=t; %initialize stats object from t-test output
%        prop.dat=sum(dat.dat'>0)'; 
%        BF_prop=estimateBF(prop,'prop'); %estimate BF
%        BF_prop_th = threshold(BF_prop, [-6 6], 'raw-outside');
%        orthviews(BF_prop_th);
% 
%
% :See also:
%
% stat_image, fmri_data

% ..
%    Programmers' notes:
%    List dates and changes here, and author of changes
%
%   Phil: created, Dec 2018
%
%   TODO: add optional input for base-rate when testing proportions
%
%   11/21/2019 : Tor Wager added documentation and explanation to help.
% ..

if isempty(varargin)
    input_type='t';
else
   input_type=varargin{1}; 
end

if ~any(strcmp(input_type,{'t','prop','r'}))
   error('Wrong input type, please specify ''t'', ''prop'', or ''r'''); 
end

if strcmp(input_type,'t')
    

    % Calculates JZS Bayes Factor for a one-sample t-test given t and sample size n.
    % The optional input r is the scale factor which defaults to 0.707.
    % This quantifies the evidence in favour of the alternative hypothesis.
    % See Rouder et al, 2009, Psychon Bull Rev for details.

    parfor i=1:length(stat_image.N)
        bf10(i,1) = t1smpbf(stat_image.dat(i),stat_image.N(i));
    end
    
elseif strcmp(input_type,'r')
    
    
    parfor i=1:length(stat_image.N)
        bf10(i,1) = corrbf(stat_image.dat(i),stat_image.N(i));
    end
    
elseif strcmp(input_type,'prop')
    
    
    parfor i=1:length(stat_image.N)
        bf10(i,1) = binombf(stat_image.dat(i),stat_image.N(i),.5); %assume baserate of .5
    end
    
else
    

    
end

BF=stat_image;
BF.dat=2*log(bf10); % scale according to kass and raftery 1995
BF.type='BF';
BF.p=[];
% BF.dat(abs(BF.dat)<2)=0; %threshold using cutoff of 2