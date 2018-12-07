function BF = estimateBayesFactor(stat_image,varargin)
% Compute voxel-wise Bayes Factors from a statistic image object.
% This code is a wrapper function for code from Sam Schwarzkopf at UCL
% (see https://figshare.com/articles/Bayes_Factors_Matlab_functions/1357917),
% it computes Bayes Factors for t-tests (Rouder et al., 2009), simple 
% correlations (Wetzels et al., 2012), and proportions (http://pcl.missouri.edu/bf-binomial)
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
% ..

if isempty(varargin)
    input_type='t';
else
   input_type=varargin{1}; 
end

if ~any(strcmp(input_type,{'t','prop','r'}));
   error('Wrong input type, please specify ''t'', ''prop'', or ''r'''); 
end

if strcmp(input_type,'t')
    
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