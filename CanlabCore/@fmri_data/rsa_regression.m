function stats = rsa_regression(obj,design,sudy,varargin)
% Representational similarity-based analysis, including inferences about a stimulus/task model 
% Constructs a rep. dissim. matrix (RDM) based on spatial covariance across images.
% Takes a stimulus/experimental design (design), which is a set of binary regressors specifying groups of images. 
% Include a second grouping variable (study) of integers coding for blocks of images to resample within.
% The function regresses the RSA matrix on the design and returns bootstrapped 
% statistics on each design regressor.  
%   This method can be used to analyze generalizability across constructs and studies. 
% e.g., see reference to Kragel et al. 2018 below. In this case, this method tests whether
% patterns of activity in an fmri data object are generalizable across
% different groupings specified in the matrix 'design'. 'study' is an integer
% valued vector grouping images according to study. This script does not
% currently implement repeated measures designs (only include one image per
% subject). See Kragel et al. (2018) Nature Neuroscience for details.
%
% Usage:
% ::
%
%    stats = rsa_regression(obj,design);
%
% ..
%     Author and copyright information:
%
%     Copyright (C) 2018 Phil Kragel
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
%   **obj:**
%        An image object with one or more images loaded
%
%   **design:**
%        A binary matrix where each column corresponds to a grouping
%        variable and each row corresponds to an observation.
%
% :Optional inputs:
%
%	**'average_Euclidean'**
%		Use Euclidean distance of the response averaged across all voxels (univariate distance)
%	**'euclidean'**
%		Use Euclidean distance computed using all voxels
%	**'cosine'**
%		Compute distance using cosine similarity computed using all voxelst
%	**'correlation'**
%		Compute distance using correlation distance (1 - Pearson R) computed using all voxels (default)
% :Outputs:
%   **stats:**
%    Structure including: 
%      gen_index  - sample estimate of generalization indices
%      bs_gen_index - bootstrap distribution for each generalization index
%      ste  - standard error of generalization index
%      Z  - Z score for each parameter
%      p  - p value for each parameter
%      sig - significance test at FDR q < .05
%
% :Examples:
% Example 1:
% ----------------------------------------------------------------------
% Download images from Kragel et al. 2018 Nature Neuroscience
% Kragel, P. A., Kano, M., Van Oudenhove, L., Ly, H. G., Dupont, P., Rubio, A., ? Wager, T. D. (2018). Generalizable representations of pain, cognitive control, and negative emotion in medial frontal cortex. Nature Neuroscience, 21(2), 283?289. doi:10.1038/s41593-017-0051-7 
% 270 subject-level images systematically sampled from 18 studies across 3 domains
% [data_obj, names] = load_image_set('kragel18_alldata');
% 
% Amygdala=select_atlas_subset(load_atlas('Canlab2018'),{'Amy'});
% masked_dat=apply_mask(data_obj,Amygdala);
% 
% dsgn = condf2indic(ceil(data_obj.metadata_table.Studynumber/6));
% study = data_obj.metadata_table.Studynumber;
% stats = rsa_regression(masked_dat,dsgn,study);
%
% ..
%    Programmers' notes:
% List dates and changes here, and author of changes
%    01/08/2018, created function, Phil Kragel
% ..


%% convert design matrix into dissimilarity matrix
modelRDM=[];
for i=1:size(design,2) %for each specified grouping
   modelRDM(:,i)=pdist(design(:,i),'seuclidean'); %#ok<*AGROW> %compute squared euclidean distance
   modelRDM(:,i)=100000*modelRDM(:,i)/sum(modelRDM(:,i)); %normalize and scale
end

%supress distance warning - 0 value distances get removed...
warning('off','stats:pdist:ConstantPoints')

if any(strcmp([varargin{:}],'average_Euclidean'))
    brainRDM=pdist(nanmean(obj.dat)'); %distance is based on region average
elseif any(strcmp([varargin{:}],'euclidean'))
    brainRDM=pdist(obj.dat'); %euclidean

elseif any(strcmp([varargin{:}],'cosine'))
    brainRDM=pdist(obj.dat','cosine'); %euclidean

else %default to correlation
    brainRDM=pdist(obj.dat','correlation');
end
brainRDM(brainRDM<.00001)=NaN;
gen_index= glmfit([ones(length(modelRDM),1) double(modelRDM)],brainRDM','normal','constant','off');

num_it=1000;
bs_gen_index=zeros(num_it,size(gen_index,1));
parfor it=1:num_it
bs_gen_index(it,:) = random_resample_within_study(modelRDM,obj,study,varargin);
end

%compute stats from bootstrap distribution (normal approx)  
b_ste = squeeze(nanstd(bs_gen_index));
b_mean = squeeze(nanmean(bs_gen_index));
b_ste(b_ste == 0) = Inf;
b_Z = b_mean ./ b_ste;
b_P = 2 * (1 - normcdf(abs(b_Z)));

% assign stats to output
stats.gen_index=gen_index;
stats.Z=b_Z;
stats.p=b_P;
stats.sig=b_P<FDR(b_P,.05);
stats.ste=b_ste;
stats.bs_gen_index=bs_gen_index;
stats.RDM=squareform(brainRDM);

end % main function

function beta = random_resample_within_study(modelRDM,obj,study,varargin)

bs_inds=zeros(size(obj.dat,2),1); %initialize

for i=1:max(study) %for each study
    
    study_inds=find(study==i); %find images that belong to this study
    bs_inds(study_inds)=study_inds(randi([1,length(study_inds)],1,length(study_inds))); %randomly resample with replacement

end
%supress distance warning - 0 value distances get removed...
warning('off','stats:pdist:ConstantPoints')

resampled_data=obj.dat(:,bs_inds)'; %random subsample and transpose for correlation


if any(strcmp([varargin{:}],'average_Euclidean'))
    brainRDM=pdist(nanmean(resampled_data')'); %distance is based on region average
elseif any(strcmp([varargin{:}],'euclidean'))
    brainRDM=pdist(resampled_data); %euclidean

elseif any(strcmp([varargin{:}],'cosine'))
    brainRDM=pdist(resampled_data,'cosine'); %euclidean

else %default to correlation
    brainRDM=pdist(resampled_data,'correlation');
end

brainRDM(brainRDM<.00001)=NaN; %exclude pairs that have exactly 0 distance


beta = glmfit([ones(length(modelRDM),1) double(modelRDM)],brainRDM','normal','constant','off'); %get parameter estimates using ols



end