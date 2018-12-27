function [fitresult, gof] = evaluate_spatial_scale(data_obj,parcel_obj,varargin)
% Evaluate information coding at multiple spatial scales, using a prespecified parcellation
% and the predict function (and associated optional inputs). This function evaluates the
% spatial scale of information coding by constructing predictive models using random draws
% of features from different parcels to test whether information is contained 1) at the
% whole-brain scale, 2) within individual parcels, or 3) across all parcels (may be similar
% to whole-brain if parcels cover most cortical and subcortical areas)
% Usage:
% ::
%
%    [fitresult, gof] = evaluate_spatial_scale(obj, parel_obj, varargin);
%
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
%   **data_obj:**
%   **parcel_obj:**
%
% :Optional inputs:
%
%   **colors:**
%       Cell array of RGB values that indicate colors for each parcel
%
%   **labels:**
%       Cell array of strings that indicate names for each parcel
%
%   **num_voxels:**
%      Vector of integers that indicate how many voxels to randomly sample
%      these should be smaller than the total number of voxels in the
%      smallest parcel
%
%
%   **do_plot:**
%      Automatically make plot of sampling curves
%
%   **render_brain**
%      Make a brain rendering of parcels on plot
%
%
% :Outputs:
%
%   **fitresult:**
%        cell array of fit objects for each parcel
%
%   **gof:**
%        cell array of with performance metrics for each parcel
%
% :Examples:
% ::
%
%  % Example 1: Load bmrk3 data and predict pain intensity sampling
%  within/across different resting-state networks
% % -------------------------------------------------------------------------
%   canlab_help_set_up_pain_prediction_walkthrough;
%   [mask_obj, networknames, imagenames] = load_image_set('bucknerlab');
%   mycolors={[120 18 134]/255 [70 130 180]/255  [0 118 14]/255  [196 58 250]/255  [220 248 164]/255  [230 148 34]/255  [205 62 78]/255 };
%   [fitresult, gof] = evaluate_spatial_scale(image_obj,mask_obj,'cv_lassopcr', 'nfolds',[rem(subject_id,2)+1],'verbose',0,'labels',networknames,'colors',mycolors);
%
%
%  % Example 1: Load bmrk3 data and predict pain intensity sampling
%  within/across different brainstem/thalamus regions
% % -------------------------------------------------------------------------
%   canlab_help_set_up_pain_prediction_walkthrough;
%   [mask_obj, networknames, imagenames] = load_atlas('brainstem');
%   [fitresult, gof] = evaluate_spatial_scale(image_obj,mask_obj,'cv_lassopcr', 'nfolds',[rem(subject_id,2)+1],'verbose',0,'labels',networknames,'colors',mycolors);
%
%
% :References:
%   Kragel, P. A., Koban, L., Barrett, L. F., & Wager, T. D. (2018). Representation, pattern information, and brain signatures: from neurons to neuroimaging. Neuron, 99(2), 257-273.
%   Kragel, P. A., Reddan, M., LaBar, K. S., & Wager, T. D. (2018). Emotion Schemas are Embedded in the Human Visual System. bioRxiv, 470237.
%
% :See also:
%
% predict, load_atlas

% ..
%    Programmers' notes:
%    List dates and changes here, and author of changes
%
%   Phil: created, Dec 2018
%
%   TODOs: - change outcome metric to optionally incude accuracy/AUC/etc.
%          - change outcome to optionally reflect % of noise ceiling
%          - account for regions that have small number of voxels
% ..

%if we have a 1d atlas object convert to matrix for
parcel_obj=resample_space(parcel_obj,data_obj,'nearest');

if size(parcel_obj.dat,2)==1
    parcel_obj.dat=condf2indic(parcel_obj.dat);
    parcel_obj.dat=parcel_obj.dat(:,2:end);

end

% defaults
num_iterations = 10;
num_parcels = size(parcel_obj.dat,2);
num_vox = [ 50 150 250 500 750 1000 1500 5000 10000];

do_plot=true;
render_brain=true;

mycolors=scn_standard_colors(num_parcels);

atlas_labels=cell(num_parcels,1);
for i=1:length(atlas_labels)
    atlas_labels{i}='';
end

% optional inputs
if any(strcmp(varargin,'colors'))
    mycolors=varargin{find(strcmp(varargin,'colors'))+1};
end

if any(strcmp(varargin,'labels'))
    atlas_labels=varargin{find(strcmp(varargin,'labels'))+1};
end

if any(strcmp(varargin,'num_voxels'))
    num_vox=varargin{find(strcmp(varargin,'num_voxels'))+1};
end

if any(strcmp(varargin,'iterations'))
    num_iterations=varargin{find(strcmp(varargin,'iterations'))+1};
end

%% adjust number of voxels to sample based on size of parcels

for n=1:num_parcels %for each parcel
    mask=parcel_obj;
    mask.dat=mask.dat(:,n);
    nf(n)=length(find(mask.dat>0));
end

num_vox(num_vox>(min(nf)))=[];
%% test with increasing number of features randomly sampled from each parcel
num_feats_within = nan(num_parcels,length(num_vox)+1);
pred_outcome_r=nan(num_iterations,num_parcels,length(num_vox)+1);
for n=1:num_parcels %for each parcel
    mask=parcel_obj;
    mask.dat=mask.dat(:,n);
    masked_data=apply_mask(data_obj,mask);

    num_feats_within(n,:)=[num_vox size(masked_data.dat,1)];
    
    
    
    for it=1:num_iterations
        temp_inds=randperm(size(masked_data.dat,1));
        cs=1;
        for nf=num_feats_within(n,:)
            try
                subsampled_masked_data=masked_data;
                
                if nf==size(subsampled_masked_data.dat,1) %if subsampling all data, do nothing
                else
                    subsampled_masked_data.dat=subsampled_masked_data.dat(temp_inds(1:nf),:); %otherwise grab first 1:nf reordered features
                end
                
                if ~any(strcmp(varargin,'nfolds')) %if cv not specified do 10-fold
                    [~, stats]=predict(subsampled_masked_data,varargin{:},'nfolds',10);
                else
                    [~, stats]=predict(subsampled_masked_data,varargin{:});
                end
                
                try
                    pred_outcome_r(it,n,cs)=stats.pred_outcome_r;
                    
                catch
                    pred_outcome_r(it,n,cs)=1-stats.cverr;
                    
                end
                cs=cs+1;
            catch
            end
            
        end
    end
end

%% random samples from entire brain

num_feats_full=[num_vox size(data_obj.dat,1)];
pred_outcome_r_full=nan(num_iterations,length(num_vox)+1);

for it=1:num_iterations
    temp_inds=randperm(size(data_obj.dat,1));
    cs=1;
    for nf=num_feats_full
        try
            subsampled_masked_data=data_obj;
            
            if nf==size(subsampled_masked_data.dat,1) %if subsampling all data, do nothing
            else
                subsampled_masked_data.dat=subsampled_masked_data.dat(temp_inds(1:nf),:); %otherwise grab first 1:nf reordered features
            end
            
            if ~any(strcmp(varargin,'nfolds')) %if cv not specified do 10-fold
                [~, stats]=predict(subsampled_masked_data,varargin{:},'nfolds',10);
            else
                [~, stats]=predict(subsampled_masked_data,varargin{:});
            end
            
            try
                pred_outcome_r_full(it,cs)=stats.pred_outcome_r;
                
            catch
                pred_outcome_r_full(it,cs)=1-stats.cverr;
                
            end
            
            
            cs=cs+1;
        catch
        end
        
    end
end


%% random samples across all parcels

mask=parcel_obj;

mask.dat=double(any(mask.dat')');

masked_data=apply_mask(data_obj,mask);
num_feats_all_parcels=[num_vox size(masked_data.dat,1)];
pred_outcome_r_all_parcels=nan(num_iterations,length(nf));

for it=1:num_iterations
    temp_inds=randperm(size(masked_data.dat,1));
    cs=1;
    for nf=num_feats_all_parcels
        try
            subsampled_masked_data=masked_data;
            
            if nf==size(subsampled_masked_data.dat,1) %if subsampling all data, do nothing
            else
                subsampled_masked_data.dat=subsampled_masked_data.dat(temp_inds(1:nf),:); %otherwise grab first 1:nf reordered features
            end
            
            if ~any(strcmp(varargin,'nfolds')) %if cv not specified do 10-fold
                [~, stats]=predict(subsampled_masked_data,varargin{:},'nfolds',10);
            else
                [~, stats]=predict(subsampled_masked_data,varargin{:});
            end
            
            
            
            try
                pred_outcome_r_all_parcels(it,cs)=stats.pred_outcome_r;
                
            catch
                pred_outcome_r_all_parcels(it,cs)=1-stats.cverr;
                
            end
            
            
            
            cs=cs+1;
        catch
        end
        
    end
end



%% do curve fitting and make plots

if do_plot
    
    [fitresult{1}, gof{1},x_mean_wb,y_mean_wb] = createFit(num_feats_full(:)', mean(pred_outcome_r_full),'color',[.5 .5 .5],'linewidth',2);
    [fitresult{2}, gof{2},x_mean_parcels,y_mean_parcels] = createFit(num_feats_all_parcels(:)', mean(pred_outcome_r_all_parcels),'color',[.1 .1 .1],'linewidth',2);
    
    for r=1:num_parcels
        
        mean_r_within=squeeze(mean(pred_outcome_r(:,r,1:size(pred_outcome_r_all_parcels,2))));
        [fitresult{2+r}, gof{2+r},x_mean_parcel{r},y_mean_parcel{r}] = createFit(num_feats_within(r,:)', mean_r_within(:),'color',mycolors{r},'linewidth',2);
        
    end
    
    
    for it=1:num_iterations
        [~, ~,x_wb(it,:),y_wb(it,:)] = createFit(num_feats_full(:)', pred_outcome_r_full(it,:),'color',[.5 .5 .5],'linewidth',2);
        [~, ~,x_all_parcels(it,:),y_all_parcels(it,:)] = createFit(num_feats_all_parcels(:)', pred_outcome_r_all_parcels(it,:),'color',[.1 .1 .1],'linewidth',2);
        
        for r=1:num_parcels
            
            r_within=squeeze((pred_outcome_r(it,r,1:size(pred_outcome_r_all_parcels,2))));
            [~, ~,x{it,r},y{it,r}] = createFit(num_feats_within(r,:)', r_within(:),'color',mycolors{r},'linewidth',2);
            
        end
    end
    
    figure;
    hold on;
    boundedline(mean(x_wb),mean(y_wb), std(y_wb),'alpha', 'cmap',[.5 .5 .5]);
    boundedline(mean(x_all_parcels), mean(y_all_parcels), std(y_all_parcels),'alpha', 'cmap',[.1 .1 .1]);
    
    for r=1:num_parcels
        clear to_plot_x to_plot_y
        for i=1:size(x,1)
            to_plot_x(i,:)=x{i,r};
            to_plot_y(i,:)=y{i,r};
        end
        boundedline(mean(to_plot_x), mean(to_plot_y), std(to_plot_y),'alpha', 'cmap',mycolors{r});
        
    end
    h=findobj(gca,'type','line');
    set(h,'linewidth',2);
    set(gca,'XScale','log')
    grid on
    legend(h(1:end),{atlas_labels{:} 'All Parcels' 'Whole-brain'});
    
    xlabel 'Number of voxels'
    ylabel 'Model performance'
    
    
    % add brain render to plot
    
    axis tight;
    if render_brain
        
        
        clear o2
        o2=fmridisplay();
        o2 = surface(o2, 'axes', [0.5 0 .5 .5], 'direction', 'surface left', 'orientation', 'lateral','trans');
        o2 = surface(o2, 'axes', [0.5 .5 .4 .4], 'direction', 'surface right', 'orientation', 'lateral');
        
        pname = 'R.pial_MSMAll_2_d41_WRN_DeDrift.32k.mat'; % from Glasser_et_al_2016_HCP
        
        try
            load(which(pname))
        catch
            error('Unable to find ''R.pial_MSMAll_2_d41_WRN_DeDrift.32k.mat'' please check your path');
        end
        
        set(o2.surface{1}.object_handle(1),'FaceVertexCData',cdata)
        set(o2.surface{2}.object_handle(1),'FaceVertexCData',cdata)
        colormap gray
        
        
        %add mask regions
        for i=1:size(parcel_obj.dat,2)
            tv=parcel_obj;
            tv.dat=tv.dat(:,i);
            
            o2=addblobs(o2,region(tv),'splitcolor',{mycolors{i} mycolors{i} mycolors{i} mycolors{i}});
            
        end
        
        %clear brainstem
        delete(o2.surface{1}.object_handle(2))
        delete(o2.surface{2}.object_handle(2))
        
        %position axes
        ha=findobj(gcf,'type','axes');
        newpos=[[.315 .1 .3 .3];[.61 .1 .3 .3]];
        for i=1:2
            set(ha(i),'Position',newpos(i,:))
        end
        set(ha([1]),'view',[-90 0]);
        set(ha([2]),'view',[90 0]);
        
        
    end
end

