function [cl] = extract_ind_peak(imnames,cl,varargin)
% function [clusters] = extract_ind_peak(imnames [can be empty],clusters,[vols])
%
% this function gets individual spatial peaks from clusters.
% input: 
%	imnames: a matrix of image names, in spm_list_files output format
%            if imnames is empty, enter full data (vols) as 3rd argument
%            or else this program will use cl.all_data to get data 
%   clusters: see tor_extract_rois
%
% 11/2/03 by Tor Wager
%
% NOTE (WARNING): WORKS ON XYZ VOXEL COORDINATES - NO TRANSFORMATION TO DIFFERENT SPACES
% see transform_coordinates.m for this.
%
% Functions called
%   C:\matlabR12\toolbox\matlab\datatypes\squeeze.m
%   c:\tor_scripts\voistatutility\nanmean.m
%   (calls other spm functions)
%
% see also cluster_manova

verbose = 1;


 	% ----------------------------------------------------------------------------------
	% load image files, if possible
	% ----------------------------------------------------------------------------------   
    if ~isempty(imnames),
        if size(imnames,1) < 60,
            % Image loading a la SPM, and manual extraction.
            V = spm_vol(imnames);
            vols = spm_read_vols(V);
        end
    elseif length(varargin) > 0
        vols = varargin{1};
    end
    
    
for i = 1:length(cl)
   
		% ----------------------------------------------------------------------------------
		% get the data for the cluster
		% ----------------------------------------------------------------------------------
		O.coords = cl(i).XYZ';
        mm = cl(i).XYZmm';
        all_data = [];
        
        if exist('vols') == 1
            for co = 1:size(O.coords,1)
                all_data(:,co) = squeeze(vols(O.coords(co,1),O.coords(co,2),O.coords(co,3),:));
            end
        else
            if i == 1,disp('Using data from clusters.all_data'),end
            all_data = cl(i).all_data;
        end
        
        for subj = 1:size(all_data,1)
            wh = find(abs(all_data(subj,:)) == max(abs(all_data(subj,:))));
            
            % randomly select if multiple peaks
            if length(wh) > 1
                warning('Multiple matches to max!')
                tmp = randn(size(wh)); wh = wh(tmp==max(tmp));
            end
            
            if isempty(wh)
                warning(['Cluster ' num2str(i) ': No individual peak data for subject ' num2str(subj) '?'])
                cl(i).peakXYZ(subj,:) = [NaN NaN NaN];
                cl(i).peakXYZmm(subj,:) = [NaN NaN NaN];
                cl(i).peakdata(subj,1) = NaN;
            else
                cl(i).peakXYZ(subj,:) = O.coords(wh,:);
                cl(i).peakXYZmm(subj,:) = mm(wh,:);
                cl(i).peakdata(subj,1) = all_data(subj,wh);
            end
        
        end
        
 end
                
        


return

