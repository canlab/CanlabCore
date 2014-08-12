% What is it?
% ------------------------------
% fmri_mask_image is a class of objects containing masks for neuroimaging
% data (not just fmri, really...)
%
% An object of this type is a structure that helps keep track of mask
% information, and what has happened to it, i.e., if it was resampled into
% the space of another image.
%
% Properties:
% ------------------------------
% dat
% dat_descrip = 'Mask data';
% 
% volInfo = an SPM-style volume info structure (see spm_vol) with a bit
% more information, from iimg_read_image.m
% volInfo_descrip = 'Volume and in-area mask info from iimg_read_img';
% 
% space_defining_image_name = char('');
% 
%history = {};
%        
% Methods include:
% ------------------------------
% resample_to_image_space
%
% Usage:
% ------------------------------
% Class constuctor method: Usage for creating objects of this type:
% - Loads first image in 4-D volume only
% - Returns values within mask (excludes zeros and NaNs)
% - entering 'implicit' as an optional keyword will use fmri_mask_image to get an implicit mask
%
% obj = fmri_mask_image(image_name, ['implicit'])
% This command loads in an image file stored in image_name and returns a 
% mask image object in obj with some fields filled in, including data.
%
%
%  Properties:
% dat: [517845x1 single]
% dat_descrip: [1x92 char]
% volInfo: [1x1 struct]
% volInfo_descrip: [1x47 char]
% space_defining_image_name: ''
% history: {}
% 
% To get a list in original full-image space of in-mask voxels, use this:
% indx = obj.volInfo.image_indx;
% indx(obj.volInfo.wh_inmask(obj.removed_voxels)) = false;



classdef fmri_mask_image < image_vector
    
    properties
        
        %dat
        %dat_descrip = 'Mask data';
        %history = {};
        
        %volInfo = struct('XYZ', [0 0 0]', 'XYZmm', [0 0 0]', 'M', zeros(4));
        volInfo_descrip = 'Volume and in-area mask info from iimg_read_img';
        
        space_defining_image_name = char('');
        
        % volInfo has this:
        %          fname
        %            mat
        %            dim
        %             dt
        %          pinfo
        %              n
        %        descrip
        %        private
        %           nvox
        %     image_indx
        
    end
    
    methods
        
        function obj = fmri_mask_image(image_name, varargin)
            
            % class constructor method
            % takes as input a mask image name
            
            % ---------------------------------
            % Create empty image_vector object, and return if no additional
            % arguments
            % ---------------------------------
            
            obj.dat = [];
            obj.volInfo = [];
            obj.image_names = [];
            obj.fullpath = char([]);
            obj.files_exist = false;
            obj.history = {};
            obj.volInfo_descrip = 'Empty';
            obj.space_defining_image_name = char('');
                
            if nargin == 0
                return
                
            elseif isa(image_name, 'image_vector')
                % map from existing image vector
                % copy legal fields
                
                N = fieldnames(image_name);
                Nmask = fieldnames(obj);
                for i = 1:length(Nmask)
                    if ~isempty(strmatch(Nmask{i}, N))
                        obj.(Nmask{i}) = image_name.(Nmask{i});
                    end
                end

                obj.dat = any(obj.dat, 2);

            elseif ischar(image_name)
                
            % initialize empty, return if no inputs
            obj.dat = [];
            obj.dat_descrip = 'Mask data';
            obj.volInfo = struct('XYZ', [0 0 0]', 'XYZmm', [0 0 0]', 'M', zeros(4));
            obj.volInfo_descrip = 'Volume and in-area mask info from iimg_read_img';
            obj.space_defining_image_name = char('');
            obj.history = {};
            
            if nargin == 0
                return
            end
            
            [obj.volInfo, obj.dat] = iimg_read_img(image_name, 2, 1, 1); % reads first vol only
            obj.dat = single(obj.dat(obj.volInfo.wh_inmask));
            
            obj.dat_descrip = image_name;
            obj.removed_voxels = false(size(obj.dat, 1), 1);
            
            doimplicit = strmatch('implicit', varargin);
            
            if doimplicit
                
                fprintf('Calculating implicit mask ');
                [dummy, dummy, nvox, is_inmask] = fmri_mask_thresh_canlab(obj);
                fprintf('%3.0f voxels in, %3.0f voxels out of mask.\n', nvox, sum(~is_inmask));
                
                obj = remove_empty(obj, ~is_inmask, []);
                
            end
            
%             if dofull
%                 maskobj.dat = ones(size(maskobj.volInfo.image_indx), 'single');
%                 maskobj.removed_voxels = [];
%                 maskobj.volInfo.wh_inmask = (1:maskobj.volInfo.nvox)';
%             end
%             

            end % data type/nargin switch
            
        end % function
        
        
    end % methods
    
end % classdef
