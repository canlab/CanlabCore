classdef fmri_timeseries < fmri_data
    
    properties
        subject = 'Subject_ID_here';
        study = 'Study_ID_here';
        
        X = [];
        X_descrip = 'time (observations) x voxels data';
        
        volInfo = struct('XYZ', [0 0 0]', 'XYZmm', [0 0 0]', 'M', zeros(4));
        volInfo_descrip = 'Volume info from iimg_read_img';
        
        image_names = char('');
        
        onsets = cell(1);
        onsets_descrip = 'Cell array with onset times for events, one cell per event type per session';
        
        condition_names = cell(1);
        condition_names_descrip = 'Cell array with names of conditions';
        
        Y = [];
        Y_descrip = 'Behavioral or outcome data matrix.';
        
        covariates = 0;
        covariates_descrip = 'Nuisance covariates associated with data';
        
        history = {'raw'};
        history_descrip = 'Cell array of names of methods applied to this data, in order';
        
    end % properties
    
    methods
        
        % -------------------------------------
        % Checks for legal data
        % -------------------------------------
% %         
% %         function obj = create(obj, varargin)
% %             
% %             N = fieldnames(obj);
% %             
% %             for i = 1:length(varargin)
% %                 if ischar(varargin{i})
% %                     
% %                     % Look for a field (attribute) with the input name
% %                     wh = strmatch(varargin{i}, N, 'exact');
% %                     
% %                     if ~isempty(wh)
% %                         
% %                         obj.(varargin{i}) = varargin{i + 1};
% %                         
% %                     end
% %                     
% %                 end
% %             end
% %    
% %         end % function
        
    end % methods
    
    
end