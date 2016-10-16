% fmridisplay: Data class for storing data matrices and information
%
% 'fmridisplay' is a data class containing information about an underlay
% and activation map(s) for creating montage plots and other types of
% plots.
%
% The default brain for overlays is based on Keuken et al. 2014
% For legacy SPM8 single subject, enter as arguments:
% 'overlay', which('SPM8_colin27T1_seg.img')
%
% Creating class instances
% -----------------------------------------------------------------------
% 
%
% Examples
% -----------------------------------------------------------------------
% obj = fmridisplay;            Create object with canonical underlay
% obj = montage(obj);           Show axial montage of underlay
% obj = addblobs(obj, cl);      Add blobs from cl clusters
% 
% obj = fmridisplay('montage', 'addblobs', cl);  Do all of the above
%
% Add colored green blobs, with smoothed edges
% obj = addblobs(obj, cl, 'trans', 'color', [0 1 0], 'smooth');
%
% Create a second montage with 4 mm spacing and add blue blobs to both:
% obj = montage(obj, 'spacing', 4);
% obj = addblobs(obj, cl, 'trans', 'color', [0 0 1], 'smooth');
%
% Add axial and parasaggital montages to the same figure:
% o2 = montage(o2, 'axial', 'slice_range', [-40 50], 'onerow', 'spacing', 6);
% axh = axes('Position', [0.03 0.45 .1 .5]);
% o2 = montage(o2, 'saggital', 'wh_slice', [0 0 0], 'existing_axes', axh);
%
% Add blue outlines only:
% obj = addblobs(obj, cl, 'outline', 'color', [0 0 1]);
%
% Sagittal images and blobs:
% o2 = fmridisplay;
% o2 = montage(o2, 'saggital', 'slice_range', [-20 20], 'onerow');
% o2 = addblobs(o2, cl);
% legend(o2, 'figure')
% o2 = addblobs(o2, cl, 'contour', 'color', [0 1 0]);
% 
% Display over your custom anatomical image:
% overlay = 'mean_T1_FSL1.nii';
% o2 = fmridisplay('overlay', overlay);
%
% Overlapping sagittal and axial images with outlines
% o2 = fmridisplay;
% o2 = montage(o2, 'saggital', 'slice_range', [-10 10], 'onerow');
% enlarge_axes(gcf, 1.2);
% o2 = montage(o2, 'axial', 'slice_range', [-40 50], 'onerow', 'spacing', 4);
% o2 = addblobs(o2, cl, 'splitcolor', {[0 0 1] [.3 0 .8] [.8 .3 0] [1 1 0]});
% o2 = addblobs(o2, cl, 'maxcolor', [1 0 0], 'mincolor', [1 .3 0], 'cmaprange', [0 .01], 'outline');
%
% Add transparent blobs, either mapped with color scale or constant
% opacity value:
% o2 = addblobs(o2, cl, 'splitcolor', {[0 0 1] [0 1 1] [1 .5 0] [1 1 0]}, 'cmaprange', [-2 2], 'trans');
% o2 = addblobs(o2, cl, 'splitcolor', {[0 0 1] [0 1 1] [1 .5 0] [1 1 0]}, 'cmaprange', [-2 2], 'transvalue', .85);
%
% o2 = montage(o2, 'saggital', 'slice_range', [-6 6], 'onerow');
% o2 = montage(o2, 'axial', 'slice_range', [-20 30], 'onerow', 'spacing', 8);
% xyz = a list of [x y z] coordinates, one coord per row
% o2 = montage(o2, 'axial', 'wh_slice', xyz, 'onerow');
% squeeze_axes(o2.montage{1}.axis_handles, 20)
% squeeze_axes(o2.montage{2}.axis_handles, 10)
%
% Plot points (i.e., coordinate locations) rather than blobs:
% o2 = addpoints(o2, DB.xyz, 'MarkerFaceColor', 'b', 'Marker', 'o', 'MarkerSize', 4);
% o2 = addpoints(o2, DB.xyz, 'text', DB.textcodes, 'condf', DB.condf, 'color', {'b' 'g'});
% o2 = removepoints(o2);
%
% Add to only 2nd montage in registered vector of montages in obj.montage 
%   obj = addblobs(obj, cl, 'which_montages', 2); 
%
% Make a montage with centers of each significant cluster and add blobs to that:
% r = region(b2);
% xyz = cat(1, r.mm_center)
% o2 = montage(o2, 'axial', 'wh_slice', xyz, 'onerow');
% o2 = addblobs(o2, region(b2));
%
% Methods
% -----------------------------------------------------------------------
% See methods(fmridisplay)
%
% See fmridisplay.montage for list of all slice display options
% See fmridisplay.render_blobs for list of all rendering options
%
% Copyright 2011 - Tor Wager

classdef fmridisplay 
    
    properties
        
        overlay
        SPACE 
        activation_maps

        montage
        surface
        orthviews
        
        history
        history_descrip = 'Cell array: names of methods applied to this data, in order';      
        additional_info = struct('');
        
    end % properties
    
    methods
        
        % Class constructor
        function obj = fmridisplay(varargin)
            
            
            % ---------------------------------
            % Create empty object, and return if no additional
            % arguments
            % ---------------------------------
            
            %obj.overlay = which('SPM8_colin27T1_seg.img');  % spm8 seg cleaned up
            obj.overlay = which('keuken_2014_enhanced_for_underlay.img');

%             obj.overlay = which('clean_Q1-Q6_RelatedParcellation210_AverageT1w_restore.nii');  % spm8 seg cleaned up
  
            if any(strcmp(varargin, 'overlay'))
                wh = find(strcmp(varargin, 'overlay'));
                wh = wh(1);
                obj.overlay = varargin{wh + 1};
                varargin([wh:wh+1]) = [];
                
                % kludge for weird full path issue 
                ov2 = which(obj.overlay);
                if ~isempty(ov2), obj.overlay = ov2; end
            end
            
            obj.SPACE = [];
            
            if ~isempty(obj.overlay)
                V = spm_vol(obj.overlay);
                V = V(1);
            else
                error('Cannot find overlay image.');
            end
            
            obj.SPACE = define_sampling_space(V);
            obj.activation_maps = {};

            obj.montage = {};
            obj.surface = {};
            obj.orthviews = {};
            
%             obj.axis_handles = struct('montage', [], 'surface', [], 'orthviews', []);
%             obj.axis_handles.montage = struct('axial', [], 'sagittal', [], 'coronal', []);
%             
            obj.history = {};
            obj.history_descrip = [];
            obj.additional_info = '';
            
            if nargin == 0
                return
            end
            

            % Now parse inputs and run methods depending on what is entered
            
            if any(strcmp(varargin, 'montage'))
                wh = strcmp(varargin, 'montage');
                varargin(wh) = [];
                obj = montage(obj, varargin{:});
            end
            
            for i = 1:length(varargin)
                if ischar(varargin{i})
                    switch varargin{i}
                                                    
                        case 'addblobs'
                            
                            cl = varargin{i + 1};
                            if ~isstruct(cl) && ~isa(cl, 'region')
                                disp('addblobs arg should be followed by clusters/region structure.');
                            end
                            
                            obj = addblobs(obj, cl, varargin);
                            
                        case {'slice_range', 'smooth', 'spacing', 'onerow'}
                            % other inputs that we should ignore here
                            % but are used in subfunctions like montage,
                            % etc.
                            
                        otherwise, warning(['Unknown input string option:' varargin{i}]);
                    end
                end
            end

            
        end % constructor function
        
    end % methods
    
    
end
