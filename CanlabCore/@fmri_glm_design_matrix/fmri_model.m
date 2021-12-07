
% fmri_model
%
% Specify and/or create an fmri_model object. This object is designed to contain
% all the information about a model necessary to create a design matrix
% and, when combined with an fmri_data object, to fit the model.
%
% The information contained in this object (its attributes) are designed to
% be very similar or identical to those defined by SPM in using
% spm_fmri_design.m. This should help with inter-operability and standards
% for storing information.
%
% See methods(fmri_model) for things you can do with this.
% See the code of this function for additional details about all the
% fields.
%
% The general way to create/initialize CANlab objects is to call them with
% no or minimal input arguments. In this case, TR is required as the first
% input.  After that, pairs of arguments can be entered in the
%'fieldname', value format standard for Matlab inputs. If the fields are
% valid attributes of the object, they will be entered.
%
% For example:
%
% my_model = fmri_model(2);
% creates an empty fmri_model object with a TR of 2 (which is used to
% create a default b-spline basis set for the HRF.
%
% my_model = fmri_model(2, 'nscan', [198 198 198]);
% creates an empty structure but assigns data to the field nscans, number
% of scans per session.
%
% my_model = fmri_model(TR, 'nscan', nscan, 'units', 'secs', 'onsets', ons_temperature, 'condition_names', names_temperature);
% does the above and also adds onsets and names to the appropriate place
% in the object structure.
% my_model.Sess(1).U(2).name : names for session 1, condition 2 
% my_model.Sess(1).U(2).ons : onsets for session 1, condition 2 in UNITS(secs or TRs)
%
% add methods:
% There are some special methods to add specific types of data:
%
% 'onsets', {cell array of onsets for conditions within sessions}
%
% 'condition_names', {cell array of condition names} -- assumed to be the same for all sessions
%
% 'pm', '*name of condition to modulate*', '*modulator name*', {cell with modulators for each session}
%
% You can also specify different basis sets for different conditions, which
% SPM will not allow you to do.  Here is an example of a process of
% creating an fmri_model object from onsets, etc., building it, and then
% replacing the basis set for one event type with another one.
%
% reportmod_model = fmri_model(TR, 'nscan', nscan, 'units', 'secs', 'onsets', ons_reportmod, ...
% 'condition_names', names_reportmod, 'pm', PM, 'pmnames', PM_names);
%
% reportmod_model = build(reportmod_model);
% plot(reportmod_model)
%
% % Generate a new basis set for the Anticipation conditions (condition 1)
% [xBF_hires, xBF] = fmri_spline_basis(2, 'length', 12, 'nbasis', 3, 'order', 3, 'plot');
% reportmod_model = replace_basis_set(reportmod_model, 1, xBF_hires);
% reportmod_model = build(reportmod_model);
% plot(reportmod_model)
%

%    SPM.Sess(s)
%            U: - Input structure array
%            C: - User specified covariate structure
%          row: - scan   indices for session s
%          col: - effect indices for session s
%           Fc: - F Contrast information for input-specific effects
%
%    SPM.xX
%            X: - design matrix
%           iH: - vector of H partition (indicator variables) indices
%           iC: - vector of C partition (covariates)          indices
%           iB: - vector of B partition (block effects)       indices
%           iG: - vector of G partition (nuisance variables)  indices
%         name: - cellstr of names for design matrix columns
%

%        3rd level
%        ------------------------------------------------------------------
%        SPM.Sess(s).U
%               dt: - time bin length {seconds}
%             name: - {1 x j} cell of names for each input or cause
%              ons: - (q x 1) onsets for q  trials {in UNITS}
%              dur: - (q x 1) durations for trials {in UNITS}
%                P: - Parameter stucture
%                u: - (t x j) inputs or stimulus function matrix
%              pst: - (1 x k) peristimulus times (seconds)
%
%
%        SPM.Sess(s).C
%
%                C: - [kx1 double] of user specified regressors
%             name: - {1xk} cellstr of regressor names
%
%
%        SPM.Sess(s).Fc
%
%                i: - F Contrast colums for input-specific effects
%             name: - F Contrast names  for input-specific effects
%
%
%            4th level
%            --------------------------------------------------------------
%            SPM.Sess(s).U(i).P(p)
%
%
%                 name: - parameter name
%                    P: - (q x 1) parameter matrix
%                    h: - order of polynomial expansion (0 = none)
%                    i: - sub-indices of U(i).u for plotting

classdef fmri_model
    
    properties
        
        % custom: not in SPM
        build_method
        history
        
        % 1st level
        % --------------------------------------------------------------------------
        % SPM.
        xY % : [1x1 struct] - data structure
        nscan % : [1xs double] - nscan(s) = number of scans in session s
        xBF % : [1x1 struct] - Basis function structure
        Sess % : [1xs struct] - Session structure array
        xX % : [1x1 struct] - Design matrix structure
        
    end % properties
    
    
    
    methods
        
        % Class constructor
        function obj = fmri_model(TR, varargin)
            %
            % [obj, cl_with_averages] = fmri_data(image_names, mask_image, varargin)
            %
            % Reads a set of image files and a mask image, and returns
            % an fmri_data object with data for all in-mask voxels.
            
            % ---------------------------------
            % Create empty fmri_data object, and return if no additional
            % arguments
            % ---------------------------------
            
            if nargin == 0
                error('Must define TR, repetition time for scans, as first input.')
            end
            
            obj.build_method = 'Separate sessions';
            
            obj.xY  = [];   % : [1x1 struct] - data structure
            obj.nscan = []; % : [1xs double] - nscan(s) = number of scans in session s
            obj.xBF = [];   % : [1x1 struct] - Basis function structure
            obj.Sess = [];  % : [1xs struct] - Session structure array
            obj.xX = [];
            
            % Default basis set
            %    ----------------------------------------------------------------------
            obj.xY.RT = TR; % : - repetition time {seconds)
            
            obj.xBF = fmri_spline_basis(TR, 0);
            % obj.xBF.name: - name of basis set
            % obj.xBF.length: - support of basis set {seconds}
            % obj.xBF.order: - order of basis set
            % obj.xBF.bf: - basis set matrix
            % obj.xBF.dt = TR ./ 16; % : - length of time bin {seconds}
            
            obj.xBF.T = 16; % : - number of time bins per scan
            obj.xBF.T0 = 1; % : - first time bin (see slice timing)
            obj.xBF.UNITS = [];  %: - 'scans'|'secs' - units in which onsets are specified
            obj.xBF.Volterra = 1; % : - 1|2 - order of [Volterra] convolution
            
            % Parameters for parametric modulation
            P = struct('name', '', 'P', [], 'h', 0, 'dur', [], 'i', []);
            
            % Onset structure
            U = struct('dt', TR ./ 16, 'name', {}, 'ons', [], 'dur', 1, 'P', P, 'u', [], 'pst', []);
            
            % Covariate structure
            C = struct('C', [], 'name', {'User-specified regressors here'});
            
            % Session structure
            S = struct('U', U, 'C', C, 'row', [], 'col', [], 'Fc', []);
            
            obj.Sess = S;
            
            % Design structure
            obj.xX = struct('X', [], 'iH', [], 'iC', [], 'iB', [], 'iG', [], 'name', {});
            
            % The code below can be generic to any class definition
            % It parses 'fieldname', value pairs of inputs
            % and returns a warning if unexpected strings are found.
            
            if nargin == 1
                return
            end
            
            % all valid fieldnames
            valid_names = fieldnames(obj);
            
            for i = 1:length(varargin)
                if ischar(varargin{i})
                    
                    % Look for a field (attribute) with the input name
                    wh = strmatch(varargin{i}, valid_names, 'exact');
                    
                    % behaviors for valid fields
                    if ~isempty(wh)
                        
                        obj.(varargin{i}) = varargin{i + 1};
                        
                        % eliminate strings to prevent warnings on char
                        % inputs
                        if ischar(varargin{i + 1})
                            varargin{i + 1} = [];
                        end
                        
                        % special methods for specific fields
                        switch varargin{i}
                            
                        end
                        
                    else
                        % Unique to this method: special subfields/other valid
                        % entries
                        % Try to add using the add method, which returns a
                        % warning if the field name is invalid.
                        
                        obj = add(obj, varargin{i}, varargin{i + 1});
                        
                        varargin{i + 1} = []; % eliminate to avoid confusing the parsing
                        
                    end % not empty
                end % string input
            end % process inputs
            
            
        end % class constructor function
    end % properties
    
    
end % classdef