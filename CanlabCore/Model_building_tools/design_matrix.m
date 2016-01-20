    % design_matrix: data class for creating a design matrix to be used with a linear model including fmri data.
    %
    %
    % Inputs:
    % ---------------------------------------------------------------------
    % dat                       : M x N numeric matrix containing Observations and Variables
    %                           If dat is a file name then will try to
    %                           import data into design_matrix object.
    %                           Make sure there is only text in header and
    %                           not anywhere else.
    %
    % varname                   : Cell array containing variable names.  Must
    %                             match number of column in data matrix
    %
    % Examples:
    % ---------------------------------------------------------------------
    % DM = design_matrix([ones(10,1), (1:10)', (1:10).^2'],{'Intercept','X','X2'})
    %
    % Original version: Copyright Luke Chang 2/2014
    
    % Notes:
    % Need to add these Methods:
    % -create regressor from stim times
    % -pca

classdef design_matrix < handle
    properties
        dat = [];
        varname = {};
        fname = '';
    end
    
    methods
        function obj = design_matrix(dat, varname)
            class constructor % Initialize instance of design_matrix
            
            if(nargin > 1)
                try
                    if(~ismatrix(dat) || ~isnumeric(dat) || iscell(dat))
                        error('Make sure input data is a matrix')
                    end
                    if(length(varname) ~= size(dat,2) || ~iscell(varname));
                        error('Make sure the number of variable names corresponds to number of data columns.')
                    end
                    obj.dat = dat;
                    obj.varname = varname;
                catch err
                    error('Make sure input variable names are in a cell array with length equal to number of data columns and data is a matrix.')
                end
                obj.varname = varname;
                
            elseif(nargin > 0)
                
                %import data into design_matrix if dat is a filename
                if ischar(dat)
                    ftest = exist(dat,'file'); %check if valid file
                    if ftest == 2
                        data = importdata(dat);
                        obj.dat = data.data;
                        obj.varname = data.colheaders;
                    end
                    
                elseif(~ismatrix(dat) || ~isnumeric(dat) || iscell(dat))
                    error('Make sure input data is a matrix')
                else
                    
                    obj.dat = dat;
                end
                
            else % if nothing initialize empty object
                return
            end
        end
        
        function names(obj)
            % names(obj)
            %
            % List variable names for each regressor
            display(obj.varname)
        end
        
        function dim = size(obj, varargin)
            % dim = size(obj, varargin)
            %
            % Return dimensions of design matrix
            % Optional Input: Indicate Dimension(row = 1 or column = 2)
            
            if nargin > 1
                dim = size(obj.dat, varargin{1});
            else
                dim = size(obj.dat);
            end
        end
        
        function plot(obj)
            % plot(obj)
            %
            % Plot design matrix
            imagesc(obj.dat)
        end
        
        function save(obj, fname)
            % save(obj, fname)
            %
            % Save Design Matrix to file
            
            save(fname, obj)
        end
        
        function obj = addvariable(obj, x, varargin)
            % obj = addvariable(obj, x, varargin)
            %
            % Add regressor to design matrix
            %
            % optional inputs
            % -------------------------------------------------------------------
            % 'Name'        : 'Name' followed by Variable name
            %                  Default is 'newVx'
            %
            % 'Order'       : 'Order' followed by where to insert new data columns location
            %                  Default is end
            
            % Check Inputs
            if size(x,1) ~= size(obj,1)
                error('Make sure new variable column is the same length as design matrix')
            end
            
            % Defaults
            for i = 1:size(x,2)
                newname{i} = ['newV' num2str(i)];
            end
            isnewname = 0;
            varorder = size(obj,2); %add to end
            
            % Parse inputs
            % -------------------------------------------------------------------
            for varg = 1:length(varargin)
                if ischar(varargin{varg})
                    % reserved keywords
                    if strcmpi('name',varargin{varg})
                        if length(varargin{varg + 1}) == size(x,2)
                            newname = varargin{varg + 1};
                            isnewname = 1;
                        else
                            error('Make sure ''Name'' is follwed by a valid variable name')
                        end
                        varargin{varg} = {}; varargin{varg + 1} = {};
                    end
                    
                    if strcmpi('order',varargin{varg})
                        if isnumeric(varargin{varg + 1}) && varargin{varg + 1} <= size(obj,2)
                            varorder = varargin{varg + 1};
                        else
                            error('Make sure ''Order'' is follwed by a valid column number')
                        end
                        varargin{varg} = {}; varargin{varg + 1} = {};
                    end
                end
            end
            % -------------------------------------------------------------------
            
            % Update Name vector if not empty
            if ~isempty(obj.varname)
                if varorder == 1 %begin
                    obj.varname = [newname, obj.varname];
                elseif varorder == size(obj,2) %end
                    obj.varname = [obj.varname, newname];
                else %Somewhere inbetween
                    obj.varname = [obj.varname(1:varorder), newname, obj.varname(varorder + 1 : end)];
                end
            end
            
            % Add new variables to design Matrix
            if varorder == 1 %begin
                obj.dat = [x, obj.dat];
            elseif varorder == size(obj,2) %end
                obj.dat = [obj.dat, x];
            else %Somewhere inbetween
                obj.dat = [obj.dat(:,1:varorder), x, obj.dat(:,varorder + 1 : end)];
            end
        end
        
        function obj = zscore(obj, varargin)
            % obj = zscore(obj, varargin)
            %
            % Standardize columns of design matrix
            %
            % optional inputs
            % -------------------------------------------------------------------
            % 'center'        : Only remove mean, don't standardize
            
            % Defaults
            center = 0; %Only remove mean
            
            % optional inputs
            % -------------------------------------------------------------------
            for varg = 1:length(varargin)
                if ischar(varargin{varg})
                    % reserved keywords
                    if strcmpi('center',varargin{varg})
                        center=1;
                        varargin{varg} = {};
                    end
                end
            end
            % -------------------------------------------------------------------
            
            if center
                obj.dat = obj.dat - repmat(mean(obj.dat),size(obj,1),1);
            else
                obj.dat = zscore(obj.dat);
            end
        end
        
        function obj = removevariable(obj, x)
            % obj = removevariable(obj, x)
            %
            % Remove columns from design_matrix
            %
            % Inputs
            % -------------------------------------------------------------------
            % x             : Input vector of columns to remove
            
            obj.dat(:,x) = [];
            obj.varname(x) = [];
        end
        
        function obj = addintercept(obj)
            % obj = addintercept(obj)
            %
            % Add intercept to design matrix
            
            %Check if Intercept exists
            if sum(strcmpi(obj.varname,'intercept')) > 0
                error('Intercept Name already included in obj.varname')
            end
            
            %Check if any variable only includes ones
            for i = 1:size(obj,2)
                if sum(obj.dat(:,i)==1) == size(obj,1)
                    error('There is already a column of ones that resembles an intercept')
                end
            end
            
            obj.dat = [obj.dat, ones(size(obj,1),1)];
            obj.varname = [obj.varname, 'Intercept'];
        end
        
        function obj = removeintercept(obj)
            %  obj = removeintercept(obj)
            %
            % Remove intercept from design matrix
            
            %Check if any variable only includes ones or if intercept is in varname
            for i = 1:size(obj,2)
                whereint(i) = sum(obj.dat(:,i)==1) == size(obj,1);
            end
            if sum(whereint) == 0
                error('There does not appear to be any column of ones that resembles an intercept')
            elseif sum(strcmpi('intercept',obj.varname)) == 0
                error('Intercept Name is not included in obj.varname')
            end
            
            % now remove intercept
            obj.dat(:,whereint) = [];
            obj.varname(whereint) = [];
        end
        
        function vif = vif(obj, varargin)
            % vif = vif(obj, varargin)
            %
            % Check for multicollinearity by getting variance inflation factors
            %
            % See original getvif.m in canlab repository - OptimizeDesign11/core_functions/getvif.m
            %
            % optional inputs
            % -------------------------------------------------------------------
            % 'nointercept'       : Remove intercept (turned off by default)
            
            % Defaults
            noint = 0;
            if strcmpi('nointercept',varargin)
                noint = 1;
            end
            
            %Remove intercept if asked
            if noint
                obj = removeintercept(obj);
            end
            
            %Calculate VIF
            for i = 1:size(obj.dat,2)
                X = obj.dat;
                y = X(:,i);
                X(:,i) = [];
                b = X\y;
                fits = X * b;
                rsquare = var(fits) / var(y);
                
                if rsquare == 1,rsquare = .9999999;end
                
                vif(i) = 1 / (1 - rsquare);
            end
        end
        
        function r = corr(obj)
            % Calculate pairwise correlation of regressors in design_matrix
            
            r = corr(obj.dat);
        end
        
        function obj = normalizedrank(obj, varargin)
            % obj = normalizedrank(obj, varargin)
            %
            % Rank each regressor and normalize between [0,1]
            %
            % See normalizedrank.m for optional inputs
            
            obj.dat = normalizedrank(obj.dat, varargin);
        end
        
        function obj = conv_hrf(obj, varargin)
            % obj = conv_hrf(obj, varargin)
            %
            % Convolve each regressors with hemodynamic response function
            % Uses spm_hrf.m
            %
            % optional inputs
            % -------------------------------------------------------------------
            % 'tr'           : Input TR to use for creating HRF
            %                  (e.g., 'tr', 3)
            %
            % 'select'       : Select Input vector of regressors to convolve
            %                  (e.g., 'select', [2,4])
            %
            % 'custom_hrf'   : Use custom HRF (e.g., 'custom_hrf', [1.00, 0.59, 0.39, 0.27]
            
            % Defaults
            include = (1:size(obj,2)); %Convolve entire Design Matrix by default
            tr = 2;
            
            % Check if spm_hrf is on path
            checkspm =  which('spm_hrf.m');
            if isempty(checkspm), error('Make sure spm is in matlab path'); end
            
            % Parse inputs
            % -------------------------------------------------------------------
            for varg = 1:length(varargin)
                if ischar(varargin{varg})
                    % reserved keywords
                    if strcmpi('tr',varargin{varg})
                        if isnumeric(varargin{varg + 1})
                            tr = varargin{varg + 1};
                            crf = spm_hrf(tr);
                        else
                            error('Make sure ''tr'' is followed by valid number')
                        end
                        varargin{varg} = {}; varargin{varg + 1} = {};
                    end
                    
                    if strcmpi('select',varargin{varg})
                        if isnumeric(varargin{varg + 1}) && all(varargin{varg + 1} <= size(obj,2))
                            include = varargin{varg + 1};
                        else
                            error('Make sure ''select'' is followed by a valid column number')
                        end
                        varargin{varg} = {}; varargin{varg + 1} = {};
                    end
                    
                    if strcmpi('custom_hrf',varargin{varg})
                        if isnumeric(varargin{varg + 1}) && all(varargin{varg + 1} <= size(obj,2))
                            crf = varargin{varg + 1};
                        else
                            error('Make sure ''custom_hrf'' is followed by a vector of a canonical response function')
                        end
                        varargin{varg} = {}; varargin{varg + 1} = {};
                    end
                end
            end
            % -------------------------------------------------------------------
            
            %Convolution of task
            for i = 1:length(include)
                convdat = conv(obj.dat(:,include(i)),crf);
                
                %Cut off extra data from convolution
                obj.dat(:,include(i)) = convdat(1:size(obj,1));
            end
        end
        
        function obj = hpfilter(obj, varargin)
            % obj = hpfilter(obj, varargin)
            %
            % Add High pass filter design matrix using spm's discrete
            % cosine Transform (spm_filter.m)
            %
            % optional inputs
            % -------------------------------------------------------------------
            % 'tr'           : Input TR to use for creating HRF
            %                  (e.g., 'tr', 3) Default = 2;
            %
            % 'duration'     : Duration of high pass filter in seconds
            %                  (e.g., 'duration', 100) Default = 180;
            
            % Defaults
            tr = 2;
            filterlength = 180;
            
            % Check if spm_hrf is on path
            checkspm =  which('spm_filter.m');
            if isempty(checkspm), error('Make sure spm is in matlab path'); end
            
            % Parse inputs
            % -------------------------------------------------------------------
            for varg = 1:length(varargin)
                if ischar(varargin{varg})
                    % reserved keywords
                    if strcmpi('tr',varargin{varg})
                        if isnumeric(varargin{varg + 1})
                            tr = varargin{varg + 1};
                        else
                            error('Make sure ''tr'' is followed by valid number')
                        end
                        varargin{varg} = {}; varargin{varg + 1} = {};
                    end
                    if strcmpi('duration',varargin{varg})
                        if isnumeric(varargin{varg + 1})
                            filterlength = varargin{varg + 1};
                        else
                            error('Make sure ''duration'' is followed by valid number')
                        end
                        varargin{varg} = {}; varargin{varg + 1} = {};
                    end
                end
            end
            % -------------------------------------------------------------------
            
            %create high pass filter
            K.RT = tr;
            K.row = 1:size(obj,1);
            K.HParam = filterlength;
            nK = spm_filter(K);
            if isempty(nK.X0), error('Check if filter duration is too long'); end
            
            %Add filter to design_matrix
            obj.dat = [obj.dat, nK.X0];
            
            %Add variable names
            for i = 1:size(nK.X0,2)
                filtname{i} = ['hpfilter' num2str(i)];
            end
            obj.varname = [obj.varname, filtname];
        end
        
        function stats = regress(obj, Y, varargin)
            % [B,BINT,R,RINT,STATS] = regress(obj, Y)
            %
            % Regress design matrix on vector Y
            % Uses matlab's regress function
            %
            % optional inputs
            % -------------------------------------------------------------------
            % 'robust'           : use robust regression
            
            % Defaults
            doRobust = 0;
            for i = 1:length(varargin)
                if strcmpi(varargin(i),'robust')
                    doRobust = 1;
                    find_robust = exist('robustfit');
                    if find_robust ~= 2
                        error('Make sure robustfit is on your path, requires the stats toolbox')
                    end
                    varargin(i) = [];
                end
            end
            
            if ~doRobust
                [stats.B, stats.BINT, stats.R, stats.RINT, stats.STATS] = regress(Y, obj.dat);
            else
                [stats.B, stats.STATS] = robustfit(obj.dat, Y,[],[],'off');
            end
        end
        
        function obj = onsettimes(obj, onset, names, tr, timing )
            % obj = onsettimes(obj, onset, names, tr, timing )
            %
            %Create stimulus regressor from onset times
            %
            % Inputs
            % -------------------------------------------------------------------
            % onset        : Input cell array of onset times for each
            %                   regressor in FSL's 3 column format (e.g., onset in sec, duration in sec, weight).
            %
            % names        : Cell array of variable names corresponding
            %                to each onset cell (e.g., {'BlueOn','RedOn'})
            %
            % tr           : Repetition time (e.g., 2)
            %
            % timing       : Timing converstion from onset array to design matrix
            %                  (e.g., 'sec2tr','tr2sec','sec2sec',or
            %                  'tr2tr'). Need to know which format each
            %                  array is in.
            
            %Convert Onset Times Into Boxcar Regressors
            r = zeros(size(obj,1),length(onset));
            for i = 1:length(onset)
                for j = 1:size(onset{i},1)
                    switch timing
                        case 'sec2tr'
                            if floor(onset{i}(j,1)/tr) == 0
                                r(1 : 1 + ceil(onset{i}(j,2)), i) = onset{i}(j,3);
                            else
                                r(floor(onset{i}(j,1) / tr) : floor(onset{i}(j,1) / tr) + ceil(onset{i}(j,2) / tr) - 1, i) = onset{i}(j,3);
                            end
                        case 'tr2sec'
                            r(floor(onset{i}(j,1) * tr) : floor(onset{i}(j,1) * tr) + ceil(onset{i}(j,2) * tr) - 1, i) = onset{i}(j,3);
                        case {'sec2sec', 'tr2tr'}
                            r(floor(onset{i}(j,1)) : floor(onset{i}(j,1)) + ceil(onset{i}(j,2)) - 1, i) = onset{i}(j,3);
                    end
                end
            end
            obj.dat = [obj.dat, r];
            
            %Add Variable names
            obj.varname = [obj.varname, names];
        end
        
        function obj = write(obj, varargin)
            % obj = write(obj, varargin)
            %
            % -------------------------------------------------------------------
            % write design_matrix object into csv file.  Will use obj.fname
            % or can specify optional name
            %
            % -------------------------------------------------------------------
            % Optional Inputs
            % -------------------------------------------------------------------
            % fname        : path and file name of csv file.
            % -------------------------------------------------------------------
            
            
            if nargin > 1 %use supplied file name
                if ischar(varargin{1})
                    hdr = sprintf('%s,',obj.varname{:});
                    hdr(end) = '';
                    dlmwrite(varargin{1}, hdr,'') %Write Header 1st
                    dlmwrite(varargin{1}, obj.dat, 'delimiter',',','-append','precision',10) %Append data
                end
                
            elseif ~isempty(obj.fname) %use obj.fname
                hdr = sprintf('%s,',obj.varname{:});
                hdr(end) = '';
                dlmwrite(obj.fname, hdr, ''); %Write Header 1st
                dlmwrite(obj.fname, obj.dat, 'delimiter',',','-append','precision',10) %Append data
            else
                error('Please supply valid file name with path to save.')
            end
            
        end
        
        function c = horzcat(varargin)
            % function c = horzcat(varargin)
            % -------------------------------------------------------------------
            % Implements the horzcat ([a b]) operator on design_matrix objects across variables.
            % Requires that each object has an equal number of rows
            % -------------------------------------------------------------------
            % Examples:
            % c = [dm1 dm2];
            % -------------------------------------------------------------------
            
            %check if number of rows is the same
                        %check if varnames are the same
            nrow = [];
            for i = 1:nargin
                nrow(i) = size(varargin{i},1);
            end
            for i = 1:nargin
                for j = 1:nargin
                    if nrow(i)~=nrow(j)
                        error('objects have a different number of rows')
                    end
                end
            end
            
            dat = [];
            varname = [];
            for i = 1:nargin
                    %Check if design_matrix object
                if ~isa(varargin{i}, 'design_matrix')
                    error('Input Data is not an design_matrix object')
                end
                
                dat = [dat, varargin{i}.dat];
                varname = [varname, varargin{i}.varname];
            end
            
            c = varargin{1};
            c.dat = dat;
            c.varname = varname;
        end
        
        function c = vertcat(varargin)
            % function c = vertcat(varargin)
            % -------------------------------------------------------------------
            % Implements the vertcat ([a b]) operator on design_matrix objects across rows.
            % Requires that each object has an equal number of columns and
            % that each varname is identical
            % -------------------------------------------------------------------
            % Examples:
            % c = [dm1; dm2];
            % -------------------------------------------------------------------
            
            %check if varnames are the same
            varname = {};
            for i = 1:nargin
                varname{i} = varargin{i}.varname;
            end
            for i = 1:nargin
                for j = 1:nargin
                    if ~strcmpi(varname{i},varname{j})
                        error('variable names do not match')
                    end
                end
            end
            dat = [];
            for i = 1:nargin
                %Check if design_matrix object
                if ~isa(varargin{i}, 'design_matrix')
                    error('Input Data is not an design_matrix object')
                end
                dat = [dat; varargin{i}.dat];
            end
            
            c = varargin{1};
            c.dat = dat;
        end
        
    end %methods
end %class
