    % comp_model: data class for creating a computational model
    %
    % Computational Model Class: comp_model
    %
    %--------------------------------------------------------------------------
    % This object is used to fit a computational model to a multi-subject
    % dataset.  The object uses the design_matrix() class to for the data set
    % and has additional fields for the model and parameters for the model
    % fitting procedure such as the parameter constraints, number of
    % iterations, and type of estimation (e.g., maximum likelihood or least
    % squares).
    %
    %--------------------------------------------------------------------------
    % Inputs:
    % ---------------------------------------------------------------------
    % dat                       : M x N numeric matrix containing Observations and Variables
    %
    % varname                   : Cell array containing variable names.  Must
    %                             match number of column in data matrix
    % model                     : name of model file (must be on matlab path)
    %
    %--------------------------------------------------------------------------
    % Current Methods for comp_model (inherits from design_matrix class too)
    %--------------------------------------------------------------------------
    %
    % avg_aic                   : display average AIC value
    % avg_bic                   : display average BIC value
    % avg_params                : display average parameter estimates
    % comp_model                : class constructor
    % fit_model                 : estimate parameters using model
    % get_aic                   : extract all subject's AIC values
    % get_bic                   : extract all subject's BIC values
    % get_params                : extract all subject's estimated
    %                             parameters
    % plot                      : plot average model predictions across subjects
    % summary                   : display summary table for model
    % save                      : save object as .mat file
    % write_tables              : write out parameter estimates and trial-to-trial
    %                             predictions to csv data frame.
    %
    %--------------------------------------------------------------------------
    % Examples:
    % ---------------------------------------------------------------------
    % m1 = comp_model([ones(10,1), (1:10)', (1:10).^2'],{'Intercept','X','X2'},'Linear_Model')
    %
    % Also see CompModel_Tutorial.m in Examples
    % -------------------------------------------------------------------------
    % Author and copyright information:
    % -------------------------------------------------------------------------
    %     Copyright (C) 2014  Luke Chang
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
    % -------------------------------------------------------------------------

classdef comp_model < handle & design_matrix
    properties
        % inherits properties from design_matrix
        model = [];
        param_min = [];
        param_max = [];
        nStart = [];
        esttype = [];
        params = [];
        trial = [];
    end
    
    methods
        function obj = comp_model(dat, varname, model, varargin)
            class constructor % Initialize instance of comp_model
            
            if(nargin > 2)
                
                %Add Data matrix
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
                
                %Check if model is a matlab function on path
                testpath = exist(model);
                if testpath ~= 2
                    error(['Make sure model : ' model ' is on your matlab path.'])
                else
                    obj.model = model;
                end
                
                %Add variable names
                obj.varname = varname;
                
            elseif(nargin > 1)
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
                if(~ismatrix(dat) || ~isnumeric(dat) || iscell(dat))
                    error('Make sure input data is a matrix')
                end
                obj.dat = dat;
            else % if nothing initialize empty object
                return
            end
            
            % Parse input
            for varg = 1:length(varargin)
                if ischar(varargin{varg})
                    if strcmpi('param_min',varargin{varg})
                        obj.param_min = varargin{varg + 1};
                        varargin{varg} = []; varargin{varg + 1} = [];
                    elseif strcmpi('param_max',varargin{varg})
                        obj.param_max = varargin{varg + 1};
                        varargin{varg} = []; varargin{varg + 1} = [];
                    elseif strcmpi('nStart',varargin{varg})
                        obj.nStart = varargin{varg + 1};
                        varargin{varg} = []; varargin{varg + 1} = [];
                    elseif strcmpi('esttype',varargin{varg})
                        obj.esttype = varargin{varg + 1};
                        varargin{varg} = []; varargin{varg + 1} = [];
                    elseif strcmpi('params',varargin{varg})
                        obj.params = varargin{varg + 1};
                        varargin{varg} = []; varargin{varg + 1} = [];
                    elseif strcmpi('trial',varargin{varg})
                        obj.trial = varargin{varg + 1};
                        varargin{varg} = []; varargin{varg + 1} = [];
                    end
                end
            end
        end
        
        function obj = save(obj,varargin)
            
            % obj = save(obj)
            %
            % -------------------------------------------------------------------------
            % This function saves comp_model class as a .mat file to fname,
            % or user specified path. fullfile(fpath,obj.model)
            %
            % -------------------------------------------------------------------------
            % OPTIONAL INPUTS:
            % -------------------------------------------------------------------------
            % fpath         Specify file name path otherwise uses obj.fname
            %
            % -------------------------------------------------------------------------
            % EXAMPLES:
            % -------------------------------------------------------------------------
            % lin_model.save('~/MATLAB')
            %
            % -------------------------------------------------------------------------
            % Author and copyright information:
            % -------------------------------------------------------------------------
            %     Copyright (C) 2014  Luke Chang
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
            % -------------------------------------------------------------------------
            
            if nargin > 1 %use supplied file name
                if ischar(varargin{1})
                    save(fullfile(varargin{1},[obj.model '.mat']), 'obj')
                end
                
            elseif ~isempty(obj.fname) %use obj.fname
                save(obj.fname, 'obj')
            else
                save([obj.model '.mat'], 'obj')
            end
        end
        
        function obj = write_tables(obj, varargin)
            
            % obj = write_tables(obj, varargin)
            %
            % -------------------------------------------------------------------------
            % This function writes out a separate table for obj.params and
            % obj.trial to csv file. File name will be fullfile(fpath, [obj.model '_Params.csv'])
            %
            % -------------------------------------------------------------------------
            % OPTIONAL INPUTS:
            % -------------------------------------------------------------------------
            % fpath         Specify file name path otherwise uses obj.fname
            %               or current working directory path
            %
            % -------------------------------------------------------------------------
            % EXAMPLES:
            % -------------------------------------------------------------------------
            % write_tables(lin_model)
            %
            % lin_model.write_tables('~/MATLAB')
            %
            % -------------------------------------------------------------------------
            % Author and copyright information:
            % -------------------------------------------------------------------------
            %     Copyright (C) 2014  Luke Chang
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
            % -------------------------------------------------------------------------
            
            if nargin > 1 %use supplied file name
                if ischar(varargin{1})
                    dlmwrite(fullfile(varargin{1},[obj.model '_Params.csv']), obj.params, 'delimiter',',','precision',10) %Params
                    dlmwrite(fullfile(varargin{1},[obj.model '_Trial.csv']), obj.trial, 'delimiter',',','precision',10) %TrialData
                end
                
            elseif ~isempty(obj.fname) %use obj.fname
                dlmwrite(fullfile(obj.fname,[obj.model '_Params.csv']), obj.params, 'delimiter',',','precision',10) %Params
                dlmwrite(fullfile(obj.fname,[obj.model '_Trial.csv']), obj.trial, 'delimiter',',','precision',10) %TrialData
            else
                dlmwrite([obj.model '_Params.csv'], obj.params, 'delimiter',',','precision',10) %Params
                dlmwrite([obj.model '_Trial.csv'], obj.trial, 'delimiter',',','precision',10) %TrialData
            end
        end
        
        function obj = fit_model(obj, varargin)
            
            % model_output = fit_model(obj, varargin)
            %
            % -------------------------------------------------------------------------
            % This function will fit a model (model) using fmincon to a dataset (data)
            % multiple times with random start values (nStart) with the parameters being
            % constrained to a lower bound (param_min) and upper bound (param_max).
            % Estimates a separate parameter to each subject (indicated in 1st column
            % of dataset).  Requires some helper functions from my github repository
            % (https://github.com/ljchang/toolbox/tree/master/Matlab).  Clone this
            % repository and add paths to Matlab.  Requires that the model be a
            % named function and that it can parse the input data.
            %
            % -------------------------------------------------------------------------
            % OPTIONAL INPUTS:
            % -------------------------------------------------------------------------
            % show_subject          Displays Subject ID for every iteration.  Helpful for
            %                       debugging.  Off by default.
            %
            % persistent_fmincon    Continues running fmincon despite error.  Good for
            %                       salvaging data if a model is having difficulty converging
            %
            % -------------------------------------------------------------------------
            % OUTPUTS:
            % -------------------------------------------------------------------------
            % obj                   com_model class instance containing all of the trial-trial data and
            %                       Parameter estimates
            %
            % -------------------------------------------------------------------------
            % EXAMPLES:
            % -------------------------------------------------------------------------
            % lin_model = fit_model(lin_model)
            % lin_model = fit_model(lin_model, 'persistent_fmincon', 'show_subject')
            %
            % -------------------------------------------------------------------------
            % Author and copyright information:
            % -------------------------------------------------------------------------
            %     Copyright (C) 2014  Luke Chang
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
            % -------------------------------------------------------------------------
            
            %--------------------------------------------------------------------------
            % Setup
            %--------------------------------------------------------------------------
            
            global trialout
            
            % Defaults
            showSubject = 0;
            persist = 0;
            
            % Parse Inputs
            for varg = 1:length(varargin)
                if ischar(varargin{varg})
                    if strcmpi(varargin(varg),'show_subject')
                        showSubject = 1;
                        varargin(varg) = [];
                    elseif strcmpi(varargin,'persistent_fmincon')
                        persist = 1;
                        varargin(varg) = [];
                    end
                end
            end
            
            %--------------------------------------------------------------------------
            % Run Model separately for every subject
            %--------------------------------------------------------------------------
            
            Subjects = unique(obj.dat(:,1));
            allout = [];
            for s = 1:length(Subjects)
                if showSubject; display(['Subject ', num2str(Subjects(s))]); end %Show Subject ID for every iteration if requested
                
                sdat = obj.dat(obj.dat(:,1)==Subjects(s),:); %Select subject's data
                
                xpar = zeros(obj.nStart,length(obj.param_min)); fval = zeros(obj.nStart,1); exitflag = zeros(obj.nStart,1); out = {obj.nStart,1};  %Initialize values to workspace
                
                for iter = 1:obj.nStart  %Loop through multiple iterations of nStart
                    
                    %generate random initial starting values for free parameters
                    for ii = 1:length(obj.param_min)
                        ipar(ii) = random('Uniform',obj.param_min(ii),obj.param_max(ii),1,1);
                    end
                    
                    if ~persist
                        [xpar(iter,1:length(obj.param_min)) fval(iter) exitflag(iter) out{iter}]=fmincon(str2func(obj.model), ipar, [], [], [], [], obj.param_min, obj.param_max, [], [], sdat);
                    else
                        try
                            [xpar(iter,1:length(obj.param_min)) fval(iter) exitflag(iter) out{iter}]=fmincon(str2func(obj.model), ipar, [], [], [], [], obj.param_min, obj.param_max, [], [], sdat);
                        catch
                            display('Fmincon Could Not Converge.  Skipping Iteration')
                            xpar(iter,1:length(obj.param_min)) = nan(1,length(obj.param_min));
                            fval(iter) = nan;
                            exitflag(iter) = nan;
                            out{iter} = nan;
                        end
                    end
                end
                
                %Find Best fitting parameter if running multiple starting parameters
                [xParMin, fvalMin] = FindParamMin(xpar, fval);
                
                %output parameters
                params(s,1) = Subjects(s);
                params(s,2:length(xParMin) + 1) = xParMin;
                if obj.esttype == 'LLE'
                    params(s,length(xParMin) + 2) = -fvalMin;
                    params(s,length(xParMin) + 3) = penalizedmodelfit(-fvalMin, size(sdat,1), length(xParMin), 'type', obj.esttype, 'metric', 'AIC');
                    params(s,length(xParMin) + 4) = penalizedmodelfit(-fvalMin, size(sdat,1), length(xParMin), 'type', obj.esttype, 'metric', 'BIC');
                elseif obj.esttype =='SSE'
                    params(s,length(xParMin) + 2) = fvalMin;
                    params(s,length(xParMin) + 3) = penalizedmodelfit(fvalMin, size(sdat,1), length(xParMin), 'type', obj.esttype, 'metric', 'AIC');
                    params(s,length(xParMin) + 4) = penalizedmodelfit(fvalMin, size(sdat,1), length(xParMin), 'type', obj.esttype, 'metric', 'BIC');
                end
                
                %aggregate trials
                allout = [allout; trialout];
            end
            
            %--------------------------------------------------------------------------
            % Collate Output
            %--------------------------------------------------------------------------
            
            obj.params = params;
            obj.trial = allout;
            
        end %Function end
        
        function aic = avg_aic(obj)
            % aic = avg_aic(obj)
            %
            % -------------------------------------------------------------------------
            % This function returns the average AIC values after using
            % model_fit()
            % -------------------------------------------------------------------------
            
            aic = nanmean(obj.params(:,end - 1));
            
        end
        
        function bic = avg_bic(obj)
            % aic = avg_aic(obj)
            %
            % -------------------------------------------------------------------------
            % This function returns the average AIC values after using
            % model_fit()
            % -------------------------------------------------------------------------
            
            bic = nanmean(obj.params(:,end));
        end
        
        function avgparams = avg_params(obj)
            % avgparams = avg_params(obj)
            %
            % -------------------------------------------------------------------------
            % This function returns the average parameter estimates across
            % subjects.  Assumes that parameters start on the second column
            % -------------------------------------------------------------------------
            
            avgparams = nanmean(obj.params(:,2:length(obj.param_min)));
        end
        
        function aic = get_aic(obj)
            % aic = get_aic(obj)
            %
            % -------------------------------------------------------------------------
            % This function returns all subject's AIC values after using
            % model_fit()
            % -------------------------------------------------------------------------
            
            aic = obj.params(:,end - 1);
            
        end
        
        function bic = get_bic(obj)
            % bic = get_bic(obj)
            %
            % -------------------------------------------------------------------------
            % This function returns all subject's BIC values after using
            % model_fit()
            % -------------------------------------------------------------------------
            
            bic = obj.params(:,end);
            
        end
        
        function params = get_params(obj)
            % aic = get_params(obj)
            %
            % -------------------------------------------------------------------------
            % This function returns all subject's estimated parameter values after using
            % model_fit()
            % -------------------------------------------------------------------------
            
            params = obj.params(:,2:length(obj.param_min)+ 1);
            
        end
                    
        function summary = summary(obj,varargin)
            % summary = summary(obj)
            %
            % -------------------------------------------------------------------------
            % This function returns the average parameter estimates across
            % subjects.  Assumes that parameters start on the second
            % column. Also returns average AIC, average BIC, average final
            % minimized value, and number of subjects.
            % -------------------------------------------------------------------------
                           
            sprintf(['Summary of Model: ' obj.model ...
                '\n-----------------------------------------' ...
                '\nAverage Parameters:\t' num2str(nanmean(obj.params(:,2:length(obj.param_min) + 1))) ...
                '\nAverage AIC:\t\t' num2str(nanmean(obj.params(:,end - 1))) ...
                '\nAverage BIC:\t\t' num2str(nanmean(obj.params(:,end))) ...
                '\nAverage ' obj.esttype ':\t\t' num2str(nanmean(obj.params(:,end-3))) ...
                '\nNumber of Subjects:\t' num2str(size(obj.params,1)) ...
                '\n-----------------------------------------'])
            
        end
        
        function f1 = plot(obj, trial, columns, varargin)
            % f1 = plot(obj)
            %
            % -------------------------------------------------------------------------
            % This function will plot the predictions of a model averaged over subjects.
            %
            % -------------------------------------------------------------------------
            % INPUTS:
            % -------------------------------------------------------------------------
            % trial                 Column of obj.trial for X axis
            % columns               Columns of obj.trial to plot
            %
            % -------------------------------------------------------------------------
            % OPTIONAL INPUTS:
            % -------------------------------------------------------------------------
            % 'title'               Followed by title of Plot
            %
            % 'xlabel'              Followed by xlabel name
            %
            % 'ylabel'              Followed by ylabel name
            %
            % 'xticklabel'          Followed by cell array of xticklabel names
            %
            % 'legend'              Followed by cell array of legend names
            %
            % -------------------------------------------------------------------------
            % OUTPUTS:
            % -------------------------------------------------------------------------
            % f1                    figure handle
            %
            % -------------------------------------------------------------------------
            % EXAMPLES:
            % -------------------------------------------------------------------------
            % plot(lin_model, 6, [3,4])
            % plot(lin_model, 6, [3,4], 'title', 'Linear Model', 'xlabel','session', 'ylabel', 'Average BDI', 'legend', {'Predicted','Observed'})
            %
            % -------------------------------------------------------------------------
            % Author and copyright information:
            % -------------------------------------------------------------------------
            %     Copyright (C) 2014  Luke Chang
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
            % -------------------------------------------------------------------------
            
            sub = unique(obj.trial(:,1));
            trials = unique(obj.trial(:,trial));
            
            if nargin < 3
                error('Please add a vector indicating which columns of trial to plot')
            end
            
            counter = 1;
            for c = columns
                for t = 1:length(trials)
                    datmn(t,counter) = nanmean(obj.trial(obj.trial(:,trial)==trials(t),c));
                    datse(t,counter) = nanstd(obj.trial(obj.trial(:,trial)==trials(t),c)) / sqrt(length(sub));
                end
                counter = counter + 1;
            end
            
            figure;
            plot(datmn,'LineWidth',2)
            
            % Parse input for plot parameters
            for varg = 1:length(varargin)
                if ischar(varargin{varg})
                    if strcmpi('title',varargin{varg})
                        eval(['title(''' varargin{varg + 1} ''')']);
                        varargin{varg} = []; varargin{varg + 1} = [];
                    elseif strcmpi('xlabel',varargin{varg})
                        eval(['xlabel(''' varargin{varg + 1} ''')']);
                        varargin{varg} = []; varargin{varg + 1} = [];
                    elseif strcmpi('ylabel',varargin{varg})
                        eval(['ylabel(''' varargin{varg + 1} ''')']);
                        varargin{varg} = []; varargin{varg + 1} = [];
                    elseif strcmpi('legend',varargin{varg})
                        l = varargin{varg+1};
                        str = 'legend({';
                        for i = 1:length(l)
                            str = [str ' ''' l{i} ''''];
                        end
                        str = [str ' })'];
                        eval(str)
                        varargin{varg} = []; varargin{varg + 1} = [];
                    elseif strcmpi('xticklabel',varargin{varg})
                        l = varargin{varg+1};
                        str = 'set(gca,''XTickLabel'' , {';
                        for i = 1:length(l)
                            str = [str ' ''' l{i} ''''];
                        end
                        str = [str ' })'];
                        eval(str)
                        varargin{varg} = []; varargin{varg + 1} = [];
                    end
                end
            end
            
        end %end plot function
        
        function obj = fit_model_xval(obj, varargin)
            
            % model_output = fit_model_xval(obj, varargin)
            %
            % -------------------------------------------------------------------------
            % This function will fit a model (model) using fmincon to a dataset (data)
            % multiple times with random start values (nStart) with the parameters being
            % constrained to a lower bound (param_min) and upper bound (param_max).
            % Estimates a separate parameter to each subject (indicated in 1st column
            % of dataset).  Requires some helper functions from my github repository
            % (https://github.com/ljchang/toolbox/tree/master/Matlab).  Clone this
            % repository and add paths to Matlab.  Requires that the model be a
            % named function and that it can parse the input data.
            %
            % -------------------------------------------------------------------------
            % OPTIONAL INPUTS:
            % -------------------------------------------------------------------------
            % kfolds                Followed by number of folds runs
            %                       cross-validated parameter estimation on
            %                       k folds.  If k=1, then estimates for
            %                       the entire sample.  If k=number of
            %                       subjects then run leave-one-subject out
            %                       (default = LOSO)
            %
            % show_subject          Displays Subject ID for every iteration.  Helpful for
            %                       debugging.  Off by default.
            %
            % persistent_fmincon    Continues running fmincon despite error.  Good for
            %                       salvaging data if a model is having difficulty converging
            %
            %
            % -------------------------------------------------------------------------
            % OUTPUTS:
            % -------------------------------------------------------------------------
            % obj                   com_model class instance containing all of the trial-trial data and
            %                       Parameter estimates
            %
            % -------------------------------------------------------------------------
            % EXAMPLES:
            % -------------------------------------------------------------------------
            % lin_model = fit_model(lin_model)
            % lin_model = fit_model(lin_model, 'persistent_fmincon', 'show_subject')
            %
            % -------------------------------------------------------------------------
            % Author and copyright information:
            % -------------------------------------------------------------------------
            %     Copyright (C) 2014  Luke Chang
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
            % -------------------------------------------------------------------------
            
            %--------------------------------------------------------------------------
            % Setup
            %--------------------------------------------------------------------------
            
            global trialout
            
            % Defaults
            showSubject = 0;
            persist = 0;
            nfolds = unique(obj.dat(:,1)); %use leave one subject out by default
            
            % Parse Inputs
            for varg = 1:length(varargin)
                if ischar(varargin{varg})
                    if strcmpi(varargin(varg),'show_subject')
                        showSubject = 1;
                        varargin(varg) = [];
                    elseif strcmpi(varargin,'persistent_fmincon')
                        persist = 1;
                        varargin(varg) = [];
                    elseif strcmpi(varargin,'kfolds')
                        nfolds = varargin{varg + 1};
                        varargin(varg) = []; varargin(varg + 1) = [];
                    end
                end
            end
            
            %--------------------------------------------------------------------------
            % create cross-validation partitions
            %--------------------------------------------------------------------------
            
            Subjects = unique(obj.dat(:,1));
            
            Ind = cvpartition(length(Subjects),'KFold',nfolds);
            
            %--------------------------------------------------------------------------
            % Estimate and test model using k-fold cross-validation
            %--------------------------------------------------------------------------
            
            for fold = 1:nfolds
                trdata = obj.dat(ismember(obj.dat(:,1),Subjects(Ind.training(s))),:);
                tedata = obj.dat(ismember(obj.dat(:,1),Subjects(Ind.test(s))),:);
                
                
                allout = [];
                
                
                %Find Best fitting parameter if running multiple starting parameters
                [xParMin, fvalMin] = FindParamMin(xpar, fval);
                
                %output parameters
                params(s,1) = Subjects(s);
                params(s,2:length(xParMin) + 1) = xParMin;
                params(s,length(xParMin) + 2) = fvalMin;
                if obj.esttype == 'LLE'
                    params(s,length(xParMin) + 3) = penalizedmodelfit(-fvalMin, size(sdat,1), length(xParMin), 'type', obj.esttype, 'metric', 'AIC');
                    params(s,length(xParMin) + 4) = penalizedmodelfit(-fvalMin, size(sdat,1), length(xParMin), 'type', obj.esttype, 'metric', 'BIC');
                elseif obj.esttype =='SSE'
                    params(s,length(xParMin) + 3) = penalizedmodelfit(fvalMin, size(sdat,1), length(xParMin), 'type', obj.esttype, 'metric', 'AIC');
                    params(s,length(xParMin) + 4) = penalizedmodelfit(fvalMin, size(sdat,1), length(xParMin), 'type', obj.esttype, 'metric', 'BIC');
                end
                
                %aggregate trials
                allout = [allout; trialout];
            end
            
            %--------------------------------------------------------------------------
            % Collate Output
            %--------------------------------------------------------------------------
            
            obj.params = params;
            obj.trial = allout;
            
            %--------------------------------------------------------------------------
            % Subfunctions
            %--------------------------------------------------------------------------
            
            function f = nested_crossvalidation(xpar, trdata, tedata)
                for s = 1:length(Subjects)
                    %                 if showSubject; display(['Subject ', num2str(Subjects(s))]); end %Show Subject ID for every iteration if requested
                    
                    sdat = obj.dat(obj.dat(:,1)==Subjects(s),:); %Select subject's data
                    
                    xpar = zeros(obj.nStart,length(obj.param_min)); fval = zeros(obj.nStart,1); exitflag = zeros(obj.nStart,1); out = {obj.nStart,1};  %Initialize values to workspace
                    
                    for iter = 1:obj.nStart  %Loop through multiple iterations of nStart
                        
                        %generate random initial starting values for free parameters
                        for ii = 1:length(obj.param_min)
                            ipar(ii) = random('Uniform',obj.param_min(ii),obj.param_max(ii),1,1);
                        end
                        
                        if ~persist
                            [xpar(iter,1:length(obj.param_min)) fval(iter) exitflag(iter) out{iter}]=fmincon(str2func(obj.model), ipar, [], [], [], [], obj.param_min, obj.param_max, [], [], sdat);
                        else
                            try
                                [xpar(iter,1:length(obj.param_min)) fval(iter) exitflag(iter) out{iter}]=fmincon(str2func(obj.model), ipar, [], [], [], [], obj.param_min, obj.param_max, [], [], sdat);
                            catch
                                display('Fmincon Could Not Converge.  Skipping Iteration')
                                xpar(iter,1:length(obj.param_min)) = nan(1,length(obj.param_min));
                                fval(iter) = nan;
                                exitflag(iter) = nan;
                                out{iter} = nan;
                            end
                        end
                    end
                end
            end
            
        end %Function end
        
    end %methods
end %class

