function [dat, datcell, wh_level, descrip] = get_var(D, varname, varargin)
%
% get Subject-level or Event-level variable from dataset D and return in rect matrix and cell array.
% Multiple variables can be requested, but all data requested must be either numeric or
% text, and not a combination of the two.
%
% Usage:
% ----------------------------------------------------------------------------------
% [dat, datcell, wh_level, descrip] = get_var(D, varname, [opt inputs])
%
% Inputs:
% ----------------------------------------------------------------------------------
% D             a canlab_dataset object
% varname       the name of a variable to get from dataset
%               - Looks for var name at either level, returns error if exists at both levels
%               - can be a cell array of multiple var names
%                  in this case, dat is a n x m matrix, where n=subjs
%                  and m=variables requested
% [Optional inputs:]
%
% a vector of 1/0 values to use as wh_keep
%
% 'conditional': to be followed by a cell array; the first cell is the name
% of the variable to be conditionally selected upon, the second cell
% contains the condition which must be met.  Example: get_var(D, 'DeltaDon', 'conditional', {'trained' 1})
%  will get DeltaDon whenever trained==1.  Currently only implemented for
%  event-level data.  Could be expanded to include multiple conditions.
%
%
%
% Outputs:
% ----------------------------------------------------------------------------------
% dat: rect matrix of subjects X events (X variables)
%       - good for plotting individuals, means/std. errors across subjects
%       - is actually a cell matrix if textual data is requested.
%
% datcell: 1 x subjects cell array, each cell containing event data for one subject
%       - good for input into some stats functions, e.g., glmfit_multilevel
%       and igls.m
%
% wh_level: 1 = 'Subject'; 2 = 'Event';
%
% descrip:  the description for this variable
%
%
% % Copyright Tor Wager, 2013

dat = [];
datcell = {};


if iscell(varname)
    [wh_level,textflag] = get_varlevel(D, varname{1});
else
    [wh_level,textflag] = get_varlevel(D, varname);
end

%varargin:  wh_keep
wh_keep = true(size(D.Subj_Level.id));
do_conditional = 0;


for i=1:length(varargin)
    if islogical(varargin{i})
        wh_keep = varargin{i};
    end
    
    if ischar(varargin{i})
        switch varargin{i}
            case 'conditional'
                do_conditional=1;
                
                switch wh_level
                    case 1
                        conditionalCol = strmatch(varargin{i+1}{1}, D.Subj_Level.names, 'exact');
                    case 2
                        conditionalCol = strmatch(varargin{i+1}{1}, D.Event_Level.names, 'exact');
                end
                
                if isempty(conditionalCol), error('Conditional variable does not exist'); end
                conditionalVal = varargin{i+1}{2};
            otherwise
                warning('Unknown varargin: %s\n',varargin{i})
        end
    end
end

% initialize dat
if textflag
    dat = {};
else
    dat = [];
end

switch wh_level
    
    case 1 % Subject-level
        if iscell(varname) %Multiple vars to collect
            for i=1:length(varname)
                wh = strmatch(varname{i}, D.Subj_Level.names, 'exact');
                
                if textflag
                    dat(:,i) = D.Subj_Level.textdata(:, wh);
                else
                    dat(:,i) = D.Subj_Level.data(:, wh);
                end
                if wh > length(D.Subj_Level.descrip)
                    descrip{i} = varname{i}; %'No description.';
                else
                    descrip{i} = D.Subj_Level.descrip(wh);
                end
            end
            
        else %Single variable to collect
            
            i=1; %only 1 variable
            wh = strmatch(varname, D.Subj_Level.names, 'exact');
            if textflag
                dat(:,i) = D.Subj_Level.textdata(:, wh);
            else
                dat(:,i) = D.Subj_Level.data(:, wh);
            end
            if wh > length(D.Subj_Level.descrip)
                descrip = varname; %'No description.';
            else
                descrip = D.Subj_Level.descrip{wh};
            end
        end
        
        if do_conditional
            if ~ischar(conditionalVal)
                whrows = D.Subj_Level.data(:, conditionalCol) == conditionalVal;
            else
                whrows = strcmp(D.Subj_Level.textdata(:, conditionalCol),conditionalVal); %if this doesn't work, need to fix below as well
            end
            dat = dat(whrows, :);
            wh_keep = wh_keep(whrows);
        end %conditional if-statement
        
    case 2 %Event-Level
        if iscell(varname)
            wh=[];
            for k=1:length(varname)
                wh(end+1)= find(strcmp(varname{k}, D.Event_Level.names));
            end
        else
            wh = find(strcmp(varname, D.Event_Level.names));
        end %variable selection
        
        if do_conditional
            d=conditionalData(D, conditionalCol, conditionalVal, wh,textflag);
        else
            my_col = @(x) x(:, wh);
            if ~textflag
                d = cellfun(my_col, D.Event_Level.data, 'UniformOutput', 0);
            else
                d = cellfun(my_col, D.Event_Level.textdata, 'UniformOutput', 0);
            end
        end %conditional if-statement
        
        datcell = d; % Events X Vars within subject cells
          
        if var(cellfun(@numel,d)) == 0 % same number of items in every cell, can concat
            if ~iscell(varname)
                dat = cat(2, d{:});  
                dat = dat';  % Subj x Events
            else %Subj X Events X Vars
                if textflag
                    dat = cell(numel(d),length(d{1}(:,1)),numel(varname));
                else
                    dat = nan(numel(d),length(d{1}(:,1)),numel(varname));
                end
                for subidx = 1:numel(d)
                    dat(subidx,:,:) = d{subidx};
                end
            end
        else % can't concat
            dat = 'cannot concat, look at datcell (2nd parameter returned from get_var)';
        end
        
        if wh > length(D.Event_Level.descrip)
            descrip = 'No description.';
        else
            descrip = D.Event_Level.descrip{wh};
        end
        
end

% if isempty(dat)
%     disp(['WARNING: ' varname ' is not a valid variable name because it does not exist in this dataset.']);
% end

%wh_keep
if length(wh_keep) ~= size(dat, 1) && length(wh_keep) ~= size(datcell, 2)
    error('wh_keep is the wrong length. Check input.');
end

if ~ischar(dat), dat = dat(wh_keep,:,:); end % dat may be a char b/c assign an error msg value to it above
if ~isempty(datcell), datcell = datcell(wh_keep); end


end % function




function [varlevel, textflag] = get_varlevel(D, varname)

[varlevel, varlevel2] = deal(0); %#ok
varlevel = any(strcmp(D.Subj_Level.names, varname));
varlevel2 = 2*any(strcmp(D.Event_Level.names, varname));

if any(varlevel & varlevel2)
    error(['Variable ' varname ' exists at both Subject and Event levels. This is not allowed.'])
end

varlevel = varlevel + varlevel2;

if any(~varlevel)
    error(['Variable ' varname ' does not exist in this dataset. Check var input names.'])
end

switch varlevel
    case 1
        
        wh = strcmp(D.Subj_Level.names, varname);
        if max(find(wh)) <= numel(D.Subj_Level.type)
            type = D.Subj_Level.type{wh};
        else
            type = 'unk';
        end
    case 2
        
        wh = strcmp(D.Event_Level.names, varname);
        if max(find(wh)) <= numel(D.Event_Level.type)
            type = D.Event_Level.type{wh};
        else
            type = 'unk';
        end
end

switch lower(type)
    case {'numeric', 'Continuous'}
        textflag = 0;
    case 'text'
        textflag = 1;
    otherwise
        textflag = 0;
        warning('Missing variable type for variable %s. Type assumed to be numeric.',varname);
end

end %get_varlevel function

function arr2 = conditionalData(D, conditionalCol, conditionalVal, outCol, textflag)

if textflag
    arr = D.Event_Level.textdata;
else
    arr = D.Event_Level.data;
end

arr2 = cell(size(arr));
%First, check for data or textdata type
if ~iscell(arr{1}) %numeric data
    if ~ischar(conditionalVal) %numeric conditional
        for i=1:length(arr) %iterate over subjects
            arr2{i} = arr{i}(arr{i}(:,conditionalCol)==conditionalVal, outCol);
        end
    else %text conditional
        for i=1:length(arr) %iterate over subjects
            arr2{i} = arr{i}(strcmp(D.Event_Level.textdata{i}(:,conditionalCol),conditionalVal), outCol);
        end
    end
else %text data
    if ~ischar(conditionalVal) %numeric conditional
        for i=1:length(arr) %iterate over subjects
            arr2{i} = arr{i}(D.Event_Level.data{i}(:,conditionalCol)==conditionalVal, outCol);
        end
    else %text conditional
        for i=1:length(arr) %iterate over subjects
            arr2{i} = arr{i}(strcmp(arr{i}(:,conditionalCol),conditionalVal), outCol);
        end
    end
end

end %conditionalData function

