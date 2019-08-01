function [dat, datcell, wh_level, descrip, wh_indx, textflag] = get_var(D, varargin)
% Get Subject-level or Event-level variable from dataset D and return in
% rect matrix and cell array. Multiple variables can be requested, but 
% all data requested must be either numeric or text, and not a combination of the two.
%
% :Usage:
% ::
%
%    [dat, datcell, wh_level, descrip, wh_indx, textflag] = get_var(D)  % to list variables
%
%    [dat, datcell, wh_level, descrip, wh_indx, textflag] = get_var(D, varname, [opt inputs])
%
% :Inputs:
%
%   **D:**
%        a canlab_dataset object
%
%   **varname:**
%        the name of a variable to get from dataset
%           - Looks for var name at either level, returns error if exists at both levels
%           - can be a cell array of multiple var names
%             in this case, dat is a n x m matrix, where n=subjs and m=variables requested
%           - If empty, list variable names
%
% :Optional inputs:
%
%   **a logical array vector of 1/0 values to use as wh_keep**
%
%   **conditional:**
%        to be followed by a cell array; the first cell is the name
%        of the variable to be conditionally selected upon, the second cell
%        contains the condition which must be met.
%        Example: get_var(D, 'DeltaDon', 'conditional', {'trained' 1})
%        will get DeltaDon whenever trained==1.  Currently only implemented for
%        event-level data.  Could be expanded to include multiple conditions.
%
%
% :Outputs:
%
%   **dat:**
%        rect matrix of subjects X events (X variables)
%       - good for plotting individuals, means/std. errors across subjects
%       - returns a cell matrix if textual data is requested.
%       - if there are different numbers of events for each subject, this
%       will return a matrix padded with NaNs to make a rectangular matrix.
%
%   **datcell:**
%        1 x subjects cell array, each cell containing event data for one subject
%       - good for input into some stats functions, e.g., glmfit_multilevel
%         and igls.m
%
%   **wh_level:**
%        1 = 'Subject'; 2 = 'Event';
%
%   **descrip:**
%        the description for this variable
%
%   **wh_indx:**
%        indices of which columns/variables in data matrix are returned
%
%   **textflag:**
%        logical indicator of whether variable is text or numeric
%
% ..
%    Copyright Tor Wager, 2013
% ..
%
% Examples:
% -----------------------------------------------------------
% if DAT is a canlab_dataset object:
% get_var(DAT) to list variable names
%

dat = [];
datcell = {};
wh_indx = [];
wh_level = [];
descrip = {};

if length(varargin) == 0
    % List variables
    disp('Subject-level variables:')
    disp(D.Subj_Level.names(:))
    
    disp('Event-level variables:')
    disp(D.Event_Level.names(:))

    return
else
    varname = varargin{1};
end

if iscell(varname)
    [wh_level,textflag] = get_varlevel(D, varname{1});
else
	[wh_level,textflag] = get_varlevel(D, varname);
end

%varargin:  wh_keep
wh_keep = true(size(D.Subj_Level.id));
do_conditional = 0;


for i = 2:length(varargin)
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
                
                if length(conditionalCol) > 1
                    warning('Multiple variables match!!! Using first one.');
                    conditionalCol = conditionalCol(1);
                end
                
                if isempty(conditionalCol), error('Conditional variable does not exist'); end
                conditionalVal = varargin{i+1}{2};
            otherwise
                % suppress warnings because other functions are passing in
                % many irrelevant args through passing their varargin
                %warning('Unknown varargin: %s\n',varargin{i})
        end
    end
end

% initialize dat
if textflag
    dat = {};
else
    dat = [];
end

nsubj = length(D.Subj_Level.id);  % for empty cells

switch wh_level
    
    case 1 % Subject-level
        if iscell(varname) %Multiple vars to collect
            
            for i=1:length(varname)
                
                wh = strmatch(varname{i}, D.Subj_Level.names, 'exact');
                
                wh_indx = [wh_indx wh];
                
                if textflag
                    % Text data
                    if isempty(D.Subj_Level.textdata)
                        
                        dat(:, i) = repmat({''}, nsubj, 1);
                        
                    else
                        dat(:,i) = D.Subj_Level.textdata(:, wh);
                    end
                    
                else % Numeric data
                    if isempty(D.Subj_Level.data)
                        
                        dat(:, i) = NaN * ones(nsubj, 1); 
                        
                    else
                        dat(:,i) = D.Subj_Level.data(:, wh);
                    end
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
                % Text data
                if isempty(D.Subj_Level.textdata)
                    
                    dat(:, i) = repmat({''}, nsubj, 1);
                        
                else
                    dat(:, i) = D.Subj_Level.textdata(:, wh);
                end
                
            else % Numeric data
                if isempty(D.Subj_Level.data)
                    
                    dat(:, i) = NaN * ones(nsubj, 1); 
                    
                else
                    if iscell(D.Subj_Level.data(:, wh)), error('Subj_Level data should not be in cells.'); end
                
                    dat(:, i) = D.Subj_Level.data(:, wh);
                end
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
        
    case 2 % Event-Level
        
        if iscell(varname)
            wh=[];
            for k=1:length(varname)
                wh(end+1)= find(strcmp(varname{k}, D.Event_Level.names));
            end
        else
            
            wh = find(strcmp(varname, D.Event_Level.names));
            
        end %variable selection
        
        wh_indx = [wh_indx wh];

        if do_conditional
            
            d = conditionalData(D, conditionalCol, conditionalVal, wh, textflag, wh_keep);
            
        elseif ~iscell(D.Event_Level.data) || isempty(D.Event_Level.data)
            % no data
            d = {};
            
        else
            my_col = @(x) x(:, wh);
            if ~textflag
                d = cellfun(my_col, D.Event_Level.data, 'UniformOutput', 0);
            else
                d = cellfun(my_col, D.Event_Level.textdata, 'UniformOutput', 0);
            end
            
        end % conditional if-statement
        
        % Pad cell array if needed
        % ----------------------------------------------------------
        if ~var(cellfun(@numel, d)) == 0
            % cells in d are different length; pad with NaN to make equal
            
             d = pad_cells_with_nan(d, textflag);
             
        end
        
        datcell = d; % Events X Vars within subject cells
        
        % Concatenate
        % ----------------------------------------------------------
        %if var(cellfun(@numel, d)) == 0 % same number of items in every cell, can concat
            
            if ~iscell(varname)
                dat = cat(2, d{:});  
                dat = dat';  % Subj x Events
                
            else %Subj X Events X Vars
                
                if textflag
                    dat = cell(numel(d),length(d{1}(:, 1)),numel(varname));
                else
                    dat = nan(numel(d),length(d{1}(:, 1)),numel(varname));
                end
                
                for subidx = 1:numel(d)
                    dat(subidx,:,:) = d{subidx};
                end
                
            end
            
            % Should not be needed
            %         elseif any(cellfun(@iscell, d))
            %             % Entries are a cell array
            %
            %             dat = 'cannot concat, look at datcell (2nd parameter returned from get_var)';
            %
            %         else % PAD with NaNs to concatenate
            %
            %             slen = max(cellfun(@length, d)); % max length for any subject
            %
            %             refvec = ones(slen, 1);
            %
            %             for i = 1:length(d)
            %                 d{i} = padwithnan(d{i}, refvec, 1);
            %             end
            %
            %             dat = cat(2, d{:});
            %             dat = dat';  % Subj x Events
            %
            %         end
        
        if wh > length(D.Event_Level.descrip)
            descrip = 'No description.';
        else
            descrip = D.Event_Level.descrip{wh};
        end
        
end

% if isempty(dat)
%     disp(['WARNING: ' varname ' is not a valid variable name because it does not exist in this dataset.']);
% end

if iscell(descrip)
    descrip = [descrip{:}];
end

if isempty(dat)
    return
end

%wh_keep
if isempty(wh_keep)
    wh_keep = ones(size(dat, 1));
    
elseif length(wh_keep) ~= size(dat, 1) && length(wh_keep) ~= size(datcell, 2)
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
    case {'numeric', 'Continuous', 'continuous'}
        textflag = 0;
    case 'text'
        textflag = 1;
    otherwise
        textflag = 0;
        warning('Missing variable type for variable %s. Type assumed to be numeric.',varname);
end

end %get_varlevel function

function arr2 = conditionalData(D, conditionalCol, conditionalVal, outCol, textflag, wh_keep)

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
            if wh_keep(i)
                arr2{i} = arr{i}(arr{i}(:,conditionalCol)==conditionalVal, outCol);
            else
                arr2{i} = NaN;
            end
        end
    else %text conditional
        for i=1:length(arr) %iterate over subjects
            if wh_keep(i)
                arr2{i} = arr{i}(strcmp(D.Event_Level.textdata{i}(:,conditionalCol),conditionalVal), outCol);
            else
                arr2{i} = NaN;
            end
        end
    end
else %text data
    if ~ischar(conditionalVal) %numeric conditional
        for i=1:length(arr) %iterate over subjects
            if wh_keep(i)
                arr2{i} = arr{i}(D.Event_Level.data{i}(:,conditionalCol)==conditionalVal, outCol);
            else
                arr2{i} = NaN;
            end
        end
    else %text conditional
        for i=1:length(arr) %iterate over subjects
            if wh_keep(i)
                arr2{i} = arr{i}(strcmp(arr{i}(:,conditionalCol),conditionalVal), outCol);
            else
                arr2{i} = NaN;
            end
        end
    end
end

end %conditionalData function




function d = pad_cells_with_nan(d, textflag)

slen = max(cellfun(@length, d)); % max length for any subject

refvec = ones(slen, 1);

for i = 1:length(d)
    
    if textflag % text var
        reftext = cell(1, abs(length(d{i}) - length(refvec)))';
        [reftext{:}] = deal('');
        
        d{i} = [d{i}; reftext];
        
    else % numeric
        d{i} = padwithnan(d{i}, refvec, 1);
    end
    
end

end

