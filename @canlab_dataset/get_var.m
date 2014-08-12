function [dat, datcell, wh_level, descrip] = get_var(D, varname, varargin)
%
% get Subject-level or Event-level variable from dataset D and return in rect matrix and cell array.
%
% Usage:
% ----------------------------------------------------------------------------------
% [dat, datcell, wh_level, descrip] = get_var(D, varname, [opt inputs])
%
% Inputs:
% ----------------------------------------------------------------------------------
% D             a canlab_dataset object
% varname       the name of a variable to get from dataset
%               - Looks for var name at either level, returns Event level if exists at both levels
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
% Outputs:
% ----------------------------------------------------------------------------------
% dat: rect matrix of subjects x events
%       - good for plotting individuals, means/std. errors across subjects
%
% datcell: 1 x subjects cell array, each cell containing event data for one subject
%       - good for input into some stats functions, e.g., glmfit_multilevel
%       and igls.m
%
% wh_level: either 'Subject' or 'Event', depending on which level the
% variable exists at.
%
% descrip:  the description for this variable
%
%
% % Copyright Tor Wager, 2013

dat = [];
datcell = {};

if iscell(varname)
    wh_level = get_varlevel(D, varname{1});
else
    wh_level = get_varlevel(D, varname);
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
        end
    end
end


switch wh_level
    case 1
        
        % Subject-level
        if iscell(varname)
            for i=1:length(varname)
                wh = strmatch(varname{i}, D.Subj_Level.names, 'exact');
                dat(:,i) = D.Subj_Level.data(:, wh);
                
                if wh > length(D.Subj_Level.descrip)
                    descrip{i} = varname; %'No description.';
                else
                    descrip{i} = D.Subj_Level.descrip(wh);
                end
            end
        else
            wh = strmatch(varname, D.Subj_Level.names, 'exact');
            dat = D.Subj_Level.data(:, wh);
            if wh > length(D.Subj_Level.descrip)
                descrip = varname; %'No description.';
            else
                descrip = D.Subj_Level.descrip{wh};
            end
        end
        
        if do_conditional
            whrows = D.Subj_Level.data(:, conditionalCol) == conditionalVal;
            dat = dat(whrows, :);
            wh_keep = wh_keep(whrows);
        end
        
    case 2
        
        if iscell(varname)
            wh=[];
            for k=1:length(varname)
                wh(end+1)= find(strcmp(varname{k}, D.Event_Level.names));
            end
        else
            wh = find(strcmp(varname, D.Event_Level.names));
        end
        
        % Event-level
        if do_conditional
            d=conditionalData(D.Event_Level.data, conditionalCol, conditionalVal, wh);
            
            
        else
            my_col = @(x) x(:, wh);
            d = cellfun(my_col, D.Event_Level.data, 'UniformOutput', 0);
            
        end
        
        datcell = d;
        
        dat = cat(2, d{:});  % d is Events x Subjects
        
        dat = dat';  % Subj x Events
        
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

dat = dat(wh_keep,:);
if ~isempty(datcell), datcell = datcell(wh_keep); end


end % function




function varlevel = get_varlevel(D, varname)

[varlevel, varlevel2] = deal(0);

varlevel = any(strcmp(D.Subj_Level.names, varname));

varlevel2 = 2*any(strcmp(D.Event_Level.names, varname));

if any(varlevel & varlevel2)
    error(['Variable ' varname ' exists at both Subject and Event levels. This is not allowed.'])
end

varlevel = varlevel + varlevel2;

if any(~varlevel)
    error(['Variable ' varname ' does not exist in this dataset. Check var input names.'])
end

end

function arr2 = conditionalData(arr, conditionalCol, conditionalVal, outCol)
arr2 = cell(size(arr));
for i=1:length(arr)
    arr2{i} = arr{i}(arr{i}(:,conditionalCol)==conditionalVal, outCol);
end
end

