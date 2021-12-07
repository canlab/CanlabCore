function obj = add(obj, varargin)

if nargin == 1
    return
end

% Identify which inputs are valid attribute names:
% Strings that are not preceded by another attribute name.

% all valid fieldnames
valid_names = fieldnames(obj);

[isattributename, isothermethod] = deal(false(1, length(varargin)));

for i = 1:length(varargin)
    if ischar(varargin{i}) && ~isempty(varargin{i}) && (i == 1 || ~isattributename(i - 1))
        
        % Look for a field (attribute) with the input name
        wh = strmatch(varargin{i}, valid_names, 'exact');
        
        % is valid (existing) field?
        if ~isempty(wh)
            isattributename(i) = true;
        else
            % must be another potential input with a special
            % method (defined below); or returns warning
            isothermethod(i) = true;
        end
    end
end

for i = find(isattributename)
    % Replace this field value with the input
    obj.(varargin{i}) = varargin{i + 1};
    
    varargin{i + 1} = []; % eliminate to avoid parsing warnings
    isothermethod(i + 1) = false;
    
    % special methods for specific fields
    switch varargin{i}
        
    end
end

% Now process arguments that may have a special meaning for this object
% class

for i = find(isothermethod)
    
    if ~isothermethod(i)
        % we have removed this because it was a string value for
        % another input
        continue
    end
    
    switch varargin{i}
        
        
        
        case 'units'
            valid_entries = {'secs', 'scans'};
            target_field = 'UNITS';  % field to add to
            
            % check for valid inputs and add if OK
            obj.xBF = add_values_string(obj.xBF, valid_entries, target_field, varargin{i}, varargin{i + 1});
            
            isothermethod(i + 1) = false;
            
        case 'onsets'
            
            % special function for checking SPM-style format rules
            [ons, nsess, nconds] = check_req_and_length(obj, varargin{i+1});
            
            indx = 1;
            for ss = 1:nsess
                for cond = 1:nconds
                    obj.Sess(ss).U(cond).ons = ons{indx};
                    indx = indx + 1;
                end
            end
            
            isothermethod(i + 1) = false;
            
        case 'durations'
            
            % special function for checking SPM-style format rules
            [durs, nsess, nconds] = check_req_and_length(obj, varargin{i+1});
            
            indx = 1;
            for ss = 1:nsess
                for cond = 1:nconds
                    obj.Sess(ss).U(cond).dur = durs{indx};
                    indx = indx + 1;
                end
            end
            
            isothermethod(i + 1) = false;
            
        case 'condition_names'
            
            % special function for checking SPM-style format rules
            [names, nsess, nconds] = check_req_and_length(obj, varargin{i+1});
            
            indx = 1;
            for ss = 1:nsess
                for condname = 1:nconds
                    obj.Sess(ss).U(condname).name = names{indx};
                    indx = indx + 1;
                end
            end
            
            isothermethod(i + 1) = false;
            
        case {'pm', 'mod', 'modulator', 'PM', 'pmname', 'pmnames', 'PMnames'}
            %
            % Format: 'pm', *number of condition to modulate*, '*modulator name*', {cell with modulators for each session}
            
            switch varargin{i}
                case {'pm', 'mod', 'modulator', 'PM'}
                    add2field = 'P';
                    fprintf('Adding modulator values for conditions\n')
                    
                case {'pmname', 'pmnames', 'PMnames'}
                    add2field = 'name';
                    fprintf('Adding modulator names for conditions\n')
                    
                otherwise
            end
            
            [pm, nsess, nconds] = check_req_and_length(obj, varargin{i + 1});
            
            
            
            indx = 1;
            for ss = 1:nsess
                for condname = 1:nconds
                    obj.Sess(ss).U(condname).P.(add2field) = pm{indx};
                    indx = indx + 1;
                end
            end
            
            isothermethod(i + 1) = false;
            
        otherwise
            warning('inputargs:BadInput', 'Unknown field: %s', varargin{i});
            
    end % other valid entries
    
end % process inputs


end





function inputstruct = add_values_string(inputstruct, valid_entries, target_field, infieldname, values)

if isempty(strmatch(values, valid_entries))
    
    fprintf('Invalid entry for %s\nValid entries are:', infieldname);
    disp(char(valid_entries{:}))
    error('Exiting.');
    
else
    
    inputstruct.(target_field) = values;
    
end

end




function [inputvalue, nsess, nconds] = check_req_and_length(obj, inputvalue)

fnames = fieldnames(obj);

if isempty(strmatch('nscan', fnames))
    error('Define number of sessions with .nscan field before adding onsets, durations, names');
end

nsess = length(obj.nscan);

if mod(length(inputvalue), nsess)
    error('Length mismatch: length (n cells) of onsets, etc. must be evenly divisible by num sessions.');
end

nconds = length(inputvalue) ./ nsess;

end

