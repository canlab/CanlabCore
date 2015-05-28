% print_summary(D, varargin) 
% prints summaries for every variable, or specified variables
%
% input:
%    - D:  dataset
%    'subj': followed by a cell array of subject level var names, to only see those vars
%    'event': followed by a cell array of event level var names, to only see those vars
%
% if either varargin is unspecified, all variables will be printed
function print_summary(D, varargin)

    fprintf('\n\n --------- DATASET VARS -------- \n\n');
  
    fprintf('%d subjects, %d subject-level vars, %d event-level vars\n', ...
        length(D.Subj_Level.id), length(D.Subj_Level.names), length(D.Event_Level.names));
    
    fprintf('\n\n --------- SUBJECT LEVEL VARS -------- \n\n');
    
    subj_varnames = D.Subj_Level.names;
    svars = find(strcmp('subj', varargin));
    if ~isempty(svars), subj_varnames = varargin{svars+1}; end

    
    for i=1:length(subj_varnames)
        
        vname = subj_varnames{i};
        [var,~,~,descrip] = get_var(D, vname);
        
        if iscell(var)
            % istext
            fprintf('%s (%s): Text. Unique values: %d\t\n', ...
            vname, descrip, length(unique(var)));
        
        else
            % isnumeric
        fprintf('%s (%s): min:%3.2f\t max:%3.2f\t mean:%3.2f\t sd:%3.2f NaNs:%d\n', ...
            vname, descrip, min(var), max(var), nanmean(var), nanstd(var), sum(isnan(var)));
        end
        
    end
    
    
    fprintf('\n\n --------- EVENT LEVEL VARS -------- \n\n');
   
    event_varnames = D.Event_Level.names;
    evars = find(strcmp('event', varargin));
    if ~isempty(evars), event_varnames = varargin{evars+1}; end
    
    for i=1:length(event_varnames)
        vname = event_varnames{i};
        [var,~,~,descrip] = get_var(D, vname);
        fprintf('%s (%s): min:%3.2f\t max:%3.2f\t mean:%3.2f NaNs:%d\n', ...
            vname, descrip, min(min(var)), max(max(var)), nanmean(nanmean(var)), sum(isnan(isnan(var))));
    end
end