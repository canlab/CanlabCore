function wh_cols = scn_spm_get_events_of_interest(SPM, varargin)
% Gets events of interest.
%
% :Usage:
% ::
%
%     wh_cols = scn_spm_get_events_of_interest(SPM, varargin)
%
% All regressors, or events only if 'events_only' is input as keyword
% 'from_multireg':  followed by an integer, to include first n columns from
% the multireg R matrix as "of interest".  only works with 'events_only'
% flag, of course.

events_only = false;  % false for all regressors, true for events only
numcols_to_add_from_multireg = 0;

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            % reserved keywords
            case 'events_only', events_only = true;
            case 'from_multireg', numcols_to_add_from_multireg = varargin{i+1};
                
            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end


if events_only
    
    nsess = length(SPM.Sess);
    
    % Tor's code to get only events
    for i = 1:nsess
        ncols_in_run(i, 1) = length(SPM.Sess(i).col);
        
        %n_ofinterest(i) = length(SPM.Sess(i).U);
        % Modified to handle parametric modulators as well. 7/2012
        n_ofinterest(i) = length(cat(2, SPM.Sess(i).U.name));
        
    end
    
    if ~all(n_ofinterest == n_ofinterest(1))
        disp('Warning! Different numbers of events of interest in each run.');
        disp('If this is what you intended, OK, but will use max val and so will include some of no interest')
    end
        
    firstcol = cumsum([1; ncols_in_run(1:end-1)]);
    
    for i = 1:max(n_ofinterest) + numcols_to_add_from_multireg
        wh_cols(:, i) = firstcol + (i - 1);
    end
    
    wh_cols = wh_cols';
    wh_cols = wh_cols(:);
    
    wh_cols = SPM.xX.iC(wh_cols);  % iC should be index of all...so redundant, but include for logical consistency
else
    wh_cols = SPM.xX.iC;
end

end % function

%%
