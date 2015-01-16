function [h, p, ci, stats] = ttest2(D, varname, wh_keep1, wh_keep2, varargin)
% Two sample ttest for two samples of one subject-level variable
%
% ttest2(D, varname, wh_keep1, wh_keep2, varargin)
%
%
% Examples:
% 
%
% Inputs:
% ----------------------------------------------------------------------------------
% D             a canlab_dataset object
% varname       the name of a valid variable to get from dataset
% wh_keep1      subjects forming first sample              
% wh_keep2      subjects forming second sample              
%
% [Optional inputs:]
%
% varargin:     passed directly to MATLAB's ttest2
%               'noverbose' will suppress print out of results and bargraph
%
% Outputs:      as from MATLAB's ttest2, 
% ----------------------------------------------------------------------------------
% Copyright Tor Wager, 2013

if any(wh_keep1 & wh_keep2), warning('YOUR SAMPLES ARE OVERLAPPING!!'); end

verbose=1;
if any(strcmp('noverbose', varargin))
    verbose=0;
    varargin(find(strcmp('noverbose', varargin))) = [];
end

x1 = get_var(D, varname, wh_keep1);
x2 = get_var(D, varname, wh_keep2);

[h, p, ci, stats] = ttest2(x1, x2, varargin{:});

if verbose, ttest2_printout(x1,x2, 1); end

end