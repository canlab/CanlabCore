function [figure_handles, figure_numbers, all_axis_han, is_valid_handle] = activate_figures(o2, varargin)
% activate_figures Activate figures associated with an fmridisplay object.
%
% Activate all figures associated with an fmridisplay object, or a subset
% for specific montage(s).
%
% :Usage:
% ::
%
%     [figure_handles, figure_numbers, all_axis_han, is_valid_handle] = activate_figures(o2)
%     [figure_handles, figure_numbers, all_axis_han, is_valid_handle] = activate_figures(o2, wh_montages)
%
% :Inputs:
%
%   **o2:**
%        An fmridisplay object with montages attached.
%
% :Optional Inputs:
%
%   **wh_montages:**
%        Vector of montage indices to activate. Default: all montages
%        (1:length(o2.montage)). A warning is issued and the function
%        returns early if any requested index exceeds the number of
%        montages.
%
% :Outputs:
%
%   **figure_handles:**
%        Cell array of handles for unique figures associated with the
%        fmridisplay object.
%
%   **figure_numbers:**
%        Vector of figure numbers for unique figures associated with
%        the fmridisplay object.
%
%   **all_axis_han:**
%        Vector of all axis handles associated with montages/figure
%        elements.
%
%   **is_valid_handle:**
%        Vector of which handles are valid (not deleted).
%
% :Limitations:
%
%   No support for surfaces yet; montages only.
%
% :See also:
%   - fmridisplay
%   - montage

% Initialize outputs
% ----------------------------------------------------
figure_handles = {};    % cell array of handles for unique figures associated with fmridisplay object
figure_numbers = [];    % vector of figure numbers for unique figures associated with fmridisplay object
all_fig_han = [];
all_axis_han = [];      % vector of all axis handles associated with montages/figure elements
is_valid_handle = [];   % vector of which handles are valid (not deleted)


% Default - all montages
% ----------------------------------------------------

n = length(o2.montage);
wh_montages = 1:n;

% Process inputs
% ----------------------------------------------------

if length(varargin) > 0 
    wh_montages = varargin{1};
    
    if any(wh_montages) > n
        warning('activate_figures: montage requested does not exist.');
        return
    end
end

% Collect axis handles
% ----------------------------------------------------

all_axis_han = {};

for i = wh_montages
    
all_axis_han{i} = o2.montage{i}.axis_handles;

end

all_axis_han = cat(2, all_axis_han{:});

is_valid_handle = ishandle(all_axis_han);

if ~any(is_valid_handle)
    warning('activate_figures: Cannot activate because all axes were deleted.');
    return
end

% Get figure handles
% ----------------------------------------------------
all_fig_han = get(all_axis_han(is_valid_handle), 'Parent');

fig_number = zeros(1, length(all_fig_han));

% fig numbers for each axis handle

for i = 1:length(all_fig_han) % cell 
   if iscell(all_fig_han)                       % all_fig_han is not always a cell - Fix - Michael Sun 10/29/2021
       fig_number(1, i) = all_fig_han{i}.Number;
   else
       fig_number(1, i) = all_fig_han.Number;
   end
    
end

[figure_numbers, wh] = unique(fig_number);  % unique fig numbers 
figure_handles = all_fig_han(wh)';          % cell array, 1 x number of unique figures

for i = 1:length(all_fig_han)
    if iscell(all_fig_han)                       % all_fig_han is not always a cell - Fix - Michael Sun 10/29/2021
        figure(all_fig_han{i});
    else
        figure(all_fig_han)
    end
end


end % main function



