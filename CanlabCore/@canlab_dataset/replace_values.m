function D = exclude_values(D, varname, function_handle)
% Replace values identified with your custom input function with NaN
%
% e.g., 
% Take the variable 'ratings' in canlab_dataset object D
% Replace all ratings <= 0 with NaN

D = omit_values(D, 'ratings', @(x) x <= 0)
