function varargout = effect_size_map(varargin)
% DEPRECATED wrapper: effect_size_map has been replaced by canlab_effect_size_map
%
% The original effect_size_map.m was a SCRIPT that analyzed between- and
% within-person effect sizes in an SPM directory and reported effect
% sizes across voxels, the optimal balance of within-person and
% between-person sample sizes, and power for a replication. It relied on
% variables (mask, groupcon, groupt, resms_images, SPM, wh_contrast) that
% had to be defined in the base workspace before running it.
%
% It has been rewritten as three documented functions:
%
%   - canlab_effect_size_map          the data-driven analysis (takes a
%                                     folder of first-level SPM analyses
%                                     or fmri_data objects)
%   - canlab_scan_time_allocation     the scanner-hours allocation model
%   - canlab_power_allocation_curves  the normative (no-data) power curves
%                                     and contour plots
%
% This file remains only so that old code that calls effect_size_map
% keeps working: any inputs are passed straight through to
% canlab_effect_size_map. Type help canlab_effect_size_map for usage.
%
% :Usage:
% ::
%
%     [results, maps] = effect_size_map(input_data, [optional inputs])   % same as canlab_effect_size_map
%
% :See also:
%   - canlab_effect_size_map, canlab_scan_time_allocation,
%     canlab_power_allocation_curves, power_from_variance
%

% ..
%    Programmers' notes:
%    2026-09: Tor Wager. Replaced the 2010 script with this deprecation
%    wrapper. See canlab_effect_size_map for the list of fixes.
% ..

warning('effect_size_map:Deprecated', ...
    ['effect_size_map is deprecated and now only forwards to canlab_effect_size_map. ' ...
    'Call canlab_effect_size_map directly (type help canlab_effect_size_map).']);

if nargin == 0
    help canlab_effect_size_map
    return
end

[varargout{1:max(nargout, 1)}] = canlab_effect_size_map(varargin{:});

end
