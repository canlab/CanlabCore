function [X, e, ons] = create_random_er_design(TR, ISI, eventduration, freqConditions, HPlength, dononlin, varargin)
% Create and plot design for one or more event types
%
% % :Usage:
% ::
%
% [X, e, ons] = create_random_er_design(TR, ISI, eventduration, freqConditions, HPlength, dononlin, varargin)
%
%
% :Inputs:
%
% Enter all values in sec:
% TR = time repetition for scans
% ISI = inter-stimulus interval
% eventduration = duration in sec of events
% freqConditions = vector of frequencies of each event type, e.g. [.2 .2 .2] for 3 events at 20% each (remainder is rest)
%   - do not have to sum to one
%
% HPlength is high-pass filter length in sec, or Inf or [] for no filter.
%
% :Outputs:
%
%   **X:**
%        design matrix, sampled at TR. [images x regressors] 
%        Intercept is added as last column.
%
%   **e:**
%        design efficiency
%
%   **ons:**
%        A cell array of onsets and durations (in seconds) for each event
%        type. ons{1} corresponds to Condition 1, ons{2} to Condition 2,
%        and so forth. ons{i} can be an [n x 2] array, where the first
%        column is onset time for each event, and the second column is the event duration
%        (in sec)
%
% Examples:
% create_figure;
% [X, e] = create_random_er_design(1, 1.3, 1, [.2 .2], 180, 0);
% axis tight

% ..
%     Author and copyright information:
%
%     Copyright (C) <year>  <name of author>
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
% ..
%

% -------------------------------------------------------------------------
% DEFAULT ARGUMENT VALUES
% -------------------------------------------------------------------------

scanLength = 200; % in sec

% -------------------------------------------------------------------------
% OPTIONAL INPUTS
% -------------------------------------------------------------------------

% This is a compact way to assign multiple variables. The input argument
% names and variable names must match, however:

allowable_inputs = {'scanLength'};

% keyword_inputs = {'scanlength' 'length'};

% optional inputs with default values - each keyword entered will create a variable of the same name

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            
            case allowable_inputs
                
                eval([varargin{i} ' = varargin{i+1}; varargin{i+1} = [];']);
                
            case 'scanlength', scanLength = varargin{i+1}; varargin{i+1} = [];
            case 'length', scanLength = varargin{i+1}; varargin{i+1} = [];
                
            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end

% -------------------------------------------------------------------------
% OTHER FIXED INPUTS
% -------------------------------------------------------------------------

conditions = 1:length(freqConditions);

if ~isempty(HPlength) && ~isinf(HPlength), dohpfilt = 1; else dohpfilt = 0; end

if dononlin
    nonlinstr = 'nonlinsaturation';
else
    nonlinstr = 'nononlin';
end

nconditions = length(conditions);
%len = ceil(scanLength ./ TR);
ons = cell(1, nconditions);

rowsz = [];
doplot = 0;
basistype = 'spm+disp';
% initalize optional variables to default values here.


% Create onsets
% ----------------------------------------------------------------

ons = create_random_onsets(scanLength, ISI, freqConditions, eventduration);


% Create contrasts
% ----------------------------------------------------------------

% there should be one column per condition in your design, *including* the
% intercept. One row per contrast.

contrasts = create_orthogonal_contrast_set(nconditions);
contrasts(:, end+1) = 0; % for intercept

% there should be one column per condition in your design, *including* the
% intercept. One row per contrast.



% Build Design Matrix
% ----------------------------------------------------------------

% X = onsets2fmridesign(ons, TR, len, 'hrf', nonlinstr);
X = onsets2fmridesign(ons, TR, scanLength, 'hrf', nonlinstr);

% high-pass filtering
if dohpfilt
    
    X(:, 1:end-1) = hpfilter(X(:, 1:end-1), TR, HPlength, round(scanLength ./ TR));
    
end

plotDesign(ons, [], TR, 'samefig', nonlinstr); % % 'durs', eventduration,   This is redundant if we pass in durs

set(gca, 'XLim', [0 scanLength], 'XTick', round(linspace(0, scanLength, 10)));

% Test efficiency
% ----------------------------------------------------------------

e = calcEfficiency(ones(1, size(contrasts, 1)), contrasts, pinv(X), []);

% related to (same as without contrasts, autocorr, filtering):
% 1 ./ diag(inv(X' * X))


end % function




