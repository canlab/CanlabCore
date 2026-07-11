function htmlpath = orthviews_niivue(obj, varargin)
% Lightweight web "orthviews" for a statistic_image: write a self-contained
% interactive NiiVue page and open it in your browser.
%
% This is a thin convenience wrapper around canlab_niivue. Where `orthviews`
% opens the SPM three-plane viewer inside MATLAB, `orthviews_niivue` writes a
% single portable .html (NiiVue, WebGL2) and opens it in your web browser —
% point-and-click slices with colormap / threshold / opacity controls, a
% crosshair MNI-coordinate + value readout, and an atlas region readout that
% names and outlines the region under the crosshair.
%
% It writes a standalone page to a fresh temporary folder by default (so it
% does not clutter your working directory) and opens it. Pass an `'outdir'`
% (and/or `'fname'`) to keep the page somewhere permanent, or `'noopen'` to
% write without launching a browser.
%
% :Usage:
% ::
%
%     htmlpath = orthviews_niivue(obj, [canlab_niivue options])
%
% For objects: Type methods(object_name) for a list of special commands
%              Type help object_name.method_name for help on specific
%              methods.
%
% ..
%     Author and copyright information:
%
%     Copyright (C) 2026  Tor Wager and CANlab
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
% :Inputs:
%
%   **obj:**
%        A statistic_image object. The thresholded map is written by default
%        (so only suprathreshold voxels show), with color limits auto-derived
%        from the suprathreshold magnitude.
%
% :Optional Inputs:
%
%   Any canlab_niivue option/value pair is forwarded as-is, e.g.:
%
%   **'outdir', [char]:**
%        Where to write the page. Default: a fresh temporary folder.
%
%   **'fname', [char]:**
%        HTML file name. Default: 'index.html'.
%
%   **'colormap' / 'cal_min' / 'cal_max' / 'opacity' / 'underlay':**
%        Overlay appearance (see help canlab_niivue).
%
%   **'atlas' / 'noatlas':**
%        Atlas region readout source (on by default). See help canlab_niivue.
%
%   **'noopen':**
%        Write the page but do not launch a browser.
%
% :Outputs:
%
%   **htmlpath:**
%        Full path to the written .html page.
%
% :Examples:
% ::
%
%    % Threshold a group t-map and pop it open in the browser:
%    dat = load_image_set('emotionreg');
%    t   = threshold(ttest(dat), .005, 'unc', 'k', 10);
%    orthviews_niivue(t);
%
%    % Keep the page in a report folder, with fixed color limits, no browser:
%    orthviews_niivue(t, 'outdir', '~/report', 'fname', 'tmap.html', ...
%        'cal_min', 2, 'cal_max', 6, 'noopen');
%
% :References:
%   NiiVue: Hanayik et al., niivue.com (BSD-2-Clause). CANlab: canlab.github.io
%
% :See also:
%   - canlab_niivue
%   - statistic_image.orthviews
%   - image_vector.surface, image_vector.montage, canlab_results_fmridisplay

% ..
%    Programmers' notes:
%    2026/06  Created (Tor Wager / CANlab). Thin wrapper over canlab_niivue;
%             defaults 'outdir' to a unique tempname folder so repeated calls
%             do not overwrite one another and the working directory stays clean.
% ..

% Default the output directory to a unique temp folder unless the caller gave
% one, so the page does not land in (or clobber files in) the working directory.
if ~any(strcmpi(varargin, 'outdir'))
    varargin = [varargin, {'outdir', tempname}];
end

% canlab_niivue defaults to standalone + open-in-browser, which is exactly the
% "lightweight, just show it" behavior we want here.
htmlpath = canlab_niivue(obj, varargin{:});

end % function
