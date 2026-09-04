function caps = rsa_toolbox_status(varargin)
% rsa_toolbox_status  Report which optional RSA dependencies are available.
%
% The CanlabCore RSA tools work standalone, but use the Kriegeskorte
% rsatoolbox (and a few File Exchange utilities) for enhanced functionality
% when they are on the MATLAB path. This prints a status report of what is
% detected and which built-in fallback will be used otherwise.
%
% Usage
% -----
%   rsa_toolbox_status            % print a report
%   caps = rsa_toolbox_status;    % also return the capability struct
%   caps = rsa_toolbox_status('quiet');   % return struct, no printing
%
% Output
% ------
%   caps  struct of logicals for each detected capability (see
%         @rsm/private/probe_rsatoolbox.m for the field list), plus:
%           .ICC          File Exchange ICC.m available
%           .rangesearch  Statistics Toolbox rangesearch (fast searchlight)
%           .fitlme       Statistics Toolbox fitlme (LME modeling)

quiet = ~isempty(varargin) && any(strcmpi(varargin, 'quiet'));

% Reuse the class-private probe by calling it through a tiny shim rsm method.
% probe_rsatoolbox lives in @rsm/private, so call it via a helper.
caps = local_probe();

% Additional non-rsatoolbox dependencies
caps.ICC         = exist('ICC', 'file') == 2;
caps.rangesearch = exist('rangesearch', 'file') == 2;
caps.fitlme      = exist('fitlme', 'file') == 2 || exist('fitlme', 'builtin') == 5;

if quiet, return; end

% =========================================================================
% Print report
% =========================================================================
fprintf('\n=== CanlabCore RSA toolbox status ===\n\n');

fprintf('Kriegeskorte rsatoolbox (optional enhancements):\n');
print_line('rankTransform / RDMcolormap (rank-transformed plots)', caps.rdm_rankTransform || caps.util_RDMcolormap, 'plain imagesc + parula');
print_line('categoricalRDM (model RDM construction)',             caps.rdm_categoricalRDM, 'built-in same-vs-different');
print_line('covdiag (Ledoit-Wolf whitening)',                     caps.stat_covdiag, 'built-in ledoit_wolf_shrinkage');
print_line('MDSConditions (MDS plots)',                           caps.fig_MDSConditions, 'built-in cmdscale');
print_line('compareRefRDM2candRDMs (formal RDM inference)',       caps.compareRefRDM2candRDMs, 'built-in rsm.compare engine');
print_line('bootstrapRDMs',                                       caps.bootstrapRDMs, 'built-in within-subject bootstrap');

fprintf('\nOther dependencies:\n');
print_line('ICC.m (File Exchange reliability)',         caps.ICC,         'built-in ICC fallback');
print_line('rangesearch (Statistics Toolbox)',          caps.rangesearch, 'O(n^2) loop (slow searchlight)');
print_line('fitlme (Statistics Toolbox)',               caps.fitlme,      'REQUIRED for rsa_lme -- no fallback');

if ~caps.fitlme
    fprintf('\n  WARNING: fitlme not found. The LME methods (rsa_lme, rsa_compare_models,\n');
    fprintf('  rsa_parcelwise ''lme'' path) require the Statistics and Machine Learning Toolbox.\n');
end

if caps.any
    fprintf('\nrsatoolbox detected: enhanced plotting/inference available.\n');
else
    fprintf('\nrsatoolbox NOT detected: all RSA tools work via built-in fallbacks.\n');
    fprintf('  To enable enhancements: addpath(genpath(''path/to/rsatoolbox''))\n');
end
fprintf('\n');

end


function print_line(label, available, fallback)
if available
    status = 'FOUND  ';
    note   = '';
else
    status = 'missing';
    note   = sprintf('  -> using: %s', fallback);
end
fprintf('  [%s] %-52s%s\n', status, label, note);
end


function caps = local_probe()
% Mirror of @rsm/private/probe_rsatoolbox.m (which is class-private and not
% directly callable from here). Kept in sync manually.
caps = struct();
caps.rdm_rankTransform        = ~isempty(which('rsa.rdm.rankTransform')) || ~isempty(which('rankTransform'));
caps.rdm_squareRDM            = ~isempty(which('rsa.rdm.squareRDM'))     || ~isempty(which('squareRDM'));
caps.rdm_categoricalRDM       = ~isempty(which('rsa.rdm.categoricalRDM'))|| ~isempty(which('categoricalRDM'));
caps.util_RDMcolormap         = ~isempty(which('rsa.util.RDMcolormap')) || ~isempty(which('RDMcolormap'));
caps.fig_MDSConditions        = ~isempty(which('rsa.fig.MDSConditions'))|| ~isempty(which('MDSConditions'));
caps.fig_dendrogramConditions = ~isempty(which('rsa.fig.dendrogramConditions')) || ~isempty(which('dendrogramConditions'));
caps.stat_covdiag             = ~isempty(which('rsa.stat.covdiag'))     || ~isempty(which('covdiag'));
caps.compareRefRDM2candRDMs   = ~isempty(which('rsa.compareRefRDM2candRDMs')) || ~isempty(which('compareRefRDM2candRDMs'));
caps.bootstrapRDMs            = ~isempty(which('rsa.bootstrapRDMs'))    || ~isempty(which('bootstrapRDMs'));
caps.any                      = any(structfun(@(x) x, caps));
end
