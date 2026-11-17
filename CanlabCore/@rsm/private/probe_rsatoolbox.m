function caps = probe_rsatoolbox()
% probe_rsatoolbox  Detect which Kriegeskorte rsatoolbox functions are available.
%
% Returns a struct of booleans for each capability the @rsm class may
% optionally use. Methods that want to take advantage of rsatoolbox features
% should call this once and branch on the returned struct.
%
% Output fields:
%   .rdm_rankTransform        rsa.rdm.rankTransform / rankTransform
%   .rdm_squareRDM            rsa.rdm.squareRDM / squareRDM
%   .rdm_categoricalRDM       rsa.rdm.categoricalRDM
%   .util_RDMcolormap         rsa.util.RDMcolormap / RDMcolormap
%   .fig_MDSConditions        rsa.fig.MDSConditions / MDSConditions
%   .fig_dendrogramConditions rsa.fig.dendrogramConditions
%   .stat_covdiag             rsa.stat.covdiag
%   .compareRefRDM2candRDMs   rsa.compareRefRDM2candRDMs
%   .bootstrapRDMs            rsa.bootstrapRDMs
%   .any                      true if any of the above were found
%
% Implementation note: tries `which` for both packaged and bare-function
% forms. The rsatoolbox package layout uses `rsa.rdm.rankTransform`, but
% adding the toolbox via `addpath(genpath(...))` exposes the underlying .m
% files at the top level (`rankTransform.m`) too.

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

caps.any = any(structfun(@(x) x, caps));

end
