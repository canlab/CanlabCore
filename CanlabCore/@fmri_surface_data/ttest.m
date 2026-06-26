function sobj = ttest(obj, varargin)
% ttest Voxel-wise (grayordinate-wise) one-sample t-test across maps.
%
% :Usage:
% ::
%     sobj = ttest(obj)
%     sobj = ttest(obj, pthresh, k)
%
% Performs a one-sample t-test at each grayordinate across the maps (columns) of
% obj.dat, mirroring fmri_data.ttest. Computed by delegating to fmri_data.ttest
% on a proxy (each grayordinate treated as a voxel), so the statistics are
% identical; the result is returned as an fmri_surface_data carrying the t-map in
% .dat and the parallel statistics in .additional_info.statistic (.t, .p, .ste,
% .sig, .dfe). Any extra arguments (e.g. threshold) are passed through.
%
% (A dedicated fmri_surface_statistic_image subclass with native threshold()
% support is planned; for now use the .additional_info.statistic fields.)
%
% :Inputs:
%   **obj:** fmri_surface_data with multiple maps (columns = observations).
%
% :Outputs:
%   **sobj:** fmri_surface_data; .dat is the t-statistic per grayordinate.
%
% :Examples:
% ::
%     all_subs = cat(sub1, sub2, sub3, ...);
%     t = ttest(all_subs);
%     surface(t, 'clim', [-5 5]);
%
% :See also: regress, cat, surface, fmri_surface_data

proxy = as_fmri_data_proxy(obj);
st = ttest(proxy, varargin{:});       % -> statistic_image

sobj = rebuild_like(obj, double(st.dat));
sobj.intent = 'dscalar';
sobj.image_names = {'t'};

stat = struct();
stat.t = double(st.dat);
for f = {'p', 'ste', 'sig', 'dfe', 'N'}
    if isprop(st, f{1}) || isfield(st, f{1})
        try, stat.(f{1}) = st.(f{1}); catch, end %#ok<NOCOM>
    end
end
sobj.additional_info.statistic = stat;
sobj.history{end+1} = sprintf('ttest across %d maps (t-map in .dat; p/ste/sig in additional_info.statistic)', size(obj.dat,2));
end
