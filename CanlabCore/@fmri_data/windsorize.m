function obj = windsorize(obj, varargin)
% windsorize Clip extreme values in an fmri_data object to a MAD-based limit.
%
% :Usage:
% ::
%
%     obj = windsorize(obj)
%     obj = windsorize(obj, madlimit)
%
% Replaces any value in obj.dat that lies further than madlimit
% Median Absolute Deviations from the global median with the
% corresponding cap value. Operates across the entire .dat matrix
% (i.e., across both rows/voxels and columns/images jointly). Prints
% a descriptives summary and the count of outlier values, and
% appends an entry to obj.history recording the MAD limit used.
%
% :Inputs:
%
%   **obj:**
%        fmri_data object whose .dat will be windsorized.
%
%   **madlimit:** *(optional)*
%        Numeric scalar; values farther than madlimit MADs from the
%        global median are clipped to median +/- madlimit*MAD.
%        Default: 5.
%
% :Outputs:
%
%   **obj:**
%        fmri_data object with .dat windsorized in place. Caps and
%        outlier counts are printed to the command window. The
%        operation is appended to obj.history as
%        'dat windsorized to <madlimit> MADs'.
%
% :Examples:
% ::
%
%     obj = canlab_get_sample_fmri_data();
%     obj = windsorize(obj);          % default 5 MADs
%     obj = windsorize(obj, 3);       % more aggressive clipping
%
% :See also:
%   - rescale, preprocess (other obj-wide value transforms)
%
% ..
%    Calculate and display descriptives, then clip in place.
% ..

madlimit = 5;
if ~isempty(varargin), madlimit = varargin{1}; end

dattmp = obj.dat(:);
med = median(dattmp);
mabsd = mad(dattmp);
climits = [med - madlimit * mabsd med + madlimit * mabsd];
wh_out = obj.dat < climits(1) | obj.dat > climits(2);
nout = sum(wh_out(:));
percout = 100 * nout ./ prod(size(obj.dat));

nout_by_case = sum(wh_out, 2);

sep = sprintf('____________________________________________\n');
fprintf('%sfmri_data structure .dat field (data)\n%s', sep, sep);
fprintf('Median: %3.3f\nMean: %3.3f\n', med, mean(dattmp));
fprintf('MAD: %3.3f\nSTD: %3.3f\n', mabsd, std(dattmp));
fprintf('Max:%3.3f\nMin:%3.3f\n', max(dattmp), min(dattmp));
fprintf('Control limits: %3.3f to %3.3f\n', climits);
fprintf('Outliers: %3.0f values, %3.2f%% of all values\n', nout, percout);
fprintf('Outliers by case (row): Max = %3.0f, Median = %3.0f, Min = %3.0f\n', max(nout_by_case), median(nout_by_case), min(nout_by_case));

datadj = obj.dat;
datadj(obj.dat < climits(1)) = climits(1);
datadj(obj.dat > climits(2)) = climits(2);

obj.dat = datadj;
obj.history{end+1} = sprintf('dat windsorized to %3.1f MADs', madlimit);


end % function
