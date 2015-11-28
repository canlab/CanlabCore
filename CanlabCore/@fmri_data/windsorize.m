function obj = windsorize(obj, varargin)
% Windsorize an fMRI data object to madlimit Median Absolute Deviations.
% Default = 5 MADs.
% Works across rows and columns.
% Registers this step in history.
%
% :Usage:
% ::
%
%     obj = windsorize(obj, [madlimit])

% ..
%    Calculate and display descriptives
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
