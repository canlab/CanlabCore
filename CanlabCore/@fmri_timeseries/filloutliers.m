function data_obj = filloutliers(data_obj)
% Find and replace outliers in each voxel of a timeseries object using a
% 60-sec moving median outlier detection with spline interpolation.

if isnan(data_obj.TR) || isempty(data_obj.TR) || data_obj.TR == 0
    error('Enter a TR in data_obj.TR')
end

% Fast, sensible, uses moving average. Use 60-sec window,
% ceil(60/TR)t elements.
win = ceil(60 ./ data_obj.TR);

[newdat, TF] = filloutliers(data_obj.dat', 'spline', 'movmedian', win); 

data_obj.dat = newdat';

n_outliers = sum(TF, 1); % mean outliers per voxel
m = mean(n_outliers);

data_obj.history{end+1} = sprintf('filloutliers: mov-median, spline interp: replaced average %3.2f outliers per voxel, %3.0f%% of values', m, 100 .* m ./ size(data_obj.dat, 2));

end

% programmers' notes:
%
% disp('Linear median')
% tic, newdat = filloutliers(data_obj.dat', 'linear'); toc
% 
% disp('pchip median') %% SOMETHING IS SERIOUSLY WRONG WITH PCHIP MEEDIAN.
% tic, newdat2 = filloutliers(data_obj.dat', 'pchip'); toc
% 
% disp('Linear GESD')
% tic, newdat3 = filloutliers(data_obj.dat', 'linear', 'gesd'); toc
% 
% disp('Linear movmedian')
% tic, newdat4 = filloutliers(data_obj.dat', 'linear', 'movmedian', 30); toc
% 
% % BEST CHOICE: Fast, sensible, uses moving average. Use 60-sec window,
% % ceil(60/TR)t elements.
% disp('spline movmedian')
% tic, newdat5 = filloutliers(data_obj.dat', 'spline', 'movmedian', 30); toc
% 
% disp('pchip GESD')
% tic, newdat6 = filloutliers(data_obj.dat', 'pchip', 'gesd'); toc
% 
% rr = corr([odat(:) newdat(:) newdat2(:) newdat3(:) newdat4(:) newdat5(:) newdat6(:)])
% 
% figure; imagesc(rr), colorbar
% 
% %%
% 
% rr = corr([odat(:) newdat(:) newdat3(:) newdat4(:) newdat5(:) newdat6(:)])
% 
% figure; imagesc(rr), colorbar
%
% piecewise demean at change points???
% 
% filloutliers(resp_trace,’linear’,’movmedian’,100))
% envelope(zscore_resp_trace, 4000,’rms’))
% RV measure, defined as the standard deviation of the treated respiratory trace within a 6-second window, was calculated following (Chang and Glover, 2009) (Matlab command: movstd(zscore_resp_trace,2400,’endpoints’,’shrink’))
% (findpeaks(zscore_resp_trace,’minpeakdistance’,800,’minpeakprominence’,.5))
% 'gesd'      - Applies the generalized extreme Studentized deviate
%                     test for outliers.
% If A is a matrix or a table, filloutliers operates on each column
%     separately.
% generalized extreme studentized deviate (GESD) test (Rosner, 1983)
% Rosner B (1983) Percentage Points for a Generalized ESD Many-Outlier Procedure. Technometrics 25: 165-172. 
% Time/computational efficiency: 'linear' mad default is about 1.8 sec per run with 520 time points. 'pchip' is not much more (2 sec). GESD takes much longer, 70 sec. 'movmedian' is fast, 3 sec.  'movmedian' 'spline' is just as fast.
% 'movmedian' - Returns all elements more than 3 local scaled MAD from
%                     the local median, over a sliding window of length WL.
% 'clip' replaces the outlier with the upper threshold value.
% 'linear'Fills using linear interpolation of neighboring, non-outlier values
% 'spline'Fills using piecewise cubic spline interpolation
% 'pchip'Fills using shape-preserving piecewise cubic spline interpolation
