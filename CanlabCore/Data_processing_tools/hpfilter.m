% :Usage:
% ::
%
%     [y, pI, S, I] = hpfilter(y, TR, HP, spersess, varargin)
%
% :Inputs:
%
%   **y:**
%        data to be filtered
%
%   **TR:**
%        TR in seconds
%
%   **HP:**
%        Highpass filter cutoff in seconds
%
%   **spersess:**
%        scans per sessions; vector with length = # of sessions
%
% :Optional inputs: in fixed order -
%
%   **pI:**
%        pseudoinverse of I matrix, for fast application to new y 
%
%   **dummy:**
%        vector of fixed time points to model with dummy regressors in each run
%        e.g., 1:2 models first two time points in each run with
%        separate dummy regressor (appended to I)
%
%   **dooutliers:**
%        flag (1/0).  If 1, imputes session mean to time points
%        with data > 4 median absolute deviations from the session median (after
%        filtering). Not part of scnlab standard preprocessing. Use with caution.
%
% :Outputs:
%
%   **y:**
%        filtered data, session intercepts removed
%
%   **pI:**
%        intercept model X*pinv, such that y - pI * y removes intercept
%
%   **I:**
%        intercept and dummy scan design matrix
%
%   **S:**
%        smoothing model, such that S * y does HP filtering
%           if sess lengths are unequal, make as long as the longest
%           session; may do funny things to shorter ones?
%
% :Examples:
% ::
%
%    % slowest, creates intercept and smoothing matrix
%    [y, I, S] = hpfilter(data, 2, 120, [169 169 169 173 173]);
%
%    % for subsequent voxels with the same session info,
%    y = hpfilter(data, [], S, [169 169 169 173 173], I);
%
%    y = clusters(1).indiv_timeseries(:, 1);
%    [y, I, S] = hpfilter(y, 2, 120, [169 169 169 173 173]);
%    y = clusters(1).indiv_timeseries(:, 1);
%    [y] = hpfilter(y, [], S, [169 169 169 173 173], I);
%
%    % Regress out average activity in first 2 scans of each session (artifacts)
%    y = hpfilter(raw, 2, 100, EXPT.FIR.nruns, [], 1:2);
%
%    % Another example set up first, then run on multi-voxel matrix:
%    [y, pI, S, I] = hpfilter(data{1}(:,1), HPDESIGN.TR, HPDESIGN.HP, HPDESIGN.spersess, [], 1);
%    tic, y = hpfilter(data{1}, [], S, HPDESIGN.spersess, pI); toc
%
%    % But if you just have one matrix, no need to set up, so this is just as fast:
%    tic , [y, pI, S, I] = hpfilter(data{1}, HPDESIGN.TR, HPDESIGN.HP,HPDESIGN.spersess, [], 1);, toc
%
% ..
%    tor wager
%    Modified 5/12/06 to include dummy images for each session
%    Modified April 07 to add sep. dummy covs for each session and add outlier option
%    Also: y can now be a matrix; hpfilter operates on columns (faster than looping)
% ..

function [y, pI, HP, incpt] = hpfilter(y, TR, HP, spersess, varargin)
    incpt = [];
    pI = [];
    dummyscans = [];
    dooutliers = 0;
    
    spersess = spersess(:)'; % enforce as row vector

    if ~isempty(varargin), pI = varargin{1}; end
    if length(varargin) > 1, dummyscans = varargin{2}; end
    if length(varargin) > 2, dooutliers = varargin{3}; end

    if isempty(pI) || ~ismatrix(pI)
        % if we haven't entered this, then compute it from spersess
        % spersess = [number of images in each run] vector

        incpt = intercept_model(spersess, dummyscans);
        pI = incpt * pinv(incpt);
    end

    n = size(y, 1);
    if size(pI, 1) > n
        disp('Warning: Intercept has more observations than data. spersess is wrong?')
        if size(incpt, 1) > n
            incpt = incpt(1:n, :);
        end
        
        pI = pI(1:n, 1:n);

        tmp = cumsum(spersess);
        while tmp(end) > length(y),
            disp('Trying removal of a session.');
            spersess = spersess(1:end-1);
            tmp = cumsum(spersess);
        end
    end

    % remove intercept
    y = y - pI * y;

    if ~ismatrix(HP) || size(HP, 1) ~= size(y, 1)
        % it's not already an SPM smoothing matrix
        len = max(spersess); % max number of images in a session (run)
        HP = use_spm_filter(TR, len, 'none', 'specify', HP);
    end

    % starting and ending images for each session
    st = cumsum([1 spersess(1:end-1)]);
    en = cumsum(spersess);

    % high-pass filter here
    for i = 1:length(st)
        y(st(i):en(i), :) = HP(1:spersess(i), 1:spersess(i)) * y(st(i):en(i), :);
    end
    
    % because intercept removal is applied before HP filtering, some residual
    % effects of intercepts may remain, so remove these:
    y = y - pI * y;

if dooutliers

    n = size(y,1);

    % identify outliers and replace with timeseries mean
    % (mean has no influence on model betas)
    mady = repmat(mad(y), n, 1);
    wh = (abs(y) > 4 * mady);

    subjectmeans = repmat(mean(y, 1), n, 1);
    y(wh) = subjectmeans(wh);
end

    
end


% ISMATRIX: Returns 1 if the input matrix is 2+ dimensional, 0 if it is a scalar 
%           or vector.
%
%     Usage ismat = ismatrix(X)
%
% RE Strauss, 5/19/00

function ismat = ismatrix(X)
  [r,c] = size(X);
  if (r>1 && c>1)
    ismat = 1;
  else
    ismat = 0;
  end

end
