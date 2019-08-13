function OUT = xcorr_multisubject(data, varargin)
% Cross-correlation and partial correlation matrices for 3-D data, i.e., a cell array of subject data matrices
%
% :Usage:
% ::
%
%     OUT = xcorr_multisubject(data, [optional inputs])
%
% ..
%     Author and copyright information:
%     -------------------------------------------------------------------------
%     Copyright (C) 2009 Tor Wager
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
%   **data:**
%        A cell array, one cell per subject/replicate, of n x k
%        data to be inter-correlated.
%
% :Optional Inputs:
%
%   **'partialr':**
%        Use partial correlations obtained via ridge regression (k = 1 fixed)
%
%   **'shift_by':**
%        Followed by integer value for max number of time points to shift
%
% :Output:
%
%   **OUT:**
%        A structure containing subject correlation matrices, the mean
%        matrix, and raw and FDR-thresholded group matrix
%
% :Examples:
% ::
%
%    % cl cell structure:
%    for i = 1:length(clpos_data), data{i} = cat(2, clpos_data{i}.timeseries); end
%
%    % cl structure:
%    for i = 1:size(cl(1).all_data, 2), for c = 1:length(cl), data{i}(:,c) = cl(c).all_data(:, i); end, end
%
%    % parcel_cl_avgs or clpos_data2 structure:
%    for i = 1:length(parcel_cl_avgs), for j = 1:N, data{j}(:,i) = parcel_cl_avgs(i).timeseries{j}; end, end
%
%    create_figure('Xcorr', 1, 3);
%    imagesc(OUT.stats.mean);
%    subplot(1, 3, 2);
%    imagesc(OUT.stats.sig);
%    subplot(1, 3, 3);
%    imagesc(OUT.stats.fdrsig);
%    colormap gray
%
% Example of MDS and plotting total (not decomposed) relationships:
% ::
%
%    OUT.stats.D = (1 - OUT.stats.mean) ./ 2;
%    [OUT.stats_mds.GroupSpace,OUT.stats_mds.obs,OUT.stats_mds.implied_dissim] = shepardplot(OUT.stats.D,[]);
%    OUT.stats_mds = nmdsfig_tools('cluster_solution',OUT.stats_mds, OUT.stats_mds.GroupSpace, 2:5, 1000, []);
%    nmdsfig(OUT.stats_mds.GroupSpace,'classes',OUT.stats_mds.ClusterSolution.classes,'names',OUT.stats_mds.names,'sig',OUT.stats.fdrsig);
%
% Example of MDS and plotting direct relationships:
% ::
%
%    OUT.ridge = matrix_direct_effects_ridge(data);
%    D = OUT.ridge.mean; D(find(eye(size(D)))) = 1;
%    D = (D' + D) ./ 2;
%    OUT.ridge.D = (1 - D) ./ 2;
%    [OUT.stats_mds.GroupSpace,OUT.stats_mds.obs,OUT.stats_mds.implied_dissim] = shepardplot(OUT.ridge.D,[]);
%    OUT.stats_mds = nmdsfig_tools('cluster_solution',OUT.stats_mds, OUT.stats_mds.GroupSpace, 2:10, 1000, []);
%    nmdsfig(OUT.stats_mds.GroupSpace,'classes',OUT.stats_mds.ClusterSolution.classes,'names',OUT.stats_mds.names,'sig',OUT.ridge.fdrsig);
%    hh = nmdsfig_fill(OUT.stats_mds);
%    axis image, axis equal
%
%    OUT.stats_mds = nmdsfig_tools('cluster_solution',OUT.stats_mds, OUT.stats_mds.GroupSpace, 2:5, 1000, []);
%
%    % data is cell, one cell per subject
%    [OUT.stats_mds.GroupSpace,OUT.stats_mds.obs,OUT.stats_mds.implied_dissim] = shepardplot(OUT.stats.D,[]);
%    OUT.stats_mds = nmdsfig_tools('cluster_solution',OUT.stats_mds, OUT.stats_mds.GroupSpace, 2:5, 1000, []);
%
% See also: 
% ttest3d, plot_correlation_matrix -- for more compact ways of estimating
% and plotting correlation matrices. See also cellfun and canlab_mat2cell
% for conversion options. 

% ..
%    Defaults
% ..
shift_by = 0;   % currently, if shift_by = 1, runs cross-correls and returns betas
robustflag = 0; % used if shift_by > 0, can do robust correlations
dopartialr = 0; % partial correlations
betaflag = 1;   % used if shift_by > 0, return betas instead of correlation
nconditions = 1; % fixed at 1 now

% optional inputs with default values
% -----------------------------------

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            
            case 'partialr', dopartialr = 1;
            case 'shift_by', shift_by = varargin{i+1}; varargin{i+1} = [];
                
                %case 'basistype', basistype = varargin{i+1}; varargin{i+1} = [];
                
                %otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end

fprintf(1,'Multi-subject cross correlation\n')

if shift_by
    fprintf(1,'Estimating time-shifted (cross-lagged) correlations\n')
    
elseif dopartialr
    fprintf(1,'Estimating partial correlations\n')
end

fprintf('\t   ')

warning off     % for robustfit iteration limit

numsub = length(data);

for i = 1:numsub            %size(DATA.dat,3);
    
    subjdat = data{i};
    m = size(subjdat, 2);
    
    fprintf(1,'\b\b\b%03d', i);
    
    warning off
    
    for n=1:nconditions
        
        if shift_by > 0
            
            if betaflag
                betastr = 'max_cross_lagged_beta'; % or corr
            else
                betastr = 'max_cross_lagged_correlation';
            end
            
            for j = 1:m - 1
                for k = j + 1 : m
                    % now can handle shift by 0 in shift_correl
                    
                    [sval,myc, mylat, myxc] = shift_correl(subjdat(:,j), subjdat(:,k), shift_by, robustflag, betaflag);
                    
                    if isempty(mylat), mylat = NaN; end
                    
                    nxl(j,k) = mylat;
                    nxc(j,k) = myxc;
                end
            end
        elseif dopartialr
            betastr = 'partial_correlation';
            nxc = calc_partial_r(subjdat);
            nxl = [];
            
        else % %if there's no shift by; fastest
            betastr = 'correlation';
            nxc = corr(subjdat);
            nxl = [];
        end
        
        warning on
        
        if shift_by > 0
            % adjust matrices
            nxc(end+1,:) = 0;
            nxc=nxc + nxc' + eye(size(nxc,1));
            
            nxl(end+1,:) = 0;
            nxl=nxl+nxl';
            
        end
        
        OUT.metric_returned = betastr;
        OUT.shift_by = shift_by;
        OUT.shift_explanation = '0 for no latency est., n for cross-correlations shifting up to n time points forward/back';
        
        OUT.pairwise_assoc{n}(:,:,i) = nxc;
        clear nxc;
        
        OUT.latency{n}(:,:,i)  = nxl;
        clear nxl;
        
    end         % condition loop
    
    
end             % subject loop

fprintf('\n');

fprintf('stats...')
[mxc,t,sig,OUT.stats] = ttest3d(OUT.pairwise_assoc{1});

t = OUT.stats.t .* abs(OUT.stats.fdrsig);
t = (t + t') ./ 2;
OUT.stats.fdr_thresholded_tvalues = t;

fprintf('\n');

end % Main function


function [b p] = calc_partial_r(X)

X = zscore(X);

p = [];

for i = 1:size(X, 2)
    y = X(:, i);
    xx = X;
    xx(:, i) = 1;  % intercept; need it, and also placeholder
    
    % ols version
    b(:, i) = xx \ y;
    b(i, i) = 1;        % use 1 for self, not intercept
    
    % ridge version
    % note: we have double intercept...
    b(:, i) = ridge(y, xx, 1);  
    b(i, i) = 1;        % use 1 for self, not intercept

    %create_figure; plot(b(:, i)); hold on; plot(b1, 'r'); drawnow; pause(.05)
    
end

end % function

