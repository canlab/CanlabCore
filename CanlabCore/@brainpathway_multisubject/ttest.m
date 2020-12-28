function [h,p,ci,stats, p_square, t_square, q_square, OUT] = ttest(bs, wh_group, varargin)
%TTEST Conducts a one- or two-sampled ttest between subjects at every node
%
% FDR correction: positive false discovery rate (pFDR) from
%   the p-values P of multiple-hypothesis testing using the procedure
%   described by Storey (2002).
%
%   wh_group: a logical array, n_subjects x 1, that defines the two groups.
%             Results show wh_group == T > wh_group == F
%             If empty, performs one-sample t-test
%   
%   n_subjects must be equal to size(bs.connectivity.regions.r, 3). For
%   now, this function defaults to testing the matrices in
%   bs.connectivity.regions.r, but could be expanded later.
%
%   Input:
%       'doplot': show a plot
%       'labels': followed by a 2-element cell array, containing labels for
%                 the top two plots
%
%   Output: 
%       [h,p,ci,stats] from matlab's ttest2
%       p_square, t_square: p and t values, in the square format
%       q_square: fdr corrected p-values, from mafdr()
%       OUT: structure of summary stats compatible with plot_correlation_matrix.m
%
% examples:
%
%  [h,p,ci,stats, p_square, t_square, q_square, OUT] = ttest(obj, []);
%  plot_correlation_matrix(OUT, 'nofigure', varargin{:});
%  title('Mean connectivity, q < 0.05 pFDR thresholded');
%        
%   Yoni Ashar, March 2020

% Programmers' notes:
% Tor Wager updated Dec 2020 to add one-sample t-test

% vectorize. mat is then subjects x nodes in lower triangle

if isempty(which('mafdr')), error('Sorry, you need the Matlab bioinformatics toolbox on your Matlab path to use the Storey 2002 mafdr function called here.'); end

if nargin < 2
    twosample = false;
else
    twosample = ~isempty(wh_group);
end

if twosample && ~islogical(wh_group), error('wh_group must be a logical array'), end

doplot = 0;
labels = {'Mean r, Group == 1', 'Mean r, Group==0'};

for i=1:length(varargin)
    
    if ischar(varargin{i})
        switch varargin{i}
            case 'doplot'
                doplot = 1;
            case 'labels'
                labels = varargin{i+1};
                
        end
    end
end

mat = bs.flatten_conn_matrices();

% test at each node

if twosample
    % -- ttest2 tests columnwise
    [h,p,ci,stats] = ttest2(mat(wh_group, :), mat(~wh_group, :));
else
    [h,p,ci,stats] = ttest(mat);
end

% revert back to square form
t_square = squareform(stats.tstat, 'tomatrix');

p_square = squareform (p, 'tomatrix');
p_square(logical(eye(size(p_square)))) = 1; % set the diagonal p values to be 1

[fdr] = mafdr(p');%, 'BHFDR', true);
q_square = squareform (fdr, 'tomatrix');
q_square(logical(eye(size(q_square)))) = 1; % set the diagonal p values to be 1

if twosample
    
    OUT = struct('descrip', 'Two-sample t-test; r is connectivity diff group 1 - 2; p field is pFDR q-values', ...
        'r', mean(bs.connectivity.regions.r(:,:,wh_group), 3) - mean(bs.connectivity.regions.r(:,:, ~wh_group), 3), ...
        'p', q_square, 'sig', q_square < 0.05, 'wh_group', wh_group);
    
else
    
    OUT = struct('descrip', 'One-sample t-test; r is mean correlation, p field is pFDR q-values', 'r', mean(bs.connectivity.regions.r, 3), 'p', q_square, 'sig', q_square < 0.05);
    
end

if doplot
    
    if twosample
        
        create_figure('brainpathways multi two-sample t-test', 2, 2)
        mat = bs.connectivity.regions.r;
        imagesc(mean(mat(:,:,wh_group), 3)); colorbar, title(labels{1})
        subplot(2,2,2)
        imagesc(mean(mat(:,:,~wh_group), 3)); colorbar, title(labels{2})
        subplot(2,2,3)
        imagesc(t_square); colorbar, title('T-statistic: group difference')
        subplot(2,2,4)
        imagesc(q_square); colorbar, title('FDR corrected p values: group difference')
        
    else
        
        plot_correlation_matrix(OUT, 'nofigure', varargin{:});
        
        title('Mean connectivity, q < 0.05 pFDR thresholded');
        
%         create_figure('brainpathways multi one-sample t-test', 1, 3)
%         mat = bs.connectivity.regions.r;
%         m = mean(mat, 3);
%         
%         imagesc(m);
%         colorbar, title('Mean correlation')
%         subplot(1, 3, 2)
%         
%         imagesc(t_square); colorbar, title('T-statistic: group')
%         
%         subplot(1, 3, 3)
%         mt = m;
%         mt(q_square >= 0.05) = 0;
%         imagesc(mt); colorbar, title('FDR corrected p values: group difference')
        
        drawnow
    end
end

end