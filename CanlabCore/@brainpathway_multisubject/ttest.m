function [h,p,ci,stats, p_square, t_square, q_square] = ttest(bs, wh_group, varargin)
%TTEST Conducts a two-sampled ttest between subjects at every node
%   wh_group: a logical array, n_subjects x 1, that defines the two groups.
%   Results show wh_group == T > wh_group == F
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
%
%   Yoni Ashar, March 2020

% vectorize. mat is then subjects x nodes in lower triangle

if ~islogical(wh_group), error('wh_group must be a logical array'), end

doplot = 0;
labels = {'Group == 1', 'Group==0'};

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

% test at each node -- ttest2 tests columnwise
[h,p,ci,stats] = ttest2(mat(wh_group, :), mat(~wh_group, :));

% revert back to square form
t_square = squareform (stats.tstat, 'tomatrix');

p_square = squareform (p, 'tomatrix');
p_square(logical(eye(size(p_square)))) = 1; % set the diagonal p values to be 1

[fdr] = mafdr(p');%, 'BHFDR', true);
q_square = squareform (fdr, 'tomatrix');
q_square(logical(eye(size(q_square)))) = 1; % set the diagonal p values to be 1

if doplot
    create_figure('brainpathways multi t-test',2,2)
    mat = bs.connectivity.regions.r;
    imagesc(mean(mat(:,:,wh_group), 3)); colorbar, title(labels{1})
    subplot(2,2,2)
    imagesc(mean(mat(:,:,~wh_group), 3)); colorbar, title(labels{2})
    subplot(2,2,3)
    imagesc(t_square); colorbar, title('T-statistic: group difference')
    subplot(2,2,4)
    imagesc(q_square); colorbar, title('FDR corrected p values: group difference')
end

end