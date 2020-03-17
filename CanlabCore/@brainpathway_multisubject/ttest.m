function [h,p,ci,stats, p_square, t_square, q_square] = ttest(bs, wh_group, varargin)
%TTEST Conducts a two-sampled ttest between subjects at every node
%   wh_group: a logical array, n_subjects x 1, that defines the two groups.
%   Results show wh_group == T > wh_group == F
%   n_subjects must be equal to size(bs.connectivity.regions.r, 3). For
%   now, this function defaults to testing the matrices in
%   bs.connectivity.regions.r, but could be expanded later.
%
%   'doplot': show a plot
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
if strcmp(varargin{1}, 'doplot'), doplot = 1; end

n_subjects = size(bs.connectivity.regions.r, 3);
n_nodes = size(bs.connectivity.regions.r, 2);

% there is probably a faster way to do this
for i=1:n_subjects
    tmp = tril(bs.connectivity.regions.r(:,:,i), -1);
    mat(i,:) = tmp(tril(true(n_nodes), -1));
end

% test at each node -- ttest2 tests columnwise
[h,p,ci,stats] = ttest2(mat(wh_group, :), mat(~wh_group, :));

% revert back to square form
t_square = tril(ones(n_nodes), -1); % sets indices for assignment in next line
t_square(t_square>0) = stats.tstat; % assign in values in correct (columnwise) order
t_square = t_square + t_square'; % flip to be mirrored on upper tri too
   
p_square = tril(ones(n_nodes), -1); % sets indices for assignment in next line
p_square(p_square>0) = p; % assign in values in correct (columnwise) order
p_square = p_square + p_square'; % flip to be mirrored on upper tri too
p_square = p_square + eye(size(p_square)); % set the diagonal p values to be 1

[fdr] = mafdr(p');%, 'BHFDR', true);
q_square = tril(ones(n_nodes), -1); % sets indices for assignment in next line
q_square(q_square>0) = fdr; % assign in values in correct (columnwise) order
q_square = q_square + q_square'; % flip to be mirrored on upper tri too
q_square = q_square + eye(size(q_square)); % set the diagonal p values to be 1

if doplot
    create_figure('brainpathways multi t-test',2,2)
    mat = bs.connectivity.regions.r;
    imagesc(mean(mat(:,:,wh_group), 3)); colorbar, title('Group == 1')
    subplot(2,2,2)
    imagesc(mean(mat(:,:,~wh_group), 3)); colorbar, title('Group == 0')
    subplot(2,2,3)
    imagesc(t_square); colorbar, title('T-statistic: group difference')
    subplot(2,2,4)
    imagesc(q_square); colorbar, title('FDR corrected p values: group difference')
end

end