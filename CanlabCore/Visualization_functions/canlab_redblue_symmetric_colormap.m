function canlab_redblue_symmetric_colormap
% This function sets the colormap for the current axis to be
% symmetrical around 0, with positive values shown in red and
% negative values in blue. It is useful for displaying images and
% matrices where the sign of the elements and the distribution of
% positive and negative values is meaningful.
% This is often the case with correlation and covariance matrices,
% and with brain component maps.
% The problem with some default colormaps, e.g., pareto in Matlab,
% is that while they look nice, zero does not reflect a
% neutral/background color (e.g., white as used here), and
% different values within either the positive or negative range
% take on different colors, turning quantitative distinctions into
% qualitative differences in perception. Put another way, the
% difference between [0.2 and 0.3] and [0.3 and 0.4] should look
% similar on the plot, but one may cross a transition from green to
% yellow (and be very salient), whereas the other may not. To
% interpret relative values when comparisons are on the same scale,
% we should use a single color dimension (here, intensity of
% redness for positive or blueness for negative values).
%
% Examples:
% -----------------------------------------------------------------
% Generate sample data
% covmtx = blkdiag(toeplitz([1 .6 .6]), toeplitz([1 .6 .6]));  % two orthogonal latent sources, 3 measures each
% x = mvnrnd(ones(1, 6), covmtx, 500);
% 
% % Image correlation matrix
% create_figure('corr'); 
% imagesc(corr(x)); colorbar; set(gca, 'YDir', 'Reverse'); axis tight;
% canlab_redblue_symmetric_colormap;
% 
% % 
% create_figure('corr2'); 
% plot_correlation_matrix(x);

cm = colormap_tor([0 0 1], [1 .2 0], [1 1 1]);
colormap(gca, cm);

% Make colormap symmetric
cl = get(gca, 'CLim');
mx = max(abs(cl));
cl = [-mx mx];
set(gca, 'CLim', cl);

end
