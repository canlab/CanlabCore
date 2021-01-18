function embedding=spectral_layout_binned(data, init, limit, n_components)
%SPECTRAL_LAYOUT_BINNED Given a large set of data, bin the samples into
% 8192 bins, then compute the adjacency graph of the bin means according to
% UMAP defaults. Then compute the spectral embedding as in
% spectral_layout.m. If the graph is not sufficiently large, this returns
% an empty array.
%
% embedding = SPECTRAL_LAYOUT_BINNED(data, init, limit)
%
% Parameters
% ----------
% data: array of shape (n_samples, n_features)
%     The source data
% 
% init: char array
%     How to initialize the low dimensional embedding.
% 
% limit: double
%     If there are fewer than "limit" data points, return an empty array.
% 
% Returns
% -------
% embedding: array of shape (8192, dim)
%     The spectral embedding of the binned data points.
%
% See also: SPECTRAL_LAYOUT
%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu>
%   Math Lead & Secondary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Provided by the Herzenberg Lab at Stanford University 
%   License: BSD 3 clause
%

embedding=[];
if ~strcmp(init, UMAP.INIT_EIGEN) && limit>=0
    sz=size(data,1);
    needLobpcg=sz>UMAP.EIGEN_LIMIT && ~exist('lobpcg.m', 'file');
    if sz>limit || needLobpcg
        probability_bins=probability_means_weights_ptrs(data);
        umap=UMAP('n_components', n_components);
        umap.init=UMAP.INIT_EIGEN;
        try
            umap.fit_transform(probability_bins.means);
            embedding=umap.embedding(probability_bins.ptrs, :);
            if needLobpcg
                showMsg(Html.WrapHr([...
                    'Larger samples are faster with lobpcg.m.<br>Download it from ' ...
                    'MathWorks File Exchange<br><br>Google "<b>matlab lobpcg</b>"']), ...
                    'MathWorks File Exchange', 'north east+', false, false, 22);
            end
        catch ex
            ex.getReport
        end
    end
end
end