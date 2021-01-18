function eigenvectors = spectral_layout(a_graph, dim)
%SPECTRAL_LAYOUT Given a graph compute the spectral embedding of the graph.
% This is simply the eigenvectors of the Laplacian of the graph. Here we
% use the normalized Laplacian.
%
% eigenvectors = SPECTRAL_LAYOUT(a_graph, dim)
%
% Parameters
% ----------
% a_graph: sparse matrix of size (n_samples, n_samples)
%     The (weighted) adjacency matrix of the graph as a sparse matrix.
% 
% dim: double
%     The dimension of the space into which to embed.
% 
% Returns
% -------
% eigenvectors: array of shape (n_samples, dim)
%     The spectral embedding of the graph.
%
%   AUTHORSHIP
%   Math Lead & Primary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Secondary Developer: Stephen Meehan <swmeehan@stanford.edu>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Provided by the Herzenberg Lab at Stanford University 
%   License: BSD 3 clause


    n_samples = size(a_graph, 1);
    
    try
        n_components = max(conncomp(graph(a_graph)));

        if n_components > 1
              warning('The adjacency graph is not connected!');
        end

        diag_data = sum(a_graph)';
        D = spdiags(ones(n_samples,1)./sqrt(diag_data), 0, n_samples, n_samples);
        L = speye(n_samples) - D * a_graph * D;
        sL = (L + L')/2;
        %mL = max(L, L');
        k = dim+1;
        if n_samples<=UMAP.EIGEN_LIMIT || ~exist('lobpcg.m', 'file')
            num_lanczos_vectors = floor(max(2*k + 1, sqrt(n_samples)));    
            opts.p = num_lanczos_vectors;
            opts.v0 = ones(n_samples, 1);
            opts.maxit = 5*n_samples;
            opts.tol = 1e-4;
            %[eigenvectors,  eigenvalues] = eigs(L, k, 'sm',opts);
            [eigenvectors,  eigenvalues] = eigs(sL, k, 'sm',opts);
            %[eigenvectors3,  eigenvalues3] = eigs(mL, k, 'sm',opts);
            %opts.issym = true;
            %f = @(x)L*x;
            %[eigenvectors2,  eigenvalues2] = eigs(f,n_samples,k,0,opts);
            eigenvalues = diag(eigenvalues);
        else
            [eigenvectors,  eigenvalues] = lobpcg(randn(n_samples, k), L, 1e-4,5*n_samples);
        end
        [~,eigenorder] = sort(eigenvalues);
        eigenorder = eigenorder(2:end);
        eigenvectors = eigenvectors(:, eigenorder);
    catch ex
        ex.getReport
        warning(['WARNING: spectral initialisation failed! The eigenvector solver '...
            'failed. This is likely due to too small an eigengap. Consider '...
            'adding some noise or jitter to your data. '...
            'Falling back to random initialisation!']);

        eigenvectors = -10 + 20*rand(n_samples, dim);
    end
end