function embedding = optimize_layout(head_embedding, tail_embedding, ...
    head, tail, n_epochs, n_vertices, epochs_per_sample, a, b, ...
    gamma, initial_alpha, negative_sample_rate, verbose)
%OPTIMIZE_LAYOUT Improve an embedding using stochastic gradient descent to
% minimize the fuzzy set cross entropy between the 1-skeletons of the high
% dimensional and low dimensional fuzzy simplicial sets. In practice this
% is done by sampling edges based on their membership strength (with the
% (1-p) terms coming from negative sampling similar to word2vec). This
% function is only called if the UMAP method is 'MATLAB'.
%
% embedding = OPTIMIZE_LAYOUT(head_embedding, tail_embedding, head, tail,
% n_epochs, n_vertices, epochs_per_sample, a, b)
%
% Parameters
% ----------
% head_embedding: array of size (n_samples, n_components)
%     The initial embedding to be improved by SGD.
% 
% tail_embedding: array of size (source_samples, n_components)
%     The reference embedding of embedded points. If not embedding new
%     previously unseen points with respect to an existing embedding this
%     is simply the head_embedding (again); otherwise it provides the
%     existing embedding to embed with respect to.
% 
% head: array of size (n_1_simplices, 1)
%     The indices of the heads of 1-simplices with non-zero membership.
% 
% tail: array of size (n_1_simplices, 1)
%     The indices of the tails of 1-simplices with non-zero membership.
% 
% n_epochs: double
%     The number of training epochs to use in optimization.
% 
% n_vertices: double
%     The number of vertices (0-simplices) in the dataset.
% 
% epochs_per_samples: array of size (n_1_simplices, 1)
%     A double value of the number of epochs per 1-simplex. 1-simplices with
%     weaker membership strength will have more epochs between being sampled.
% 
% a: double
%     Parameter of differentiable approximation of right adjoint functor.
% 
% b: double
%     Parameter of differentiable approximation of right adjoint functor.
% 
% gamma: double (optional, default 1)
%     Weight to apply to negative samples.
% 
% initial_alpha: double (optional, default 1)
%     Initial learning rate for the SGD.
% 
% negative_sample_rate: double (optional, default 5)
%     Number of negative samples to use per positive sample.
% 
% verbose: boolean (optional, default false)
%     Whether to report information on the current progress of the algorithm.
% 
% Returns
% -------
% embedding: array of size (n_samples, n_components)
%     The optimized embedding.
%
% See also: OPTIMIZE_LAYOUT2
%
%   AUTHORSHIP
%   Math Lead & Primary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Secondary Developer: Stephen Meehan <swmeehan@stanford.edu>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Provided by the Herzenberg Lab at Stanford University 
%   License: BSD 3 clause

    if nargin < 13
        verbose = false;
        if nargin < 12
            negative_sample_rate = 5;
            if nargin < 11
                initial_alpha = 1;
                if nargin < 10
                    gamma = 1;
                end
            end
        end
    end
    
    dim = size(head_embedding, 2);
    n_1_simplices = size(epochs_per_sample, 1);
    same_embedding = isequal(head_embedding, tail_embedding);
    alpha = initial_alpha;
    ONES=ones(1, dim);
    BG2S=2*gamma* b*ONES;
    FOURS=4*ONES;
    ABNEG2=-2.0*a*b;
    BNEG1=b-1;
    
    epochs_per_negative_sample = epochs_per_sample / single(negative_sample_rate);
    epoch_of_next_negative_sample = epochs_per_negative_sample;
    epoch_of_next_sample = epochs_per_sample;
    if verbose
        fprintf('\t0/%d epochs done\n', int32(n_epochs));
    end
    for n = 1:n_epochs
        for i = 1:n_1_simplices 
            if epoch_of_next_sample(i) <= n
                j = head(i);
                k = tail(i);

                current = head_embedding(j,:);
                other = tail_embedding(k,:);

                dist_squared = norm(current - other).^2;

                grad_coeff = (dist_squared > 0)*(ABNEG2*(dist_squared).^BNEG1)./ (a*dist_squared.^b + 1);
                grad_coeff(isnan(grad_coeff)) = 0;

                grad= max(-4, min(4, grad_coeff .* (current - other)));
                current = current + grad * alpha;

                epoch_of_next_sample(i) = epoch_of_next_sample(i) + epochs_per_sample(i);

                n_neg_samples = floor((single(n) - epoch_of_next_negative_sample(i)) / epochs_per_negative_sample(i));
                
                if same_embedding
                    other = other - grad * alpha;
                    head_embedding(k,:) = other;
                    tail_embedding(k,:) = other;
                end

                for p = 1:n_neg_samples
                    k = int32(randi(n_vertices));

                    other = tail_embedding(k,:);

                    dist_squared = norm(current - other).^2;
                    
                    if same_embedding && j == k
                        continue
                    end

                    grad_coeff = (dist_squared > 0)*BG2S./(0.001 + dist_squared)./(a*dist_squared.^b + 1);
                    grad_coeff(isnan(grad_coeff)) = 0;
                    grad=(grad_coeff > 0).*max(-4, min(4, grad_coeff .* (current - other))) + ~(grad_coeff > 0).*FOURS;
                    current = current + grad * alpha;
                end
                
                head_embedding(j,:) = current;
                if same_embedding
                    tail_embedding(j,:) = current;
                end

                epoch_of_next_negative_sample(i) = epoch_of_next_negative_sample(i)+(n_neg_samples * epochs_per_negative_sample(i));
            end

        end
        alpha = initial_alpha * (1 - single(n)/single(n_epochs));
        
        progress_checkpoint = min(floor(n_epochs / 10), 50);
        
        if verbose && mod(n, progress_checkpoint) == 0
            fprintf('\t%d/%d epochs done\n', int32(n), int32(n_epochs));
        end
    end
    
    embedding = head_embedding;
 end
    
