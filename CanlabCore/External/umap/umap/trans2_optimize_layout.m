function embedding = trans2_optimize_layout(head_embedding, tail_embedding, ...
    head, tail, n_epochs, n_vertices, epochs_per_sample, a, b,...
    gamma, initial_alpha, negative_sample_rate, verbose, min_dist, move_point)
%OPTIMIZE_LAYOUT2 Another version of optimize_layout.m that slightly
% optimizes the speed for MATLAB; in particular, several quantities are
% vectorized. This function is only called if the UMAP method is 'MATLAB
% vectorized'.
%
% embedding = OPTIMIZE_LAYOUT2(head_embedding, tail_embedding, head, tail,
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
% See also: OPTIMIZE_LAYOUT
%
%   AUTHORSHIP
%   Math Lead & Primary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Secondary Developer: Stephen Meehan <swmeehan@stanford.edu>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Provided by the Herzenberg Lab at Stanford University 
%   License: BSD 3 clause
%
    if nargin < 15
        move_point = vertcat(true(size(head_embedding, 1), 1), false(size(tail_embedding, 1)-size(head_embedding, 1), 1));
        if nargin < 14
            min_dist = 0.3;
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
        end
    end
    dim = size(head_embedding, 2);
    ONES=ones(1, dim);
    BG2S=2*gamma* b*ONES;
    FOURS=4*ONES;
    ABNEG2=-2.0*a*b;
    BNEG1=b-1;
    
    %same_embedding = isequal(head_embedding, tail_embedding);
    alpha = initial_alpha;
    
    epochs_per_negative_sample = epochs_per_sample / single(negative_sample_rate);
    epoch_of_next_negative_sample = epochs_per_negative_sample;
    epoch_of_next_sample = epochs_per_sample;
    N=size(epochs_per_sample, 1);
    weights = ones(N,1)./epochs_per_sample; %We probably should have passed in weights to this instead...
    %tc=tic;
    
            dists = sqrt(sum((head_embedding(head,:) - tail_embedding(tail,:)).^2, 2));
            %spread = 1; %Should pass in spread from UMAP object.

            result = cross_entropy(double(head_embedding), double(tail_embedding), double(head), double(tail), double(weights), a, b, same_embedding);

            disp(['The current cross entropy is ' num2str(result)]);
    
    for n = 1:n_epochs
        l=epoch_of_next_sample <= n;
        idxs=int32(find(l));
        N=int32(sum(l));
        nNegSamples = floor((single(n) - epoch_of_next_negative_sample(l)) ./ epochs_per_negative_sample(l));
        mxNeg=max(nNegSamples);
        randy=int32(randi(n_vertices, N, mxNeg));
        next=epoch_of_next_sample(l) + epochs_per_sample(l);
        next_neg=epoch_of_next_negative_sample(l)+(nNegSamples .* epochs_per_negative_sample(l));        
        %Can NOT vectorize dist_squared OUTSIDE next for loop since head_embedding and tail_embedding
        %change WITHIN next for loop since head and tail contain repeating values
        %   dist_squareds=sqrt(sum((head_embedding(head(l),:)-tail_embedding(tail(l),:)).^2, 2)).^2;
            
        for lIdx = 1:N
            i=idxs(lIdx);
            
            j = head(i);
            k = tail(i);
            
            current = head_embedding(j,:);
            other = tail_embedding(k,:);
            
            dist_squared = norm(current - other).^2;
            %Can NOT vectorize dist_squared OUTSIDE this loop
            %  assert(dist_squareds(li)==dist_squared);
            grad_coeff=(dist_squared > 0)*(ABNEG2*(dist_squared).^BNEG1)./ (a*dist_squared.^b + 1);
            grad_coeff(isnan(grad_coeff)) = 0;
            
            %optimize clip()->grad = clip(grad_coeff .* (current - other));
            grad= max(-4, min(4, grad_coeff .* (current - other)));
            current = current + grad * alpha;
            %if same_embedding
            if move_point(k)
                other = other - grad * alpha;
                head_embedding(k,:) = other;
                tail_embedding(k,:) = other;
            end
            
            %optimize by vectorizing->epoch_of_next_sample(i) = epoch_of_next_sample(i) + epochs_per_sample(i);
            
            %optimize by vectorizing->n_neg_samples = floor((n - epoch_of_next_negative_sample(i)) / epochs_per_negative_sample(i));
            
            n_neg_samples=nNegSamples(lIdx);
            for p = 1:n_neg_samples
                %k = randi(n_vertices);
                k=randy(lIdx,p);
                if move_point(k) && j == k
                    continue
                end
                
                other = tail_embedding(k,:);
                
                dist_squared = norm(current - other).^2;
                
                %optimize->grad_coeff = (dist_squared > 0)*2*gamma* b*ones(1, dim)./(0.001 + dist_squared)./(a*dist_squared.^b + 1);
                grad_coeff=(dist_squared > 0)*BG2S./(0.001 + dist_squared)./(a*dist_squared.^b + 1);
                grad_coeff(isnan(grad_coeff)) = 0;
                
                %optimize-->grad = (grad_coeff > 0).*clip(grad_coeff .* (current - other)) + ~(grad_coeff > 0)*4.*ones(1, dim);
                grad=(grad_coeff > 0).*max(-4, min(4, grad_coeff .* (current - other))) + ~(grad_coeff > 0).*FOURS;
                
                %the ABOVE statement remains faster than this optimization
                %if any(grad_coeff <= 0)
                %    grad=(grad_coeff > 0).*max(-4, min(4, grad_coeff .* (current - other))) + ~(grad_coeff > 0).*FOURS;
                %else
                %    grad=max(-4, min(4, grad_coeff .* (current - other)));
                %end   
                current = current + grad * alpha;
            end
            
            head_embedding(j,:) = current;
            if move_point(j)
                tail_embedding(j,:) = current;
            end
            
            %optimize->epoch_of_next_negative_sample(i) = epoch_of_next_negative_sample(i)+(n_neg_samples * epochs_per_negative_sample(i));
            
        end
        epoch_of_next_sample(l)=next;
        epoch_of_next_negative_sample(l)=next_neg;
        alpha = initial_alpha * (1 - single(n)/single(n_epochs));
        
        progress_checkpoint = min(floor(n_epochs / 10), 50);
        if verbose && mod(n, progress_checkpoint) == 0
            fprintf('\t%d/%d epochs done\n', int32(n), int32(n_epochs));
            
%              dists = sqrt(sum((head_embedding(head,:) - tail_embedding(tail,:)).^2, 2));
%             %spread = 1; %Should pass in spread from UMAP object.
% 
%             result = cross_entropy(double(head_embedding), double(tail_embedding), double(head), double(tail), double(weights), a, b, same_embedding);
% 
%             disp(['The current cross entropy is ' num2str(result)]);
        end
    end
    
    embedding = head_embedding;
end    
    
