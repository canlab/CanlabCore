function K = spikernel( k, d1, d2, ind1, ind2, kerparam)

% K = spikernel( k, d1, d2, ind1, ind2, kerparam)
% computes the spikernel between two sequences of spiking activity
%
% kerparam{1} = N; max subsequence lengths (the kernel compares subsequences of lengths 1 to N and returns a sum of kernels)
% kerparam{2} = lam; \lambda parameter from article
% kerparam{3} = mu; \mu parameter from article
% kerparam{4} = p; the q parameter in the article - the parameter that is used in the sum of kernels to weigh them differently
% kerparam{5} = bins; the number of bins used;          we need this to reshape the data correctly
% kerparam{6} = nneurons; the number of neurons;    we need this to reshape the data correctly
%
%   This code uses the spikernel function, defined in:
%   Spikernels: Embedding Spiking Neurons in Inner-Product Spaces 
%   Lavi Shpigelman, Yoram Singer, Rony Paz and Eilon Vaadia 
%   Advances in Neural Information Processing Systems (NIPS) 15 
%   MIT Press, Cambridge, MA, 2003. 
%
% for further details see Fspikernel.m

bins = kerparam{5};
nneurons = kerparam{6}; 

xx = get_x( d1, ind1);
yy = get_x( d2, ind2);
K = zeros( size( yy,1), size( xx,1));

for i = 1:size( xx, 1)
    for j = 1:size( yy,1)
        K( j, i) = Fspikernel( reshape( yy( j, :), bins, nneurons), reshape( xx( i,:), bins, nneurons), kerparam{1}, kerparam{2}, kerparam{3}, kerparam{4});
%        K( j, i) = K( i, j);
%       if bins >= 100 disp(sprintf('we''re at: (%d, %d)', i,j)); end   %<-- if it takes longer, print out, where we are
    end
end