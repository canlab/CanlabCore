function res = Fspikernel(seqS,seqT,N,lam,mu,p)

%spikernel kernel function
%
%   This code implements the spikernel function, defined in:
%   Spikernels: Embedding Spiking Neurons in Inner-Product Spaces
%   Lavi Shpigelman, Yoram Singer, Rony Paz and Eilon Vaadia
%   Advances in Neural Information Processing Systems (NIPS) 15
%   MIT Press, Cambridge, MA, 2003.
%
%   res= Fspikernel (seqS,seqT,N,lam,mu,p)
%
%   returns res = K(seqS,seqT) = sum_{i=1}^N p^i K_i(seqS,seqT) where the kernel parameters are:
%
%   Inputs
%       seqS and seqT are the two sequences to be compared.
%           Each row of these sequences is one 'letter' in the sequence (must be of same length)
%           Their column lengths are the sequence lengths (can be of different length).
%       N = max subsequence lengths (the kernel compares subsequences of lengths 1 to N and returns a sum of kernels)
%       lam=\lambda parameter from article
%       mu=\mu parameter from article
%       p= the q prarameter in the article - the parameter that is used in the sum of kernels to weigh them differently
%
%    This implementation assumes that the function d(x,y) is the squared \ell_2 norm
%

% ------------- parameter check -------------------------
lens=size(seqS,1);
lent=size(seqT,1);
q=size(seqS,2);

if (q~=size(seqT,2))
    error('the dimesions of the sequence letters in seqS and seqT must be the same')
end
if (length(N)~=1)
    error('N should be a scalar value');
end
if (length(lam)~=1)
    error('lam should be a scalar value');
end
if (length(mu)~=1)
    error('mu should be a scalar value');
end
if (length(p)~=1)
    error('p should be a scalar value');
end
% --------------------------------------------------------

dmu = mu.^(0.5*( repmat( sum( seqS.^2, 2), 1, lent) + ...
              ( repmat( sum(seqT.^2, 2), 1, lens))' - 2*seqS*seqT'))';

res = fspike( N, p, lam, dmu);
