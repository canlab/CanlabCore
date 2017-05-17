function K = emp_feat_map( kern, dat1, dat2, ind1, ind2, kerparam)

% K = emp_feat_map( kern, dat1, dat2, ind1, ind2, kerparam)
% 
% Implements the use of an empirical feature map in a kernel. In a first
% step, the feature vectors \Phi(x) = k_base(x, prototype) are computed.
% Then a second kernel is applied to these feature vectors.
% 
% The kernel parameters contain the data- and the kernel-objects
% needed.
%
% Kernelparameters (with defaults):
% kerparam = { prototypes, kernel_base, kernel_2nd }
% 
% prototypes                -- data object that contains the prototypes
% kernel_base               -- kernel object to compute feature vectors
%                              (needs not to be a proper kernel, btw)
% kernel_2nd = 'linear'     -- kernel object applied to these feature vectors


prototypes = kerparam{ 1};
kernel_base = kerparam{ 2};

if length( kerparam) < 3
    kerparam{ 3} = [];
end
if isempty( kerparam{ 3})
    kernel_2nd = kernel;        % assume linear kernel
else
    kernel_2nd = kerparam{ 3};
end


x2 = get( dat2, ind2);
x1 = get( dat1, ind1);

feat_vec1 = calc( kernel_base, prototypes, x1);
feat_vec2 = calc( kernel_base, prototypes, x2);


K = calc( kernel_2nd, data( feat_vec1), data( feat_vec2));