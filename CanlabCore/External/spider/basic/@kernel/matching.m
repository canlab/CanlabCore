function K = matching( k, d1, d2, ind1, ind2, par)

%
%  Matching type similarity function (actually not a kernel) as described in the paper:
%  "Recognition with local features: the kernel recipe"
%  by C. Wallraven, B.Caputo and A.Graf
%
%  takes 2 parameters: { minorKernelName, minorKernelParameter}
%
%  implemented minor kernel names are:
%   'linear', 'rbf' ....  (standard kernels)
%   'normCC' ...........  normalized cross-correlation kernel
%   'spider' ...........  uses a spider kernel object (passed in minorKernelParameter),
%                         (is slower due to object-handling overhead)
%  
  
  global VERBOSITY
  
  x1 = get_x( d1, ind1);
  x2 = get_x( d2, ind2);
  sz1 = size( x1, 1);
  sz2 = size( x2, 1);
  K = zeros( sz2, sz1);
  nTot = sz1*sz2;

  for ii = 1:sz1
    for jj = 1:sz2
      K( jj, ii) = Fkernel( double( x1{ ii}), double( x2{ jj}), par);
    end
    if VERBOSITY Fprogress( (sz2*(ii-1)+jj)/nTot); end
  end

  K = 0.5*( K + K');


function k = Fkernel( a, b, par)

  switch par{ 1}
   case 'linear'
    M = b*a';
    
   case 'rbf'
    M = rbf( a, b, par{ 2});
    
   case 'normCC'
    M = normCC( a, b, par{ 2});
   
   case 'spider'
    M = calc( par{ 2}, data( a), data( b));
    
   otherwise
    error( 'Kernel name unknown or not implemented: %s', par{ 1});
  end

  k = mean( max( M));
  

  
function K = rbf( a, b, sigma)
% exponential kernel 
  sigma2 = sigma * sigma;
  K = b*a';
  K = K + K; % *2
  K1 = sum( a.^2, 2);
  K2 = sum( b.^2, 2);
  K = ones( length( K2), 1)*K1' + K2*ones( 1, length( K1)) - K;
  K = exp( -K ./ ( 2*sigma2));
  
  
function K = normCC( dat1, dat2, rho)
% normalized cross-correlation
  
  dat1 = dat1 - repmat( mean( dat1), size( dat1, 1), 1);
  dat2 = dat2 - repmat( mean( dat2), size( dat2, 1), 1);
  dat1 = dat1 ./ repmat( sqrt( diag( dat1*dat1')), 1, size( dat1, 2));
  dat2 = dat2 ./ repmat( sqrt( diag( dat2*dat2')), 1, size( dat2, 2));
  
  K = exp( -rho*( 1 - dat1*dat2'));
