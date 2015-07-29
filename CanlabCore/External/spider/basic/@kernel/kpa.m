function K = kpa( k, d1, d2, ind1, ind2, par)

%
%  Set-Kernel based on (K)ernel (P)rincipal (A)ngles
%  as described in the paper:  
%  "Learning over Sets using Kernel Principal Angles", JMLR, page 921
%  by L.Wolf and A.Shashua
%
%  takes 4 parameters: {logScale, symmetryFlag, minorKernelName, minorKernelParameter}
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

  if par{ 2} && (sz1 == sz2)   % compute only upper half of the matrix
    nTot = 0.5*sz1*( sz1 + 1);
    for ii = 1:sz1
      for jj = 1:ii
        K( jj, ii) = Fkernel( double( x1{ ii}), double( x2{ jj}), par);
      end
      if VERBOSITY Fprogress( (0.5*ii*(ii-1)+jj)/nTot); end
    end
    K = K + K' - diag( diag( K));  % make a full matrix out of the upper half
    
  else  
    nTot = sz1*sz2;
    for ii = 1:sz1
      for jj = 1:sz2
        K( jj, ii) = Fkernel( double( x1{ ii}), double( x2{ jj}), par);
      end
      if VERBOSITY Fprogress( (sz2*(ii-1)+jj)/nTot); end
    end
  end  
  K = exp( K.*par{ 1});
  
  
function k = Fkernel( X1, X2, par)
  
  switch par{ 3}
   case 'linear'
    K1 = X1*X1';
    K2 = X2*X2';
    K3 = X1*X2';
    
   case 'rbf'
    K1 = rbf( X1, X1, par{ 4});
    K2 = rbf( X2, X2, par{ 4});
    K3 = rbf( X2, X1, par{ 4});
    
   case 'normCC'
    K1 = normCC( X1, X1, par{ 4});
    K2 = normCC( X2, X2, par{ 4});
    K3 = normCC( X2, X1, par{ 4});
   
   case 'spider'
    K1 = calc( par{ 4}, data( X1));
    K2 = calc( par{ 4}, data( X2));
    K3 = calc( par{ 4}, data( X2), data( X1));
    
   otherwise
    error( 'Kernel name unknown or not implemented: %s', par{ 3});
  end
  
  [ U1, D1] = eig( K1);
  [ U2, D2] = eig( K2);
  
  U1 = U1 .* ( ones( length( U1), 1) * ( diag( D1).^-.5)' );
  U2 = U2 .* ( ones( length( U2), 1) * ( diag( D2).^-.5)' );
  
  [ U, S, V] = svd( U1'*K3*U2, 0);
%  keyboard
  k = sum( log( diag( S)));
%  k = prod( diag( S))^2;
  


function K = rbf( a, b, sigma)
% exponential kernel 
  sigma2 = sigma * sigma;
  K = b*a';
  K = K + K; % *2
  K1 = sum( a.^2, 2);
  K2 = sum( b.^2, 2);
  K = ones( length( K2), 1)*K1' + K2*ones( 1, length( K1)) - K;
  K = exp( -K ./ ( 2*sigma2));
  
  
function K = normCC( dat1, dat2, sigma)
% normalized cross-correlation
  rho = 0.5/sigma*sigma;    % to have same scale as rbf
  
  dat1 = dat1 - repmat( mean( dat1), size( dat1, 1), 1);
  dat2 = dat2 - repmat( mean( dat2), size( dat2, 1), 1);
  dat1 = dat1 ./ repmat( sqrt( diag( dat1*dat1')), 1, size( dat1, 2));
  dat2 = dat2 ./ repmat( sqrt( diag( dat2*dat2')), 1, size( dat2, 2));
  
  K = exp( -rho*( 1 - dat1*dat2'));
