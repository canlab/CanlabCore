function K = kondor( k, d1, d2, ind1, ind2, kpara)

%
%  Set-Kernel as described in the paper:
%  "Kernels between sets of vectors"
%  by R.Kondor and T.Jebara
%  
%  takes 4 parameters: { eta (ridge), symmetryFlag, minorKernelName, minorKernelParameter}
%
%  implemented minor kernel names are:
%   'linear', 'rbf' ...... (standard kernels)
%   'linear_direct' ...... no cholesky decomposition, works directly on the input vectors
%   'normCC' ............. normalized cross-correlation kernel
%   'spider' ............. uses a spider kernel object (passed in minorKernelParameter),
%                          (is slower due to object-handling overhead)
%  
%  if the symmetryFlag is set, the kernel matrix is supposed to be symmetric and only the upper half
%  is actually computed
  
  global VERBOSITY

  x1 = get_x( d1, ind1);
  x2 = get_x( d2, ind2);
  sz1 = size( x1, 1);
  sz2 = size( x2, 1);
  K = zeros( sz2, sz1);

  if kpara{ 2} && (sz1 == sz2)   % compute only upper half of the matrix
    nTot = 0.5*sz1*( sz1 + 1);
    for ii = 1:sz1
      for jj = 1:ii
        K( jj, ii) = Fkernel( double( x1{ ii}), double( x2{ jj}), kpara);
      end
      if VERBOSITY Fprogress( (0.5*ii*(ii-1)+jj)/nTot); end
    end
    K = K + K' - diag( diag( K));  % make a full matrix out of the upper half
    
  else  
    nTot = sz1*sz2;
    for ii = 1:sz1
      for jj = 1:sz2
        K( jj, ii) = Fkernel( double( x1{ ii}), double( x2{ jj}), kpara);
      end
      if VERBOSITY Fprogress( (sz2*(ii-1)+jj)/nTot); end
    end
  end
  
  K = exp( K);



  
function k = Fkernel( a, b, par)
  
  eta = par{ 1};
  [ mA, nA] = size( a);
  [ mB, nB] = size( b);
  
  switch par{ 3}
   case 'linear_direct'
    phiA = a;
    phiB = b;
    
   case 'linear'
    K = [ a; b]*[ a; b]';
    reg = max( max( abs( K)));
    K = K + (reg*1e-13).*eye( length( K));
    R = chol( K')';  % gives us one representation in feature space
    phiA = R( 1:mA, :);
    phiB = R( mA+1:end, :);
    
   case 'rbf'
    K = rbf( [ a; b], [ a; b], par{ 4});
    reg = max( max( abs( K)));
    K = K + (reg*1e-13).*eye( length( K));
    R = chol( K')';  % gives us one representation in feature space
    phiA = R( 1:mA, :);
    phiB = R( mA+1:end, :);
    
   case 'normCC'
    K = normCC( [ a; b], [ a; b], par{ 4});
    reg = max( max( abs( K)));
    K = K + (reg*1e-13).*eye( length( K));
    R = chol( K')';  % gives us one representation in feature space
    phiA = R( 1:mA, :);
    phiB = R( mA+1:end, :);
   
   case 'spider'
    K = calc( par{ 4}, data( [ a; b]));
    reg = max( max( abs( K)));
    K = K + (reg*1e-13).*eye( length( K));
    R = chol( K')';  % gives us one representation in feature space
    phiA = R( 1:mA, :);
    phiB = R( mA+1:end, :);
    
   otherwise
    error( 'Kernel name unknown or not implemented: %s', par{3});
  end

  muA = mean( phiA);
  muB = mean( phiB);
  
  Ca = cov( phiA, 1);
  Cb = cov( phiB, 1);

  Ca = Ca + trace( Ca)*eta .* eye( size( Ca));
  Cb = Cb + trace( Cb)*eta .* eye( size( Cb));
  %  Ca = (Ca / trace(Ca)) + (eta / ( mA+mB) .* eye( mA+mB));
  %  Ca = Ca + ( eta / ( mA) .* eye( mA+mB)) * trace( Ca);
  %  Cb = (Cb / trace(Cb)) + (eta / ( mA+mB) .* eye( mA+mB));
  %  Cb = Cb + ( eta / ( mB) .* eye( mA+mB)) * trace( Cb);

 
  [ iCa, logDetA] = inv_logdet_pd( Ca);
  [ iCb, logDetB] = inv_logdet_pd( Cb);
  
  [ Cp, logDetP] = inv_logdet_pd( 0.5*( iCa + iCb));
  muP = 0.5*( iCa*muA' + iCb*muB')';
  
  logk = -0.25*( logDetA + logDetB) - 0.5*logDetP;
  logk = logk -0.25*( muA*iCa*muA' + muB*iCb*muB') + 0.5*muP*Cp*muP';
  k = logk;
  

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
