function K = chisquared(k,d1,d2,ind1,ind2,kerparam),

% DISTANCE 
% 
% chi^2 - distance (I used it for histograms)
%
% d( q, v) = Sum_k (q_k - v_k)^2 / (q_k + v_k) 
%
% no parameters

  K = zeros( length( ind2), length( ind1));
  x1 = get_x( d1, ind1);
  x2 = get_x( d2, ind2);

  %% in C++
  % K = fchisq( x1, x2);
  
  for ii = 1:length( ind2)
    for jj = 1:length( ind1)
      q = x1( jj, :) - x2( ii, :);
      v = x1( jj, :) + x2( ii, :);
      mask = ( v ~= 0);
      K( ii, jj) = (q( mask)./v( mask))*q( mask)';
%      fprintf( 'we''re at: i = %d,  j = %d \r', ii, jj);
    end
  end
  