function [results,a] =  training(a,d)
       
  % [results,algorithm] =  training(algorithm,data,loss)

  disp(['training ' get_name(a) '.... '])

  %calculate kernel
  K=calc(a.child,d,[]);
  a.Kt=K;
  
  %center kernel in feature space
  if (a.center_data)
    %  this procedure is faster and numerically more stable
    %  than the matrix version (which is cubic):  K=(I-O)*K*(I-O)
    MM = mean( mean( K, 1));
    A = repmat( mean( K, 1), size( K, 1), 1) + repmat( mean( K, 2), 1, size( K, 2));
    K = K - A + MM;
    
    % I=eye(length(K));
    % O=ones(length(K))/length(K);
    % K=(I-O)*K*(I-O);
    
  end
  
  % don't know how many features? compute rank first
  if (a.feat == 0)
    a.feat=rank(K);
  end

  % compute eigensystem% 
  if ( a.feat < size( K, 1) - 1)
    opts.disp=0;
    [ a.e_vec, a.e_val] = eigs( K, a.feat, 'LM', opts);
  else    
    [a.e_vec, a.e_val]= eig(K);
  end
  a.e_val = real( diag( a.e_val));
  %  a.feat=sum(a.e_val>1e-10);

  %sort eigenvalues and eigenvector according to absolute size of eigenvals
  [vals ind] = sort( -abs(a.e_val));
  a.e_val = a.e_val( ind( 1:a.feat));
  a.e_vec = a.e_vec( :, ind( 1:a.feat));

  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% new section to treat degenerate subspaces %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  vals = Ftrunc( a.e_val( 1:a.feat), 8);   % round to 8 significant digits
  if length( unique( vals)) < a.feat    % is there a degenerated subspace ??
    disp(['degenerate subspace encountered ... orthogonalizing']);
    for vv = reshape( unique( vals), 1,length(unique( vals)))
      subind = ( vals == vv);
      if sum( subind) > 1             % do orthonormalization of the degenerate subspace
        subsp = a.e_vec( :, subind);
        a.e_vec( :, subind) = orth( subsp);  % * sqrt( abs( vv));   normalisation is done below  
      end
    end
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  a.e_vec = a.e_vec ./ ( ones( size( a.e_vec, 1), 1) * sqrt( a.e_val)');
  
  a.dat=d;
  results = test(a,d);
