%% SOBOL_TEST01 tests EXOR.
%
fprintf ( 1, '\n' );
fprintf ( 1, 'SOBOL_TEST01\n' );
fprintf ( 1, '  EXOR returns the exclusive OR of two integers.\n' );

fprintf ( 1, '\n' );
fprintf ( 1, '     I     J     EXOR(I,J)\n' );
fprintf ( 1, '\n' );

seed = get_seed ( 0 );

for ( test = 1 : 10 )

  [ i, seed ] = i_random ( 0, 100, seed );
  [ j, seed ] = i_random ( 0, 100, seed );
  k = exor ( i, j );

  fprintf ( 1, '%6d %6d %6d\n', i, j, k );

end 
