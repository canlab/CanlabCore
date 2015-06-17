%% SOBOL_TEST02 tests BIT_HI1_BASE_2.
%
fprintf ( 1, '\n' );
fprintf ( 1, 'SOBOL_TEST02\n' );
fprintf ( 1, '  BIT_HI1_BASE_2 returns the location of the high 1 bit.\n' );

fprintf ( 1, '\n' );
fprintf ( 1, '     I     BIT_HI1_BASE_2\n' );
fprintf ( 1, '\n' );

seed = get_seed ( 0 );

for ( test = 1 : 10 )

  [ i, seed ] = i_random ( 0, 100, seed );

  j = bit_hi1_base_2 ( i );

  fprintf ( 1, '%6d %6d %6d\n', i, j );

end 
