%% SOBOL_TEST05 tests SOBOL.
%
global SOBOL_seed;

fprintf ( 1, '\n' );
fprintf ( 1, 'SOBOL_TEST05\n' );
fprintf ( 1, '  SOBOL computes the next element of a Sobol sequence.\n' );
fprintf ( 1, '\n' );
fprintf ( 1, '  In this test, we demonstrate how the SEED can be\n' );
fprintf ( 1, '  manipulated to skip ahead in the sequence, or\n' );
fprintf ( 1, '  to come back to any part of the sequence.\n' );
fprintf ( 1, '\n' );

dim_num = 3;
SOBOL_seed = 0;

fprintf ( 1, '\n' );
fprintf ( 1, '  Using dimension DIM_NUM =   %d\n', dim_num );

seed = 0;

fprintf ( 1, '\n' );
fprintf ( 1, '  Seed  Seed   Sobol\n' );
fprintf ( 1, '  In    Out\n' );
fprintf ( 1, '\n' );

for ( i = 0 : 10 )
  [ r, seed_out ] = sobol ( dim_num, seed );
  fprintf ( 1, '%6d %6d  ', seed, seed_out );
  for ( j = 1 : dim_num )
    fprintf ( 1, '%10f  ', r(j) );
  end
  fprintf ( 1, '\n' );
  seed = seed_out;
end

fprintf ( 1, '\n' );
fprintf ( 1, '  Jump ahead by increasing SEED:\n' );
fprintf ( 1, '\n' );

seed = 100;

fprintf ( 1, '\n' );
fprintf ( 1, '  Seed  Seed   Sobol\n' );
fprintf ( 1, '  In    Out\n' );
fprintf ( 1, '\n' );

for ( i = 1 : 5 )
  [ r, seed_out ] = sobol ( dim_num, seed );
  fprintf ( 1, '%6d %6d  ', seed, seed_out );
  for ( j = 1 : dim_num )
    fprintf ( 1, '%10f  ', r(j) );
  end
  fprintf ( 1, '\n' );
  seed = seed_out;
end

fprintf ( 1, '\n' );
fprintf ( 1, '  Jump back by decreasing SEED:\n' );
fprintf ( 1, '\n' );

seed = 3;

fprintf ( 1, '\n' );
fprintf ( 1, '  Seed  Seed   Sobol\n' );
fprintf ( 1, '  In    Out\n' );
fprintf ( 1, '\n' );

for ( i = 0 : 10 )
  [ r, seed_out ] = sobol ( dim_num, seed );
  fprintf ( 1, '%6d %6d  ', seed, seed_out );
  for ( j = 1 : dim_num )
    fprintf ( 1, '%10f  ', r(j) );
  end
  fprintf ( 1, '\n' );
  seed = seed_out;
end

fprintf ( 1, '\n' );
fprintf ( 1, '  Jump back by decreasing SEED:\n' );
fprintf ( 1, '\n' );

seed = 98;

fprintf ( 1, '\n' );
fprintf ( 1, '  Seed  Seed   Sobol\n' );
fprintf ( 1, '  In    Out\n' );
fprintf ( 1, '\n' );

for ( i = 1 : 5 )
  [ r, seed_out ] = sobol ( dim_num, seed );
  fprintf ( 1, '%6d %6d  ', seed, seed_out );
  for ( j = 1 : dim_num )
    fprintf ( 1, '%10f  ', r(j) );
  end
  fprintf ( 1, '\n' );
  seed = seed_out;
end


