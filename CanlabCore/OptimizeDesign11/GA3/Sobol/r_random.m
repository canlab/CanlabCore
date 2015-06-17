function [ r, seed_new ] =  r_random ( rlo, rhi, seed )

%% R_RANDOM returns a random real in a given range.
%
%  Modified:
%
%    21 March 2003
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, real RLO, RHI, the minimum and maximum values.
%
%    Input, integer SEED, a seed for UNIFORM, the random number generator.
%
%    Output, real R, the randomly chosen value.
%
%    Output, integer SEED_NEW, the updated random number seed.
%

%
%  Pick T, a random number in (0,1).
%
  [ t, seed_new ] = uniform ( seed );
%
%  Set R in ( RLO, RHI ).
%
  r = ( 1.0E+00 - t ) * rlo + t * rhi;

