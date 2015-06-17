function [ i, seed_new ] = i_random ( ilo, ihi, seed )

%% I_RANDOM returns a random integer in a given range.
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
%    Input, integer ILO, IHI, the minimum and maximum acceptable values.
%
%    Input, integer SEED, a seed for UNIFORM, the random number generator.
%
%    Output, integer I, the randomly chosen integer.
%
%    Output, integer SEED_NEW, the updated random number seed.
%
ilo = round ( ilo ) - 0.5;
ihi = round ( ihi ) + 0.5;

[ r, seed_new ] = r_random ( ilo, ihi, seed );

i = round ( r );

i = max ( i, ilo );
i = min ( i, ihi );

