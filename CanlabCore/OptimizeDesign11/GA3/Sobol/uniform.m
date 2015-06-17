function [ r, seed_new ] = uniform ( seed )

%% UNIFORM is a uniform random number generator.
%
%  Discussion:
%
%    This routine implements the recursion
%
%      seed_new = 16807 * seed mod ( 2**31 - 1 )
%      uniform = seed_new / ( 2**31 - 1 )
%
%    The integer arithmetic never requires more than 32 bits,
%    including a sign bit.
%
%  Note:
%
%    Development of this routine has been hampered by the fact that
%    MATLAB does not allow a function argument to be modified,
%    and does not make it easy for a function to retain local memory.
%    While later research may turn up methods of dealing with these
%    difficulties, for now I've decided to botch the job and move on!
%
%    Contrary to common practice, the function has two output arguments, 
%    the second of which is the updated seed.  Contrary to common practice,
%    the values of certain internal "constants" are precomputed and
%    the literal values inserted here by hand, rather than allowing the
%    code to compute them on first call, and reuse them thereafter.
%
%    John Burkardt
%
%  Modified:
%
%    07 March 2003
%
%  Reference:
%
%    Paul Bratley, Bennett Fox, L E Schrage,
%    A Guide to Simulation,
%    Springer Verlag, pages 201-202, 1983.
%
%    Bennett Fox,
%    Algorithm 647:
%    Implementation and Relative Efficiency of Quasirandom
%    Sequence Generators,
%    ACM Transactions on Mathematical Software,
%    Volume 12, Number 4, pages 362-376, 1986.
%
%  Parameters:
%
%    Input, integer SEED, the integer "seed" used to generate
%    the output random number.
%
%    Output, real R, a random value between 0 and 1.
%
%    Output, integer SEED_NEW, the updated seed.  This would
%    normally be used as the input seed on the next call.
%
  seed = floor ( seed );

  seed = mod ( seed, 2147483647 );

  if ( seed < 0 ) 
    seed = seed + 2147483647;
  end 

  k = floor ( seed / 127773 );

  seed_new = 16807 * ( seed - k * 127773 ) - k * 2836;

  if ( seed_new < 0 )
    seed_new = seed_new + 2147483647;
  end

  r = seed_new * 4.656612875E-10;
