function seed = get_seed ( dummy )

%% GET_SEED returns a random seed for the random number generator.
%
%  Modified:
%
%    04 March 2003
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, integer dummy, a dummy input value, since MATLAB
%    will not allow functions with no arguments.
%
%    Output, integer seed, a random seed value.
%
  I_MAX = 2147483647;

  time_array = clock;

  hour = time_array(4);
  minute = time_array(5);
  second = time_array(6);

  seed = second + 60 * ( minute + 60 * hour );
%
%  We want values in [1,43200], not [0,43199].
%
  seed = seed + 1;
%
%  Remap SEED from [1,43200] to [1,I_MAX].
%
  seed = I_MAX * ( seed  / ( 60.0 * 60.0 * 12.0 ) );

  seed = floor ( seed );
%
%  Never use a seed of 0.
%
  if ( seed == 0 )
    seed = 1;
  end
