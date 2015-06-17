function k = exor ( i, j )

%% EXOR calculates the exclusive OR of two integers.
%
%  Modified:
%
%    31 March 2003
%
%  Author:
%
%   John Burkardt
%
%  Reference:
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
%    Input, integer I, J, two values whose exclusive OR is needed.
%
%    Output, integer K, the exclusive OR of I and J.
%
  k = 0;
  l = 1;
%
  i = floor ( i );
  j = floor ( j );

  while ( i ~= 0 | j ~= 0 )
%
%  Check the current right-hand bits of I and J.
%  If they differ, set the appropriate bit of K.
%
    i2 = floor ( i / 2 );
    j2 = floor ( j / 2 );

    if ( ...
      ( ( i == 2 * i2 ) & ( j ~= 2 * j2 ) ) | ...
      ( ( i ~= 2 * i2 ) & ( j == 2 * j2 ) ) )
      k = k + l;
    end

    i = i2;
    j = j2;
    l = 2 * l;

  end
