function bit = bit_hi1_base_2 ( n )

%% BIT_HI1_BASE_2 returns the position of the high 1 bit base 2 in an integer.
%
%  Example:
%
%       N    Binary     BIT
%    ----    --------  ----
%       0           0     0
%       1           1     1
%       2          10     2
%       3          11     2 
%       4         100     3
%       5         101     3
%       6         110     3
%       7         111     3
%       8        1000     4
%       9        1001     4
%      10        1010     4
%      11        1011     4
%      12        1100     4
%      13        1101     4
%      14        1110     4
%      15        1111     4
%      16       10000     5
%      17       10001     5
%    1023  1111111111    10
%    1024 10000000000    11
%    1025 10000000001    11
%
%  Modified:
%
%    28 March 2003
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, integer N, the integer to be measured.
%    N should be nonnegative.  If N is nonpositive, BIT will always be 0.
%
%    Output, integer BIT, the number of bits base 2.
%
  i = floor ( n );
  bit = 0;

  while ( 1 )

    if ( i <= 0 )
      break;
    end

    bit = bit + 1;
    i = i / 2;

  end
