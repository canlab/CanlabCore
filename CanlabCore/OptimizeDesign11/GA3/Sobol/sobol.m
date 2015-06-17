function [ quasi, seed_new ] = sobol ( dim_num, seed )

%% SOBOL generates a new quasirandom Sobol vector with each call.
%
%  Discussion:
%
%    The routine adapts the ideas of Antonov and Saleev.
%
%  Modified:
%
%    30 March 2003
%
%  Reference:
%
%    Antonov and Saleev,
%    USSR Computational Mathematics and Mathematical Physics,
%    Volume 19, 1980, pages 252 - 256.
%
%    Paul Bratley and Bennett Fox,
%    Algorithm 659:
%    Implementing Sobol's Quasirandom Sequence Generator,
%    ACM Transactions on Mathematical Software,
%    Volume 14, Number 1, pages 88-100, 1988.
%
%    Bennett Fox,
%    Algorithm 647:
%    Implementation and Relative Efficiency of Quasirandom 
%    Sequence Generators,
%    ACM Transactions on Mathematical Software,
%    Volume 12, Number 4, pages 362-376, 1986.
%
%    I Sobol,
%    USSR Computational Mathematics and Mathematical Physics,
%    Volume 16, pages 236-242, 1977.
%
%    I Sobol and Levitan, 
%    The Production of Points Uniformly Distributed in a Multidimensional 
%    Cube (in Russian),
%    Preprint IPM Akad. Nauk SSSR, 
%    Number 40, Moscow 1976.
%
%  Parameters:
%
%    Input, integer DIM_NUM, the number of spatial dimensions.
%    DIM_NUM must satisfy 2 <= DIM_NUM <= 40.
%
%    Input/output, integer SEED, the "seed" for the sequence.
%    This is essentially the index in the sequence of the quasirandom
%    value to be generated.  On output, SEED has been set to the
%    appropriate next value, usually simply SEED+1.
%    If SEED is less than 0 on input, it is treated as though it were 0.
%    An input value of 0 requests the first (0-th) element of the sequence.
%
%    Output, real QUASI(DIM_NUM), the next quasirandom vector.
%
  global SOBOL_lastq;
  global SOBOL_seed;

  dim_max = 40;
%
%  Initialize (part of) V.
%
    v(1:40,1:30) = zeros(40,30);

    v(1:40,1) = [ ...
      1, 1, 1, 1, 1, 1, 1, 1, 1, 1, ...
      1, 1, 1, 1, 1, 1, 1, 1, 1, 1, ...
      1, 1, 1, 1, 1, 1, 1, 1, 1, 1, ...
      1, 1, 1, 1, 1, 1, 1, 1, 1, 1 ]';

    v(3:40,2) = [ ...
            1, 3, 1, 3, 1, 3, 3, 1, ...
      3, 1, 3, 1, 3, 1, 1, 3, 1, 3, ...
      1, 3, 1, 3, 3, 1, 3, 1, 3, 1, ...
      3, 1, 1, 3, 1, 3, 1, 3, 1, 3 ]';

    v(4:40,3) = [ ...
               7, 5, 1, 3, 3, 7, 5, ...
      5, 7, 7, 1, 3, 3, 7, 5, 1, 1, ...
      5, 3, 3, 1, 7, 5, 1, 3, 3, 7, ...
      5, 1, 1, 5, 7, 7, 5, 1, 3, 3 ]';

    v(6:40,4) = [ ...
                     1, 7, 9,13,11, ...
      1, 3, 7, 9, 5,13,13,11, 3,15, ...
      5, 3,15, 7, 9,13, 9, 1,11, 7, ...
      5,15, 1,15,11, 5, 3, 1, 7, 9 ]';
  
    v(8:40,5) = [ ...
                           9, 3,27, ...
     15,29,21,23,19,11,25, 7,13,17, ...
      1,25,29, 3,31,11, 5,23,27,19, ...
     21, 5, 1,17,13, 7,15, 9,31, 9 ]';

    v(14:40,6) = [ ...
              37,33, 7, 5,11,39,63, ...
     27,17,15,23,29, 3,21,13,31,25, ...
      9,49,33,19,29,11,19,27,15,25 ]';

    v(20:40,7) = [ ...
                                         13, ...
     33,115, 41, 79, 17, 29,119, 75, 73,105, ...
      7, 59, 65, 21,  3,113, 61, 89, 45,107 ]';

    v(38:40,8) = [ ...
                                  7, 23, 39 ]';
%
%  Set POLY.
%
    poly(1:40)= [ ...
        1,   3,   7,  11,  13,  19,  25,  37,  59,  47, ...
       61,  55,  41,  67,  97,  91, 109, 103, 115, 131, ...
      193, 137, 145, 143, 241, 157, 185, 167, 229, 171, ...
      213, 191, 253, 203, 211, 239, 247, 285, 369, 299 ];
%
%  Check parameters.
%
    if ( dim_num < 2 | dim_max < dim_num )
      fprintf ( 1, '\n' );
      fprintf ( 1, 'SOBOL - Fatal error!\n' );
      fprintf ( 1, '  The spatial dimension DIM_NUM should satisfy:\n' );
      fprintf ( 1, '    2 <= DIM_NUM <= %d\n', dim_max );
      fprintf ( 1, '  But this input value is DIM_NUM = %d\n', dim_num );
      return
    end

    atmost = 2^30 - 1;
%
%  Find the number of bits in ATMOST.
%
    maxcol = bit_hi1_base_2 ( atmost );
%
%  Initialize row 1 of V.
%
    v(1,1:maxcol) = 1;
%
%  Initialize the remaining rows of V.
%
    for ( i = 2 : dim_num )
%
%  The bit pattern of the integer POLY(I) gives the form
%  of polynomial I.
%
%  Find the degree of polynomial I from binary encoding.
%
      j = poly(i);
      m = 0;

      while ( 1 )

        j = floor ( j / 2 );

        if ( j <= 0 )
          break;
        end

        m = m + 1;

      end
%
%  We expand this bit pattern to separate components of the logical array INCLUD.
%
      j = poly(i);
      for ( k = m : -1 : 1 )
        j2 = floor ( j / 2 );
        includ(k) = ( j ~= 2 * j2 );
        j = j2;
      end
%
%  Calculate the remaining elements of row I as explained
%  in Bratley and Fox, section 2.
%
      for ( j = m + 1 : maxcol )

        newv = v(i,j-m);
        l = 1;

        for ( k = 1 : m )

          l = 2 * l;

          if ( includ(k) )
            newv = exor ( newv, l * v(i,j-k) );
          end

        end

        v(i,j) = newv;

      end

    end
%
%  Multiply columns of V by appropriate power of 2.
%
    l = 1;
    for ( j = maxcol-1 : -1 : 1 )
      l = 2 * l;
      v(1:dim_num,j) = v(1:dim_num,j) * l;
    end
%
%  RECIPD is 1/(common denominator of the elements in V).
%
    recipd = 1.0E+00 / ( 2 * l );

  seed = floor ( seed );

  if ( seed < 0 )
    seed = 0;
  end

  if ( seed == 0 )

    l = 1;
    SOBOL_lastq(1:dim_num) = 0;

  elseif ( seed == SOBOL_seed + 1 )
%
%  Find the position of the right-hand zero in SEED.
%
    l = bit_lo0_base_2 ( seed );

  elseif ( seed <= SOBOL_seed )

    SOBOL_seed = 0;
    l = 1;
    SOBOL_lastq(1:dim_num) = 0;

    for ( seed_temp = SOBOL_seed : seed-1 )

      l = bit_lo0_base_2 ( seed_temp );

      for ( i = 1 : dim_num )
        SOBOL_lastq(i) = exor ( SOBOL_lastq(i), v(i,l) );
      end

    end

    l = bit_lo0_base_2 ( seed );

  elseif ( SOBOL_seed+1 < seed )

    for ( seed_temp = SOBOL_seed+1 : seed-1 )

      l = bit_lo0_base_2 ( seed_temp );

      for ( i = 1 : dim_num )
        SOBOL_lastq(i) = exor ( SOBOL_lastq(i), v(i,l) );
      end

    end

    l = bit_lo0_base_2 ( seed );

  end
%
%  Check that the user is not calling too many times!
%
  if ( maxcol < l )
    fprintf ( 1, '\n' );
    fprintf ( 1, 'SOBOL - Fatal error!\n' );
    fprintf ( 1, '  Too many calls!\n' );
    fprintf ( 1, '  MAXCOL = %d\n', maxcol );
    fprintf ( 1, '  L =      %d\n', l );
    return
  end
%
%  Calculate the new components of QUASI.
%
  for ( i = 1 : dim_num )

    quasi(i) = SOBOL_lastq(i) * recipd;

    SOBOL_lastq(i) = exor ( SOBOL_lastq(i), v(i,l) );

  end

  SOBOL_seed = seed;
  seed_new = seed + 1;
