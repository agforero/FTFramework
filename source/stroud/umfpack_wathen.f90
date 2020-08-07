program main

!*****************************************************************************80
!
!  Purpose:
!
!    MAIN is the main program for UMFPACK_WATHEN.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 July 2014
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Timothy Davis,
!    UMFPACK User Guide,
!    Version 5.6.2, 25 April 2013
!    http://suitesparse.com
!
  implicit none
 
  real ( kind = 8 ), allocatable :: acc(:)
  real ( kind = 8 ), allocatable :: ast(:)
  real ( kind = 8 ), allocatable :: b(:)
  integer ( kind = 4 ), allocatable :: ccc(:)
  real ( kind = 8 ) control(20)
  integer ( kind = 4 ) filenum
  integer ( kind = 4 ) i
  integer ( kind = 4 ), allocatable :: icc(:)
  integer ( kind = 4 ), allocatable :: ist(:)
  integer ( kind = 4 ), allocatable :: jst(:)
  real ( kind = 8 ) info(90)
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) ncc
  integer ( kind = 4 ) nst
  integer ( kind = 8 ) numeric
  integer ( kind = 4 ) nx
  integer ( kind = 4 ) ny
  real ( kind = 8 ) r
  real ( kind = 8 ) r8vec_diff_norm
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) status
  integer ( kind = 8 ) symbolic
  integer ( kind = 4 ) sys
  real ( kind = 8 ), allocatable :: x1(:)
  real ( kind = 8 ), allocatable :: x2(:)

  call timestamp ( );
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'UMFPACK_WATHEN:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Use UMFPACK for the sparse linear system A*x=b.'
!
!  Set the default control parameters.
!
  call umf4def ( control )
!
!  Get the size of the ST matrix.
!
  nx = 3
  ny = 3
  call wathen_st_size ( nx, ny, nst )

  write ( *, '(a)' ) ''
  write ( *, '(a,i4)' ) '  Number of ST values = ', nst
!
!  Set the formal matrix size
!
  m = 3 * nx * ny + 2 * nx + 2 * ny + 1
  n = m
!
!  Set a random vector.
!
  seed = 123456789
  allocate ( x1(1:n) )
  call r8vec_uniform_01 ( n, seed, x1 )
!
!  Allocate space.
!
  allocate ( ast(1:nst) )
  allocate ( ist(1:nst) )
  allocate ( jst(1:nst) )
!
!  Create the ST matrix.
!
  seed = 123456789
  call wathen_st ( nx, ny, nst, seed, ist, jst, ast )
!
!  Compute B = AST * X1
!
  allocate ( b(1:n) )
  call st_mv ( m, n, nst, ist, jst, ast, x1, b )
!
!  Get the CC size.
!
  call st_to_cc_size ( nst, ist, jst, ncc )

  write ( *, '(a,i4)' ) '  Number of CC values = ', ncc
!
!  Create the CC indices.
!
  allocate ( icc(1:ncc) )
  allocate ( ccc(1:n+1) )
  allocate ( acc(1:ncc) )

  call st_to_cc_index ( nst, ist, jst, ncc, n, icc, ccc )
!
!  Create the CC values.
!
  call st_to_cc_values ( nst, ist, jst, ast, ncc, n, icc, ccc, acc )
!
!  Decrement the row and column indices.
!
  call i4vec_dec ( ncc, icc )
  call i4vec_dec ( n + 1, ccc )
!
!  Print the matrix.
!
  call cc_print ( m, n, ncc, icc, ccc, acc, '  The CC matrix:' )
!
!  From the matrix data, create the symbolic factorization information.
!
  call umf4sym ( n, n, ccc, icc, acc, symbolic, control, info )

  if ( info(1) < 0.0D+00 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'UMFPACK_WATHEN - Fatal error!'
    write ( *, '(a,g14.6)' ) '  UMF4SYM returns INFO(1) = ', info(1)
    stop 1
  end if
!
!  From the symbolic factorization information, carry out the numeric factorization.
!
  call umf4num ( ccc, icc, acc, symbolic, numeric, control, info )

  if ( info(1) < 0.0D+00 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'UMFPACK_WATHEN - Fatal error!'
    write ( *, '(a,g14.6)' ) '  UMF4NUM returns INFO(1) = ', info(1)
    stop 1
  end if
!
!  Free the memory associated with the symbolic factorization.
!
  call umf4fsym ( symbolic )
!
!  Solve the linear system.
!
  sys = 0
  allocate ( x2(1:n) )
  call umf4sol ( sys, x2, b, numeric, control, info )

  if ( info(1) < 0.0D+00 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'UMFPACK_WATHEN - Fatal error!'
    write ( *, '(a,g14.6)' ) '  UMF4SOL returns INFO(1) = ', info(1)
    stop 1
  end if
!
!  Free the memory associated with the numeric factorization.
!
  call umf4fnum ( numeric )
!
!  Print the error.
!
  r = r8vec_diff_norm ( n, x1, x2 )
  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6)' ) '  L2 error ||X1 - X2|| = ', r
!
!  Free memory.
!
  deallocate ( acc )
  deallocate ( ast )
  deallocate ( b )
  deallocate ( ccc )
  deallocate ( icc )
  deallocate ( ist )
  deallocate ( jst )
  deallocate ( x1 )
  deallocate ( x2 )
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'UMFPACK_WATHEN:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ''
  call timestamp ( )

  stop
end
subroutine cc_print ( m, n, ncc, icc, ccc, acc, title )

!*****************************************************************************80
!
!! CC_PRINT prints a sparse matrix in CC format.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 July 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in the matrix.
!
!    Input, integer ( kind = 4 ) N, the number of columns in the matrix.
!
!    Input, integer ( kind = 4 ) NCC, the number of CC elements.
!
!    Input, integer ( kind = 4 ) ICC(NCC), the CC rows.
!
!    Input, integer ( kind = 4 ) CCC(N+1), the compressed CC columns.
!
!    Input, real ( kind = 8 ) ACC(NCC), the CC values.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) ncc

  real ( kind = 8 ) acc(ncc)
  integer ( kind = 4 ) ccc(ncc)
  integer ( kind = 4 ) icc(ncc)
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  character ( len = * ) title

  call cc_print_some ( 1, m, 1, n, ncc, n, icc, ccc, acc, title )

  return
end
subroutine cc_print_some ( i_min, i_max, j_min, j_max, ncc, n, icc, ccc, acc, &
  title )

!*****************************************************************************80
!
!! CC_PRINT_SOME prints some of a sparse matrix in CC format.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 July 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I_MIN, IMAX, the first and last rows to print.
!
!    Input, integer ( kind = 4 ) J_MIN, J_MAX, the first and last columns 
!    to print.
!
!    Input, integer ( kind = 4 ) NCC, the number of CC elements.
!
!    Input, integer ( kind = 4 ) N, the number of columns.
!
!    Input, integer ( kind = 4 ) ICC(NCC), the CC rows.
!
!    Input, integer ( kind = 4 ) CCC(N+1), the compressed CC columns.
!
!    Input, real ( kind = 8 ) ACC(NCC), the CC values.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) ncc

  real ( kind = 8 ) acc(ncc)
  integer ( kind = 4 ) ccc(n+1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i_max
  integer ( kind = 4 ) i_min
  integer ( kind = 4 ) icc(ncc)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j_max
  integer ( kind = 4 ) j_min
  integer ( kind = 4 ) jnext
  integer ( kind = 4 ) k
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) '     #     I     J       A'
  write ( *, '(a)' ) '  ----  ----  ----  --------------'
  write ( *, '(a)' ) ' '

  j = 1
  jnext = ccc(2)

  do k = 1, ncc

    i = icc(k)
    do while ( jnext <= k )
      j = j + 1
      jnext = ccc(j+1)
    end do
 
    if ( i_min <= i .and. i <= i_max .and. &
         j_min <= j .and. j <= j_max ) then
      write ( *, '(2x,i4,2x,i4,2x,i4,2x,g16.8)' ) k, i, j, acc(k)
    end if
  end do

  return
end
subroutine i4vec_copy ( n, a1, a2 )

!*****************************************************************************80
!
!! I4VEC_COPY copies an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the length of the vectors.
!
!    Input, integer ( kind = 4 ) A1(N), the vector to be copied.
!
!    Output, integer ( kind = 4 ) A2(N), a copy of A1.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a1(n)
  integer ( kind = 4 ) a2(n)

  a2(1:n) = a1(1:n)

  return
end
subroutine i4vec_dec ( n, a )

!*****************************************************************************80
!
!! I4VEC_DEC decrements an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 July 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input/output, integer ( kind = 4 ) A(N), the vector to be decremented.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)

  a(1:n) = a(1:n) - 1

  return
end
subroutine i4vec2_compare ( n, a1, a2, i, j, isgn )

!*****************************************************************************80
!
!! I4VEC2_COMPARE compares entries of an I4VEC2.
!
!  Discussion:
!
!    An I4VEC2 is a pair of I4VEC's.
!
!    An I4VEC is a vector of I4's.
!
!    Entry K of an I4VEC2 is the pair of values located
!    at the K-th entries of the two I4VEC's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 October 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of data items.
!
!    Input, integer ( kind = 4 ) A1(N), A2(N), contain the two components
!    of each item.
!
!    Input, integer ( kind = 4 ) I, J, the items to be compared.
!
!    Output, integer ( kind = 4 ) ISGN, the results of the comparison:
!    -1, item I < item J,
!     0, item I = item J,
!    +1, item J < item I.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a1(n)
  integer ( kind = 4 ) a2(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j

  isgn = 0

       if ( a1(i) < a1(j) ) then

    isgn = -1

  else if ( a1(i) == a1(j) ) then

         if ( a2(i) < a2(j) ) then
      isgn = -1
    else if ( a2(i) < a2(j) ) then
      isgn = 0
    else if ( a2(j) < a2(i) ) then
      isgn = +1
    end if

  else if ( a1(j) < a1(i) ) then

    isgn = +1

  end if

  return
end
subroutine i4vec2_sort_a ( n, a1, a2 )

!*****************************************************************************80
!
!! I4VEC2_SORT_A ascending sorts a vector of pairs of integers.
!
!  Discussion:
!
!    An I4VEC2 is a pair of I4VEC's.
!
!    An I4VEC is a vector of I4's.
!
!    Entry K of an I4VEC2 is the pair of values located
!    at the K-th entries of the two I4VEC's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 September 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of items of data.
!
!    Input/output, integer ( kind = 4 ) A1(N), A2(N), the data to be sorted.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a1(n)
  integer ( kind = 4 ) a2(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j
  integer ( kind = 4 ) t

  if ( n <= 1 ) then
    return
  end if
!
!  Initialize.
!
  i = 0
  indx = 0
  isgn = 0
  j = 0
!
!  Call the external heap sorter.
!
  do

    call sort_heap_external ( n, indx, i, j, isgn )
!
!  Interchange the I and J objects.
!
    if ( 0 < indx ) then

      t     = a1(i)
      a1(i) = a1(j)
      a1(j) = t

      t     = a2(i)
      a2(i) = a2(j)
      a2(j) = t
!
!  Compare the I and J objects.
!
    else if ( indx < 0 ) then

      call i4vec2_compare ( n, a1, a2, i, j, isgn )

    else if ( indx == 0 ) then

      exit

    end if

  end do

  return
end
subroutine i4vec2_sorted_unique_count ( n, a1, a2, unique_num )

!*****************************************************************************80
!
!! I4VEC2_SORTED_UNIQUE_COUNT counts unique elements in a sorted I4VEC2.
!
!  Discussion:
!
!    An I4VEC2 is a pair of I4VEC's.
!
!    An I4VEC is a vector of I4's.
!
!    Entry K of an I4VEC2 is the pair of values located
!    at the K-th entries of the two I4VEC's.
!
!    Item I is stored as the pair A1(I), A2(I).
!
!    The items must have been sorted, or at least it must be the
!    case that equal items are stored in adjacent vector locations.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 July 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of items.
!
!    Input, integer ( kind = 4 ) A1(N), A2(N), the items.
!
!    Output, integer ( kind = 4 ) UNIQUE_NUM, the number of unique items.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a1(n)
  integer ( kind = 4 ) a2(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iu
  integer ( kind = 4 ) unique_num

  if ( n <= 0 ) then
    unique_num = 0
    return
  end if

  iu = 1
  unique_num = 1

  do i = 2, n

    if ( a1(i) /= a1(iu) .or. a2(i) /= a2(iu) ) then
      unique_num = unique_num + 1
      iu = i
    end if

  end do

  return
end
subroutine i4vec2_sorted_uniquely ( n1, a1, b1, n2, a2, b2 )

!*****************************************************************************80
!
!! I4VEC2_SORTED_UNIQUELY copies unique elements from a sorted I4VEC2.
!
!  Discussion:
!
!    An I4VEC2 is a pair of I4VEC's.
!
!    An I4VEC is a vector of I4's.
!
!    Entry K of an I4VEC2 is the pair of values located
!    at the K-th entries of the two I4VEC's.
!
!    Item I is stored as the pair A1(I), A2(I).
!
!    The items must have been sorted, or at least it must be the
!    case that equal items are stored in adjacent vector locations.
!
!    If the items were not sorted, then this routine will only
!    replace a string of equal values by a single representative.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 July 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N1, the number of items.
!
!    Input, integer ( kind = 4 ) A1(N1), B1(N1), the array of items.
!
!    Input, integer ( kind = 4 ) N2, the number of unique items.
!
!    Output, integer ( kind = 4 ) A2(N2), B2(N2), the array of unique items.
!
  implicit none

  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2

  integer ( kind = 4 ) a1(n1)
  integer ( kind = 4 ) a2(n2)
  integer ( kind = 4 ) b1(n1)
  integer ( kind = 4 ) b2(n2)
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2

  i1 = 1
  i2 = 1
  a2(i2) = a1(i1)
  b2(i2) = b1(i1)

  do i1 = 2, n1

    if ( a1(i1) /= a2(i2) .or. b1(i1) /= b2(i2) ) then

      i2 = i2 + 1

      a2(i2) = a1(i1)
      b2(i2) = b1(i1)

    end if

  end do

  return
end
function r8_uniform_01 ( seed )

!*****************************************************************************80
!
!! R8_UNIFORM_01 returns a unit pseudorandom R8.
!
!  Discussion:
!
!    An R8 is a real ( kind = 8 ) value.
!
!    For now, the input quantity SEED is an integer variable.
!
!    This routine implements the recursion
!
!      seed = 16807 * seed mod ( 2^31 - 1 )
!      r8_uniform_01 = seed / ( 2^31 - 1 )
!
!    The integer arithmetic never requires more than 32 bits,
!    including a sign bit.
!
!    If the initial seed is 12345, then the first three computations are
!
!      Input     Output      R8_UNIFORM_01
!      SEED      SEED
!
!         12345   207482415  0.096616
!     207482415  1790989824  0.833995
!    1790989824  2035175616  0.947702
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Springer Verlag, pages 201-202, 1983.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley Interscience, page 95, 1998.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!    Peter Lewis, Allen Goodman, James Miller
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which should
!    NOT be 0. On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R8_UNIFORM_01, a new pseudorandom variate,
!    strictly between 0 and 1.
!
  implicit none

  integer ( kind = 4 ), parameter :: i4_huge = 2147483647
  integer ( kind = 4 ) k
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_UNIFORM_01 - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop 1
  end if

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + i4_huge
  end if

  r8_uniform_01 = real ( seed, kind = 8 ) * 4.656612875D-10

  return
end
function r8vec_diff_norm ( n, a, b )

!*****************************************************************************80
!
!! R8VEC_DIFF_NORM returns the L2 norm of the difference of R8VEC's.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    The vector L2 norm is defined as:
!
!      R8VEC_NORM_L2 = sqrt ( sum ( 1 <= I <= N ) A(I)^2 ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 April 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in A.
!
!    Input, real ( kind = 8 ) A(N), B(N), the vectors
!
!    Output, real ( kind = 8 ) R8VEC_DIFF_NORM, the L2 norm of A - B.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) b(n)
  real ( kind = 8 ) r8vec_diff_norm

  r8vec_diff_norm = sqrt ( sum ( ( a(1:n) - b(1:n) )**2 ) )

  return
end
subroutine r8vec_uniform_01 ( n, seed, r )

!*****************************************************************************80
!
!! R8VEC_UNIFORM_01 returns a unit pseudorandom R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Springer Verlag, pages 201-202, 1983.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!    Peter Lewis, Allen Goodman, James Miller
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R(N), the vector of pseudorandom values.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  integer ( kind = 4 ) seed
  real ( kind = 8 ) r(n)

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8VEC_UNIFORM_01 - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop 1
  end if

  do i = 1, n

    k = seed / 127773

    seed = 16807 * ( seed - k * 127773 ) - k * 2836

    if ( seed < 0 ) then
      seed = seed + 2147483647
    end if

    r(i) = real ( seed, kind = 8 ) * 4.656612875D-10

  end do

  return
end
subroutine sort_heap_external ( n, indx, i, j, isgn )

!*****************************************************************************80
!
!! SORT_HEAP_EXTERNAL externally sorts a list of items into ascending order.
!
!  Discussion:
!
!    The actual list of data is not passed to the routine.  Hence this
!    routine may be used to sort integers, reals, numbers, names,
!    dates, shoe sizes, and so on.  After each call, the routine asks
!    the user to compare or interchange two items, until a special
!    return value signals that the sorting is completed.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 February 2004
!
!  Author:
!
!    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms for Computers and Calculators,
!    Academic Press, 1978,
!    ISBN: 0-12-519260-6,
!    LC: QA164.N54.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of items to be sorted.
!
!    Input/output, integer ( kind = 4 ) INDX, the main communication signal.
!    The user must set INDX to 0 before the first call.
!    Thereafter, the user should not change the value of INDX until
!    the sorting is done.
!    On return, if INDX is
!    *greater than 0,
!    ...interchange items I and J;
!    ...call again.
!    *less than 0,
!    ...compare items I and J;
!    ...set ISGN = -1 if I < J, ISGN = +1 if J < I;
!    ...call again.
!    * equal to 0, 
!    ...the sorting is done.
!
!    Output, integer ( kind = 4 ) I, J, the indices of two items.
!    On return with INDX positive, elements I and J should be interchanged.
!    On return with INDX negative, elements I and J should be compared, and
!    the result reported in ISGN on the next call.
!
!    Input, integer ( kind = 4 ) ISGN, results of comparison of elements
!    I and J. (Used only when the previous call returned INDX less than 0).
!    ISGN <= 0 means I is less than or equal to J;
!    0 <= ISGN means I is greater than or equal to J.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ), save :: i_save = 0
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j
  integer ( kind = 4 ), save :: j_save = 0
  integer ( kind = 4 ), save :: k = 0
  integer ( kind = 4 ), save :: k1 = 0
  integer ( kind = 4 ) n
  integer ( kind = 4 ), save :: n1 = 0
!
!  INDX = 0: This is the first call.
!
  if ( indx == 0 ) then

    i_save = 0
    j_save = 0
    k = n / 2
    k1 = k
    n1 = n
!
!  INDX < 0: The user is returning the results of a comparison.
!
  else if ( indx < 0 ) then

    if ( indx == -2 ) then

      if ( isgn < 0 ) then
        i_save = i_save + 1
      end if

      j_save = k1
      k1 = i_save
      indx = -1
      i = i_save
      j = j_save
      return

    end if

    if ( 0 < isgn ) then
      indx = 2
      i = i_save
      j = j_save
      return
    end if

    if ( k <= 1 ) then

      if ( n1 == 1 ) then
        i_save = 0
        j_save = 0
        indx = 0
      else
        i_save = n1
        n1 = n1 - 1
        j_save = 1
        indx = 1
      end if

      i = i_save
      j = j_save
      return

    end if

    k = k - 1
    k1 = k
!
!  0 < INDX, the user was asked to make an interchange.
!
  else if ( indx == 1 ) then

    k1 = k

  end if

  do

    i_save = 2 * k1

    if ( i_save == n1 ) then
      j_save = k1
      k1 = i_save
      indx = -1
      i = i_save
      j = j_save
      return
    else if ( i_save <= n1 ) then
      j_save = i_save + 1
      indx = -2
      i = i_save
      j = j_save
      return
    end if

    if ( k <= 1 ) then
      exit
    end if

    k = k - 1
    k1 = k

  end do

  if ( n1 == 1 ) then
    i_save = 0
    j_save = 0
    indx = 0
    i = i_save
    j = j_save
  else
    i_save = n1
    n1 = n1 - 1
    j_save = 1
    indx = 1
    i = i_save
    j = j_save
  end if

  return
end
subroutine st_mv ( m, n, nst, ist, jst, ast, x, b )

!*****************************************************************************80
!
!! ST_MV multiplies an R8SP matrix by an R8VEC.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 July 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns of 
!    the matrix.
!
!    Input, integer ( kind = 4 ) NST, the number of nonzero elements in
!    the matrix.
!
!    Input, integer ( kind = 4 ) IST(NST), JST(NST), the row and 
!    column indices of the nonzero elements.
!
!    Input, real ( kind = 8 ) AST(NST), the nonzero elements of the matrix.
!
!    Input, real ( kind = 8 ) X(N), the vector to be multiplied by A.
!
!    Output, real ( kind = 8 ) B(M), the product vector A*X.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nst

  real ( kind = 8 ) ast(nst)
  real ( kind = 8 ) b(m)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ist(nst)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jst(nst)
  integer ( kind = 4 ) k
  real ( kind = 8 ) x(n)

  b(1:m) = 0.0D+00

  do k = 1, nst
    i = ist(k)
    j = jst(k)
    b(i) = b(i) + ast(k) * x(j)
  end do

  return
end
subroutine st_to_cc_index ( nst, ist, jst, ncc, n, icc, ccc )

!*****************************************************************************80
!
!! ST_TO_CC_INDEX creates CC indices from ST data.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 July 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NST, the number of ST elements.
!
!    Input, integer ( kind = 4 ) IST(NST), JST(NST), the ST rows and columns.
!
!    Input, integer ( kind = 4 ) NCC, the number of CC elements.
!
!    Input, integer ( kind = 4 ) N, the number of columns in the matrix.
!
!    Output, integer ( kind = 4 ) ICC(NCC), the CC rows.
!
!    Output, integer ( kind = 4 ) CCC(N+1), the compressed CC columns.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) ncc
  integer ( kind = 4 ) nst

  integer ( kind = 4 ) ccc(n+1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icc(ncc)
  integer ( kind = 4 ) ist(nst)
  integer ( kind = 4 ) ist2(nst)
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  integer ( kind = 4 ) jcc(ncc)
  integer ( kind = 4 ) jst(nst)
  integer ( kind = 4 ) jst2(nst)
!
!  Make copies so the sorting doesn't confuse the user.
!
  call i4vec_copy ( nst, ist, ist2 )
  call i4vec_copy ( nst, jst, jst2 )
!
!  Sort the elements.
!
  call i4vec2_sort_a ( nst, jst2, ist2 )
!
!  Get the unique elements.
!
  call i4vec2_sorted_uniquely ( nst, jst2, ist2, ncc, jcc, icc )
!
!  Compress the column index.
!
  ccc(1) = 1
  jlo = 1
  do i = 1, ncc
    jhi = jcc(i)
    if ( jhi /= jlo ) then
      ccc(jlo+1:jhi) = i
      jlo = jhi
    end if
  end do
  jhi = n + 1
  ccc(jlo+1:jhi) = ncc + 1

  return
end
subroutine st_to_cc_size ( nst, ist, jst, ncc )

!*****************************************************************************80
!
!! ST_TO_CC_SIZE sizes CC indexes based on ST data.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 July 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NST, the number of ST elements.
!
!    Input, integer ( kind = 4 ) IST(NST), JST(NST), the ST rows and columns.
!
!    Output, integer ( kind = 4 ) NCC, the number of CC elements.
!
  implicit none

  integer ( kind = 4 ) nst

  integer ( kind = 4 ) ist(nst)
  integer ( kind = 4 ) ist2(nst)
  integer ( kind = 4 ) jst2(nst)
  integer ( kind = 4 ) jst(nst)
  integer ( kind = 4 ) ncc
!
!  Make copies so the sorting doesn't confuse the user.
!
  call i4vec_copy ( nst, ist, ist2 )
  call i4vec_copy ( nst, jst, jst2 )
!
!  Sort by column first, then row.
!
  call i4vec2_sort_a ( nst, jst2, ist2 )
!
!  Count the unique pairs.
!
  call i4vec2_sorted_unique_count ( nst, jst2, ist2, ncc )

  return
end
subroutine st_to_cc_values ( nst, ist, jst, ast, ncc, n, icc, ccc, acc )

!*****************************************************************************80
!
!! ST_TO_CC_VALUES creates CC values from ST data.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 July 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NST, the number of ST elements.
!
!    Input, integer ( kind = 4 ) IST(NST), JST(NST), the ST rows and columns.
!
!    Input, real ( kind = 8 ) AST(NST), the ST values.
!
!    Input, integer ( kind = 4 ) NCC, the number of CC elements.
!
!    Input, integer ( kind = 4 ) N, the number of columns.
!
!    Input, integer ( kind = 4 ) ICC(NCC), the CC rows.
!
!    Input, integer ( kind = 4 ) CCC(N+1), the CC compressed columns.
!
!    Output, real ( kind = 8 ) ACC(NCC), the CC values.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) ncc
  integer ( kind = 4 ) nst

  real ( kind = 8 ) ast(nst)
  real ( kind = 8 ) acc(ncc)
  integer ( kind = 4 ) ccc(n+1)
  integer ( kind = 4 ) chi
  integer ( kind = 4 ) clo
  logical fail
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icc(ncc)
  integer ( kind = 4 ) ist(nst)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jst(nst)
  integer ( kind = 4 ) kcc
  integer ( kind = 4 ) kst

  acc(1:ncc) = 0.0D+00

  do kst = 1, nst

    i = ist(kst)
    j = jst(kst)

    clo = ccc(j)
    chi = ccc(j+1)

    fail = .true.

    do kcc = clo, chi - 1
      if ( icc(kcc) == i ) then
        acc(kcc) = acc(kcc) + ast(kst)
        fail = .false.
        exit
      end if                  
    end do

    if ( fail ) then
      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'ST_TO_CC_VALUES - Fatal error!'
      write ( *, '(a)' ) '  ST entry cannot be located in CC array.'
      write ( *, '(a,i4)' ) '  ST index KST    = ', kst
      write ( *, '(a,i4)' ) '  ST row IST(KST) = ', ist(kst)
      write ( *, '(a,i4)' ) '  ST col JST(KST) = ', jst(kst)
      write ( *, '(a,g14.6)' ) '  ST val AST(KST) = ', ast(kst)
      stop 1
    end if

  end do

  return
end
subroutine timestamp ( )

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    31 May 2001   9:45:54.872 AM
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 May 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none

  character ( len = 8 ) ampm
  integer ( kind = 4 ) d
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  integer ( kind = 4 ) values(8)
  integer ( kind = 4 ) y

  call date_and_time ( values = values )

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  if ( h < 12 ) then
    ampm = 'AM'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
      ampm = 'Noon'
    else
      ampm = 'PM'
    end if
  else
    h = h - 12
    if ( h < 12 ) then
      ampm = 'PM'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'Midnight'
      else
        ampm = 'AM'
      end if
    end if
  end if

  write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
subroutine wathen_st ( nx, ny, nz_num, seed, row, col, a )

!*****************************************************************************80
!
!! WATHEN_ST: Wathen matrix stored in sparse triplet (ST) format.
!
!  Discussion:
!
!    When dealing with sparse matrices in MATLAB, it can be much more efficient
!    to work first with a triple of I, J, and X vectors, and only once
!    they are complete, convert to MATLAB's sparse format.
!
!    The Wathen matrix is a finite element matrix which is sparse.
!
!    The entries of the matrix depend in part on a physical quantity
!    related to density.  That density is here assigned random values between
!    0 and 100.
!
!    The matrix order N is determined by the input quantities NX and NY,
!    which would usually be the number of elements in the X and Y directions.
!
!    The value of N is
!
!      N = 3*NX*NY + 2*NX + 2*NY + 1,
!
!    The matrix is the consistent mass matrix for a regular NX by NY grid
!    of 8 node serendipity elements.
!
!    The local element numbering is
!
!      3--2--1
!      |     |
!      4     8
!      |     |
!      5--6--7
!
!    Here is an illustration for NX = 3, NY = 2:
!
!     23-24-25-26-27-28-29
!      |     |     |     |
!     19    20    21    22
!      |     |     |     |
!     12-13-14-15-16-17-18
!      |     |     |     |
!      8     9    10    11
!      |     |     |     |
!      1--2--3--4--5--6--7
!
!    For this example, the total number of nodes is, as expected,
!
!      N = 3 * 3 * 2 + 2 * 2 + 2 * 3 + 1 = 29
!
!    The matrix is symmetric positive definite for any positive values of the
!    density RHO(X,Y).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 July 2014
!
!  Author:
!
!    John Burkardt.
!
!  Reference:
!
!    Nicholas Higham,
!    Algorithm 694: A Collection of Test Matrices in MATLAB,
!    ACM Transactions on Mathematical Software,
!    Volume 17, Number 3, September 1991, pages 289-305.
!
!    Andrew Wathen,
!    Realistic eigenvalue bounds for the Galerkin mass matrix,
!    IMA Journal of Numerical Analysis,
!    Volume 7, Number 4, October 1987, pages 449-457.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NX, NY, values which determine the size of 
!    the matrix.
!
!    Input, integer ( kind = 4 ) NZ_NUM, the number of values used to 
!    describe the matrix.
!
!    Input/output, integer ( kind = 4 ) SEED, the random number seed.
!
!    Output, integer ( kind = 4 ) ROW(NZ_NUM), COL(NZ_NUM), the row and 
!    column indices of the nonzero entries.
!
!    Output, real ( kind = 8 ) A(NZ_NUM), the nonzero entries of the matrix.
!
  implicit none

  integer ( kind = 4 ) nx
  integer ( kind = 4 ) ny
  integer ( kind = 4 ) nz_num

  real ( kind = 8 ) a(nz_num)
  integer ( kind = 4 ) col(nz_num)
  real ( kind = 8 ), dimension ( 8, 8 ), save :: em =  reshape ( (/ &
     6.0, -6.0,  2.0, -8.0,  3.0, -8.0,  2.0, -6.0, &
    -6.0, 32.0, -6.0, 20.0, -8.0, 16.0, -8.0, 20.0, &
     2.0, -6.0,  6.0, -6.0,  2.0, -8.0,  3.0, -8.0, &
    -8.0, 20.0, -6.0, 32.0, -6.0, 20.0, -8.0, 16.0, &
     3.0, -8.0,  2.0, -6.0,  6.0, -6.0,  2.0, -8.0, &
    -8.0, 16.0, -8.0, 20.0, -6.0, 32.0, -6.0, 20.0, &
     2.0, -8.0,  3.0, -8.0,  2.0, -6.0,  6.0, -6.0, &
    -6.0, 20.0, -8.0, 16.0, -8.0, 20.0, -6.0, 32.0 /), &
    (/ 8, 8 /) )
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kcol
  integer ( kind = 4 ) krow
  integer ( kind = 4 ) node(8)
  real ( kind = 8 ) r8_uniform_01
  real ( kind = 8 ) rho
  integer ( kind = 4 ) row(nz_num)
  integer ( kind = 4 ) seed

  row(1:nz_num) = 0
  col(1:nz_num) = 0
  a(1:nz_num) = 0.0D+00
  
  k = 0

  do j = 1, ny
    do i = 1, nx

      node(1) = 3 * j * nx + 2 * j + 2 * i + 1
      node(2) = node(1) - 1
      node(3) = node(1) - 2
      node(4) = ( 3 * j - 1 ) * nx + 2 * j + i - 1
      node(5) = ( 3 * j - 3 ) * nx + 2 * j + 2 * i - 3
      node(6) = node(5) + 1
      node(7) = node(5) + 2
      node(8) = node(4) + 1

      rho = 100.0D+00 * r8_uniform_01 ( seed )

      do krow = 1, 8
        do kcol = 1, 8
          k = k + 1
          row(k) = node(krow)
          col(k) = node(kcol)
          a(k) = rho * em(krow,kcol)
        end do
      end do

    end do
  end do

  return
end
subroutine wathen_st_size ( nx, ny, nz_num )

!*****************************************************************************80
!
!! WATHEN_ST_SIZE: Size of Wathen matrix stored in sparse triplet format.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 June 2014
!
!  Author:
!
!    John Burkardt.
!
!  Reference:
!
!    Nicholas Higham,
!    Algorithm 694: A Collection of Test Matrices in MATLAB,
!    ACM Transactions on Mathematical Software,
!    Volume 17, Number 3, September 1991, pages 289-305.
!
!    Andrew Wathen,
!    Realistic eigenvalue bounds for the Galerkin mass matrix,
!    IMA Journal of Numerical Analysis,
!    Volume 7, Number 4, October 1987, pages 449-457.
!
!  Parameters:
!
!    Input, integer NX, NY, values which determine the size of the matrix.
!
!    Output, integer NZ_NUM, the number of items of data used to describe
!    the matrix.
!
  implicit none

  integer ( kind = 4 ) nx
  integer ( kind = 4 ) ny
  integer ( kind = 4 ) nz_num

  nz_num = nx * ny * 64

  return
end