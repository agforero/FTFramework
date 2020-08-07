program main

!*****************************************************************************80
!
!! MAIN is the main program for PARTIAL_DIGEST_TEST.
!
!  Discussion:
!
!    PARTIAL_DIGEST_TEST tests the PARTIAL_DIGEST library.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 January 2018
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PARTIAL_DIGEST_TEST:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the PARTIAL_DIGEST library.'

  call find_distances_test ( )
  call i4_uniform_ab_test ( )
  call i4vec_max_last_test ( )
  call i4vec_print_test ( )
  call partial_digest_recur_test01 ( )
  call partial_digest_recur_test02 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PARTIAL_DIGEST_TEST:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine find_distances ( l_length, l, x_length, x, y, success )

!*****************************************************************************80
!
!! FIND_DISTANCES determines if the "free" distances include every ||X(I)-Y||.
!
!  Discussion:
!
!    This routine is given a candidate point Y, a set of placed points
!    X(1:X_LENGTH), and a list of unused or "free" distances in
!    L(1:L_LENGTH).  The routine seeks to find in L a copy of the
!    distance from Y to each X.
!
!    If so, then the L array is reordered so that entries
!    L(L_LENGTH-X_LENGTH+1:L_LENGTH) contain theses distances.
!
!    In other words, Y can be added into X, and L_LENGTH reduced to
!    L_LENGTH-X_LENGTH.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 November 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Pavel Pevzner,
!    Computational Molecular Biology,
!    MIT Press, 2000,
!    ISBN: 0-262-16197-4,
!    LC: QH506.P47.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) L_LENGTH, the length of the array.
!
!    Input/output, integer ( kind = 4 ) L(L_LENGTH), the array.  On output,
!    some entries have been shuffled.  In particular, if SUCCESS is TRUE,
!    the entries L(L_LENGTH-X_LENGTH+1:L_LENGTH) contain the distances
!    of X(1:X_LENGTH) to Y.
!
!    Input, integer ( kind = 4 ) X_LENGTH, the number of entries in X.
!
!    Input, integer ( kind = 4 ) X(X_LENGTH), the number of points
!    already accepted.
!
!    Input, integer ( kind = 4 ) Y, a new point that we are considering.
!
!    Output, logical SUCCESS, is TRUE if the entries of L included
!    the values of the distance of Y to each entry of X.
!
  implicit none

  integer ( kind = 4 ) l_length
  integer ( kind = 4 ) x_length

  integer ( kind = 4 ) d
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) l(l_length)
  integer ( kind = 4 ) l2_length
  logical success
  integer ( kind = 4 ) x(x_length)
  integer ( kind = 4 ) y

  l2_length = l_length

  do i = 1, x_length

    d = abs ( x(i) - y )

    success = .false.

    do j = 1, l2_length

      if ( l(j) == d ) then
        l(j) = l(l2_length)
        l(l2_length) = d
        l2_length = l2_length - 1
        success = .true.
        exit
      end if

    end do

    if ( .not. success ) then
      return
    end if

  end do

  success = .true.

  return
end
subroutine find_distances_test ( ) 

!*****************************************************************************80
!
!! FIND_DISTANCES_TEST tests FIND_DISTANCES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 January 2018
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  integer ( kind = 4 ) i4vec_max_last
  integer ( kind = 4 ) l(n*(n-1)/2)
  integer ( kind = 4 ) l_length
  integer ( kind = 4 ) l_max
  logical success
  integer ( kind = 4 ) x(n)
  integer ( kind = 4 ) x_length
  integer ( kind = 4 ) y

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'FIND_DISTANCES_TEST:'
  write ( *, '(a)' ) '  FIND_DISTANCES takes a candidate location Y'
  write ( *, '(a)' ) '  and determines whether its distance to each point'
  write ( *, '(a)' ) '  in the X array is listed in the L array.'

  l_length = n * ( n - 1 ) / 2
  l = (/  13, 15, 38, 90, 2, 25, 77, 23, 75, 52 /) 
  call i4vec_print ( l_length, l, '  Initial L array:' )

  l_max = i4vec_max_last ( l_length, l )
  l_length = l_length - 1

  x(1) = 0
  x(2) = l_max
  x_length = 2
!
!  Solution is X = (/ 0, 13, 15, 38, 90 /) or (/ 0, 52, 75, 77, 90 /)
!  So Y = 13, 15, 38, 52, 75 or 77 will be acceptable.
!
  y = i4vec_max_last ( l_length, l )
  call find_distances ( l_length, l, x_length, x, y, success )

  write ( *, '(a)' ) ''
  write ( *, '(a,i4)' ) '  Consider Y = ', y
  write ( *, '(a)' ) ''
  if ( success ) then
    write ( *, '(a)' ) '  This Y is acceptable'
    l_length = l_length - x_length
    x_length = x_length + 1
    x(x_length) = y
    call i4vec_print ( x_length, x, '  New X array:' )
    call i4vec_print ( l_length, l, '  New L array:' )
  else
    write ( *, '(a)' ) '  This Y is not acceptable'
  end if

  y = 35
  call find_distances ( l_length, l, x_length, x, y, success )

  write ( *, '(a)' ) ''
  write ( *, '(a,i4)' ) '  Consider Y = ', y
  write ( *, '(a)' ) ''
  if ( success ) then
    write ( *, '(a)' ) '  This Y is acceptable'
    l_length = l_length - x_length
    x_length = x_length + 1
    x(x_length) = y
    call i4vec_print ( x_length, x, '  New X array:' )
    call i4vec_print ( l_length, l, '  New L array:' )
  else
    write ( *, '(a)' ) '  This Y is not acceptable'
  end if

  return
end
function i4_uniform_ab ( a, b, seed )

!*****************************************************************************80
!
!! I4_UNIFORM_AB returns a scaled pseudorandom I4 between A and B.
!
!  Discussion:
!
!    An I4 is an integer ( kind = 4 ) value.
!
!    The pseudorandom number will be scaled to be uniformly distributed
!    between A and B.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 October 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Second Edition,
!    Springer, 1987,
!    ISBN: 0387964673,
!    LC: QA76.9.C65.B73.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, December 1986, pages 362-376.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley, 1998,
!    ISBN: 0471134031,
!    LC: T57.62.H37.
!
!    Peter Lewis, Allen Goodman, James Miller,
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, Number 2, 1969, pages 136-143.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) A, B, the limits of the interval.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, integer ( kind = 4 ) I4_UNIFORM_AB, a number between A and B.
!
  implicit none

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  integer ( kind = 4 ), parameter :: i4_huge = 2147483647
  integer ( kind = 4 ) i4_uniform_ab
  integer ( kind = 4 ) k
  real ( kind = 4 ) r
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) value

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4_UNIFORM_AB - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop 1
  end if

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + i4_huge
  end if

  r = real ( seed, kind = 4 ) * 4.656612875E-10
!
!  Scale R to lie between A-0.5 and B+0.5.
!
  r = ( 1.0E+00 - r ) * ( real ( min ( a, b ), kind = 4 ) - 0.5E+00 ) & 
    +             r   * ( real ( max ( a, b ), kind = 4 ) + 0.5E+00 )
!
!  Use rounding to convert R to an integer between A and B.
!
  value = nint ( r, kind = 4 )

  value = max ( value, min ( a, b ) )
  value = min ( value, max ( a, b ) )

  i4_uniform_ab = value

  return
end
subroutine i4_uniform_ab_test ( )

!*****************************************************************************80
!
!! I4_UNIFORM_AB_TEST tests I4_UNIFORM_AB.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 October 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: a = -100
  integer ( kind = 4 ), parameter :: b = 200
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_uniform_ab
  integer ( kind = 4 ) j
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'I4_UNIFORM_AB_TEST'
  write ( *, '(a)' ) '  I4_UNIFORM_AB computes pseudorandom values '
  write ( *, '(a)' ) '  in an interval [A,B].'

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  The lower endpoint A = ', a
  write ( *, '(a,i12)' ) '  The upper endpoint B = ', b
  write ( *, '(a,i12)' ) '  The initial seed is ', seed
  write ( *, '(a)' ) ' '

  do i = 1, 20

    j = i4_uniform_ab ( a, b, seed )

    write ( *, '(2x,i8,2x,i8)' ) i, j

  end do

  return
end
function i4vec_max_last ( l_length, l )

!*****************************************************************************80
!
!! I4VEC_MAX_LAST moves the maximum entry of an I4VEC to the last position.
!
!  Discussion:
!
!    This routine finds the largest entry in an array and moves
!    it to the end of the array.
!
!    If we ignore this last array entry, then the effect is the same
!    as "deleting" the maximum entry from the array.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 January 2018
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Pavel Pevzner,
!    Computational Molecular Biology,
!    MIT Press, 2000,
!    ISBN: 0-262-16197-4,
!    LC: QH506.P47.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) L_LENGTH, the length of the array.
!
!    Input, integer ( kind = 4 ) L(L_LENGTH), the array.  On output,
!    the maximum entry has been "deleted", that is, the array has
!    been shifted so that this entry occurs at the end.
!
!    Output, integer ( kind = 4 ) I4VEC_MAX_LAST, the maximum entry in the
!    input array.
!
  implicit none

  integer ( kind = 4 ) l_length

  integer ( kind = 4 ) i4vec_max_last
  integer ( kind = 4 ) l(l_length)
  integer ( kind = 4 ) max_index(1)
  integer ( kind = 4 ) value

  max_index = maxloc ( l(1:l_length) )
  value = l(max_index(1))
  l(max_index(1):l_length-1) = l(max_index(1)+1:l_length)
  l(l_length) = value

  i4vec_max_last = value

  return
end
subroutine i4vec_max_last_test ( )

!*****************************************************************************80
!
!! I4VEC_MAX_LAST_TEST tests I4VEC_MAX_LAST.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 January 2018
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_uniform_ab
  integer ( kind = 4 ) i4vec_max_last
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) x(n)
  integer ( kind = 4 ) x_max

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'I4VEC_MAX_LAST_TEST'
  write ( *, '(a)' ) '  I4VEC_MAX_LAST identifies the largest element in an'
  write ( *, '(a)' ) '  I4VEC, and moves it to the final entry.'

  seed = 123456789

  do i = 1, n
    x(i) = i4_uniform_ab ( 1, 30, seed )
  end do

  call i4vec_print ( n, x, '  Input vector:' )

  x_max = i4vec_max_last ( n, x )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Maximum:                  ', x_max

  call i4vec_print ( n, x, '  Output vector:' )

  return
end
subroutine i4vec_print ( n, a, title )

!*****************************************************************************80
!
!! I4VEC_PRINT prints an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of integer values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components of the vector.
!
!    Input, integer ( kind = 4 ) A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i8,2x,i12)' ) i, a(i)
  end do

  return
end
subroutine i4vec_print_test ( )

!*****************************************************************************80
!
!! I4VEC_PRINT_TEST tests I4VEC_PRINT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 November 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 4

  integer ( kind = 4 ), dimension ( n ) :: a = (/ &
    91, 92, 93, 94 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'I4VEC_PRINT_TEST'
  write ( *, '(a)' ) '  I4VEC_PRINT prints an I4VEC'

  call i4vec_print ( n, a, '  The I4VEC:' )

  return
end
subroutine partial_digest_recur ( n, l )

!*****************************************************************************80
!
!! PARTIAL_DIGEST_RECUR uses recursion on the partial digest problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 November 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Pavel Pevzner,
!    Computational Molecular Biology,
!    MIT Press, 2000,
!    ISBN: 0-262-16197-4,
!    LC: QH506.P47.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of nodes.
!
!    Input, integer ( kind = 4 ) L((N*(N-1))/2), the distances between all pairs
!    of distinct nodes.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i4vec_max_last
  integer ( kind = 4 ) l((n*(n-1))/2)
  integer ( kind = 4 ) l_length
  integer ( kind = 4 ) width
  integer ( kind = 4 ) x(n)
  integer ( kind = 4 ) x_length
!
!  How long is L?
!
  l_length = ( n * ( n - 1 ) ) / 2
!
!  Find WIDTH, the largest element of L, and move it to the last position.
!
  width = i4vec_max_last ( l_length, l )
!
!  Think of L as being 1 entry shorter.
!
  l_length = l_length - 1
!
!  Using WIDTH, set the first two entries of X.
!
  x(1) = 0
  x(2) = width
  x_length = 2
!
!  Begin recursive operation.
!
  call place ( l_length, l, x_length, x )

  return
end
subroutine partial_digest_recur_test01 ( )

!*****************************************************************************80
!
!! PARTIAL_DIGEST_RECUR_TEST01 tests PARTIAL_DIGEST_RECUR.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 November 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5
  integer ( kind = 4 ), parameter :: nn2 = ( n * ( n - 1 ) ) / 2
!
!  Set the distance array.
!
  integer ( kind = 4 ), dimension ( ((n-1)*n)/2 ) :: dist = (/ &
    2, 2, 3, 3, 4, 5, 6, 7, 8, 10 /)
  integer ( kind = 4 ) i

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'PARTIAL_DIGEST_RECUR_TEST01'
  write ( *, '(a)' ) '  PARTIAL_DIGEST_RECUR generates solutions to the partial'
  write ( *, '(a)' ) '  digest problem, using recursion'

  write ( *, '(a)' ) ''
  write ( *, '(a,i8)' ) '  The number of objects to place is N = ', n
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  The original placement was 0,3,6,8,10.'
  write ( *, '(a)' ) '  These placements generate the following distances:'

  call i4vec_print ( nn2, dist, '  Distance array:' )

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  PARTIAL_DIGEST_RECUR may recover the original placements'
  write ( *, '(a)' ) '  from the pairwise distances.  It may also find other'
  write ( *, '(a)' ) '  placements that have the same distance array.'

  call partial_digest_recur ( n, dist )

  return
end
subroutine partial_digest_recur_test02 ( )

!*****************************************************************************80
!
!! PARTIAL_DIGEST_RECUR_TEST02 considers tests from a library.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 January 2018
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), allocatable :: d(:)
  integer ( kind = 4 ) dmax
  integer ( kind = 4 ) k
  integer ( kind = 4 ), allocatable :: locate ( : )
  integer ( kind = 4 ) seed

  call timestamp ( )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'PARTIAL_DIGEST_RECUR_TEST02:'
  write ( *, '(a)' ) '  PARTIAL_DIGEST_RECUR generates solutions to the partial'
  write ( *, '(a)' ) '  digest problem, using recursion'
  write ( *, '(a)' ) '  TEST_PARTIAL_DIGEST creates test problems for the'
  write ( *, '(a)' ) '  partial digest problem.'
!
!  Request a sample problem.
!
  k = 6
  dmax = 20
  seed = 123456789

  write ( *, '(a)' ) ''
  write ( *, '(a,i4)' ) '  Number of nodes = ', k
  write ( *, '(a,i4)' ) '  Maximum distance = ', dmax

  allocate ( locate(1:k) )
  allocate ( d(1:k*(k-1)/2) )
  
  call test_partial_digest ( k, dmax, seed, locate, d )
!
!  Sort the data.
!
  call i4vec_sort_heap_a ( k, locate )
  call i4vec_sort_heap_a ( k*(k-1)/2, d )
!
!  Print the data.
!
  call i4vec_print ( k, locate, '  Locations:' )
  call i4vec_print ( k * ( k - 1 ) / 2, d, '  Distances:' )
!
!  Solve the problem.
!
  call partial_digest_recur ( k, d )
!
!  Free memory.
!
  deallocate ( d )
  deallocate ( locate )
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'PARTIAL_DIGEST_RECUR_TEST02:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ''
  call timestamp ( )

  return
end
recursive subroutine place ( l_length, l, x_length, x )

!*****************************************************************************80
!
!! PLACE tries to place the next point for the partial digest problem.
!
!  Discussion:
!
!    Note that this is a recursive subroutine.  A solution to the
!    partial digest problem is sought by calling this routine repeatedly.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 November 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Pavel Pevzner,
!    Computational Molecular Biology,
!    MIT Press, 2000,
!    ISBN: 0-262-16197-4,
!    LC: QH506.P47.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) L_LENGTH, the number of entries in L.
!
!    Input/output, integer ( kind = 4 ) L(L_LENGTH), the array of distances.

!    Input/output, integer ( kind = 4 ) X_LENGTH, the number of entries in X.
!
!    Input/output, integer ( kind = 4 ) X(X_LENGTH), the current partial solution.
!
  implicit none

  integer ( kind = 4 ) l_length
  integer ( kind = 4 ) x_length

  integer ( kind = 4 ) i4vec_max_last
  integer ( kind = 4 ) l(l_length)
  integer ( kind = 4 ) l_length2
  logical success
  integer ( kind = 4 ) x(x_length)
  integer ( kind = 4 ) y
!
!  Are we done?
!
  if ( l_length <= 0 ) then
    call i4vec_print ( x_length, x, '  Solution:' )
    return
  end if
!
!  Find the maximum remaining distance.
!
  y = i4vec_max_last ( l_length, l )
!
!  We can add a point at Y if L contains all the distances from Y to
!  the current X's.
!
  call find_distances ( l_length, l, x_length, x, y, success )

  if ( success ) then
    l_length2 = l_length - x_length
    x_length = x_length + 1
    x(x_length) = y
    call place ( l_length2, l, x_length, x )
    x_length = x_length - 1
  end if
!
!  We must also consider the case where Y represents the distance
!  to X(2), not X(1).
!
  y = x(2) - y

  call find_distances ( l_length, l, x_length, x, y, success )

  if ( success ) then
    l_length2 = l_length - x_length
    x_length = x_length + 1
    x(x_length) = y
    call place ( l_length2, l, x_length, x )
    x_length = x_length - 1
  end if

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
