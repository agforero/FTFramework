program main

!*****************************************************************************80
!
!! D_SAMPLE calls SUPERLU to solve a double precision real system.
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
  implicit none

  integer ( kind = 4 ), parameter :: n = 5
  integer ( kind = 4 ), parameter :: m = n
  integer ( kind = 4 ), parameter :: ncc = 12

  real ( kind = 8 ), dimension ( ncc ) :: acc = (/ &
    19.0, 12.0, 12.0, &
    21.0, 12.0, 12.0, &
    21.0, 16.0, &
    21.0,  5.0, &
    21.0, 18.0 /)
  real ( kind = 8 ) b(n)
  real ( kind = 8 ) b2(n)
  integer ( kind = 4 ), dimension ( n + 1 ) :: ccc = (/ &
    1, 4, 7, 9, 11, 13 /)
  integer ( kind = 4 ) factors(8)
  integer ( kind = 4 ) i
  integer ( kind = 4 ), dimension ( ncc ) :: icc = (/ &
    1, 2, 5, &
    2, 3, 5, &
    1, 3, &
    1, 4, &
    4, 5 /)
  integer ( kind = 4 ) info
  integer ( kind = 4 ) iopt
  integer ( kind = 4 ) ldb
  integer ( kind = 4 ) nrhs

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'D_SAMPLE'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  DGSSV factors and solves a linear system'
  write ( *, '(a)' ) '  using double precision real arithmetic.'

  write ( *, '(a)' ) ''
  write ( *, '(a,i6)' ) '  Matrix order N = ', n
  write ( *, '(a,i6)' ) '  Matrix nonzeros NCC = ', ncc
!
!  Print the matrix.
!
  call cc_print ( m, n, ncc, icc, ccc, acc, '  CC matrix:' )

  nrhs = 1
  ldb = n
  do i = 1, n
    b(i) = 1.0D+00
  end do
!
!  Factor the matrix.
!
  iopt = 1
  call c_fortran_dgssv ( iopt, n, ncc, nrhs, acc, icc, &
    ccc, b, ldb, factors, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'D_SAMPLE - Fatal error!'
    write ( *, '(a)' ) '  Factorization failed'
    write ( *, '(a,i4)' ) '  INFO = ', info
    stop 1
  end if

  write ( *, '(a)' ) '  Factorization succeeded.'
!
!  Solve the factored system.
!
  iopt = 2
  call c_fortran_dgssv ( iopt, n, ncc, nrhs, acc, icc, &
    ccc, b, ldb, factors, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'D_SAMPLE - Fatal error!'
    write ( *, '(a)' ) '  Backsolve failed'
    write ( *, '(a,i4)' ) '  INFO = ', info
    stop 1
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Computed solution:'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(g14.6)' ) b(i)
  end do
!
!  B now contains the solution X.
!  Set B2 = A * X.
!
  call cc_mv ( m, n, ncc, icc, ccc, acc, b, b2 )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Product A*X:'
  write ( *, '(a)' ) ''
  do i = 1, n
    write ( *, '(g14.6)' ) b2(i)
  end do
!
!  Free memory.
!
  iopt = 3
  call c_fortran_dgssv ( iopt, n, ncc, nrhs, acc, icc, &
    ccc, b, ldb, factors, info )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'D_SAMPLE:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ''
  call timestamp ( )

  stop 0
end
subroutine cc_mv ( m, n, ncc, icc, ccc, acc, x, b )

!*****************************************************************************80
!
!! CC_MV multiplies a CC matrix by a vector
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
!  Reference:
!
!    Iain Duff, Roger Grimes, John Lewis,
!    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
!    October 1992
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows.
!
!    Input, integer ( kind = 4 ) N, the number of columns.
!
!    Input, integer ( kind = 4 ) NCC, the number of CC values.
!
!    Input, integer ( kind = 4 ) ICC(NCC), the CC rows.
!
!    Input, integer ( kind = 4 ) CCC(N+1), the compressed CC columns
!
!    Input, real ( kind = 8 ) ACC(NCC), the CC values.
!
!    Input, real ( kind = 8 ) X(N), the vector to be multiplied.
!
!    Output, real ( kind = 8 ) B(M), the product A*X.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) ncc

  real ( kind = 8 ) acc(ncc)
  real ( kind = 8 ) b(m)
  integer ( kind = 4 ) ccc(n+1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icc(ncc)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) x(n)

  b(1:m) = 0.0D+00

  do j = 1, n
    do k = ccc(j), ccc(j+1) - 1
      i = icc(k)
      b(i) = b(i) + acc(k) * x(j)
    end do
  end do

  return
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

  integer ( kind = 4 ) n
  integer ( kind = 4 ) ncc

  real ( kind = 8 ) acc(ncc)
  integer ( kind = 4 ) ccc(n+1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icc(ncc)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jnext
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(2x,i4,a,i4,a)' ) m, ' rows by ', n, ' columns.'
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
 
    write ( *, '(2x,i4,2x,i4,2x,i4,2x,g16.8)' ) k, i, j, acc(k)

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
