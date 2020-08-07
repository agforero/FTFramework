subroutine i4_to_digits_binary ( i, n, c )

!*****************************************************************************80
!
!! I4_TO_DIGITS_BINARY produces the binary digits of an I4.
!
!  Discussion:
!
!    An I4 is an integer.
!
!  Example:
!
!     I    N     C               Binary
!    --  ---   ---         ------------
!     0    1   0                      0
!     0    2   0, 0                  00
!     1    3   1, 0, 0              100
!     2    3   0, 1, 0              010
!     3    3   1, 1, 0              011
!     4    3   0, 0, 1              100
!     8    3   0, 0, 0           (1)000
!     8    5   0, 0, 0, 1, 0      01000
!    -8    5   0, 0, 0, 1, 0  (-) 01000
!
!     0    3   0, 0, 0
!     1    3   1, 0, 0
!     2    3   0, 1, 0
!     3    3   1, 1, 0
!     4    3   0, 0, 1
!     5    3   1, 0, 1
!     6    3   0, 1, 1
!     7    3   1, 1, 1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 December 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, an integer to be represented.
!
!    Input, integer ( kind = 4 ) N, the number of binary digits to produce.
!
!    Output, integer ( kind = 4 ) C(N), the first N binary digits of I,
!    with C(1) being the units digit.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) c(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i_copy
  integer ( kind = 4 ) j

  i_copy = abs ( i )

  do j = 1, n

    c(j) = mod ( i_copy, 2 )
    i_copy = i_copy / 2

  end do

  return
end
subroutine subset_sum_serial ( n, weight, target, choice )

!*****************************************************************************80
!
!! SUBSET_SUM_SERIAL seeks a subset of a set that has a given sum.
!
!  Discussion:
!
!    This function tries to compute a target value as the sum of
!    a selected subset of a given set of weights.
!
!    This function works by brute force, that is, it tries every
!    possible subset to see if it sums to the desired value.
!
!    Given N weights, every possible selection can be described by 
!    one of the N-digit binary numbers from 0 to 2^N-1.
!
!    It is possible that there may be multiple solutions of the problem.  
!    This function will only return the first solution found.
!
!  Example:
!
!    n = 6
!    target = 22
!    w = (/ 1, 2, 4, 8, 16, 32 /)
!
!    choice = (/ 0, 1, 1, 0, 1, 0 /)
!    w(choice) = 2 + 4 + 16 = 22
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 July 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of weights.
!
!    Input, integer ( kind = 4 ) WEIGHT(N), the weights.
!
!    Input, integer ( kind = 4 ) TARGET, the target value.
!
!    Output, integer ( kind = 4 ) CHOICE(N), contains a 1 for each
!    weight that is chosen.  If no solution was found, all entries
!    are returned as -1.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) choice(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i_max
  integer ( kind = 4 ) target
  integer ( kind = 4 ) w_sum
  integer ( kind = 4 ) weight(n)

  i_max = ( 2 ** n ) - 1

  do i = 0, i_max
!
!  Convert I to a string of binary digits.
!
    call i4_to_digits_binary ( i, n, choice )
!
!  Combine the weights whose binary digit is 1.
!
    w_sum = dot_product ( choice(1:n), weight(1:n) )
!
!  Return if we matched our target sum.
!
    if ( w_sum == target ) then
      return
    end if
    
  end do

  choice(1:n) = -1

  return
end
function c = i4_to_digits_binary ( i, n )

%*****************************************************************************80
%
%% I4_TO_DIGITS_BINARY produces the binary digits of an I4.
%
%  Discussion:
%
%    An I4 is an integer.
%
%  Example:
%
%     I    N     C               Binary
%    --  ---   ---         ------------
%     0    1   0                      0
%     0    2   0, 0                  00
%     1    3   1, 0, 0              100
%     2    3   0, 1, 0              010
%     3    3   1, 1, 0              011
%     4    3   0, 0, 1              100
%     8    3   0, 0, 0           (1)000
%     8    5   0, 0, 0, 1, 0      01000
%    -8    5   0, 0, 0, 1, 0  (-) 01000
%
%     0    3   0, 0, 0
%     1    3   1, 0, 0
%     2    3   0, 1, 0
%     3    3   1, 1, 0
%     4    3   0, 0, 1
%     5    3   1, 0, 1
%     6    3   0, 1, 1
%     7    3   1, 1, 1
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    19 December 2011
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, integer I, an integer to be represented.
%
%    Input, integer N, the number of binary digits to produce.
%
%    Output, integer C(N), the first N binary digits of I,
%    with C(1) being the units digit.
%
  i_copy = floor ( abs ( i ) );

  c = zeros ( n, 1 );

  for j = 1 : n

    c(j) = mod ( i_copy, 2 );
    i_copy = floor ( i_copy / 2 );

  end

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

