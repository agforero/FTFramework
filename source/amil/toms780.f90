FUNCTION REXPU() RESULT(fn_val)
!
!     Random-number generator for the exponential distribution
!     Algorithm EA from J. H. Ahrens and U. Dieter,
!     Communications of the ACM, 31 (1988) 1330--1337.
!     Coded by K. G. Hamilton, December 1996, with corrections.
!

IMPLICIT NONE
REAL :: fn_val

REAL, PARAMETER :: alog2 = 0.6931471805599453,  a =  5.7133631526454228,  &
                       b = 3.4142135623730950,  c = -1.6734053240284925,  &
                       p = 0.9802581434685472, aa =  5.6005707569738080,  &
                      bb = 3.3468106480569850, hh =  0.0026106723602095,  &
                      dd = 0.0857864376269050
REAL :: g, u, up, y
!
!       Comment out the following lines if your RNG can never return 0.0
DO
  call random_number(u)
  if (u > 0.0) EXIT               ! Zero-protector
END DO

g = c
DO
  u = u + u
  if (u >= 1.0) EXIT
  g = g + alog2
END DO

u = u - 1.0
if (u <= p) THEN
  fn_val = g + aa/(bb-u)
  return
END IF
DO
  call random_number(u)
  y = a/(b-u)
  call random_number(up)
  if ((up*hh + dd)*(b-u)**2 <= exp(-(y+c))) EXIT
END DO

fn_val = g + y

return
end FUNCTION REXPU


FUNCTION REXPS() RESULT(fn_val)
!
!     Random-number generator for the exponential distribution
!     Algorithm EA from J. H. Ahrens and U. Dieter,
!     Communications of the ACM, 31 (1988) 1330--1337.
!     Coded by K. G. Hamilton, December 1996, with corrections.
!
IMPLICIT NONE
REAL :: fn_val

REAL, PARAMETER :: alog2 = 0.6931471805599453,  a =  5.7133631526454228,  &
                       b = 3.4142135623730950,  c = -1.6734053240284925,  &
                       p = 0.9802581434685472, aa =  5.6005707569738080,  &
                      bb = 3.3468106480569850, hh =  0.0026106723602095,  &
                      dd = 0.0857864376269050
REAL :: g, u, up, y
!
!       Comment out the following lines if your RNG can never return 0.0
DO
  call random_number(u)
  if (u > 0.0) EXIT               ! Zero-protector
END DO

g = c
u = u + u
do while (u < 1.0)
   g = g + alog2
   u = u + u
end do
u = u - 1.0
if (u <= p) then
  fn_val = g + aa/(bb-u)
  return
end if

do
  call random_number(u)
  y = a/(b-u)
  call random_number(up)
  if ((up*hh + dd)*(b-u)**2 <= exp(-(y+c))) then
    fn_val = g + y
    return
  end if
end do

RETURN
end FUNCTION REXPS
