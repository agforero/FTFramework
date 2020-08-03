SUBROUTINE Easter(year, day, month)
! U.S. Naval Observatory Astronomical Applications Department

!            The Date of Easter (Sunday)

! The algorithm is due to J.-M. Oudin (1940) and is reprinted in the Explanatory
! Supplement to the Astronomical Almanac, ed. P. K. Seidelmann (1992).
! See Chapter 12, "Calendars", by L. E. Doggett.

! The following are dates of Easter from 1980 to 2024:

! 1980  April  6        1995  April 16        2010  April  4
! 1981  April 19        1996  April  7        2011  April 24
! 1982  April 11        1997  March 30        2012  April  8
! 1983  April  3        1998  April 12        2013  March 31
! 1984  April 22        1999  April  4        2014  April 20
! 1985  April  7        2000  April 23        2015  April  5
! 1986  March 30        2001  April 15        2016  March 27
! 1987  April 19        2002  March 31        2017  April 16
! 1988  April  3        2003  April 20        2018  April  1
! 1989  March 26        2004  April 11        2019  April 21
! 1990  April 15        2005  March 27        2020  April 12
! 1991  March 31        2006  April 16        2021  April  4
! 1992  April 19        2007  April  8        2022  April 17
! 1993  April 11        2008  March 23        2023  April  9
! 1994  April  3        2009  April 12        2024  March 31

! N.B. The date of Easter for the Eastern Orthodox Church may be different.

! This code assembled by Alan Miller
! Reference web site:
! http://aa.usno.navy.mil/faq/docs/easter.html
! Latest revision 8 April 2002

IMPLICIT NONE
INTEGER, INTENT(IN)   :: year
INTEGER, INTENT(OUT)  :: day, month

! Local variables
INTEGER  :: c, i, j, k, l, n

c = year / 100
n = year - 19 * ( year / 19 )
k = ( c - 17 ) / 25
i = c - c / 4 - ( c - k ) / 3 + 19 * n + 15
i = i - 30 * ( i / 30 )
i = i - (i / 28) * (1 - (i / 28) * (29 / (i + 1 )) * ( (21 - n) / 11) )
j = year + year / 4 + i + 2 - c + c / 4
j = j - 7 * ( j / 7 )
l = i - j
month = 3 + ( l + 40 ) / 44
day = l + 28 - 31 * ( month / 4 )

RETURN
END SUBROUTINE Easter



PROGRAM test_Easter
! Generate the table contained in the subroutine.

IMPLICIT NONE

INTERFACE
  SUBROUTINE Easter(year, day, month)
    IMPLICIT NONE
    INTEGER, INTENT(IN)   :: year
    INTEGER, INTENT(OUT)  :: day, month
  END SUBROUTINE Easter
END INTERFACE

INTEGER            :: day(1980:2024), month(1980:2024), year, y2, y3
CHARACTER (LEN=5)  :: mon(3:4) = (/ 'March', 'April' /)

DO year = 1980, 2024
  CALL Easter(year, day(year), month(year))
END DO

DO year = 1980, 1994
  y2 = year + 15
  y3 = y2 + 15
  WRITE(*, '(3("  ", i4, "  ", a5, i3, "      "))')  &
           year, mon(month(year)), day(year),  &
           y2, mon(month(y2)), day(y2),  &
           y3, mon(month(y3)), day(y3)
END DO
STOP
END PROGRAM test_Easter
