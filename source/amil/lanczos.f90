FUNCTION lngamma(z) RESULT(lanczos)

!  Uses Lanczos-type approximation to ln(gamma) for z > 0.
!  Reference:
!       Lanczos, C. 'A precision approximation of the gamma
!               function', J. SIAM Numer. Anal., B, 1, 86-96, 1964.
!  Accuracy: About 14 significant digits except for small regions
!            in the vicinity of 1 and 2.

!  Programmer: Alan Miller
!              1 Creswick Street, Brighton, Vic. 3187, Australia
!  e-mail: amiller @ bigpond.net.au
!  Latest revision - 14 October 1996

IMPLICIT NONE
INTEGER, PARAMETER          :: doub_prec = SELECTED_REAL_KIND(15, 60)
REAL(doub_prec), INTENT(IN) :: z
REAL(doub_prec)             :: lanczos

! Local variables

REAL(doub_prec)  :: a(9) = (/ 0.9999999999995183D0, 676.5203681218835D0, &
                              -1259.139216722289D0, 771.3234287757674D0, &
                              -176.6150291498386D0, 12.50734324009056D0, &
                              -0.1385710331296526D0, 0.9934937113930748D-05, &
                               0.1659470187408462D-06 /), zero = 0.D0,   &
                               one = 1.d0, lnsqrt2pi = 0.9189385332046727D0, &
                               half = 0.5d0, sixpt5 = 6.5d0, seven = 7.d0, tmp
INTEGER          :: j

IF (z <= zero) THEN
  WRITE(*, *)'Error: zero or -ve argument for lngamma'
  RETURN
END IF

lanczos = zero
tmp = z + seven
DO j = 9, 2, -1
  lanczos = lanczos + a(j)/tmp
  tmp = tmp - one
END DO
lanczos = lanczos + a(1)
lanczos = LOG(lanczos) + lnsqrt2pi - (z + sixpt5) + (z - half)*LOG(z + sixpt5)
RETURN

END FUNCTION lngamma


