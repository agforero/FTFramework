SUBROUTINE diptst(x, n, dip, xl, xu, ifault, gcm, lcm, mn, mj)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2001-03-13  Time: 22:28:38

!     ALGORITHM AS 217 APPL. STATIST. (1985) VOL.34, NO.3

!     Does the dip calculation for an ordered vector X using the
!     greatest convex minorant and the least concave majorant, skipping
!     through the data using the change points of these distributions.
!     It returns the dip statistic 'DIP' and the modal interval
!     (XL, XU).

IMPLICIT NONE

REAL, INTENT(IN)      :: x(:)
INTEGER, INTENT(IN)   :: n
REAL, INTENT(OUT)     :: dip
REAL, INTENT(OUT)     :: xl
REAL, INTENT(OUT)     :: xu
INTEGER, INTENT(OUT)  :: ifault
INTEGER, INTENT(OUT)  :: gcm(:)
INTEGER, INTENT(OUT)  :: lcm(:)
INTEGER, INTENT(OUT)  :: mn(:)
INTEGER, INTENT(OUT)  :: mj(:)

INTEGER          :: high, ic, icv, icva, icx, icxa, ig, igcm, igcm1, igcmx, &
                    ih, iv, ix, j, jb, je, jk, jr, k, kb, ke, kr,  &
                    lcm1, lcmiv, lcmiv1, low, mjk, mjmjk, mnj, mnmnj, na
REAL             :: a, b, const, d, dipnew, dl, du, dx, fn, t, temp
REAL, PARAMETER  :: zero = 0.0, half = 0.5, one = 1.0

ifault = 1
IF (n <= 0) RETURN
ifault = 0

!     Check if N = 1

IF (n /= 1) THEN
  
!     Check that X is sorted
  
  ifault = 2
  DO  k = 2, n
    IF (x(k) < x(k-1)) RETURN
  END DO
  ifault = 0
  
!     Check for all values of X identical,
!     and for 1 < N < 4.
  
  IF (x(n) > x(1) .AND. n >= 4) GO TO 20
END IF
xl = x(1)
xu = x(n)
dip = zero
RETURN

!     LOW contains the index of the current estimate of the lower end
!     of the modal interval, HIGH contains the index for the upper end.

20 fn = n
low = 1
high = n
dip = one / fn
xl = x(low)
xu = x(high)

!     Establish the indices over which combination is necessary for the
!     convex minorant fit.

mn(1) = 1
DO  j = 2, n
  mn(j) = j - 1
  30 mnj = mn(j)
  mnmnj = mn(mnj)
  a = mnj - mnmnj
  b = j - mnj
  IF (mnj /= 1 .AND. (x(j)-x(mnj))*a >= (x(mnj)-x(mnmnj))*b) THEN
    mn(j) = mnmnj
    GO TO 30
  END IF
END DO

!     Establish the indices over which combination is necessary for the
!     concave majorant fit.

mj(n) = n
na = n - 1
DO  jk = 1, na
  k = n - jk
  mj(k) = k + 1
  50 mjk = mj(k)
  mjmjk = mj(mjk)
  a = mjk - mjmjk
  b = k - mjk
  IF (mjk /= n .AND. (x(k)-x(mjk))*a >= (x(mjk)-x(mjmjk))*b) THEN
    mj(k) = mjmjk
    GO TO 50
  END IF
END DO

!     Start the cycling.
!     Collect the change points for the GCM from HIGH to LOW.

70 ic = 1
gcm(1) = high

80 igcm1 = gcm(ic)
ic = ic + 1
gcm(ic) = mn(igcm1)
IF (gcm(ic) > low) GO TO 80
icx = ic

!     Collect the change points for the LCM from LOW to HIGH.

ic = 1
lcm(1) = low
90 lcm1 = lcm(ic)
ic = ic + 1
lcm(ic) = mj(lcm1)
IF (lcm(ic) < high) GO TO 90
icv = ic

!     ICX, IX, IG are counters for the convex minorant,
!     ICV, IV, IH are counters for the concave majorant.

ig = icx
ih = icv

!     Find the largest distance greater than 'DIP' between the GCM and
!     the LCM from LOW to HIGH.

ix = icx - 1
iv = 2
d = zero
IF (icx == 2 .AND. icv == 2) THEN
  d = one / fn
  GO TO 120
END IF

100 igcmx = gcm(ix)
lcmiv = lcm(iv)
IF (igcmx <= lcmiv) THEN
  
!     If the next point of either the GCM or LCM is from the LCM,
!     calculate the distance here.
  
  lcmiv1 = lcm(iv-1)
  a = lcmiv - lcmiv1
  b = igcmx - lcmiv1 - 1
  dx = (x(igcmx)-x(lcmiv1)*a) / (fn*(x(lcmiv)-x(lcmiv1))) - b / fn
  ix = ix - 1
  IF (dx < d) GO TO 110
  d = dx
  ig = ix + 1
  ih = iv
ELSE
  
!     If the next point of either the GCM or LCM is from the GCM,
!     calculate the distance here.
  
  lcmiv = lcm(iv)
  igcm = gcm(ix)
  igcm1 = gcm(ix+1)
  a = lcmiv - igcm1 + 1
  b = igcm - igcm1
  dx = a / fn - ((x(lcmiv)-x(igcm1))*b) / (fn*(x(igcm)-x(igcm1)))
  iv = iv + 1
  IF (dx >= d) THEN
    d = dx
    ig = ix + 1
    ih = iv - 1
  END IF
END IF

110 IF (ix < 1) ix = 1
IF (iv > icv) iv = icv
IF (gcm(ix) /= lcm(iv)) GO TO 100

120 IF (d >= dip) THEN
  
!     Calculate the DIPs for the current LOW and HIGH.
  
!     The DIP for the convex minorant.
  
  dl = zero
  IF (ig /= icx) THEN
    icxa = icx - 1
    DO  j = ig, icxa
      temp = one / fn
      jb = gcm(j+1)
      je = gcm(j)
      IF (je-jb > 1) THEN
        IF (x(je) /= x(jb)) THEN
          a = je - jb
          const = a / (fn*(x(je)-x(jb)))
          DO  jr = jb, je
            b = jr - jb + 1
            t = b / fn - (x(jr)-x(jb)) * const
            IF (t > temp) temp = t
          END DO
        END IF
      END IF
      IF (dl < temp) dl = temp
    END DO
  END IF
  
!     The DIP for the concave majorant.
  
  du = zero
  IF (ih /= icv) THEN
    icva = icv - 1
    DO  k = ih, icva
      temp = one / fn
      kb = lcm(k)
      ke = lcm(k+1)
      IF (ke-kb > 1) THEN
        IF (x(ke) /= x(kb)) THEN
          a = ke - kb
          const = a / (fn*(x(ke)-x(kb)))
          DO  kr = kb, ke
            b = kr - kb - 1
            t = (x(kr)-x(kb)) * const - b / fn
            IF (t > temp) temp = t
          END DO
        END IF
      END IF
      IF (du < temp) du = temp
    END DO
  END IF
  
!     Determine the current maximum.
  
  dipnew = dl
  IF (du > dl) dipnew = du
  IF (dip < dipnew) dip = dipnew
  low = gcm(ig)
  high = lcm(ih)
  
!     Recycle
  
  GO TO 70
END IF

dip = half * dip
xl = x(low)
xu = x(high)

RETURN
END SUBROUTINE diptst
