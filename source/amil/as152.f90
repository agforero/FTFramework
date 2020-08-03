FUNCTION chyper(point, kk, ll, mm, nn, ifault) RESULT(fn_val)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2002-09-24  Time: 22:57:18

!     Cumulative hypergeometric probabilities

!     ALGORITHM AS R77  APPL. STATIST. (1989), VOL.38, NO.1
!     Replaces AS 59 and AS 152
!     Incorporates AS R86 from vol.40(2)

!     Auxiliary routines required: ALNFAC (AS 245), ALNORM (AS 66)

IMPLICIT NONE
INTEGER, PARAMETER    :: dp = SELECTED_REAL_KIND(15, 100)

LOGICAL, INTENT(IN)   :: point
INTEGER, INTENT(IN)   :: kk
INTEGER, INTENT(IN)   :: ll
INTEGER, INTENT(IN)   :: mm
INTEGER, INTENT(IN)   :: nn
INTEGER, INTENT(OUT)  :: ifault
REAL(dp)              :: fn_val

INTEGER   :: k, l, m, n, i, j, nl, kl, mnkl
REAL(dp)  :: p, pt,  alnfac,  mean, sig, alnorm, arg
LOGICAL   :: dir
REAL(dp), PARAMETER  :: zero = 0.0_dp, half = 0.5_dp, one = 1.0_dp
INTEGER, PARAMETER   :: mvbig = 1000, mbig = 600
REAL(dp), PARAMETER  :: elimit = -88.0_dp, sxteen = 16.0_dp, scale = 1.0D35, &
                      rootpi = 2.506628274631001_dp, hundrd = 100.0_dp

k = kk + 1
l = ll + 1
m = mm + 1
n = nn + 1
dir = .true.

!     Check arguments are within permitted limits

ifault = 1
fn_val = zero
IF (n < 1 .OR. m < n .OR. k < 1 .OR. k > m) RETURN

ifault = 2
IF (l < 1 .OR. k-l > m-n) RETURN
IF (.NOT. point) fn_val = one
IF (l > n .OR. l > k) RETURN
ifault = 0
fn_val = one
IF (k == 1 .OR. k == m .OR. n == 1 .OR. n == m) RETURN
IF (.NOT. point .AND. ll == MIN(kk, nn)) RETURN

p = REAL(nn, KIND=dp) / REAL(mm - nn, KIND=dp)
IF (REAL(MIN(kk, mm-kk)) > sxteen * MAX(p, one/p) .AND.  &
      mm > mvbig .AND. elimit > -hundrd) THEN
  
!     Use a normal approximation
  
  mean = REAL(kk) * REAL(nn) / REAL(mm)
  sig = SQRT(mean * (REAL(mm-nn) / REAL(mm)) * (REAL(mm-kk) / (REAL(mm-1))))
  IF (point) THEN
    arg = -half * (((REAL(ll) - mean) / sig)**2)
    fn_val = zero
    IF (arg >= elimit) fn_val = EXP(arg) / (sig * rootpi)
  ELSE
    fn_val = alnorm((REAL(ll) + half - mean) / sig, .false.)
  END IF
  
ELSE
  
!     Calculate exact hypergeometric probabilities.
!     Interchange K and N if this saves calculations.
  
  IF (MIN(k-1, m-k) > MIN(n-1, m-n)) THEN
    i = k
    k = n
    n = i
  END IF
  IF (m-k < k-1) THEN
    dir = .NOT. dir
    l = n - l + 1
    k = m - k + 1
  END IF
  
  IF (mm > mbig) THEN
    
!     Take logarithms of factorials.
    
    p = alnfac(nn) - alnfac(mm) + alnfac(mm-kk) + alnfac(kk) +  &
        alnfac(mm-nn) - alnfac(ll) - alnfac(nn-ll) - alnfac(kk-ll)  &
        - alnfac(mm-nn-kk+ll)
    fn_val = zero
    IF (p >= elimit) fn_val = EXP(p)
  ELSE
    
!     Use Freeman/Lund algorithm
    
    DO  i = 1, l-1
      fn_val = fn_val * REAL(k-i) * REAL(n-i) / (REAL(l-i) * REAL(m-i))
    END DO
    IF (l /= k) THEN
      j = m - n + l
      DO  i = l, k-1
        fn_val = fn_val * REAL(j-i) / REAL(m-i)
      END DO
    END IF
    
  END IF
  
  IF (point) RETURN
  IF (fn_val == zero) THEN
    
!     We must recompute the point probability since it has underflowed.
    
    IF (mm <= mbig) p = alnfac(nn) - alnfac(mm) + alnfac(kk) +  &
        alnfac(mm-nn) - alnfac(ll) - alnfac(nn-ll) - alnfac(kk-ll) -  &
        alnfac(mm-nn-kk+ll) + alnfac(mm-kk)
    p = p + LOG(scale)
    IF (p < elimit) THEN
      ifault = 3
      IF (ll > REAL(nn*kk + nn + kk +1)/(mm +2)) fn_val = one
      RETURN
    ELSE
      p = EXP(p)
    END IF
  ELSE
    
!     Scale up at this point.
    
    p = fn_val * scale
  END IF
  
  pt = zero
  nl = n - l
  kl = k - l
  mnkl = m - n - kl + 1
  IF (l <= kl) THEN
    DO  i = 1, l-1
      p = p * REAL(l-i) * REAL(mnkl-i) / (REAL(nl+i) * REAL(kl+i))
      pt = pt + p
    END DO
    IF (p == zero) ifault = 3
  ELSE
    dir = .NOT. dir
    DO  j = 0, kl-1
      p = p * REAL(nl-j) * REAL(kl-j) / (REAL(l+j) * REAL(mnkl+j))
      pt = pt + p
    END DO
    IF (p == zero) ifault = 3
  END IF
  
  IF (dir) THEN
    fn_val = fn_val + (pt / scale)
  ELSE
    fn_val = one - (pt / scale)
  END IF
  
END IF

END FUNCTION chyper
