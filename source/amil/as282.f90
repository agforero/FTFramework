MODULE robust_regression

IMPLICIT NONE
INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(12, 60)

CONTAINS


SUBROUTINE mvelms (x, y, ncas, npre, ioptn, maxtry, nclow, nchigh, coeffs,  &
                   eprmin, resid, robdsq, cvemin, ifault)

!     ALGORITHM AS 282 APPL.STATIST. (1993), VOL.42, NO.2

!     High breakdown regression and multivariate estimation

! ROUTINE TO CALCULATE LEAST MEDIAN OF SQUARES REGRESSION, MINIMUM VOLUME
! ELLIPSOID, AND ASSOCIATED STATISTICS

! N.B. This version has fewer arguments than AS 282.

! Converted to be compatible with ELF90 by Alan Miller
! amiller @ bigpond.net.au
! http://users.bigpond.net.au/amiller

! Latest revision - 11 July 1997

REAL (dp), INTENT(IN)  :: x(:,:), y(:)
REAL (dp), INTENT(OUT) :: coeffs(:,:), eprmin(:), resid(:,:), cvemin(:),  &
                          robdsq(:,:)
INTEGER, INTENT(IN)    :: ncas, npre, ioptn, maxtry, nclow, nchigh
INTEGER, INTENT(OUT)   :: ifault

!   OPTIONAL ADDITIONAL ARGUMENTS FOR MORE SPECIALIZED USE

!   DIMENSION XINV(NVDIM,NVDIM,*), ELSCAL(*), LMSBAS(NVDIM,*), MVEBAS(NVDIM,*)

! Local variables
INTEGER                :: itemp, nvar, int, j2, j3, j4, i2, i3, i4, i5, nvr2, &
                          lfresh, ncount, iptr, ncrang, i, ixlo, j, ii,  &
                          inxptr, ncol, kx, kxsto, icount, kkx, jj, k
REAL (dp)              :: rpp1, addcon, deter, powmed, ymed, ymad, detadj,  &
                          target, offs, short, criter, basvol, volu
LOGICAL                :: lms, mve, isint, exh, exact
REAL (dp), PARAMETER   :: toler = 1.e-9_dp, big = 1.e30_dp, half = 0.5_dp,  &
                          one = 1._dp, ten = 10._dp, zero = 0._dp
INTEGER, PARAMETER     :: nfresh = 50
INTEGER, ALLOCATABLE   :: iwork(:)
REAL (dp), ALLOCATABLE :: data(:,:), work(:)

!     EXTRACT THE OPTIONS REQUESTED

itemp = ioptn
lms = MOD(itemp,2) == 0
itemp = itemp / 2
mve = MOD(itemp,2) == 0
itemp = itemp / 2
isint = MOD(itemp,2) == 0
itemp = itemp / 2
exh = MOD(itemp,2) == 0
IF (isint) THEN
  nvar = npre + 1
  addcon = one / DBLE(nvar)
  rpp1 = SQRT(addcon)
  int = 1
ELSE
  nvar = npre
  addcon = zero
  rpp1 = one
  int = 0
END IF

! CHECK FAULT PARAMETERS

ifault = 0
IF (nclow < nvar .OR. nclow > nchigh .OR. nchigh > ncas) ifault = 2
IF (ncas < nvar) ifault = 1
IF (isint .AND. nvar == 1 .AND. mve) ifault = ifault+4
IF (ifault > 0) GO TO 490

! Allocate arrays IWORK, DATA & WORK

ALLOCATE( iwork(4*nvar+ncas), data(0:nvar,nvar+ncas), work(3*ncas+nvar) )

exact = .false.
j2 = ncas
j3 = 2*ncas
j4 = 3*ncas
i2 = nvar
i3 = 2*nvar
i4 = 3*nvar
i5 = i4 + ncas
nvr2 = 2*nvar
deter = rpp1
powmed = DBLE(npre)*half
DO i = 1, nvar
  iwork(i) = ncas + i - nvar
END DO
iptr = 1
ncrang = nchigh - nclow + 1
DO i = 1, ncrang
  cvemin(i) = big
  eprmin(i) = big
END DO
lfresh = 0
ncount = 0

! TARGET VARIABLE, IF PRESENT, IS TREATED IN STANDARDIZED FORM SO THAT EXACT
! FIT CAN BE DETECTED MORE EASILY

IF (lms) THEN
  ixlo = 0
  DO i = 1, ncas
    work(i) = y(i)
  END DO
  
! THE SUBROUTINE SORTSUB(RA,N,ILOW) SORTS ENTRIES ILOW+1 THROUGH ILOW+N
! IN ASCENDING ORDER FROM A VECTOR RA
  
  CALL sortsub(work, ncas, 0)
  ymed = (work(ncas/2+1) + work((ncas+1)/2))*half
  DO i = 1, ncas
    work(i) = ABS(y(i) - ymed)
  END DO
  CALL sortsub(work, ncas, 0)
  ymad = (work(ncas/2+1) + work((ncas+1)/2))*half
  IF (ymad == zero) ymad = one
ELSE
  ixlo = 1
END IF
detadj = zero
DO j = 1, npre
  DO i = 1, ncas
    work(i) = ABS(x(i,j))
  END DO
  CALL sortsub(work, ncas, 0)
  work(j4+j+int) = (work(ncas/2+1) + work((ncas+1)/2))*half
  IF (work(j4+j+int) == zero) work(j4+j+int) = work(ncas)
  detadj = detadj + LOG10(work(j4+j+int))
END DO
IF (isint) work(j4+1) = one

! DATA IS TRANSFERED TO WORKAREA; INITIAL SIMPLEX TABLEAU IS SET UP

CALL refresh(data, ixlo, x, ncas, npre, y, ymad, nvar, work,  &
             j4, iwork, i2, i3, i5, 0, isint, int, lms, deter, lfresh, toler)

! INITIAL BASIS IS SET UP. CHECKS ARE MADE THAT INITIAL BASIS MEMBERS ARE IN
! GENERAL POSITION

IF (exh) THEN
  iwork(i2+1) = 1
  DO i = 1, nvar
    70 j = iwork(i2+i)
    DO ii = i, nvar
      IF (ABS(data(ii,j)) >= toler) THEN
        CALL pivot(data, ixlo, nvar, ncas, iwork, i5, ii, j, deter)
        iwork(i3+i) = j
        IF (ii /= i) CALL swap(data, nvar, ncas, ii, i)
        GO TO 90
      END IF
    END DO
    IF (j == ncas) THEN
      ifault = ifault + 8
      IF (ifault >= 24) ifault = ifault - 16
      GO TO 490
    ELSE
      ifault = 16
      iwork(i2+i) = iwork(i2+i) + 1
      GO TO 70
    END IF
    90 IF (i < nvar) iwork(i2+i+1) = iwork(i2+i) + 1
  END DO
  iptr = nvar
  GO TO 230
ELSE
  DO i = 1, ncas
    iwork(i4+i) = i
  END DO
  CALL perm(iwork, i4, ncas, 0)
  inxptr = 0
  DO i = 1, nvar
    120 inxptr = inxptr+1
    j = iwork(i4+inxptr)
    DO ii = i, nvar
      IF (ABS(DATA(ii,j)) >= toler) THEN
        CALL pivot(DATA, ixlo, nvar, ncas, iwork, i5, ii, j, deter)
        iwork(i3+i) = j
        IF (ii /= i) CALL swap(DATA, nvar, ncas, ii, i)
        iwork(i2+i) = j
        CYCLE
      END IF
    END DO
    IF (inxptr == ncas) THEN
      ifault = ifault+8
      IF (ifault >= 24) ifault = ifault - 16
      GO TO 490
    ELSE
      ifault = 16
      GO TO 120
    END IF
  END DO
  GO TO 230
END IF


! MAIN ANALYSIS LOOP. GENERATE ALL SUBSETS (IF EXH) OR A SUBSET (IF NOT EXH)

150 IF (exh) THEN
  
! IF EXHAUSTIVE, SUCCESSIVE BASES ARE CONSIDERED
  
  160 iwork(i2+iptr) = iwork(i2+iptr) + 1
  IF (iwork(i2+iptr) > iwork(iptr)) THEN
    iptr = iptr - 1
    IF (iptr == 0) THEN
      GO TO 440
    ELSE
      GO TO 160
    END IF
  ELSE
    
! IF DATA IS NOT IN GENERAL POSITION, NEXT POSITION IN LIST MUST BE FOUND
    
    IF (ABS(DATA(iptr,iwork(i2+iptr))) >= toler) GO TO 210
    IF (lfresh > nvr2) THEN
      deter = rpp1
      CALL refresh(DATA, ixlo, x, ncas, npre, y, ymad, nvar, work, j4, iwork, &
                   i2, i3, i5, iptr-1, isint, INT, lms, deter, lfresh, toler)
      IF (ABS(DATA(iptr,iwork(i2+iptr))) >= toler) GO TO 210
      ifault = 16
      IF (iptr == nvar) THEN
        GO TO 160
      ELSE
        GO TO 170
      END IF
    END IF
    ifault = 16
    IF (iptr == nvar) GO TO 160
    deter = rpp1
    CALL refresh(DATA, ixlo, x, ncas, npre, y, ymad, nvar, work, j4, iwork, &
                 i2, i3, i5, iptr-1, isint, INT, lms, deter, lfresh, toler)
    170 j = iwork(i2+iptr)
    DO ii = iptr, nvar
      IF (ABS(DATA(ii,j)) >= toler) THEN
        CALL pivot(DATA, ixlo, nvar, ncas, iwork, i5, ii, j, deter)
        iwork(i3+iptr) = j
        IF (ii /= iptr) CALL swap(DATA, nvar, ncas, ii, iptr)
        GO TO 220
      END IF
    END DO
    iwork(i2+iptr) = iwork(i2+iptr)+1
    IF (iwork(i2+iptr) > iwork(iptr)) THEN
      iptr = iptr-1
      IF (iptr == 0) THEN
        GO TO 440
      ELSE
        iwork(i2+iptr) = iwork(i2+iptr) + 1
      END IF
    END IF
    GO TO 170
  END IF
ELSE
  
! IF NOT EXHAUSTIVE, THE NEXT ENTRY IN RANDOM PERMUTATION VECTOR IS ENTERED
! INTO THE BASIS
  
  IF (ncount > maxtry) GO TO 440
  iptr = MOD(iptr,nvar) + 1
  190 inxptr = inxptr+1
  IF (inxptr > ncas) THEN
    
! IF COME TO THE END OF THE PERMUTATION VECTOR, GENERATE A NEW ONE
    
    DO i = 1, nvar
      iwork(i4+ncas-nvar+i) = iwork(i4+i)
      iwork(i4+i) = iwork(i2+i)
    END DO
    CALL perm(iwork, i4, ncas, nvar)
    inxptr = nvar + 1
  END IF
  ncol = iwork(i4+inxptr)
  IF (ABS(DATA(iptr,ncol)) < toler) THEN
    IF (lfresh > nvr2) THEN
      deter = rpp1
      CALL refresh(DATA, ixlo, x, ncas, npre, y, ymad, nvar, work, j4, iwork, &
                   i2, i3, i5, iptr-1, isint, INT, lms, deter, lfresh, toler)
      IF (ABS(DATA(iptr,ncol)) < toler) THEN
        ifault = 16
        GO TO 190
      END IF
    ELSE
      ifault = 16
      GO TO 190
    END IF
  END IF
  iwork(i2+iptr) = ncol
END IF
210 ncount = ncount + 1
lfresh = lfresh + 1

! CARRY OUT NEXT PIVOT, THEREBY CREATING A NEW BASIS

CALL pivot (DATA, ixlo, nvar, ncas, iwork, i3, iptr, iwork(i2+iptr), deter)

! RECOMPUTE THE INVERSE BASIS EVERY NFRESH SIMPLEX PIVOTS

IF (lfresh > nfresh) THEN
  deter = rpp1
  CALL refresh(DATA, ixlo, x, ncas, npre, y, ymad, nvar, work, j4, iwork,  &
               i2, i3, i5, nvar, isint, INT, lms, deter, lfresh, toler)
END IF
220 IF (exh .AND. iptr < nvar) THEN
  iptr = iptr + 1
  iwork(i2+iptr) = iwork(i2+iptr-1)
  GO TO 150
END IF
230 IF (lms) THEN
  
! CHECK TO SEE IF THIS BASIS HAS A SMALLER LMS CRITERION VALUE THAN PREVIOUS
! BASES. IF THERE IS NO INTERCEPT, A PRELIMINARY COMPARISON IS MADE TO SEE IF
! THIS BASIS COULD BE OPTIMAL.
  
  IF (exact) GO TO 340
  IF (isint) THEN
    DO j = 1, ncas
      work(j) = DATA(0,j)
    END DO
  ELSE
    DO j = 1, ncas
      work(j) = ABS(DATA(0,j))
    END DO
    DO kx = nclow, nchigh
      kxsto = kx - nclow + 1
      target = eprmin(kxsto)
      icount = 0
      DO j = 1, ncas
        IF (work(j) <= target) icount = icount+1
      END DO
      IF (icount >= kx) GO TO 280
    END DO
    GO TO 340
  END IF
  280 CALL sortsub(work, ncas, 0)
  DO kx = nclow, nchigh
    kxsto = kx - nclow + 1
    offs = zero
    IF (isint) THEN
      
! CALCULATE PROPER OFFSET FOR INTERCEPT TERM, IF REQUIRED, AND DETERMINE
! CRITERION FOR THIS BASIS
      
      short = big
      DO kkx = kx, ncas
        IF (work(kkx)-work(kkx-kx+1) < short) THEN
          short = work(kkx) - work(kkx-kx+1)
          criter = short*half
          offs = (work(kkx) + work(kkx-kx+1))*half
        END IF
      END DO
    ELSE
      criter = work(kx)
    END IF
    
! UPDATE SUMMARY STATISTICS IF THIS BASIS IS CURRENTLY OPTIMAL
    
    IF (criter < eprmin(kxsto)) THEN
      eprmin(kxsto) = criter
!         DO 300 J = 1, NVAR
!300 LMSBAS(J,KXSTO) = IWORK(I3+J)
      DO j = 1, nvar
        coeffs(j,kxsto) = -DATA(0, ncas+j)*ymad
      END DO
      IF (isint) coeffs(1, kxsto) = coeffs(1,kxsto) + offs*ymad
      DO j = 1, ncas
        resid(j,kxsto) = (DATA(0,j) - offs)*ymad
      END DO
      
! A CHECK FOR EXACT FIT IS MADE
      
      IF ((criter < toler).AND.(kx == nchigh)) THEN
        exact = .true.
        IF (.NOT.mve) GO TO 460
      END IF
    END IF
  END DO
  
! IF THIS IS A UNIVARIATE LOCATION PROBLEM, NO FURTHER ANALYSIS IS NEEDED
  
  IF ((isint).AND.(nvar == 1)) GO TO 460
END IF
340 IF (mve) THEN
  
! CHECK TO SEE IF THIS BASIS HAS A SMALLER MVE CRITERION VALUE THAN PREVIOUS
! BASES.  A PRELIMINARY COMPARISON IS MADE TO SEE IF THIS BASIS COULD BE
! OPTIMAL.
  
  basvol = LOG10(ABS(deter))
  DO j = 1, ncas
    work(j+j2) = -addcon
    DO i = 1, nvar
      work(j+j2) = work(j+j2) + DATA(i,j)**2
    END DO
    IF (work(j+j2) > zero) THEN
      work(j+j2) = LOG10(work(j+j2))
    ELSE
      work(j+j2) = -big
    END IF
    work(j+j3) = work(j+j2)
  END DO
  DO kx = nclow, nchigh
    kxsto = kx-nclow+1
    target = (cvemin(kxsto)-basvol)/powmed
    icount = 0
    DO j = 1, ncas
      IF (work(j+j3) <= target) icount = icount+1
    END DO
    IF (icount >= kx) GO TO 390
  END DO
  GO TO 150
  390 CALL sortsub(work, ncas, j3)
  DO kx = nclow, nchigh
    kxsto = kx-nclow+1
    volu = basvol+powmed*work(kx+j3)
    IF (volu < cvemin(kxsto)) THEN
      
! UPDATE SUMMARY STATISTICS IF THIS BASIS IS CURRENTLY OPTIMAL
      
      cvemin(kxsto) = volu
!         ELSCAL(KXSTO) = WORK(KX+J3)
!         DO 400 J = 1, NVAR
!400 MVEBAS(J,KXSTO) = IWORK(I3+J)
      DO jj = 1, ncas
        robdsq(jj,kxsto) = work(jj+j2)-work(kx+j3)
      END DO
!         DO 420 J = 1, NVAR
!         DO 420 JJ = 1, NVAR
!420 XINV(JJ,J,KXSTO) = DATA(J,JJ+NCAS)/WORK(J4+JJ)
    END IF
  END DO
END IF

! END OF MAIN PROCESSING LOOP. GO TO BEGINNING AND GENERATE ANOTHER BASIS

GO TO 150
440 IF (mve) THEN
  
! CONVERT ROBUST DISTANCES AND SCALING BACK TO ORIGINAL SCALE
  
  DO k = 1, ncrang
!       ELSCAL(K) = TEN**ELSCAL(K)
    cvemin(k) = cvemin(k) + detadj
    DO j = 1, ncas
      IF (robdsq(j,k) > (-big)) THEN
        robdsq(j,k) = ten**robdsq(j,k)
      ELSE
        robdsq(j,k) = zero
      END IF
    END DO
  END DO
END IF
460 IF (lms) THEN
  DO k = 1, ncrang
    eprmin(k) = eprmin(k)*ymad
    DO i = 1, npre
      coeffs(i+int,k) = coeffs(i+int,k) / work(i+j4+int)
    END DO
  END DO
END IF

490 IF (ALLOCATED(iwork)) DEALLOCATE( iwork, data, work )
RETURN
END SUBROUTINE mvelms



SUBROUTINE refresh(data, ixlo, x, ncas, npre, y, ymad,  &
                   nvar, work, j4, iwork, i2, i3, i5, iup, isint, int, lms,  &
                   deter, lfresh, toler)

! SUBROUTINE TO REFRESH FIRST IUP ENTRIES OF SIMPLEX BASIS

INTEGER, INTENT(IN)       :: ixlo, ncas, npre, nvar, j4, i2, i3, i5, iup, int
INTEGER, INTENT(IN OUT)   :: iwork(:)
INTEGER, INTENT(OUT)      :: lfresh
REAL (dp), INTENT(IN)     :: x(:,:), y(:), ymad, work(:), toler
REAL (dp), INTENT(IN OUT) :: data(0:,:), deter
LOGICAL, INTENT(IN)       :: isint, lms

! Local variables
REAL (dp), PARAMETER :: zero = 0._dp, one = 1._dp
INTEGER              :: i, j, ii

lfresh = 0
IF (isint) THEN
  DATA(1,1:ncas) = one
END IF

DO i = 1, ncas
  DATA(int+1:int+npre,i) = x(i,1:npre) / work(j4+int+1:j4+int+npre)
END DO

DO i = 1, nvar
  data(1:nvar,ncas+i) = zero
  data(i,ncas+i) = one
END DO

IF (lms) THEN
  data(0,1:ncas) = y(1:ncas) / ymad
  data(0,ncas+1:ncas+nvar) = zero
END IF

DO i = 1, iup
  j = iwork(i2+i)
  DO ii = i, nvar
    IF (ABS(data(ii,j)) >= toler) THEN
      CALL pivot(data, ixlo, nvar, ncas, iwork, i5, ii, j, deter)
      IF (ii /= i) CALL swap(data, nvar, ncas, ii, i)
      iwork(i3+i) = j
      EXIT
    END IF
  END DO
END DO

RETURN
END SUBROUTINE refresh



SUBROUTINE pivot(x, ixlo, nord, ncas, iwork, i3, nrow, ncol, deter)

! SUBROUTINE TO PIVOT ENTRY CORRESPONDING TO (NROW,NCOL) INTO SIMPLEX TABLEAU

INTEGER, INTENT(IN)       :: ixlo, nord, ncas, i3, nrow, ncol
REAL (dp), INTENT(IN OUT) :: x(:,:), deter
INTEGER, INTENT(IN OUT)   :: iwork(:)

! Local variables
INTEGER              :: nhigh, j, i
REAL (dp)            :: fmult, pivt
REAL (dp), PARAMETER :: zero = 0._dp, one = 1._dp

pivt = x(nrow,ncol)
deter = deter*pivt
nhigh = nord + ncas
DO j = 1, nhigh
  IF (j /= ncol) THEN
    fmult = x(nrow,j)/pivt
    DO i = ixlo, nord
      IF (i /= nrow) x(i,j) = x(i,j) - fmult*x(i,ncol)
    END DO
    x(nrow,j) = fmult
  END IF
END DO
x(ixlo:nord,ncol) = zero
x(nrow,ncol) = one
iwork(i3+nrow) = ncol

RETURN
END SUBROUTINE pivot



SUBROUTINE perm(index, iaa, n, iab)

! SUBROUTINE TO RETURN PSEUDORANDOM PERMUTATION OF N - IAB ENTRIES IN
! VECTOR INDEX STARTING AT IAA+1

INTEGER, INTENT(IN OUT) :: index(:)
INTEGER, INTENT(IN)     :: iaa, n, iab

! Local variables
INTEGER :: m, j, itemp

m = n-iab

! GENERATE A RANDOM DIGIT FROM 1 TO M USING THE PSEUDORANDOM U(0,1) VARIATE
! RANDOM() (SUCH AS FROM ALGORITHM AS 183)

DO
  j = INT(random()*DBLE(m))+1

! SWAP ENTRIES IN INDEX CORRESPONDING TO J AND M, OFFSET BY IAB

  IF (j /= m) THEN
    itemp = index(iaa+j+iab)
    index(iaa+j+iab) = index(iaa+m+iab)
    index(iaa+m+iab) = itemp
  END IF
  m = m-1
  IF (m <= 1) EXIT
END DO

RETURN
END SUBROUTINE perm



SUBROUTINE swap(data, nvar, ncas, ir1, ir2)

! SUBROUTINE TO SWAP ROWS IR1 AND IR2 OF MATRIX DATA

REAL (dp), INTENT(IN OUT) :: data(0:,:)
INTEGER, INTENT(IN)       :: nvar, ncas, ir1, ir2

! Local variables
INTEGER   :: nhigh, i
REAL (dp) :: temp

nhigh = nvar + ncas
DO i = 1, nhigh
  temp = data(ir1,i)
  data(ir1,i) = data(ir2,i)
  data(ir2,i) = temp
END DO

RETURN
END SUBROUTINE swap



FUNCTION random() RESULT(r)
REAL (dp) :: r

CALL RANDOM_NUMBER(r)
RETURN
END FUNCTION random




SUBROUTINE sortsub(ra, n, ilow)

! SUBROUTINE TO SORT ENTRIES ILOW+1 THROUGH ILOW+N IN ASCENDING ORDER FROM A
! VECTOR RA, USING HEAPSORT

REAL (dp), INTENT(IN OUT) :: ra(:)
INTEGER, INTENT(IN)       :: n, ilow

! Local variables
INTEGER   :: l, ir, i, j
REAL (dp) :: rra

l = n/2+1
ir = n
10 IF (l > 1) THEN
  l = l-1
  rra = ra(l+ilow)
ELSE
  rra = ra(ir+ilow)
  ra(ir+ilow) = ra(1+ilow)
  ir = ir-1
  IF (ir == 1) THEN
    ra(1+ilow) = rra
    RETURN
  END IF
END IF

i = l
j = l+l
20 IF (j <= ir) THEN
  IF (j < ir) THEN
    IF (ra(j+ilow) < ra(j+1+ilow)) j = j+1
  END IF
  IF (rra < ra(j+ilow)) THEN
    ra(i+ilow) = ra(j+ilow)
    i = j
    j = j+j
  ELSE
    j = ir+1
  END IF
  GO TO 20
END IF
ra(i+ilow) = rra
GO TO 10

END SUBROUTINE sortsub

END MODULE robust_regression
