module str2int_mod ! from https://stackoverflow.com/questions/24071722/
  contains 
        elemental subroutine str2int(str,int,stat)
        implicit none
        ! Arguments
        character(len=*),intent(in) :: str
        integer,intent(out)         :: int
        integer,intent(out)         :: stat
  
        read(str,*,iostat=stat)  int
  end subroutine str2int
end module

PROGRAM stream
use str2int_mod
!     IMPLICIT NONE
!     .. Parameters ..
INTEGER n,offset,ndim,ntimes,maxthreads
PARAMETER (offset=0,ntimes=10) !
!     removed:
!     n = 1024000 and
!     ndim = n+offset

!     ..
!     .. Local Scalars ..
DOUBLE PRECISION scalar,t
INTEGER j,k,nbpw,quantum
!     ..
!     .. Local Arrays ..
DOUBLE PRECISION maxtime(4),mintime(4),avgtime(4),times(4,ntimes)
INTEGER bytes(4)
CHARACTER label(4)*11
!     ..
!     .. External Functions ..
DOUBLE PRECISION mysecond
INTEGER checktick,realsize
EXTERNAL mysecond,checktick,realsize
!$    INTEGER omp_get_num_threads
!$    EXTERNAL omp_get_num_threads
!     ..
!     .. Intrinsic Functions ..
!
INTRINSIC dble,max,min,nint,sqrt
!     ..
!     .. Arrays in Common ..
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: a,b,c
!     ..

!     .. Framework for command-line arguments ..
INTEGER status
CHARACTER*32 args(2)

!     .. Common blocks ..
!     COMMON a,b,c
!     ..
!     .. Data statements ..
DATA avgtime/4*0.0D0/,mintime/4*1.0D+36/,maxtime/4*0.0D0/
!DATA label/'Copy:      ','Scale:     ','Add:       ','Triad:     '/
DATA bytes/2,2,3,3/
!     ..

!       --- SETUP --- determine precision and check timing ---
!     allocate depending on arguments
!     first, we gotta get n, and calcualte ndim.
CALL getarg(1, args(1))

!     using str2int() to assign string version of args(1) to n
CALL str2int(args(1), n, status)

!     exit if invalid input. just to be safe.
if (status /= 0) then
  print *, "Usage: mystream.F90 nsize maxthreads"
  CALL EXIT(1)
end if

!     defining ndim
ndim = n + offset

!     allocating the arrays
ALLOCATE(a(ndim))
ALLOCATE(b(ndim))
ALLOCATE(c(ndim))

!     then, find number of cores.
CALL getarg(2, args(2))

!     using str2int() to assign string version of args(2) to maxthreads
CALL str2int(args(2), maxthreads, status)

!     exit if invalid input. just to be safe.
if (status /= 0) then
  print *, "Usage: mystream.F90 nsize maxthreads"
  CALL EXIT(1)
end if

CALL OMP_SET_NUM_THREADS(maxthreads) ! the user should enter the highest amount of threads possible; implied to be $(nproc) from generate.sh

nbpw = realsize()

!PRINT *,'----------------------------------------------'
!PRINT *,'STREAM Version $Revision: 5.6 $'
!PRINT *,'----------------------------------------------'
!WRITE (*,FMT=9010) 'Array size = ',n
!WRITE (*,FMT=9010) 'Offset     = ',offset
PRINT *, 3*nbpw*n/ (1024*1024)
!WRITE (*,FMT=9030) 'You are running each test ',ntimes,' times'
!WRITE (*,FMT=9030) '--'
!WRITE (*,FMT=9030) 'The *best! time for each test is used'
!WRITE (*,FMT=9030) '*EXCLUDING! the first and last iterations'

!$OMP PARALLEL
!$OMP MASTER
!PRINT *,'----------------------------------------------'
!PRINT *,'Number of Threads = ',OMP_GET_NUM_THREADS()
!$OMP END MASTER
!$OMP END PARALLEL

!PRINT *,'----------------------------------------------'
!$OMP PARALLEL
!PRINT *,'Printing one line per active thread....'
!$OMP END PARALLEL

!$OMP PARALLEL DO
DO 10 j = 1,n
a(j) = 2.0d0
b(j) = 0.5D0
c(j) = 0.0D0
10 CONTINUE
t = mysecond()
!$OMP PARALLEL DO
DO 20 j = 1,n
a(j) = 0.5d0*a(j)
20 CONTINUE
t = mysecond() - t
!PRINT *,'----------------------------------------------------'
quantum = checktick()
!WRITE (*,FMT=9000)'Your clock granularity/precision appears to be ',quantum,' microseconds'
!PRINT *,'----------------------------------------------------'

!       --- MAIN LOOP --- repeat test cases NTIMES times ---
scalar = 0.5d0*a(1)
DO 70 k = 1,ntimes

t = mysecond()
a(1) = a(1) + t
!$OMP PARALLEL DO
DO 30 j = 1,n
    c(j) = a(j)
30     CONTINUE
t = mysecond() - t
c(n) = c(n) + t
times(1,k) = t

t = mysecond()
c(1) = c(1) + t
!$OMP PARALLEL DO
DO 40 j = 1,n
    b(j) = scalar*c(j)
40     CONTINUE
t = mysecond() - t
b(n) = b(n) + t
times(2,k) = t

t = mysecond()
a(1) = a(1) + t
!$OMP PARALLEL DO
DO 50 j = 1,n
    c(j) = a(j) + b(j)
50     CONTINUE
t = mysecond() - t
c(n) = c(n) + t
times(3,k) = t

t = mysecond()
b(1) = b(1) + t
!$OMP PARALLEL DO
DO 60 j = 1,n
    a(j) = b(j) + scalar*c(j)
60     CONTINUE
t = mysecond() - t
a(n) = a(n) + t
times(4,k) = t
70 CONTINUE

!       --- SUMMARY ---
DO 90 k = 2,ntimes
DO 80 j = 1,4
    avgtime(j) = avgtime(j) + times(j,k)
    mintime(j) = min(mintime(j),times(j,k))
    maxtime(j) = max(maxtime(j),times(j,k))
80     CONTINUE
90 CONTINUE
!WRITE (*,FMT=9040)
DO 100 j = 1,4
avgtime(j) = avgtime(j)/dble(ntimes-1)
WRITE (*,FMT=9050) n*bytes(j)*nbpw/mintime(j)/1.0D6,avgtime(j),mintime(j),maxtime(j)
100 CONTINUE
!PRINT *,'----------------------------------------------------'
CALL checksums (a,b,c,n,ntimes)
!PRINT *,'----------------------------------------------------'

9000 FORMAT (1x,a,i6,a)
9010 FORMAT (1x,a,i10)
9020 FORMAT (1x,a,i4,a)
9030 FORMAT (1x,a,i3,a,a)
!9040 FORMAT ('Function',5x,'Rate (MB/s)  Avg time   Min time  Max time')
9050 FORMAT (4 (f12.4,2x))
END

!-------------------------------------
! INTEGER FUNCTION dblesize()
!
! A semi-portable way to determine the precision of DOUBLE PRECISION
! in Fortran.
! Here used to guess how many bytes of storage a DOUBLE PRECISION
! number occupies.
!
INTEGER FUNCTION realsize()
!     IMPLICIT NONE

!     .. Local Scalars ..
DOUBLE PRECISION result,test
INTEGER j,ndigits
!     ..
!     .. Local Arrays ..
DOUBLE PRECISION ref(30)
!     ..
!     .. External Subroutines ..
EXTERNAL confuse
!     ..
!     .. Intrinsic Functions ..
INTRINSIC abs,acos,log10,sqrt
!     ..

!       Test #1 - compare single(1.0d0+delta) to 1.0d0

10 DO 20 j = 1,30
ref(j) = 1.0d0 + 10.0d0** (-j)
20 CONTINUE

DO 30 j = 1,30
test = ref(j)
ndigits = j
CALL confuse(test,result)
IF (test.EQ.1.0D0) THEN
    GO TO 40
END IF
30 CONTINUE
GO TO 50

40 IF (ndigits.LE.8) THEN
realsize = 4
ELSE
realsize = 8
END IF
RETURN

50 PRINT *,'Hmmmm.  I am unable to determine the size.'
PRINT *,'Please enter the number of Bytes per DOUBLE PRECISION',' number : '
READ (*,FMT=*) realsize
IF (realsize.NE.4 .AND. realsize.NE.8) THEN
PRINT *,'Your answer ',realsize,' does not make sense.'
PRINT *,'Try again.'
PRINT *,'Please enter the number of Bytes per ','DOUBLE PRECISION number : '
READ (*,FMT=*) realsize
END IF
PRINT *,'You have manually entered a size of ',realsize,' bytes per DOUBLE PRECISION number'
WRITE (*,FMT='(a)') '----------------------------------------------'
END

SUBROUTINE confuse(q,r)
!     IMPLICIT NONE
!     .. Scalar Arguments ..
DOUBLE PRECISION q,r
!     ..
!     .. Intrinsic Functions ..
INTRINSIC cos
!     ..
r = cos(q)
RETURN
END

! A semi-portable way to determine the clock granularity
! Adapted from a code by John Henning of Digital Equipment Corporation
!
INTEGER FUNCTION checktick()
!     IMPLICIT NONE

!     .. Parameters ..
INTEGER n
PARAMETER (n=20)
!     ..
!     .. Local Scalars ..
DOUBLE PRECISION t1,t2
INTEGER i,j,jmin
!     ..
!     .. Local Arrays ..
DOUBLE PRECISION timesfound(n)
!     ..
!     .. External Functions ..
DOUBLE PRECISION mysecond
EXTERNAL mysecond
!     ..
!     .. Intrinsic Functions ..
INTRINSIC max,min,nint
!     ..
i = 0

10 t2 = mysecond()
IF (t2.EQ.t1) GO TO 10

t1 = t2
i = i + 1
timesfound(i) = t1
IF (i.LT.n) GO TO 10

jmin = 1000000
DO 20 i = 2,n
j = nint((timesfound(i)-timesfound(i-1))*1d6)
jmin = min(jmin,max(j,0))
20 CONTINUE

IF (jmin.GT.0) THEN
checktick = jmin
ELSE
!  PRINT *,'Your clock granularity appears to be less ','than one microsecond'
checktick = 1
END IF
RETURN

!      PRINT 14, timesfound(1)*1d6
!      DO 20 i=2,n
!         PRINT 14, timesfound(i)*1d6,
!     &       nint((timesfound(i)-timesfound(i-1))*1d6)
!   14    FORMAT (1X, F18.4, 1X, i8)
!   20 CONTINUE

END




SUBROUTINE checksums(a,b,c,n,ntimes)
!     IMPLICIT NONE
!     ..
!     .. Arguments ..
DOUBLE PRECISION a(*),b(*),c(*)
INTEGER n,ntimes
!     ..
!     .. Local Scalars ..
DOUBLE PRECISION aa,bb,cc,scalar,suma,sumb,sumc,epsilon
INTEGER k
!     ..

!     Repeat the main loop, but with scalars only.
!     This is done to check the sum & make sure all
!     iterations have been executed correctly.

aa = 2.0D0
bb = 0.5D0
cc = 0.0D0
aa = 0.5D0*aa
scalar = 0.5d0*aa
DO k = 1,ntimes
cc = aa
bb = scalar*cc
cc = aa + bb
aa = bb + scalar*cc
END DO
aa = aa*DBLE(n-2)
bb = bb*DBLE(n-2)
cc = cc*DBLE(n-2)

!     Now sum up the arrays, excluding the first and last
!     elements, which are modified using the timing results
!     to confuse aggressive optimizers.

suma = 0.0d0
sumb = 0.0d0
sumc = 0.0d0
!$OMP PARALLEL DO REDUCTION(+:suma,sumb,sumc)
DO 110 j = 2,n-1
suma = suma + a(j)
sumb = sumb + b(j)
sumc = sumc + c(j)
110 CONTINUE

epsilon = 1.D-6

IF (ABS(suma-aa)/suma .GT. epsilon) THEN
PRINT *,'Failed Validation on array a()'
PRINT *,'Target   Sum of a is = ',aa
PRINT *,'Computed Sum of a is = ',suma
ELSEIF (ABS(sumb-bb)/sumb .GT. epsilon) THEN
PRINT *,'Failed Validation on array b()'
PRINT *,'Target   Sum of b is = ',bb
PRINT *,'Computed Sum of b is = ',sumb
ELSEIF (ABS(sumc-cc)/sumc .GT. epsilon) THEN
PRINT *,'Failed Validation on array c()'
PRINT *,'Target   Sum of c is = ',cc
PRINT *,'Computed Sum of c is = ',sumc
!ELSE
!  PRINT *,'Solution Validates!'
ENDIF

END

