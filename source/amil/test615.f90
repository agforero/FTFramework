PROGRAM subl1
!****************************************************************
!      *this code solves the l1 best subset problem*
!      *submitted to acm trans.spring 1982*
!       *last revised mar. 12, 1984*
!****************************************************************

USE toms615
IMPLICIT NONE

INTEGER   :: i, idex(210), ifault, istat(20), iter, j, jmin, k, l,  &
             minin, mmax, n, nmax, nprob
REAL (dp) :: bval(210), popt, x(300,20), y(300), zl(20)

INTEGER   :: unitgo, unitin

!-----unitin    is the input unit number
!-----unitgo    is the output unit number

CHARACTER (LEN=72) :: imt
CHARACTER (LEN=40) :: name

!-----imt       is used to store the input format
!-----name      saves the problem name

unitin = 5
unitgo = 6
nmax = 300
mmax = 30
nprob = 0

!-----k        stores the number of parameters in the full model

10 READ (unitin, 70, END=60) k, name
WRITE (unitgo, 70) k, name
IF (k == 0) STOP
WRITE (unitgo, 80) name

!-----n          stores the number of observations
!-----minin      stores the minimum number of parameters considered
!-----popt       deviation from optimality allowed

READ (unitin, 140) n, minin, popt
20   FORMAT ('    Number of observation =', i5/,  &
             '    Number of parameters =', i5/  &
             '    Minimum number of parameters considered =', i5/  &
             '    Percentage deviation from optimality allowed =', f6.2/  &
             '        ** best subset lav program')
WRITE (unitgo, 20) n, k, minin, popt
READ (unitin, 150) imt
DO i = 1, n
  READ (unitin, FMT=imt) y(i), (x(i,j),j=1,k)
  
!        write (unitgo, FMT=imt) y(i), (x(i,j),j=1,k)
  
END DO

READ (unitin, 90) istat(1:k)

iter = 0

!******call the routine to find the best subsets


CALL kbest (x, y, k, n, iter, ifault, popt, minin, nmax, mmax,  &
            bval, idex, istat, zl)

WRITE (unitgo, 100) ifault

!     write the final best subset solution

j = k*(k+1)/2
jmin = (minin-1)*(minin)/2
j = j - jmin
DO i = minin,k
  WRITE (unitgo, 130) i
  WRITE (unitgo, 120) zl(i)
  DO l = 1,i
    WRITE (unitgo, 110) idex(j), bval(j)
    j = j - 1
  END DO
END DO

nprob = nprob + 1
WRITE (unitgo,160) iter
GO TO 10

60 STOP

70   FORMAT (i5, A)
80   FORMAT ('        ','problem title ', '  ', A)
90   FORMAT (10I2)
100  FORMAT (' ifault=', i3)
110  FORMAT (/' beta(', i3, ')', f15.3)
120  FORMAT ('    ', 14('*'), 'sum of absolute values = ', f15.3)
130  FORMAT (//'   best results for lav subset of size = ', i3)
140  FORMAT (i5, i2, f6.2)
150  FORMAT (A)
160  FORMAT (//'     iteration count = ', i7)

END PROGRAM subl1
