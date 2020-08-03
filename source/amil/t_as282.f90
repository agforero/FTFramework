!----------------------------------------------------------------------------

! Calling program for MVE/LMS subroutine. The input data has predictors first
! as columns, with target variable following, if this is an LMS problem.  The
! first column can be a case number if desired.   Current maximum dimensions
! are 200 observations and 20 variables, but these can be increased by
! changing the dimensions in lines 1-4.  The data are assumed to be in free
! format.

PROGRAM test_as282

USE robust_regression

IMPLICIT NONE
REAL (dp)          :: x(300,20), y(300), coeffs(20,100), eprmin(100),  &
                      resid(300,100), cvemin(100), robdsq(300,100), scalf
INTEGER            :: i, j, ncas, npre, ioptn, caseno,  &
                      int, nvar, maxtry, nclow, nchigh, ifault, ii
LOGICAL            :: mve, lms, exh, isint, case
CHARACTER (LEN=60) :: infile,outfile

WRITE(*,*) 'Input data file name?'
READ(*,10) infile
10 FORMAT(a)
WRITE(*,*) 'Does the data include a case number? Yes (T) or No',' (F)'
READ(*,*) case
OPEN (20,FILE=infile)
WRITE(*,*) 'Output file name?'
READ(*,10) outfile
OPEN(30,FILE=outfile)
WRITE(*,*) 'Sample size?'
READ(*,*) ncas
WRITE(*,*) '# predictors?'
READ(*,*) npre
WRITE(*,*) 'Enter scaling factor'
READ(*,*) scalf
ioptn = 0
WRITE(*,*) 'Do LMS? enter T or F'
READ(*,*) lms
IF (.NOT. lms) ioptn = ioptn + 1
WRITE(*,*) 'Do MVE? enter T or F'
READ(*,*) mve
IF (.NOT. mve) ioptn = ioptn + 2
WRITE(*,*) 'Intercept (constant) term? enter T or F'
READ(*,*) isint
IF (isint) THEN
  int=1
ELSE
  int=0
  ioptn = ioptn + 4
END IF
nvar=npre+int
WRITE(*,*) 'Exhaustive search? enter T or F'
READ(*,*) exh
IF (.NOT.exh) THEN
  ioptn = ioptn + 8
  WRITE(*,*) 'Max number of subsamples to check?'
  READ(*,*) maxtry
END IF
WRITE(*,*) 'Low, high coverages to check?'
READ(*,*) nclow,nchigh
DO i=1, ncas
  y(i) = 0
  IF (case) THEN
    IF (lms) THEN
      READ(20,*) caseno, (x(i,j),j=1,npre),y(i)
    ELSE
      READ(20,*) caseno, (x(i,j),j=1,npre)
    END IF
  ELSE
    IF (lms) THEN
      READ(20,*) (x(i,j),j=1,npre), y(i)
    ELSE
      READ(20,*) (x(i,j),j=1,npre)
    END IF
  END IF
  y(i) = y(i) * scalf
  x(i,1:npre) = scalf * x(i,1:npre)
END DO

IF (case) WRITE(*, *) 'Last case number was: ', caseno

CALL mvelms (x, y, ncas, npre, ioptn, maxtry, nclow, nchigh, coeffs,  &
             eprmin, resid, robdsq, cvemin, ifault)
WRITE(30,*) 'Fault parameter ', ifault
DO i=1, nchigh-nclow+1
  WRITE(30,*) 'Coverage level ', i+nclow-1
  IF (lms) THEN
    WRITE(30,*) 'LMS results'
    WRITE(30,*) 'Coefficients of LMS'
    WRITE(30,101) (coeffs(j,i),j=1,nvar)
    WRITE(30,*) 'Criterion ',eprmin(i)
    WRITE(30,*) 'Data, residuals'
    DO ii=1, ncas
      WRITE(30,103) ii, (x(ii,j),j=1,npre), y(ii), resid(ii,i)
    END DO
  END IF
  IF (mve) THEN
    WRITE(30,*) 'MVE results'
    WRITE(30,*) 'Criterion ', cvemin(i)
    WRITE(30,*) 'Data, robust d-squareds'
    DO ii=1, ncas
      WRITE(30,103) ii, (x(ii,j),j=1,npre), robdsq(ii,i)
    END DO
  END IF
END DO
101 FORMAT(5G15.6)
103 FORMAT(i3, (t5, 5G15.6))
STOP
END PROGRAM test_as282
