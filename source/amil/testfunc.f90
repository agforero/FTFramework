PROGRAM testfunktionen
! Testing the global optimization program of Tibor Csendes using
! the 30 test functions from:
!   http://solon.mat.univie.ac.at/~vpk/math/funcs.html

! ** This is not complete yet  **
! It takes a long while to program all 30 functions + data + starting values

USE global_minimum
IMPLICIT NONE
INTEGER, PARAMETER    :: dp = SELECTED_REAL_KIND(12, 60)

INTEGER, ALLOCATABLE  :: seed(:)       ! For the random number generator

INTEGER, PARAMETER    :: ipr = 9

REAL (dp)             :: amin(30), amax(30), x0(30,20), f00(20)
INTEGER               :: i, k, nparm, m, nsampl, ng, nsig, nc
CHARACTER (LEN=30)    :: heading

OPEN(UNIT=ipr, FILE='testfunc.out')

! Set the random number seed

CALL RANDOM_SEED(size=k)
ALLOCATE (seed(k))
WRITE(*, '(1x, a, i4, a)') 'Enter ', k, ' integers as random number seeds: '
READ(*, *) seed
CALL RANDOM_SEED(put=seed)
WRITE(ipr, '(a / (7i11) )') ' Random number seed(s): ', seed
WRITE(ipr, * )

DO m = 1, 30
  CALL start_values()
  WRITE(ipr, '(/" ", a)') 'Test problem: ' // heading
  nsampl = MIN(100*nparm, 1000)
  ng = 5
  nsig = 6
  CALL global(amin, amax, nparm, m, nsampl, ng, ipr, nsig, x0, nc, f00)
  DO i = 1, nc
    WRITE(ipr, '(a, g13.5)') ' Function value = ', f00(i)
    WRITE(ipr, '(a / (6g13.5))') ' Parameter values = ', x0(1:nparm, i)
  END DO
  WRITE(ipr, '(//)')
END DO

STOP


CONTAINS


SUBROUTINE start_values()

SELECT CASE ( m )
  CASE ( 1 )
    heading = 'Sphere'
    nparm = 3
    amin(1:nparm) = -5.12_dp
    amax(1:nparm) =  5.12_dp
  CASE ( 2 )
    heading = 'Rosenbrock'
    nparm = 10
    amin(1:nparm) = -2.048_dp
    amax(1:nparm) =  2.048_dp
  CASE ( 3 )
    heading = 'Step function'
    nparm = 5
    amin(1:nparm) = -5.12_dp
    amax(1:nparm) =  5.12_dp
  CASE ( 4 )
    heading = 'Hyper-Ellipsoid '
    nparm = 30
    amin(1:nparm) = -1.0_dp
    amax(1:nparm) =  1.0_dp
  CASE ( 5 )
    heading = 'Neumaier nr.3'
    nparm = 30
    amin(1:nparm) = -nparm**2
    amax(1:nparm) =  nparm**2
  CASE ( 6 )
    heading = 'Schaffer nr.2'
    nparm = 2
    amin(1:nparm) = -100._dp
    amax(1:nparm) =  100._dp
  CASE ( 7 )
    heading = 'Shekel'
    nparm = 10
    amin(1:nparm) =   0.0_dp
    amax(1:nparm) =  10.0_dp
  CASE ( 8 )
    heading = 'Hartman'
    nparm = 6
    amin(1:nparm) =  0.0_dp
    amax(1:nparm) =  1.0_dp
  CASE ( 9 )
    heading = 'Goldstein-Price'
    nparm = 2
    amin(1:nparm) =  0.0_dp
    amax(1:nparm) =  1.0_dp
  CASE ( 10 )
    heading = 'Branin'
    nparm = 2
    amin(1) = -5.0_dp
    amax(1) = 10.0_dp
    amin(2) =  0.0_dp
    amax(2) = 15.0_dp
  CASE ( 11 )
    heading = 'Six-hump camel back'
    nparm = 2
    amin(1:nparm) =  0.0_dp
    amax(1:nparm) =  1.0_dp
  CASE ( 12 )
    heading = 'Shubert'
    nparm = 2
    amin(1:nparm) = -10._dp
    amax(1:nparm) =  10._dp
  CASE ( 13 )
    heading = 'Shekels Foxholes'
    nparm = 2
    amin(1:nparm) = -65.536_dp
    amax(1:nparm) =  65.536_dp
  CASE ( 14 )
    heading = 'Modifizierte Shekels Foxholes'
    nparm = 10
    amin(1:nparm) =  0._dp
    amax(1:nparm) = 10._dp
  CASE ( 15 )
    heading = 'Coranas Parabel'
    nparm = 4
    amin(1:nparm) = -1000._dp
    amax(1:nparm) =  1000._dp
  CASE ( 16 )
    heading = 'Griewanke'
    nparm = 10
    amin(1:nparm) = -400._dp
    amax(1:nparm) =  400._dp
  CASE ( 17 )
    heading = 'Zimmermann'
    nparm = 2
    amin(1:nparm) =   0._dp
    amax(1:nparm) = 100._dp
  CASE ( 18 )
    heading = 'Katsuuras'
    nparm = 10
    amin(1:nparm) =    0._dp
    amax(1:nparm) = 1000._dp
  CASE ( 19 )
    heading = 'Rastringins'
    nparm = 20
    amin(1:nparm) = -600._dp
    amax(1:nparm) =  600._dp
  CASE ( 20 )
    heading = 'Ackleys'
    nparm = 30
    amin(1:nparm) = -30._dp
    amax(1:nparm) =  30._dp
  CASE ( 21 )
    heading = 'Schaffer nr.1'
    nparm = 2
    amin(1:nparm) = -100._dp
    amax(1:nparm) =  100._dp
  CASE ( 22 )
    heading = 'Easom'
    nparm = 2
    amin(1:nparm) = -100._dp
    amax(1:nparm) =  100._dp
  CASE ( 23 )
    heading = 'Bohachevsky nr.1'
    STOP
  CASE ( 24 )
    heading = 'Bohachevsky nr.2'
  CASE ( 25 )
    heading = 'Colville'
  CASE ( 26 )
    heading = 'Neumaier nr.2'
  CASE ( 27 )
    heading = 'Modifizierte Langerman'
  CASE ( 28 )
    heading = 'Epistatic Michalewicz'
  CASE ( 29 )
    heading = 'Odd Square'
  CASE ( 30 )
    heading = 'Chebychev-Polynom'
END SELECT

RETURN
END SUBROUTINE start_values

END PROGRAM testfunktionen



SUBROUTINE funct(x, value, nparm, m)

IMPLICIT NONE
INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(12, 60)

REAL (dp), INTENT(IN)   :: x(:)
REAL (dp), INTENT(OUT)  :: value
INTEGER, INTENT(IN)     :: nparm
INTEGER, INTENT(IN)     :: m          ! this is the function number

INTEGER   :: i, j
REAL (dp) :: temp, z
REAL (dp), PARAMETER  :: zero = 0.0_dp, one = 1.0_dp,  &
                         pi = 3.141592653589793846_dp

! Data for Shekel test function
REAL (dp), PARAMETER  :: c(10) = (/ 0.1_dp, 0.2_dp, 0.2_dp, 0.4_dp, 0.4_dp, &
                                    0.6_dp, 0.3_dp, 0.7_dp, 0.5_dp, 0.5_dp /)
REAL (dp), PARAMETER  :: a(10,4) = RESHAPE(  &
  (/ 4._dp, 1._dp, 8._dp, 6._dp, 3._dp, 2._dp, 5._dp, 8._dp, 6._dp, 7._dp, &
     4._dp, 1._dp, 8._dp, 6._dp, 7._dp, 9._dp, 5._dp, 1._dp, 2._dp, 3.6_dp, &
     4._dp, 1._dp, 8._dp, 6._dp, 3._dp, 2._dp, 3._dp, 8._dp, 6._dp, 7._dp, &
     4._dp, 1._dp, 8._dp, 6._dp, 7._dp, 9._dp, 3._dp, 1._dp, 2._dp, 3.6_dp /), &
     (/ 10, 4 /) )

! Data for Hartman test function
REAL (dp), PARAMETER  :: ha(6,4) = RESHAPE(  &
  (/ 10._dp,   3._dp, 17._dp, 3.5_dp, 1.7_dp, 8._dp,  &
     0.05_dp, 10._dp, 17._dp, 0.1_dp, 8._dp, 14._dp,  &
     3._dp,  3.5_dp, 1.7_dp,  10._dp, 17._dp, 8._dp,  &
     17._dp,  8._dp, 0.05_dp, 10._dp, 0.1_dp, 14._dp /), (/ 6, 4 /) )
REAL (dp), PARAMETER  :: hc(4) = (/ 1.0_dp, 1.2_dp, 3.0_dp, 3.2_dp /)
REAL (dp), PARAMETER  :: hp(6,4) = RESHAPE(  &
  (/ 0.1312_dp, 0.1696_dp, 0.5569_dp, 0.0124_dp, 0.8283_dp, 0.5886_dp,  &
     0.2329_dp, 0.4135_dp, 0.8307_dp, 0.3736_dp, 0.1004_dp, 0.9991_dp,  &
     0.2348_dp, 0.1451_dp, 0.3522_dp, 0.2883_dp, 0.3047_dp, 0.6650_dp,  &
     0.4047_dp, 0.8828_dp, 0.8732_dp, 0.5743_dp, 0.1091_dp, 0.0381_dp /), &
     (/ 6, 4 /) )

! Data for Shekels Foxholes
REAL (dp), PARAMETER  :: sh(25,2) = RESHAPE(  &
  (/ -32._dp, -16._dp, 0._dp, 16._dp, 32._dp,  &
     -32._dp, -16._dp, 0._dp, 16._dp, 32._dp,  &
     -32._dp, -16._dp, 0._dp, 16._dp, 32._dp,  &
     -32._dp, -16._dp, 0._dp, 16._dp, 32._dp,  &
     -32._dp, -16._dp, 0._dp, 16._dp, 32._dp,  &
     -32._dp, -32._dp, -32._dp, -32._dp, -32._dp,  &
     -16._dp, -16._dp, -16._dp, -16._dp, -16._dp,  &
       0._dp,   0._dp,   0._dp,   0._dp,   0._dp,  &
     +16._dp, +16._dp, +16._dp, +16._dp, +16._dp,  &
     +32._dp, +32._dp, +32._dp, +32._dp, +32._dp /), (/ 25, 2 /) )

! Data for Coranas Parabel
REAL (dp), PARAMETER  :: d(4) = (/ 1._dp, 1000._dp, 10._dp, 100._dp /)

SELECT CASE (m )
  CASE ( 1 )                 ! Sphere
    value = SUM( x(1:nparm)**2 )

  CASE ( 2 )                 ! Rosenbrock
    value = zero
    DO i = 1, nparm-1
      value = value + 100._dp*(x(i+1)**2 - x(i))**2 + (one - x(i))**2
    END DO

  CASE ( 3 )                 ! Step function
    value = SUM( FLOOR(x(1:nparm)) )

  CASE ( 4 )                 ! Hyper-Ellipsoid
    value = zero
    DO i = 1, nparm
      value = value + (i*x(i))**2
    END DO

  CASE ( 5 )                 ! Neumaier nr.3
    value = SUM( (x(1:nparm) - one)**2 )
    DO i = 2, nparm
      value = value - x(i)*x(i-1)
    END DO

  CASE ( 6 )                 ! Schaffer nr.2
    temp = x(1)**2 + x(2)**2
    value = temp**0.25 * (2500._dp*temp**0.2_dp + one)

  CASE ( 7 )                 ! Shekel
    value = zero
    DO i = 1, nparm
      value = value - one / (c(i) + SUM( (x(1:4) - a(i,1:4))**2 ))
    END DO

  CASE ( 8 )                 ! Hartman (n = 6)
    value = zero
    DO i = 1, 4
      value = value - hc(i) *   &
              EXP(-DOT_PRODUCT( ha(1:nparm,i), (x(1:nparm)-hp(1:nparm,i))**2 ))
    END DO

  CASE ( 9 )                 ! Goldstein-Price
    value = (one + (x(1)+x(2)+one)**2 *  &
                   (19-14*x(1)+3*x(1)**2-14*x(2)+6*x(1)*x(2)+3*x(2)**2)) *  &
            (30 + (2*x(1)-3*x(2))**2 *  &
                  (18-32*x(1)+12*x(1)**2+48*x(2)-36*x(1)*x(2)+27*x(2)**2))

  CASE ( 10 )                ! Branin
    value = (x(2) - 1.25*(x(1)/pi)**2 + 5*x(1)/pi - 6)**2 +  &
            10*(one - 0.125_dp/pi)*COS(x(1)) + 10._dp

  CASE ( 11 )                ! Six-hump camel back
    value = (4 - 2.1_dp*x(1)**2 + x(1)**4/3)*x(1)**2 + x(1)*x(2) +  &
            4*(x(2)**2 - one)*x(2)**2

  CASE ( 12 )                ! Shubert
    value = one
    DO i = 1, 2
      temp = zero
      DO j = 1, 5
        temp = temp + j*COS((j+1)*x(i) + j)
      END DO
      value = value * temp
    END DO

  CASE ( 13 )                ! Shekels Foxholes
    temp = zero
    DO i = 1, 25
      temp = temp + one / (i + (x(1)-sh(i,1))**6 + (x(2)-sh(i,2))**6)
    END DO
    value = one / (0.002_dp + temp)

  CASE ( 14 )                ! Modifizierte Shekels Foxholes
                             ! Using same pattern for A-matrix as for F7,
                             ! and 3 repeats of the c-vector.
    value = zero
    DO i = 1, 30
      temp = c( MOD(i-1,10) + 1 )
      DO j = 1, nparm
        temp = temp + (x(j) - a(j, MOD(i-1,4) + 1))**2
      END DO
      value = value - one / temp
    END DO

  CASE ( 15 )                ! Coranas Parabel
    value = zero
    DO i = 1, nparm
      z = FLOOR( ABS(5*x(i)) + 0.49999_dp )
      z = SIGN(z, x(i))*0.2_dp
      IF ( ABS(x(i)-z) < 0.05_dp ) THEN
        value = value + 0.15_dp * (z - SIGN(0.05_dp, z))**2 * d(i)
      ELSE
        value = value + d(i)*x(i)**2
      END IF
    END DO

  CASE ( 16 )                ! Griewanke
    value = one + SUM( x(1:nparm)**2 ) / 4000.
    temp = one
    DO i = 1, nparm
      temp = temp * COS(x(i) / SQRT(DBLE(i)))
    END DO
    value = value - temp

  CASE ( 17 )                ! Zimmermann
    value = MAX( h1(x), SIGN( p(h2(x)), h2(x) ), SIGN( p(h3(x)), h3(x) ),  &
                 SIGN( p(-x(1)), x(1) ), SIGN( p(-x(2)), x(2) ) )

  CASE ( 18 )                ! Katsuuras
    value = one
    DO i = 1, nparm
      temp = zero
      DO j = 1, 32
        temp = temp + abstoint(x(i), j)
      END DO
      value = value * (one + i*temp)
    END DO

  CASE ( 19 )                ! Rastringins
    value = 10*nparm
    DO i = 1, nparm
      value = value + x(i)**2 - 10*COS(2*pi*x(i))
    END DO

  CASE ( 20 )                ! Ackleys
    value = -20*EXP( -0.02_dp*SQRT( SUM( x(1:nparm)**2 ) / nparm ) ) -  &
            EXP( SUM( COS(2*pi*x(1:nparm)) ) / nparm ) + 20._dp

  CASE ( 21 )                ! Schaffer nr.1
    temp = x(1)**2 + x(2)**2
    value = 0.5 + ( SIN(SQRT(temp))**2 - 0.5 ) / (one + 0.001_dp*temp)**2

  CASE ( 22 )                ! Easom
    value = -COS(x(1)) * COS(x(2)) * EXP( -(x(1)-pi)**2 - (x(2)-pi)**2 )

  CASE ( 23 )                ! Bohachevsky nr.1

  CASE ( 24 )                ! Bohachevsky nr.2

  CASE ( 25 )                ! Colville

  CASE ( 26 )                ! Neumaier nr.2

  CASE ( 27 )                ! Modifizierte Langerman

  CASE ( 28 )                ! Epistatic Michalewicz

  CASE ( 29 )                ! Odd Square

  CASE ( 30 )                ! Chebychev-Polynom
END SELECT

RETURN


CONTAINS


FUNCTION h1(x) RESULT(fn_val)
! Function required for Zimmermann's Problem (F17)

REAL (dp), INTENT(IN)  :: x(:)
REAL (dp)              :: fn_val

fn_val = 9._dp - x(1) - x(2)

RETURN
END FUNCTION h1



FUNCTION h2(x) RESULT(fn_val)
! Function required for Zimmermann's Problem (F17)

REAL (dp), INTENT(IN)  :: x(:)
REAL (dp)              :: fn_val

fn_val = (x(1) - 3._dp)**2 + (x(2) - 2._dp)**2 - 16._dp

RETURN
END FUNCTION h2



FUNCTION h3(x) RESULT(fn_val)
! Function required for Zimmermann's Problem (F17)

REAL (dp), INTENT(IN)  :: x(:)
REAL (dp)              :: fn_val

fn_val = x(1)*x(2) - 14._dp

RETURN
END FUNCTION h3



FUNCTION p(delta) RESULT(fn_val)
! Function required for Zimmermann's Problem (F17)

REAL (dp), INTENT(IN)  :: delta
REAL (dp)              :: fn_val

fn_val = 100*(one + delta)

RETURN
END FUNCTION p



FUNCTION abstoint(x, j) RESULT(fn_val)

REAL (dp), INTENT(IN)  :: x
INTEGER, INTENT(IN)    :: j
REAL (dp)              :: fn_val

INTEGER, PARAMETER  :: large = HUGE(1)

fn_val = INT( MIN( SCALE(x, j), DBLE(large) ) )
fn_val = SCALE(fn_val, -j)

RETURN
END FUNCTION abstoint

END SUBROUTINE funct
