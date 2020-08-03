FUNCTION dsinint(xvalue) RESULT(fn_val)

!   DEFINITION:

!   This program calculates the value of the sine-integral defined as

!       DSININT(x) = Integral (0 to x) sin(t)/t  dt

!   The program uses rational approximations with the coefficients
!   given to 20sf. accuracy.


!   INPUT PARAMETER:

!     XVALUE - DOUBLE PRECISION - The argument to the function


!   MACHINE-DEPENDENT PARAMETERS:

!     XLOW - DOUBLE PRECISION - The absolute value below which
!                                    DSININT( x ) = x,
!                               to machine precision.
!                               The recommended value is
!                                    SQRT(18*EPSNEG)

!     XHIGH1 - DOUBLE PRECISION - The value above which
!                                DSININT(x) = pi/2 - cos(x)/x -sin(x)/x*x
!                                 to machine precision.
!                                 The recommended value is
!                                   SQRT(6/EPSNEG)

!     XHIGH2 - DOUBLE PRECISION - The value above which
!                                     DSININT(x) = pi/2
!                                 to machine precision.
!                                 The recommended value is
!                                     2 / min(EPS,EPSNEG)

!     XHIGH3 - DOUBLE PRECISION - The value above which it is not sensible
!                                 to compute COS or SIN. The recommended
!                                 value is     pi/EPS


!      Values of EPS and EPSNEG for certain machine/compiler
!      combinations can be found in the paper

!      W.J. CODY  Algorithm 665: MACHAR: A subroutine to dynamically
!      determine machine parameters, ACM Trans. Math. Soft. 14 (1988) 303-311.

!      The current code gives numerical values for XLOW,XHIGH1,XHIGH2,XHIGH3,
!      suitable for machines whose arithmetic conforms to the IEEE
!      standard. The codes will probably work on other machines
!      but might overflow or underflow for certain arguments.

!   AUTHOR: Allan MacLeod
!           Dept. of Mathematics and Statistics
!           University of Paisley
!           Scotland
!           (e-mail: macl_ms0@paisley.ac.uk)

IMPLICIT NONE
INTEGER, PARAMETER    :: dp = SELECTED_REAL_KIND(14, 60)
REAL (dp), INTENT(IN) :: xvalue
REAL (dp)             :: fn_val

! Local variables

INTEGER   :: i, indsgn
REAL (dp) :: cx, fival, gival, sumden, sumnum, sx, x, xhigh, xsq

!  DATA VALUES

REAL (dp), PARAMETER :: zero = 0.0_dp, one = 1.0_dp, six = 6.0_dp, &
                        twelve = 12.0_dp
REAL (dp), PARAMETER :: piby2 = 1.5707963267948966192_dp

!  MACHINE-DEPENDENT PARAMETERS (SUITABLE FOR IEEE MACHINES)

REAL (dp), PARAMETER :: xlow = 4.47E-8_dp, xhigh1 = 2.32472E8_dp
REAL (dp), PARAMETER :: xhigh2 = 9.0072E15_dp, xhigh3 = 1.4148475E16_dp

!  VALUES FOR SINE-INTEGRAL FOR 0 <= |X| <= 6

REAL (dp), PARAMETER :: asintn(0:7) = (/ 1.0_dp,  &
          -0.44663998931312457298E-1_dp, 0.11209146443112369449E-2_dp,  &
          -0.13276124407928422367E-4_dp, 0.85118014179823463879E-7_dp,  &
          -0.29989314303147656479E-9_dp, 0.55401971660186204711E-12_dp, &
          -0.42406353433133212926E-15_dp /)
REAL (dp), PARAMETER :: asintd(0:7) = (/ 1.0_dp,  &
           0.10891556624243098264E-1_dp, 0.59334456769186835896E-4_dp,  &
           0.21231112954641805908E-6_dp, 0.54747121846510390750E-9_dp,  &
           0.10378561511331814674E-11_dp, 0.13754880327250272679E-14_dp,&
           0.10223981202236205703E-17_dp /)

!  VALUES FOR FI(X) FOR 6 <= X <= 12

REAL (dp), PARAMETER :: afn1(0:7) = (/ 0.99999999962173909991_dp,   &
          0.36451060338631902917E3_dp, 0.44218548041288440874E5_dp, &
          0.22467569405961151887E7_dp, 0.49315316723035561922E8_dp, &
          0.43186795279670283193E9_dp, 0.11847992519956804350E10_dp,&
          0.45573267593795103181E9_dp /)
REAL (dp), PARAMETER :: afd1(0:7) = (/ 1.0_dp, 0.36651060273229347594E3_dp,  &
                  0.44927569814970692777E5_dp, 0.23285354882204041700E7_dp,  &
                  0.53117852017228262911E8_dp, 0.50335310667241870372E9_dp,  &
                  0.16575285015623175410E10_dp, 0.11746532837038341076E10_dp /)

!   VALUES OF GI(X) FOR 6 <= X <=12

REAL (dp), PARAMETER :: agn1(0:8) = (/ 0.99999999920484901956_dp,     &
          0.51385504875307321394E3_dp, 0.92293483452013810811E5_dp,   &
          0.74071341863359841727E7_dp, 0.28142356162841356551E9_dp,   &
          0.49280890357734623984E10_dp, 0.35524762685554302472E11_dp, &
          0.79194271662085049376E11_dp, 0.17942522624413898907E11_dp /)
REAL (dp), PARAMETER :: agd1(0:8) = (/ 1.0_dp, 0.51985504708814870209E3_dp,  &
                  0.95292615508125947321E5_dp, 0.79215459679762667578E7_dp,  &
                  0.31977567790733781460E9_dp, 0.62273134702439012114E10_dp,  &
                  0.54570971054996441467E11_dp, 0.18241750166645704670E12_dp,  &
                  0.15407148148861454434E12_dp /)

!   VALUES FOR FI(X) FOR X > 12

REAL (dp), PARAMETER :: afn2(0:7) = (/ 0.19999999999999978257E1_dp,   &
          0.22206119380434958727E4_dp, 0.84749007623988236808E6_dp,   &
          0.13959267954823943232E9_dp, 0.10197205463267975592E11_dp,  &
          0.30229865264524075951E12_dp, 0.27504053804288471142E13_dp, &
          0.21818989704686874983E13_dp /)
REAL (dp), PARAMETER :: afd2(0:7) = (/ 1.0_dp, 0.11223059690217167788E4_dp,  &
                  0.43685270974851313242E6_dp, 0.74654702140658116258E8_dp,  &
                  0.58580034751805687471E10_dp, 0.20157980379272098841E12_dp,  &
                  0.26229141857684496445E13_dp, 0.87852907334918467516E13_dp /)

!   VALUES FOR GI(X) FOR X > 12

REAL (dp), PARAMETER :: agn2(0:8) = (/ 0.59999999999999993089E1_dp,   &
          0.96527746044997139158E4_dp, 0.56077626996568834185E7_dp,   &
          0.15022667718927317198E10_dp, 0.19644271064733088465E12_dp, &
          0.12191368281163225043E14_dp, 0.31924389898645609533E15_dp, &
          0.25876053010027485934E16_dp, 0.12754978896268878403E16_dp /)
REAL (dp), PARAMETER :: agd2(0:8) = (/ 1.0_dp, 0.16287957674166143196E4_dp,  &
                  0.96636303195787870963E6_dp, 0.26839734750950667021E9_dp,  &
                  0.37388510548029219241E11_dp, 0.26028585666152144496E13_dp,  &
                  0.85134283716950697226E14_dp, 0.11304079361627952930E16_dp,  &
                  0.42519841479489798424E16_dp /)

!   START COMPUTATION

x = xvalue
indsgn = 1
IF ( x < zero ) THEN
  x = -x
  indsgn = -1
END IF

!   CODE FOR 0 <= |X| <= 6

IF ( x <= six ) THEN
  IF ( x < xlow ) THEN
    fn_val = x
  ELSE
    sumnum = zero
    sumden = zero
    xsq = x * x
    DO i = 7 , 0 , -1
      sumnum = sumnum * xsq + asintn(i)
      sumden = sumden * xsq + asintd(i)
    END DO
    fn_val = x * sumnum / sumden
  END IF
END IF

!   CODE FOR 6 < |X| <= 12

IF ( x > six .AND. x <= twelve ) THEN
  sumnum = zero
  sumden = zero
  xsq = one / ( x * x )
  DO i = 7 , 0 , -1
    sumnum = sumnum * xsq + afn1(i)
    sumden = sumden * xsq + afd1(i)
  END DO
  fival = sumnum / ( x * sumden )
  sumnum = zero
  sumden = zero
  DO i = 8 , 0 , -1
    sumnum = sumnum * xsq + agn1(i)
    sumden = sumden * xsq + agd1(i)
  END DO
  gival = xsq * sumnum / sumden
  fn_val = piby2 - fival * COS(x) - gival * SIN(x)
END IF

!   CODE FOR |X| > 12

IF ( x > twelve ) THEN
  xhigh = MIN(xhigh2, xhigh3)
  IF ( x > xhigh ) THEN
    fn_val = piby2
  ELSE
    cx = COS(x)
    sx = SIN(x)
    xsq = one / ( x * x )
    IF ( x > xhigh1 ) THEN
      fn_val = piby2 - cx / x - sx * xsq
    ELSE
      sumnum = zero
      sumden = zero
      DO i = 7 , 0 , -1
        sumnum = sumnum * xsq + afn2(i)
        sumden = sumden * xsq + afd2(i)
      END DO
      fival =  ( one - xsq * sumnum / sumden ) / x
      sumnum = zero
      sumden = zero
      DO i = 8 , 0 , -1
        sumnum = sumnum * xsq + agn2(i)
        sumden = sumden * xsq + agd2(i)
      END DO
      gival =  ( one - xsq * sumnum / sumden ) * xsq
      fn_val = piby2 - cx * fival - sx * gival
    END IF
  END IF
END IF
IF ( indsgn == -1 ) fn_val = -fn_val
RETURN
END FUNCTION dsinint
