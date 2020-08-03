FUNCTION dcosint(xvalue) RESULT(fn_val)

!   This program calculates the value of the cosine-integral
!   defined as

!   DCOSINT(x) = Gamma + Ln(x) + Integral (0 to x) [cos(t)-1]/t  dt

!                where Gamma is Euler's constant.

!   The code uses rational approximations with a maximum accuracy of 20sf.


!   INPUT PARAMETER:

!     XVALUE - DOUBLE PRECISION - The argument to the function


!   MACHINE-DEPENDENT PARAMETERS:

!     XLOW - DOUBLE PRECISION - The absolute value below which
!                                   DCOSINT( x ) = gamma + LN(x) ,
!                               to machine precision.
!                               The recommended value is SQRT(2*EPSNEG)

!     XHIGH1 - DOUBLE PRECISION - The value above which
!                                    DCOSINT(x) = sin(x)/x - cos(x)/x^2
!                                 to machine precision.
!                                 The recommended value is SQRT(6/EPSNEG)

!     XHIGH2 - DOUBLE PRECISION - The value above which the trig. functions
!                                 cannot be accurately determined.
!                                 The value of the function is
!                                        DCOSINT(x) = 0.0
!                                 The recommended value is pi/EPS.

!      Values of EPS and EPSNEG for certain machine/compiler
!      combinations can be found in the paper

!      W.J. CODY  Algorithm 665: MACHAR: A subroutine to dynamically
!      determine machine parameters, ACM Trans. Math. Soft. 14 (1988) 303-311.

!      The current code gives numerical values for XLOW,XHIGH1,XHIGH2
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
INTEGER   :: i
REAL (dp) :: cx, dif, fival, gival, logval, root, sum, sumden, sumnum, &
             sx, x, xsq

!   DATA VALUES

REAL (dp), PARAMETER :: zero = 0.0_dp, one = 1.0_dp, three = 3.0_dp, &
                        six = 6.0_dp, twelve = 12.0_dp
REAL (dp), PARAMETER :: gamma = 0.57721566490153286060_dp
REAL (dp), PARAMETER :: logl1 = 0.046875_dp, logl2 = 0.25_dp

!   MACHINE-DEPENDENT PARAMETERS (SUITABLE FOR IEEE MACHINES)

REAL (dp), PARAMETER :: xlow = 1.48996E-8_dp, xhigh1 = 2.324953E8_dp,  &
                        xhigh2 = 1.4148475E16_dp

!  VALUES FOR COS-INTEGRAL FOR 0 < X <= 3

REAL (dp), PARAMETER :: ac1n(0:5) = (/ -0.24607411378767540707_dp,   &
         0.72113492241301534559E-2_dp, -0.11867127836204767056E-3_dp,  &
         0.90542655466969866243E-6_dp, -0.34322242412444409037E-8_dp,  &
         0.51950683460656886834E-11_dp /)
REAL (dp), PARAMETER :: ac1d(0:5) = (/ 1.0_dp, 0.12670095552700637845E-1_dp, &
                 0.78168450570724148921E-4_dp, 0.29959200177005821677E-6_dp, &
                 0.73191677761328838216E-9_dp, 0.94351174530907529061E-12_dp /)

!  VALUES FOR COS-INTEGRAL FOR 3 < X <= 6

REAL (dp), PARAMETER :: ac2n(0:7) = (/ -0.15684781827145408780_dp,     &
        0.66253165609605468916E-2_dp,  -0.12822297297864512864E-3_dp,  &
        0.12360964097729408891E-5_dp,  -0.66450975112876224532E-8_dp,  &
        0.20326936466803159446E-10_dp, -0.33590883135343844613E-13_dp, &
        0.23686934961435015119E-16_dp /)
REAL (dp), PARAMETER :: ac2d(0:6) = (/ 1.0_dp, 0.96166044388828741188E-2_dp,  &
                 0.45257514591257035006E-4_dp, 0.13544922659627723233E-6_dp,  &
                 0.27715365686570002081E-9_dp, 0.37718676301688932926E-12_dp, &
                 0.27706844497155995398E-15_dp /)

!  VALUES FOR FI(X) FOR 6 <= X <= 12

REAL (dp), PARAMETER :: afn1(0:7) = (/ 0.99999999962173909991E0_dp,  &
          0.36451060338631902917E3_dp, 0.44218548041288440874E5_dp,  &
          0.22467569405961151887E7_dp, 0.49315316723035561922E8_dp,  &
          0.43186795279670283193E9_dp, 0.11847992519956804350E10_dp, &
          0.45573267593795103181E9_dp /)
REAL (dp), PARAMETER :: afd1(0:7) = (/ 1.0_dp, 0.36651060273229347594E3_dp,  &
                  0.44927569814970692777E5_dp, 0.23285354882204041700E7_dp,  &
                  0.53117852017228262911E8_dp, 0.50335310667241870372E9_dp,  &
                  0.16575285015623175410E10_dp, 0.11746532837038341076E10_dp /)

!   VALUES OF GI(X) FOR 6 <= X <=12

REAL (dp), PARAMETER :: agn1(0:8) = (/ 0.99999999920484901956E0_dp,  &
         0.51385504875307321394E3_dp,  0.92293483452013810811E5_dp,  &
         0.74071341863359841727E7_dp,  0.28142356162841356551E9_dp,  &
         0.49280890357734623984E10_dp, 0.35524762685554302472E11_dp, &
         0.79194271662085049376E11_dp, 0.17942522624413898907E11_dp /)
REAL (dp), PARAMETER :: agd1(0:8) = (/ 1.0_dp,  0.51985504708814870209E3_dp,  &
                  0.95292615508125947321E5_dp,  0.79215459679762667578E7_dp,  &
                  0.31977567790733781460E9_dp,  0.62273134702439012114E10_dp, &
                  0.54570971054996441467E11_dp, 0.18241750166645704670E12_dp, &
                  0.15407148148861454434E12_dp /)

!   VALUES FOR FI(X) FOR X > 12

REAL (dp), PARAMETER :: afn2(0:7) = (/ 0.19999999999999978257E1_dp,   &
          0.22206119380434958727E4_dp, 0.84749007623988236808E6_dp,   &
          0.13959267954823943232E9_dp, 0.10197205463267975592E11_dp,  &
          0.30229865264524075951E12_dp, 0.27504053804288471142E13_dp, &
          0.21818989704686874983E13_dp /)
REAL (dp), PARAMETER :: afd2(0:7) = (/ 1.0_dp,  0.11223059690217167788E4_dp,  &
                  0.43685270974851313242E6_dp,  0.74654702140658116258E8_dp,  &
                  0.58580034751805687471E10_dp, 0.20157980379272098841E12_dp, &
                  0.26229141857684496445E13_dp, 0.87852907334918467516E13_dp /)

!   VALUES FOR GI(X) FOR X > 12

REAL (dp), PARAMETER :: agn2(0:8) = (/  0.59999999999999993089E1_dp,  &
          0.96527746044997139158E4_dp,  0.56077626996568834185E7_dp,  &
          0.15022667718927317198E10_dp, 0.19644271064733088465E12_dp, &
          0.12191368281163225043E14_dp, 0.31924389898645609533E15_dp, &
          0.25876053010027485934E16_dp, 0.12754978896268878403E16_dp /)
REAL (dp), PARAMETER :: agd2(0:8) = (/ 1.0_dp,  0.16287957674166143196E4_dp,  &
                  0.96636303195787870963E6_dp,  0.26839734750950667021E9_dp,  &
                  0.37388510548029219241E11_dp, 0.26028585666152144496E13_dp,  &
                  0.85134283716950697226E14_dp, 0.11304079361627952930E16_dp,  &
                  0.42519841479489798424E16_dp /)

!   VALUES FOR AN APPROXIMATION TO LN(X/ROOT)

REAL (dp), PARAMETER :: p(0:2) = (/   0.83930008362695945726E1_dp,  &
        -0.65306663899493304675E1_dp, 0.569155722227490223_dp /)
REAL (dp), PARAMETER :: q(0:1) = (/ 0.41965004181347972847E1_dp,  &
                                   -0.46641666676862479585E1_dp /)

!   VALUES OF THE FIRST TWO ROOTS OF THE COSINE-INTEGRAL

REAL (dp), PARAMETER :: rt1n = 631.0_dp, rt1d = 1024.0_dp, &
                        rt1r = 0.29454812071623379711E-3_dp
REAL (dp), PARAMETER :: rt2n = 3465.0_dp, rt2d = 1024.0_dp, &
                        rt2r = 0.39136005118642639785E-3_dp

!   START COMPUTATION

x = xvalue
IF ( x <= zero ) THEN
  fn_val = zero
  RETURN
END IF
IF ( x <= six ) THEN
  
!   CODE FOR 3 < X < =  6
  
  IF ( x > three ) THEN
    sumnum = zero
    sumden = zero
    xsq = x * x
    DO i = 7 , 0 , -1
      sumnum = sumnum * xsq + ac2n( i )
    END DO
    DO i = 6 , 0 , -1
      sumden = sumden * xsq + ac2d( i )
    END DO
    root = rt2n / rt2d
    dif = ( x - root ) - rt2r
    sum = root + rt2r
    IF ( ABS(dif) < logl2 ) THEN
      cx = dif / ( sum + x )
      xsq = cx * cx
      sx = p(0) + xsq * ( p(1) + xsq * p(2) )
      sx = sx / ( q(0) + xsq * ( q(1) + xsq ) )
      logval = cx * sx
    ELSE
      logval = LOG( x / sum )
    END IF
    fn_val = logval + dif * ( x + sum ) * sumnum / sumden
  ELSE
    
!   CODE FOR 0 < X < =  3
    
    IF ( x > xlow ) THEN
      sumnum = zero
      sumden = zero
      xsq = x * x
      DO i = 5 , 0 , -1
        sumnum = sumnum * xsq + ac1n( i )
        sumden = sumden * xsq + ac1d( i )
      END DO
      root = rt1n / rt1d
      dif = ( x - root ) - rt1r
      sum = root + rt1r
      IF ( ABS(dif) < logl1 ) THEN
        cx = dif / ( sum + x )
        xsq = cx * cx
        sx = p(0) + xsq * ( p(1) + xsq * p(2) )
        sx = sx / ( q(0) + xsq * ( q(1) + xsq ) )
        logval = cx * sx
      ELSE
        logval = LOG( x / sum )
      END IF
      fn_val = logval + dif * ( x + sum ) * sumnum / sumden
    ELSE
      fn_val = gamma + LOG( x )
    END IF
  END IF
END IF

!   CODE FOR 6 < X < =  12

IF ( x > six .AND. x <= twelve ) THEN
  sumnum = zero
  sumden = zero
  xsq = one / ( x * x )
  DO i = 7 , 0 , -1
    sumnum = sumnum * xsq + afn1( i )
    sumden = sumden * xsq + afd1( i )
  END DO
  fival = sumnum / ( x * sumden )
  sumnum = zero
  sumden = zero
  DO i = 8 , 0 , -1
    sumnum = sumnum * xsq + agn1( i )
    sumden = sumden * xsq + agd1( i )
  END DO
  gival = xsq * sumnum / sumden
  fn_val = fival * SIN( x ) - gival * COS( x )
END IF

!   CODE FOR X > 12

IF ( x > twelve ) THEN
  IF ( x > xhigh2 ) THEN
    fn_val = zero
  ELSE
    cx = COS( x )
    sx = SIN( x )
    xsq = one / ( x * x )
    IF ( x > xhigh1 ) THEN
      fn_val = sx / x - cx * xsq
    ELSE
      sumnum = zero
      sumden = zero
      DO i = 7 , 0 , -1
        sumnum = sumnum * xsq + afn2( i )
        sumden = sumden * xsq + afd2( i )
      END DO
      fival = ( one - xsq * sumnum / sumden ) / x
      sumnum = zero
      sumden = zero
      DO i = 8 , 0 , -1
        sumnum = sumnum * xsq + agn2( i )
        sumden = sumden * xsq + agd2( i )
      END DO
      gival = ( one - xsq * sumnum / sumden ) * xsq
      fn_val = sx * fival - cx * gival
    END IF
  END IF
END IF

RETURN
END FUNCTION dcosint
