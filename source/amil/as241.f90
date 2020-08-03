PROGRAM test_as241
IMPLICIT NONE
INTEGER, PARAMETER :: sp = SELECTED_REAL_KIND(6, 30),  &
                      dp = SELECTED_REAL_KIND(12, 60)

REAL (sp)  :: p_single, z_single
REAL (dp)  :: p_dble, z_dble
INTEGER    :: ifault

INTERFACE
  SUBROUTINE ppnd7 (p, normal_dev, ifault)
    IMPLICIT NONE
    INTEGER, PARAMETER      :: sp = SELECTED_REAL_KIND(6, 30)
    REAL (sp), INTENT(IN)   :: p
    INTEGER, INTENT(OUT)    :: ifault
    REAL (sp), INTENT(OUT)  :: normal_dev
  END SUBROUTINE ppnd7
  SUBROUTINE ppnd16 (p, normal_dev, ifault)
    IMPLICIT NONE
    INTEGER, PARAMETER      :: dp = SELECTED_REAL_KIND(12, 60)
    REAL (dp), INTENT(IN)   :: p
    INTEGER, INTENT(OUT)    :: ifault
    REAL (dp), INTENT(OUT)  :: normal_dev
  END SUBROUTINE ppnd16
END INTERFACE

DO
  WRITE(*, *)'Enter area under normal curve: '
  READ(*, *) p_dble
  p_single = p_dble
  CALL ppnd7(p_single, z_single, ifault)
  IF (ifault /= 0) THEN
    WRITE(*, '(1x, "IFAULT =", i3)') ifault
    CYCLE
  END IF
  CALL ppnd16(p_dble, z_dble, ifault)
  IF (ifault /= 0) THEN
    WRITE(*, '(1x, "IFAULT =", i3)') ifault
    CYCLE
  END IF
  WRITE(*, '(1x, a, f11.6, 2x, a, f16.11)')   &
        'Low prec. result =',  z_single, ' High prec. result =', z_dble
END DO
STOP

END PROGRAM test_as241



SUBROUTINE ppnd7 (p, normal_dev, ifault)

! ALGORITHM AS241  APPL. STATIST. (1988) VOL. 37, NO. 3, 477- 484.

! Produces the normal deviate Z corresponding to a given lower tail area of P;
! Z is accurate to about 1 part in 10**7.

! The hash sums below are the sums of the mantissas of the coefficients.
! They are included for use in checking transcription.

! This ELF90-compatible version by Alan Miller - 20 August 1996
! N.B. The original algorithm is as a function; this is a subroutine

IMPLICIT NONE

INTEGER, PARAMETER      :: sp = SELECTED_REAL_KIND(6, 30)
REAL (sp), INTENT(IN)   :: p
INTEGER, INTENT(OUT)    :: ifault
REAL (sp), INTENT(OUT)  :: normal_dev

! Local variables

REAL (sp) :: zero = 0.0, one = 1.0, half = 0.5, split1 = 0.425,  &
             split2 = 5.0, const1 = 0.180625, const2 = 1.6, q, r

! Coefficients for P close to 0.5

REAL (sp) :: a0 = 3.3871327179E+00, a1 = 5.0434271938E+01, &
                   a2 = 1.5929113202E+02, a3 = 5.9109374720E+01, &
                   b1 = 1.7895169469E+01, b2 = 7.8757757664E+01, &
                   b3 = 6.7187563600E+01
! HASH SUM AB          32.3184577772

! Coefficients for P not close to 0, 0.5 or 1.

REAL (sp) :: c0 = 1.4234372777E+00, c1 = 2.7568153900E+00, &
             c2 = 1.3067284816E+00, c3 = 1.7023821103E-01, &
             d1 = 7.3700164250E-01, d2 = 1.2021132975E-01
! HASH SUM CD    15.7614929821

! Coefficients for P near 0 or 1.

REAL (sp) :: e0 = 6.6579051150E+00, e1 = 3.0812263860E+00, &
             e2 = 4.2868294337E-01, e3 = 1.7337203997E-02, &
             f1 = 2.4197894225E-01, f2 = 1.2258202635E-02
! HASH SUM EF    19.4052910204

ifault = 0
q = p - half
IF (ABS(q) <= split1) THEN
  r = const1 - q * q
  normal_dev = q * (((a3 * r + a2) * r + a1) * r + a0) / &
               (((b3 * r + b2) * r + b1) * r + one)
  RETURN
ELSE
  IF (q < zero) THEN
    r = p
  ELSE
    r = one - p
  END IF
  IF (r <= zero) THEN
    ifault = 1
    normal_dev = zero
    RETURN
  END IF
  r = SQRT(-LOG(r))
  IF (r <= split2) THEN
    r = r - const2
    normal_dev = (((c3 * r + c2) * r + c1) * r + c0) / ((d2 * r + d1) * r + one)
  ELSE
    r = r - split2
    normal_dev = (((e3 * r + e2) * r + e1) * r + e0) / ((f2 * r + f1) * r + one)
  END IF
  IF (q < zero) normal_dev = - normal_dev
  RETURN
END IF
END SUBROUTINE ppnd7



SUBROUTINE ppnd16 (p, normal_dev, ifault)

! ALGORITHM AS241  APPL. STATIST. (1988) VOL. 37, NO. 3

! Produces the normal deviate Z corresponding to a given lower
! tail area of P; Z is accurate to about 1 part in 10**16.

! The hash sums below are the sums of the mantissas of the
! coefficients.   They are included for use in checking
! transcription.

! This ELF90-compatible version by Alan Miller - 20 August 1996
! N.B. The original algorithm is as a function; this is a subroutine

IMPLICIT NONE

INTEGER, PARAMETER      :: dp = SELECTED_REAL_KIND(12, 60)
REAL (dp), INTENT(IN)   :: p
INTEGER, INTENT(OUT)    :: ifault
REAL (dp), INTENT(OUT)  :: normal_dev

! Local variables

REAL (dp) :: zero = 0.d0, one = 1.d0, half = 0.5d0,  &
             split1 = 0.425d0, split2 = 5.d0, const1 = 0.180625d0, &
             const2 = 1.6d0, q, r

! Coefficients for P close to 0.5

REAL (dp) :: a0 = 3.3871328727963666080D0, &
             a1 = 1.3314166789178437745D+2, &
             a2 = 1.9715909503065514427D+3, &
             a3 = 1.3731693765509461125D+4, &
             a4 = 4.5921953931549871457D+4, &
             a5 = 6.7265770927008700853D+4, &
             a6 = 3.3430575583588128105D+4, &
             a7 = 2.5090809287301226727D+3, &
             b1 = 4.2313330701600911252D+1, &
             b2 = 6.8718700749205790830D+2, &
             b3 = 5.3941960214247511077D+3, &
             b4 = 2.1213794301586595867D+4, &
             b5 = 3.9307895800092710610D+4, &
             b6 = 2.8729085735721942674D+4, &
             b7 = 5.2264952788528545610D+3
! HASH SUM AB    55.8831928806149014439

! Coefficients for P not close to 0, 0.5 or 1.

REAL (dp) :: c0 = 1.42343711074968357734D0, &
             c1 = 4.63033784615654529590D0, &
             c2 = 5.76949722146069140550D0, &
             c3 = 3.64784832476320460504D0, &
             c4 = 1.27045825245236838258D0, &
             c5 = 2.41780725177450611770D-1, &
             c6 = 2.27238449892691845833D-2, &
             c7 = 7.74545014278341407640D-4, &
             d1 = 2.05319162663775882187D0, &
             d2 = 1.67638483018380384940D0, &
             d3 = 6.89767334985100004550D-1, &
             d4 = 1.48103976427480074590D-1, &
             d5 = 1.51986665636164571966D-2, &
             d6 = 5.47593808499534494600D-4, &
             d7 = 1.05075007164441684324D-9
! HASH SUM CD    49.33206503301610289036

! Coefficients for P near 0 or 1.

REAL (dp) :: e0 = 6.65790464350110377720D0, &
             e1 = 5.46378491116411436990D0, &
             e2 = 1.78482653991729133580D0, &
             e3 = 2.96560571828504891230D-1, &
             e4 = 2.65321895265761230930D-2, &
             e5 = 1.24266094738807843860D-3, &
             e6 = 2.71155556874348757815D-5, &
             e7 = 2.01033439929228813265D-7, &
             f1 = 5.99832206555887937690D-1, &
             f2 = 1.36929880922735805310D-1, &
             f3 = 1.48753612908506148525D-2, &
             f4 = 7.86869131145613259100D-4, &
             f5 = 1.84631831751005468180D-5, &
             f6 = 1.42151175831644588870D-7, &
             f7 = 2.04426310338993978564D-15
! HASH SUM EF    47.52583317549289671629

ifault = 0
q = p - half
IF (ABS(q) <= split1) THEN
  r = const1 - q * q
  normal_dev = q * (((((((a7*r + a6)*r + a5)*r + a4)*r + a3)*r + a2)*r + a1)*r + a0) / &
           (((((((b7*r + b6)*r + b5)*r + b4)*r + b3)*r + b2)*r + b1)*r + one)
  RETURN
ELSE
  IF (q < zero) THEN
    r = p
  ELSE
    r = one - p
  END IF
  IF (r <= zero) THEN
    ifault = 1
    normal_dev = zero
    RETURN
  END IF
  r = SQRT(-LOG(r))
  IF (r <= split2) THEN
    r = r - const2
    normal_dev = (((((((c7*r + c6)*r + c5)*r + c4)*r + c3)*r + c2)*r + c1)*r + c0) / &
             (((((((d7*r + d6)*r + d5)*r + d4)*r + d3)*r + d2)*r + d1)*r + one)
  ELSE
    r = r - split2
    normal_dev = (((((((e7*r + e6)*r + e5)*r + e4)*r + e3)*r + e2)*r + e1)*r + e0) / &
             (((((((f7*r + f6)*r + f5)*r + f4)*r + f3)*r + f2)*r + f1)*r + one)
  END IF
  IF (q < zero) normal_dev = - normal_dev
  RETURN
END IF
END SUBROUTINE ppnd16
