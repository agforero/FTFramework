C$Procedure ZZDSIN ( SGP4 deep space initialization )

      SUBROUTINE ZZDSIN (GEOPHS,
     .                   COSIM,  EMSQ,  ARGPO, S1,    S2,    S3,
     .                   S4,     S5,    SINIM, SS1,   SS2,   SS3,
     .                   SS4,    SS5,   SZ1,   SZ3,   SZ11,  SZ13,
     .                   SZ21,   SZ23,  SZ31,  SZ33,  T,     TC,
     .                   GSTO,   MO,    MDOT,  NO,    NODEO, NODEDOT,
     .                   XPIDOT, Z1,    Z3,    Z11,   Z13,   Z21,
     .                   Z23,    Z31,   Z33,   ECCO,  ECCSQ,
     .                   ECCM,   ARGPM, INCLM, MM,    XN,    NODEM,
     .                   IREZ,   ATIME, D2201, D2211, D3210, D3222,
     .                   D4410,  D4422, D5220, D5232, D5421, D5433,
     .                   DEDT,   DIDT,  DMDT,  DNDT,  DNODT, DOMDT,
     .                   DEL1,   DEL2,  DEL3,  XFACT, XLAMO, XLI,
     .                   XNI )

C$ Abstract
C
C     This Subroutine provides Deep Space contributions to Mean
C     Motion Dot due to geopotential resonance with half day and one
C     day orbits.
C
C$ Disclaimer
C
C     THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE
C     CALIFORNIA INSTITUTE OF TECHNOLOGY (CALTECH) UNDER A U.S.
C     GOVERNMENT CONTRACT WITH THE NATIONAL AERONAUTICS AND SPACE
C     ADMINISTRATION (NASA). THE SOFTWARE IS TECHNOLOGY AND SOFTWARE
C     PUBLICLY AVAILABLE UNDER U.S. EXPORT LAWS AND IS PROVIDED "AS-IS"
C     TO THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING ANY
C     WARRANTIES OF PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR A
C     PARTICULAR USE OR PURPOSE (AS SET FORTH IN UNITED STATES UCC
C     SECTIONS 2312-2313) OR FOR ANY PURPOSE WHATSOEVER, FOR THE
C     SOFTWARE AND RELATED MATERIALS, HOWEVER USED.
C
C     IN NO EVENT SHALL CALTECH, ITS JET PROPULSION LABORATORY, OR NASA
C     BE LIABLE FOR ANY DAMAGES AND/OR COSTS, INCLUDING, BUT NOT
C     LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND,
C     INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY AND LOST PROFITS,
C     REGARDLESS OF WHETHER CALTECH, JPL, OR NASA BE ADVISED, HAVE
C     REASON TO KNOW, OR, IN FACT, SHALL KNOW OF THE POSSIBILITY.
C
C     RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF
C     THE SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY
C     CALTECH AND NASA FOR ALL THIRD-PARTY CLAIMS RESULTING FROM THE
C     ACTIONS OF RECIPIENT IN THE USE OF THE SOFTWARE.
C
C$ Required_Reading
C
C     None.
C
C$ Keywords
C
C     None.
C
C$ Declarations

      IMPLICIT NONE

      INCLUDE 'zzsgp4.inc'

      DOUBLE PRECISION         ARGPM
      DOUBLE PRECISION         ARGPO
      DOUBLE PRECISION         ATIME
      DOUBLE PRECISION         COSIM
      DOUBLE PRECISION         D2201
      DOUBLE PRECISION         D2211
      DOUBLE PRECISION         D3210
      DOUBLE PRECISION         D3222
      DOUBLE PRECISION         D4410
      DOUBLE PRECISION         D4422
      DOUBLE PRECISION         D5220
      DOUBLE PRECISION         D5232
      DOUBLE PRECISION         D5421
      DOUBLE PRECISION         D5433
      DOUBLE PRECISION         DEDT
      DOUBLE PRECISION         DEL1
      DOUBLE PRECISION         DEL2
      DOUBLE PRECISION         DEL3
      DOUBLE PRECISION         DIDT
      DOUBLE PRECISION         DMDT
      DOUBLE PRECISION         DNDT
      DOUBLE PRECISION         DNODT
      DOUBLE PRECISION         DOMDT
      DOUBLE PRECISION         ECCM
      DOUBLE PRECISION         ECCO
      DOUBLE PRECISION         ECCSQ
      DOUBLE PRECISION         EMSQ
      DOUBLE PRECISION         GEOPHS    ( 8 )
      DOUBLE PRECISION         GSTO
      DOUBLE PRECISION         INCLM
      DOUBLE PRECISION         MDOT
      DOUBLE PRECISION         MM
      DOUBLE PRECISION         MO
      DOUBLE PRECISION         NO
      DOUBLE PRECISION         NODEDOT
      DOUBLE PRECISION         NODEM
      DOUBLE PRECISION         NODEO
      DOUBLE PRECISION         S1
      DOUBLE PRECISION         S2
      DOUBLE PRECISION         S3
      DOUBLE PRECISION         S4
      DOUBLE PRECISION         S5
      DOUBLE PRECISION         SINIM
      DOUBLE PRECISION         SS1
      DOUBLE PRECISION         SS2
      DOUBLE PRECISION         SS3
      DOUBLE PRECISION         SS4
      DOUBLE PRECISION         SS5
      DOUBLE PRECISION         SZ1
      DOUBLE PRECISION         SZ11
      DOUBLE PRECISION         SZ13
      DOUBLE PRECISION         SZ21
      DOUBLE PRECISION         SZ23
      DOUBLE PRECISION         SZ3
      DOUBLE PRECISION         SZ31
      DOUBLE PRECISION         SZ33
      DOUBLE PRECISION         T
      DOUBLE PRECISION         TC
      DOUBLE PRECISION         XFACT
      DOUBLE PRECISION         XLAMO
      DOUBLE PRECISION         XLI
      DOUBLE PRECISION         XN
      DOUBLE PRECISION         XNI
      DOUBLE PRECISION         XPIDOT
      DOUBLE PRECISION         Z1
      DOUBLE PRECISION         Z11
      DOUBLE PRECISION         Z13
      DOUBLE PRECISION         Z21
      DOUBLE PRECISION         Z23
      DOUBLE PRECISION         Z3
      DOUBLE PRECISION         Z31
      DOUBLE PRECISION         Z33
    
      INTEGER                  IREZ

C$ Brief_I/O
C
C    See Detailed_input and Detailed_Output.
C
C$ Detailed_Input
C
C    COSIM       COS of mean inclination.
C
C    SINIM       SIN of mean inclination.   
C
C    EMSQ        Eccentricity squared
C
C    ARGPO       Argument of Perigee
C
C    S1          S coefficients
C
C    S2              ...
C
C    S3              ...
C
C    S4              ...
C
C    S5              ...
C
C    SS1             ...
C
C    SS2             ...
C
C    SS3             ...
C
C    SS4             ...
C
C    SS5             ...
C
C    SZ1             ...
C
C    SZ3             ...
C
C    SZ11            ...
C
C    SZ13            ...
C
C    SZ21            ...
C
C    SZ23            ...
C
C    SZ31            ...
C
C    SZ33            ...
C
C    T           Time
C
C    TC
C
C    GSTO        Greenwich sidereal time in radians
C
C    MO          Mean anomaly
C
C    MDOT        Mean anomaly rate
C
C    NO          Mean motion
C
C    NODEO       Right ascension of ascending node
C
C    NODEDOT     Right ascension of ascending node rate
C
C    XPIDOT
C
C    Z1          Z coefficients
C
C    Z3               ...
C
C    Z11              ...
C
C    Z13              ...
C
C    Z21              ...
C
C    Z23              ...
C
C    Z31              ...
C
C    Z33              ...
C
C    ECCM        Mean eccentricity
C
C    ARGPM       Mean argument of perigee
C
C    INCLM       Mean inclination
C
C    MM          Mean anomaly
C
C    XN          Mean motion
C
C    NODEM       Mean right ascension of ascending node
C
C$ Detailed_Output
C
C    ECCM        Mean eccentricity
C
C    ARGPM       Mean argument of perigee
C
C    INCLM       Mean inclination
C
C    MM          Mean anomaly
C
C    XN          Mean motion
C
C    NODEM       Right ascension of ascending node
C
C    IREZ        Resonance flags: 0-none, 1-one day, 2-half day
C
C    ATIME       Internal SGD4 parameter.
C
C    D2201       D COEFFCIENTS
C
C    D2211           ...
C
C    D3210           ...
C
C    D3222           ...
C
C    D4410           ...
C
C    D4422           ...
C
C    D5220           ...
C
C    D5232           ...
C
C    D5421           ...
C
C    D5433           ...
C
C    DEDT        Internal SGD4 parameter.
C
C    DIDT        Internal SGD4 parameter.
C
C    DMDT        Internal SGD4 parameter.
C
C    DNDT        Internal SGD4 parameter.
C
C    DNODT       Internal SGD4 parameter.
C
C    DOMDT       Internal SGD4 parameter.
C
C    DEL1        Internal SGD4 parameter.
C
C    DEL2        Internal SGD4 parameter.
C
C    DEL3        Internal SGD4 parameter.
C
C    SES         Internal SGD4 parameter.
C
C    SGHL        Internal SGD4 parameter.
C
C    SGHS        Internal SGD4 parameter.
C
C    SGS         Internal SGD4 parameter.
C
C    SHL         Internal SGD4 parameter.
C
C    SHS         Internal SGD4 parameter.
C
C    SIS         Internal SGD4 parameter.
C
C    SLS         Internal SGD4 parameter.
C
C    THETA       Internal SGD4 parameter.
C
C    XFACT       Internal SGD4 parameter.
C
C    XLAMO       Internal SGD4 parameter.
C
C    XLI         Internal SGD4 parameter.
C
C    XNI         Internal SGD4 parameter.
C
C$ Parameters
C
C     None.
C
C$ Exceptions
C
C     None.
C
C$ Files
C
C     None.
C
C$ Particulars
C
C     This routine is based on the DSINIT code by David Vallado
C     corresponding to "Revisiting Spacetrack Report #3" [4].
C     The intent is to maintain the original Vallado algorithm,
C     changing code only to meet NAIF format standards and to
C     integrate with SPICELIB.
C
C        Removed getgravconst call, replaced with GEOPHS array.
C
C        Capitalize all variables.
C
C        ENDIF replaced with END IF.
C
C        whichconst    eliminated, function provided by GEOPHS
C
C$ Examples
C
C     None.
C
C$ Restrictions
C
C     None.
C
C$ Literature_References
C
C   [1] Hoots, F. R., and Roehrich, R. L. 1980. "Models for
C       Propagation of the NORAD Element Sets." Spacetrack Report #3.
C       U.S. Air Force: Aerospace Defense Command.
C
C   [2] Hoots, Felix R. "Spacetrack Report #6: Models for Propagation
C       of Space Command Element Sets." Space Command,
C       U. S. Air Force, CO.
C
C   [3] Hoots, Felix R., P. W. Schumacher, and R. A. Glover. 2004.
C       History of Analytical Orbit Modeling in the U. S. Space
C       Surveillance System. Journal of Guidance, Control, and
C       Dynamics. 27(2):174-185.
C
C   [4] Vallado, David, Crawford, Paul, Hujsak, Richard,
C       and Kelso, T.S. 2006. Revisiting Spacetrack Report #3. Paper
C       AIAA 2006-6753 presented at the AIAA/AAS Astrodynamics
C       Specialist Conference, August 21-24, 2006. Keystone, CO.
C
C$ Author_and_Institution
C
C     David Vallado   (AGI)
C     E. D. Wright    (JPL)
C
C$ Version
C
C-    SPICELIB Version 1.0.0, 03-NOV-2014 (EDW)
C
C        Based on routine DSINIT, 28-JUN-2005, Vallado 2006 [4].
C
C-&

C$ Index_Entries
C
C   SGP4
C
C-&

C
C     Local Variables
C
      DOUBLE PRECISION         AINV2
      DOUBLE PRECISION         AONV
      DOUBLE PRECISION         COSISQ
      DOUBLE PRECISION         EMO
      DOUBLE PRECISION         EMSQO
      DOUBLE PRECISION         EOC
      DOUBLE PRECISION         F220
      DOUBLE PRECISION         F221
      DOUBLE PRECISION         F311
      DOUBLE PRECISION         F321
      DOUBLE PRECISION         F322
      DOUBLE PRECISION         F330
      DOUBLE PRECISION         F441
      DOUBLE PRECISION         F442
      DOUBLE PRECISION         F522
      DOUBLE PRECISION         F523
      DOUBLE PRECISION         F542
      DOUBLE PRECISION         F543
      DOUBLE PRECISION         G200
      DOUBLE PRECISION         G201
      DOUBLE PRECISION         G211
      DOUBLE PRECISION         G300
      DOUBLE PRECISION         G310
      DOUBLE PRECISION         G322
      DOUBLE PRECISION         G410
      DOUBLE PRECISION         G422
      DOUBLE PRECISION         G520
      DOUBLE PRECISION         G521
      DOUBLE PRECISION         G532
      DOUBLE PRECISION         G533
      DOUBLE PRECISION         Q22
      DOUBLE PRECISION         Q31
      DOUBLE PRECISION         Q33
      DOUBLE PRECISION         ROOT22
      DOUBLE PRECISION         ROOT32
      DOUBLE PRECISION         ROOT44
      DOUBLE PRECISION         ROOT52
      DOUBLE PRECISION         ROOT54
      DOUBLE PRECISION         RPTIM
      DOUBLE PRECISION         SES
      DOUBLE PRECISION         SGHL
      DOUBLE PRECISION         SGHS
      DOUBLE PRECISION         SGS
      DOUBLE PRECISION         SHL
      DOUBLE PRECISION         SHS
      DOUBLE PRECISION         SINI2
      DOUBLE PRECISION         SIS
      DOUBLE PRECISION         SLS
      DOUBLE PRECISION         TEMP
      DOUBLE PRECISION         TEMP1
      DOUBLE PRECISION         THETA
      DOUBLE PRECISION         X2O3
      DOUBLE PRECISION         XKE
      DOUBLE PRECISION         XNO2
      DOUBLE PRECISION         ZNL
      DOUBLE PRECISION         ZNS



C
C     SPICELIB routines.
C
      DOUBLE PRECISION         PI
      DOUBLE PRECISION         TWOPI
      LOGICAL                  RETURN


C
C     Standard SPICE error handling.
C
      IF ( RETURN() ) THEN
         RETURN
      END IF

      CALL CHKIN ( 'ZZDSIN' )

C
C     Constants
C
      Q22    = 1.7891679D-6
      Q31    = 2.1460748D-6
      Q33    = 2.2123015D-7
      ROOT22 = 1.7891679D-6
      ROOT44 = 7.3636953D-9
      ROOT54 = 2.1765803D-9

C
C     This equates to 7.29211514668855e-5 rad/sec
C
      RPTim  = 4.37526908801129966D-3

      Root32 = 3.7393792D-7
      Root52 = 1.1428639D-7
      X2o3   = 2.0D0 / 3.0D0
      ZNL    = 1.5835218D-4
      ZNS    = 1.19459D-5

C
C     This code block replaces the call:
C
C     sgp4fix identify constants and allow alternate values.
C
C     CALL getgravconst( whichconst, tumin,
C     .                  mu, radiusearthkm, xke,
C     .                  j2, j3, j4, j3oj2 )
C
      XKE   = GEOPHS(K_KE)

C
C     DEEP SPACE INITIALIZATION
C
      IREZ = 0
      IF ( (XN.LT.0.0052359877D0)  .AND.
     .     (XN.GT.0.0034906585D0) ) THEN
         IREZ = 1
      END IF

      IF ( (XN   .GE. 8.26D-3) .AND.
     .     (XN   .LE. 9.24D-3) .AND.
     .     (ECCM .GE. 0.5D0)  )THEN
         IREZ = 2
      END IF

C
C     DO SOLAR TERMS
C
      SES  =  SS1*ZNS*SS5
      SIS  =  SS2*ZNS*(SZ11 + SZ13)
      SLS  = -ZNS*SS3*(SZ1 + SZ3 - 14.0D0 - 6.0D0*EMSQ)
      SGHS =  SS4*ZNS*(SZ31 + SZ33 - 6.0D0)
      SHS  = -ZNS*SS2*(SZ21 + SZ23)

C
C       sgp4fix for 180 deg incl
C
        IF ( (INCLM .LT. 5.2359877D-2) .OR.
     .       (INCLM .GT. PI()-5.2359877D-2) ) THEN
            SHS = 0.0D0
        END IF

        IF (SINIM.NE.0.0D0) THEN
            SHS = SHS/SINIM
        END IF

        SGS  = SGHS - COSIM*SHS

C
C       DO LUNAR TERMS
C
        DEDT = SES + S1*ZNL*S5
        DIDT = SIS + S2*ZNL*(Z11 + Z13)
        DMDT = SLS - ZNL*S3*(Z1 + Z3 - 14.0D0 - 6.0D0*EMSQ)
        SGHL = S4*ZNL*(Z31 + Z33 - 6.0D0)
        SHL  = -ZNL*S2*(Z21 + Z23)

C
C       sgp4fix for 180 deg incl
C

        IF ( (INCLM .LT. 5.2359877D-2 ) .OR.
     .       (INCLM .GT. PI()-5.2359877D-2)) THEN
            SHL = 0.0D0
        END IF

        DOMDT= SGS+SGHL
        DNODT= SHS

        IF (SINIM .NE. 0.0D0) THEN
            DOMDT = DOMDT-COSIM/SINIM*SHL
            DNODT = DNODT+SHL/SINIM
        END IF

C
C       CALCULATE DEEP SPACE RESONANCE EFFECTS
C
        DNDT  = 0.0D0
        THETA = DMOD(GSTO + TC*RPTIM, TWOPI() )
        ECCM  = ECCM + DEDT*T
        EMSQ  = ECCM**2
        INCLM = INCLM + DIDT*T
        ARGPM = ARGPM + DOMDT*T
        NODEM = NODEM + DNODT*T
        MM    = MM + DMDT*T

C
C   sgp4fix for negative inclinations
C   the following if statement should be commented out
C
C           IF(Inclm .lt. 0.0D0) THEN
C             Inclm  = -Inclm
C             Argpm  = Argpm-PI
C             nodem = nodem+PI
C           END IF
C

C
C       Initialize the resonance terms
C
         IF (IREZ .NE. 0) THEN

            AONV = (XN/XKE)**X2O3

C
C           GEOPOTENTIAL RESONANCE FOR 12 HOUR ORBITS
C
            IF (IREZ .EQ. 2) THEN

               COSISQ = COSIM*COSIM
               EMO    = ECCM
               EMSQO  = EMSQ
               ECCM   = ECCO
               EMSQ   = ECCSQ
               EOC    = ECCM*EMSQ
               G201   = -0.306D0-(ECCM-0.64D0)*0.440D0

               IF (ECCM.LE.0.65D0) THEN

                  G211 =   3.616D0 -  13.2470D0*ECCM +
     .                        16.2900D0*EMSQ
                  G310 = -19.302D0 + 117.3900D0*ECCM -
     .                        228.4190D0*EMSQ + 156.591D0*EOC
                  G322 = -18.9068D0+ 109.7927D0*ECCM -
     .                        214.6334D0*EMSQ + 146.5816D0*EOC
                  G410 = -41.122D0 + 242.6940D0*ECCM -
     .                        471.0940D0*EMSQ + 313.953D0*EOC
                  G422 =-146.407D0 + 841.8800D0*ECCM -
     .                        1629.014D0*EMSQ + 1083.435D0*EOC
                  G520 =-532.114D0 + 3017.977D0*ECCM -
     .                        5740.032D0*EMSQ + 3708.276D0*EOC

               ELSE

                  G211 =  -72.099D0 +  331.819D0*ECCM -
     .                        508.738D0*EMSQ + 266.724D0*EOC
                  G310 = -346.844D0 + 1582.851D0*ECCM -
     .                        2415.925D0*EMSQ + 1246.113D0*EOC
                  G322 = -342.585D0 + 1554.908D0*ECCM -
     .                        2366.899D0*EMSQ + 1215.972D0*EOC
                  G410 =-1052.797D0 + 4758.686D0*ECCM -
     .                        7193.992D0*EMSQ + 3651.957D0*EOC
                  G422 =-3581.690D0 + 16178.11D0*ECCM -
     .                        24462.77D0*EMSQ + 12422.52D0*EOC

                  IF (ECCM.GT.0.715D0) THEN
                     G520 =-5149.66D0 + 29936.92D0*ECCM -
     .                        54087.36D0*EMSQ + 31324.56D0*EOC
                  ELSE
                     G520 = 1464.74D0 -  4664.75D0*ECCM +
     .                        3763.64D0*EMSQ
                  END IF

               END IF

            IF (ECCM.LT.0.7D0) THEN

                G533 = -919.22770D0 + 4988.6100D0*ECCM-9064.7700D0*EMSQ
     .               + 5542.21D0*EOC
                G521 = -822.71072D0 + 4568.6173D0*ECCM-8491.4146D0*EMSQ
     .               + 5337.524D0*EOC
                G532 = -853.66600D0 + 4690.2500D0*ECCM-8624.7700D0*EMSQ
     .               + 5341.4D0*EOC

            ELSE

                G533 =-37995.780D0 + 161616.52D0*ECCM-229838.20D0*EMSQ+
     .              109377.94D0*EOC
                G521 =-51752.104D0 + 218913.95D0*ECCM-309468.16D0*EMSQ+
     .              146349.42D0*EOC
                G532 =-40023.880D0 + 170470.89D0*ECCM-242699.48D0*EMSQ+
     .              115605.82D0*EOC

            END IF

            SINI2 =  SINIM*SINIM
            F220  =  0.75D0* (1.0D0+2.0D0*COSIM+COSISQ)
            F221  =  1.5D0*SINI2
            F321  =  1.875D0*SINIM * (1.0D0-2.0D0*COSIM-3.0D0*COSISQ)
            F322  = -1.875D0*SINIM * (1.0D0+2.0D0*COSIM-3.0D0*COSISQ)
            F441  = 35.0D0*SINI2*F220
            F442  = 39.3750D0*SINI2*SINI2
            F522  =  9.84375D0*SINIM * (SINI2* (1.0D0-2.0D0*COSIM-
     .               5.0D0*COSISQ)+0.33333333D0 * (-2.0D0+4.0D0*COSIM+
     .               6.0D0*COSISQ) )
            F523  =  SINIM * (4.92187512D0*SINI2 * (-2.0D0-4.0D0*COSIM+
     .               10.0D0*COSISQ) + 6.56250012D0*
     .               (1.0D0+2.0D0*COSIM-3.0D0*COSISQ))
            F542  =  29.53125D0*SINIM * (2.0D0-8.0D0*COSIM+COSISQ*
     .               (-12.0D0+8.0D0*COSIM+10.0D0*COSISQ) )
            F543  = 29.53125D0*SINIM * (-2.0D0-8.0D0*COSIM+COSISQ*
     .               (12.0D0+8.0D0*COSIM-10.0D0*COSISQ) )

            XNO2   =  XN * XN
            AINV2  =  AONV * AONV
            TEMP1  =  3.0D0*XNO2*AINV2
            TEMP   =  TEMP1*ROOT22
            D2201  =  TEMP*F220*G201
            D2211  =  TEMP*F221*G211
            TEMP1  =  TEMP1*AONV
            TEMP   =  TEMP1*ROOT32
            D3210  =  TEMP*F321*G310
            D3222  =  TEMP*F322*G322
            TEMP1  =  TEMP1*AONV
            TEMP   =  2.0D0*TEMP1*ROOT44
            D4410  =  TEMP*F441*G410
            D4422  =  TEMP*F442*G422
            TEMP1  =  TEMP1*AONV
            TEMP   =  TEMP1*ROOT52
            D5220  =  TEMP*F522*G520
            D5232  =  TEMP*F523*G532
            TEMP   =  2.0D0*TEMP1*ROOT54
            D5421  =  TEMP*F542*G521
            D5433  =  TEMP*F543*G533
            XLAMO  =  DMOD(MO+NODEO+NODEO-THETA-THETA, TWOPI())
            XFACT  = MDOT + DMDT + 2.0D0 * (NODEDOT+DNODT-RPTIM) - NO

            ECCM = EMO
            EMSQ = EMSQO

         END IF

         IF (IREZ .EQ. 1) THEN

C
C           SYNCHRONOUS RESONANCE TERMS
C
            G200  = 1.0D0 + EMSQ * (-2.5D0+0.8125D0*EMSQ)
            G310  = 1.0D0 + 2.0D0*EMSQ
            G300  = 1.0D0 + EMSQ * (-6.0D0+6.60937D0*EMSQ)
            F220  = 0.75D0 * (1.0D0+COSIM) * (1.0D0+COSIM)
            F311  = 0.9375D0*SINIM*SINIM*
     .               (1.0D0+3.0D0*COSIM) - 0.75D0*(1.0D0+COSIM)
            F330  = 1.0D0+COSIM
            F330  = 1.875D0*F330*F330*F330
            DEL1  = 3.0D0*XN*XN*AONV*AONV
            DEL2  = 2.0D0*DEL1*F220*G200*Q22
            DEL3  = 3.0D0*DEL1*F330*G300*Q33*AONV
            DEL1  = DEL1*F311*G310*Q31*AONV
            XLAMO = DMOD(MO+NODEO+ARGPO-THETA, TWOPI())
            XFACT = MDOT + XPIDOT - RPTIM + DMDT + DOMDT + DNODT - NO

         END IF

C
C        FOR SGP4, INITIALIZE THE INTEGRATOR
C
         XLI   = XLAMO
         XNI   = NO
         ATIME = 0.0D0
         XN    = NO + DNDT

      END IF

      CALL CHKOUT ( 'ZZDSIN' )

      RETURN
      END
